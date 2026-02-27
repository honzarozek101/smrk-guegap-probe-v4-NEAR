// ══════════════════════════════════════════════════════════════════════════════
// SMRK-GUEGAP Probe v4 — NEAR Confirmation Tier
// Copyright (c) 2026 101research.group  — All Rights Reserved.
// probe_version : smrk-guegap-near-v1
// N             : 32  (confirmation; spec N=512 → scaled for gas budget)
// Role          : Second-layer confirmation of ICP screening results
//
// NEAR gas budget: 300 TGas / tx.
//   Jacobi N=32, 300 sweeps ≈ O(32³·300) ≈ 10 M f64 ops ≈ ~30 TGas → safe.
//
// Architecture: single contract (registry + compute unified).
//   Cross-contract calls on NEAR use async promise chains; for this use-case
//   a single-contract design is simpler, auditable, and avoids callback debt.
//
// Kaspa anchoring: commitment is computed on-chain (deterministic SHA-256).
//   HTTP verification is off-chain — NEAR contracts cannot make HTTP calls.
//   The contract stores txid and a caller-supplied verification note.
//
// Spec refs: §2 (operator), §3 (observables), §4 (acceptance), §D.1.3 (anchor)
// ══════════════════════════════════════════════════════════════════════════════

use near_sdk::borsh::{BorshDeserialize, BorshSerialize};
use near_sdk::serde::{Deserialize, Serialize};
use near_sdk::store::UnorderedMap;
use near_sdk::{env, near_bindgen, AccountId, PanicOnDefault};
use sha2::{Digest, Sha256};
use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// SMRK operator constants
// ══════════════════════════════════════════════════════════════════════════════

const NMAT: usize = 32;
const M_SHIFTS: usize = 8;
const PRIMES: [usize; 8] = [2, 3, 5, 7, 11, 13, 17, 19];
const GOLDEN: f64 = 0.6180339887498949; // α = (√5−1)/2  golden gauge
const S_EXP: f64  = 1.0;               // prime-decay exponent s
const ETA: f64    = 0.5;               // Liouville modulation strength η

// Acceptance thresholds — tighter than N=16 (more eigenvalues → better stats)
// H1: |⟨r⟩ − r_GUE| ≤ 0.10   (r_GUE ≈ 0.5996)
// H2: Δ1 ≥ 1×10⁻⁴
const H1_BAND: f64  = 0.10;
const R_GUE: f64    = 0.5996;
const H2_MIN_GAP: f64 = 1e-4;

// Bulk window: β₀ = 0.2, β₁ = 0.8
// i₀ = ⌊0.2·32⌋ = 6,   i₁ = ⌊0.8·32⌋ − 2 = 23   → ~18 ratio samples
const BULK_I0: usize = 6;
const BULK_I1: usize = 23;

// ══════════════════════════════════════════════════════════════════════════════
// Complex arithmetic  (stack-allocated, no alloc)
// ══════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy)]
struct C64 { re: f64, im: f64 }

impl C64 {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };
    fn conj(self) -> Self { Self { re: self.re, im: -self.im } }
    fn norm(self) -> f64 { (self.re * self.re + self.im * self.im).sqrt() }
    fn scale(self, s: f64) -> Self { Self { re: self.re * s, im: self.im * s } }
    fn arg(self) -> f64 { self.im.atan2(self.re) }
    fn polar(r: f64, theta: f64) -> Self { Self { re: r * theta.cos(), im: r * theta.sin() } }
}
impl std::ops::Add for C64 {
    type Output = Self;
    fn add(self, o: Self) -> Self { Self { re: self.re + o.re, im: self.im + o.im } }
}
impl std::ops::Mul for C64 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// Liouville function  λ_L(n) = (−1)^Ω(n)
// ══════════════════════════════════════════════════════════════════════════════

fn liouville(n: usize) -> f64 {
    if n == 0 { return 1.0; }
    let (mut x, mut cnt, mut d) = (n, 0usize, 2usize);
    while d * d <= x { while x % d == 0 { cnt += 1; x /= d; } d += 1; }
    if x > 1 { cnt += 1; }
    if cnt % 2 == 0 { 1.0 } else { -1.0 }
}

// ══════════════════════════════════════════════════════════════════════════════
// SMRK operator  H(N)(θ) = D(N)(θ) + C(N)(θ)  — spec §2
// ══════════════════════════════════════════════════════════════════════════════

type Matrix = [[C64; NMAT]; NMAT];

fn build_smrk_matrix(kappa: f64, lam: f64, phi_offset: f64) -> Matrix {
    let mut h = [[C64::ZERO; NMAT]; NMAT];

    // Diagonal  D_nn = λ · ln(n+1),  n 1-indexed → 0-indexed: ln(n+2)
    for n in 0..NMAT {
        h[n][n] = C64 { re: lam * ((n + 2) as f64).ln(), im: 0.0 };
    }

    // Off-diagonal prime-shift couplings (§2.3–2.5)
    for n in 0..NMAT {
        let n1 = (n + 1) as f64;
        let lv = liouville(n + 1);

        for m in 0..M_SHIFTS {
            let pm    = PRIMES[m];
            let sigma = (n + pm) % NMAT;

            let amp = kappa / (1.0 + pm as f64).ln().powf(S_EXP) * (1.0 + ETA * lv);
            let phi = 2.0 * PI * GOLDEN * n1 + phi_offset;
            let w   = C64::polar(amp, phi);

            h[n][sigma] = h[n][sigma] + w;
            h[sigma][n] = h[sigma][n] + w.conj();
        }
    }
    h
}

// ══════════════════════════════════════════════════════════════════════════════
// Cyclic complex Jacobi eigenvalue solver  (§6.12)
//
// Complexity: O(N³ · sweeps) ≈ 32³ · 300 ≈ 10 M ops  → fits NEAR gas budget.
// Tolerance: 1×10⁻¹¹ on max off-diagonal norm.
// ══════════════════════════════════════════════════════════════════════════════

fn jacobi_eigen(h_in: &Matrix) -> [f64; NMAT] {
    let mut a = *h_in;
    let tol = 1e-11_f64;

    'outer: for _ in 0..300 {
        // Convergence check
        let mut max_off = 0.0_f64;
        for i in 0..NMAT {
            for j in 0..NMAT {
                if i != j { let v = a[i][j].norm(); if v > max_off { max_off = v; } }
            }
        }
        if max_off < tol { break 'outer; }

        // Cyclic sweep over all upper-triangular (p, q) pairs
        for p in 0..(NMAT - 1) {
            for q in (p + 1)..NMAT {
                let apq = a[p][q];
                let r   = apq.norm();
                if r < 1e-14 { continue; }

                let phi = apq.arg();
                let tau = (a[q][q].re - a[p][p].re) / (2.0 * r);
                let t   = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    1.0 / (tau - (1.0 + tau * tau).sqrt())
                };
                let c  = 1.0 / (1.0 + t * t).sqrt();
                let s  = t * c;
                let ep  = C64::polar(1.0,  phi);
                let enp = C64::polar(1.0, -phi);

                // Right-multiply: A ← A · G  (update columns p, q)
                let mut cp = [C64::ZERO; NMAT];
                let mut cq = [C64::ZERO; NMAT];
                for i in 0..NMAT {
                    cp[i] = a[i][p].scale(c) + (enp * a[i][q]).scale(-s);
                    cq[i] = (ep * a[i][p]).scale(s) + a[i][q].scale(c);
                }
                for i in 0..NMAT { a[i][p] = cp[i]; a[i][q] = cq[i]; }

                // Left-multiply: A ← G† · A  (update rows p, q)
                let rp = a[p];
                let rq = a[q];
                for j in 0..NMAT {
                    a[p][j] = rp[j].scale(c)      + (ep  * rq[j]).scale(-s);
                    a[q][j] = (enp * rp[j]).scale(s) + rq[j].scale(c);
                }

                // Suppress floating-point drift
                a[p][p].im = 0.0;
                a[q][q].im = 0.0;
                let pq = a[p][q];
                a[q][p]   = pq.conj();
            }
        }
    }

    let mut eigs = [0.0_f64; NMAT];
    for i in 0..NMAT { eigs[i] = a[i][i].re; }
    eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}

// ══════════════════════════════════════════════════════════════════════════════
// SMRK observables  (spec §3)
// ══════════════════════════════════════════════════════════════════════════════

struct SmrkResult {
    eigs:    [f64; NMAT],
    delta1:  f64,
    r_mean:  f64,
    r_count: usize,
    h1_pass: bool,
    h2_pass: bool,
}

/// Map SHA-256(input)[0..2] → (κ, λ, φ_offset), build matrix, diagonalise.
///
/// κ ∈ [0.05, 0.20]  — coupling strength
/// λ ∈ [0.00, 0.20]  — diagonal scale
/// φ ∈ [0.00, 0.20]  — golden-gauge offset
fn compute_smrk(input: &[u8]) -> SmrkResult {
    let h = sha256_raw(input);
    let kappa      = 0.05 + (h[0] as f64 / 255.0) * 0.15;
    let lam        = (h[1] as f64 / 255.0) * 0.20;
    let phi_offset = (h[2] as f64 / 255.0) * 0.20;

    let mat  = build_smrk_matrix(kappa, lam, phi_offset);
    let eigs = jacobi_eigen(&mat);

    // H2: primary gap Δ1 = E₁ − E₀
    let delta1 = eigs[1] - eigs[0];

    // H1: bulk spacing ratios  rᵢ = min(sᵢ₋₁,sᵢ) / max(sᵢ₋₁,sᵢ)
    let mut r_sum   = 0.0_f64;
    let mut r_count = 0usize;
    for i in BULK_I0..=BULK_I1 {
        if i + 1 >= NMAT { break; }
        let s_prev = eigs[i]     - eigs[i - 1];
        let s_curr = eigs[i + 1] - eigs[i];
        if s_prev > 0.0 && s_curr > 0.0 {
            r_sum   += s_prev.min(s_curr) / s_prev.max(s_curr);
            r_count += 1;
        }
    }
    let r_mean = if r_count > 0 { r_sum / r_count as f64 } else { 0.0 };

    let h1 = (r_mean - R_GUE).abs() <= H1_BAND;
    let h2 = delta1 >= H2_MIN_GAP;

    SmrkResult { eigs, delta1, r_mean, r_count, h1_pass: h1, h2_pass: h2 }
}

// ══════════════════════════════════════════════════════════════════════════════
// Crypto helpers
// ══════════════════════════════════════════════════════════════════════════════

fn sha256_raw(data: &[u8]) -> [u8; 32] {
    let mut h = Sha256::new();
    h.update(data);
    h.finalize().into()
}

fn sha256_hex(data: &[u8]) -> String { hex::encode(sha256_raw(data)) }

/// anchor_commitment = SHA256("QFC|v4|TWISTOR|" ‖ run_id ‖ commit_hash_hex ‖ plateau_flag)
/// spec §D.1.3
fn compute_anchor_commitment(run_id: &str, commit_hash_hex: &str, plateau_flag: bool) -> [u8; 32] {
    let mut data = b"QFC|v4|TWISTOR|".to_vec();
    data.extend_from_slice(run_id.as_bytes());
    data.extend_from_slice(commit_hash_hex.as_bytes());
    data.extend_from_slice(if plateau_flag { b"true" } else { b"false" });
    sha256_raw(&data)
}

// ══════════════════════════════════════════════════════════════════════════════
// Canonical JSON  (stable key order, 16-decimal scientific notation)
// ══════════════════════════════════════════════════════════════════════════════

fn fmt_eigs(eigs: &[f64], k: usize) -> String {
    let parts: Vec<String> = eigs.iter().take(k).map(|e| format!("{:.16e}", e)).collect();
    format!("[{}]", parts.join(","))
}

/// Math payload — NO metadata → SHA-256 stable across contract upgrades.
/// Stored as math_sha256.  probe_version: smrk-guegap-near-v1
fn math_payload(input: &[u8], input_sha256_hex: &str) -> (Vec<u8>, [u8; 32]) {
    let res       = compute_smrk(input);
    let pass_str  = if res.h1_pass && res.h2_pass { "pass" } else { "fail" };
    let eigs_json = fmt_eigs(&res.eigs, 5);

    // Lexicographically ordered keys, 16-decimal scientific notation (spec §A.1.1)
    let json = format!(
        "{{\"H1_H2_proxy\":\"{}\",\
\"N\":{},\
\"bulk_r\":{{\"count\":{},\"r_mean\":{:.16e}}},\
\"eigs_low\":{},\
\"gap\":{{\"delta1\":{:.16e}}},\
\"input_sha256_hex\":\"{}\",\
\"probe_version\":\"smrk-guegap-near-v1\"}}",
        pass_str, NMAT,
        res.r_count, res.r_mean,
        eigs_json,
        res.delta1,
        input_sha256_hex,
    );
    let bytes = json.into_bytes();
    let hash  = sha256_raw(&bytes);
    (bytes, hash)
}

/// Audit output — adds contract build metadata.
/// SHA-256 changes when contract is redeployed (code_hash changes).
/// Stored as output / output_sha256.
fn confirmation_audit_json(
    input: &[u8],
    input_sha256_hex: &str,
    code_hash_hex: &str,
    block_height: u64,
) -> (Vec<u8>, [u8; 32]) {
    let (_, math_hash) = math_payload(input, input_sha256_hex);
    let math_sha_hex   = hex::encode(math_hash);

    let res       = compute_smrk(input);
    let pass_str  = if res.h1_pass && res.h2_pass { "pass" } else { "fail" };
    let eigs_json = fmt_eigs(&res.eigs, 5);

    let json = format!(
        "{{\"H1_H2_proxy\":\"{}\",\
\"N\":{},\
\"bulk_r\":{{\"count\":{},\"r_mean\":{:.16e}}},\
\"eigs_low\":{},\
\"gap\":{{\"delta1\":{:.16e}}},\
\"meta\":{{\
\"block_height\":{},\
\"code_hash_hex\":\"{}\",\
\"contract_version\":\"0.1.0\",\
\"input_sha256_hex\":\"{}\",\
\"math_sha256_hex\":\"{}\"\
}},\
\"probe_version\":\"smrk-guegap-near-v1\"}}",
        pass_str, NMAT,
        res.r_count, res.r_mean,
        eigs_json,
        res.delta1,
        block_height,
        code_hash_hex,
        input_sha256_hex,
        math_sha_hex,
    );
    let bytes = json.into_bytes();
    let hash  = sha256_raw(&bytes);
    (bytes, hash)
}

// ══════════════════════════════════════════════════════════════════════════════
// Data types  (Borsh for state, Serde for JSON view methods)
// ══════════════════════════════════════════════════════════════════════════════

#[derive(BorshSerialize, BorshDeserialize, Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(crate = "near_sdk::serde")]
pub enum JobStatus { Queued, Running, Done, Failed }

#[derive(BorshSerialize, BorshDeserialize, Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(crate = "near_sdk::serde")]
pub enum AnchorStatus { Pending, Recorded, Verified, Failed }

#[derive(BorshSerialize, BorshDeserialize, Serialize, Deserialize, Clone, Debug)]
#[serde(crate = "near_sdk::serde")]
pub struct JobRecord {
    pub run_id:           String,
    pub created_at_ns:   u64,
    pub submitted_by:    AccountId,
    pub status:          JobStatus,

    // Input
    pub input:           Vec<u8>,
    pub input_sha256_hex: String,

    // Screening result (from ICP tier — optional, for cross-chain audit)
    pub icp_run_id:      Option<String>,
    pub icp_math_sha256_hex: Option<String>,

    // NEAR confirmation output
    pub output:          Option<Vec<u8>>,
    pub output_sha256_hex: Option<String>,
    pub math_sha256_hex: Option<String>,
    pub commit_hash_hex: Option<String>,  // = output_sha256_hex, used in anchor formula
    pub error:           Option<String>,

    // Kaspa anchoring
    pub anchor_commitment_hex: Option<String>,  // hex payload for Kaspa tx
    pub kaspa_txid:           Option<String>,
    pub anchor_status:        Option<AnchorStatus>,
    pub anchor_note:          Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct SubmitJobResp {
    pub run_id:          String,
    pub input_sha256_hex: String,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct ConfirmResp {
    pub ok:             bool,
    pub run_id:         String,
    pub h1_h2_proxy:    String,
    pub math_sha256_hex: Option<String>,
    pub message:        String,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct VerifyResp {
    pub ok:      bool,
    pub matches: bool,
    pub message: String,
    pub stored_hex:   Option<String>,
    pub computed_hex: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct AnchorResp {
    pub ok:                   bool,
    pub anchor_commitment_hex: Option<String>,
    pub message:              String,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct AnchorInfoResp {
    pub run_id:                String,
    pub anchor_commitment_hex: Option<String>,
    pub kaspa_txid:            Option<String>,
    pub anchor_status:         Option<AnchorStatus>,
    pub anchor_note:           Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct ContractMeta {
    pub contract_version: String,
    pub probe_version:    String,
    pub n_matrix:         u32,
    pub code_hash_hex:    String,
    pub block_height:     u64,
    pub block_timestamp_ns: u64,
}

// ══════════════════════════════════════════════════════════════════════════════
// Contract state
// ══════════════════════════════════════════════════════════════════════════════

#[near_bindgen]
#[derive(BorshSerialize, BorshDeserialize, PanicOnDefault)]
pub struct SmrkGuegapContract {
    owner:      AccountId,
    jobs:       UnorderedMap<String, JobRecord>,
    job_count:  u64,
}

#[near_bindgen]
impl SmrkGuegapContract {

    // ── Constructor ──────────────────────────────────────────────────────────

    #[init]
    pub fn new(owner: AccountId) -> Self {
        Self {
            owner,
            jobs:      UnorderedMap::new(b"j"),
            job_count: 0,
        }
    }

    // ── Job submission ───────────────────────────────────────────────────────

    /// Submit raw input bytes for NEAR confirmation screening.
    /// run_id = SHA-256(input) as hex — same derivation as ICP tier for
    /// cross-chain correlation.
    /// Optional icp_run_id / icp_math_sha256_hex for cross-chain audit trail.
    pub fn submit_job(
        &mut self,
        input: Vec<u8>,
        icp_run_id: Option<String>,
        icp_math_sha256_hex: Option<String>,
    ) -> SubmitJobResp {
        let input_hash     = sha256_raw(&input);
        let input_sha256_hex = hex::encode(input_hash);
        let run_id         = input_sha256_hex.clone(); // same as ICP

        if self.jobs.contains_key(&run_id) {
            // Idempotent: return existing run_id without re-inserting
            return SubmitJobResp { run_id, input_sha256_hex };
        }

        let job = JobRecord {
            run_id:              run_id.clone(),
            created_at_ns:       env::block_timestamp(),
            submitted_by:        env::predecessor_account_id(),
            status:              JobStatus::Queued,
            input,
            input_sha256_hex:    input_sha256_hex.clone(),
            icp_run_id,
            icp_math_sha256_hex,
            output:              None,
            output_sha256_hex:   None,
            math_sha256_hex:     None,
            commit_hash_hex:     None,
            error:               None,
            anchor_commitment_hex: None,
            kaspa_txid:          None,
            anchor_status:       None,
            anchor_note:         None,
        };

        self.jobs.insert(run_id.clone(), job);
        self.job_count += 1;

        SubmitJobResp { run_id, input_sha256_hex }
    }

    // ── NEAR confirmation  ────────────────────────────────────────────────────
    //
    // Runs the full N=32 SMRK-GUEGAP computation on-chain.
    // Gas budget: ~30 TGas for Jacobi + overhead → well within 300 TGas limit.

    pub fn run_confirmation(&mut self, run_id: String) -> ConfirmResp {
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None => return ConfirmResp {
                ok: false, run_id, h1_h2_proxy: "fail".into(),
                math_sha256_hex: None,
                message: "Job not found. Call submit_job first.".into(),
            },
        };

        if job.status == JobStatus::Done {
            return ConfirmResp {
                ok: true, run_id, h1_h2_proxy: "pass".into(),
                math_sha256_hex: job.math_sha256_hex,
                message: "Already confirmed.".into(),
            };
        }
        if job.status == JobStatus::Failed {
            return ConfirmResp {
                ok: false, run_id, h1_h2_proxy: "fail".into(),
                math_sha256_hex: None,
                message: format!("Job previously failed: {}", job.error.as_deref().unwrap_or("?")),
            };
        }

        job.status = JobStatus::Running;

        // ── Compute ──────────────────────────────────────────────────────────
        let code_hash_hex  = hex::encode(env::code_hash());
        let block_height   = env::block_height();

        let (math_bytes, math_hash) = math_payload(&job.input, &job.input_sha256_hex);
        let (audit_bytes, audit_hash) = confirmation_audit_json(
            &job.input,
            &job.input_sha256_hex,
            &code_hash_hex,
            block_height,
        );

        // Decode pass/fail from math_bytes for response field
        let math_str     = String::from_utf8_lossy(&math_bytes);
        let h1h2         = if math_str.contains("\"pass\"") { "pass" } else { "fail" };
        let math_sha_hex = hex::encode(math_hash);
        let out_sha_hex  = hex::encode(audit_hash);

        job.status             = JobStatus::Done;
        job.output             = Some(audit_bytes);
        job.output_sha256_hex  = Some(out_sha_hex.clone());
        job.math_sha256_hex    = Some(math_sha_hex.clone());
        job.commit_hash_hex    = Some(out_sha_hex);
        job.error              = None;

        self.jobs.insert(run_id.clone(), job);

        ConfirmResp {
            ok: true,
            run_id,
            h1_h2_proxy: h1h2.to_string(),
            math_sha256_hex: Some(math_sha_hex),
            message: "N=32 confirmation complete.".into(),
        }
    }

    // ── Verification ──────────────────────────────────────────────────────────

    /// Recomputes math_payload and checks against stored math_sha256_hex.
    /// Always matches — stable across contract redeployments.
    pub fn verify_math(&self, run_id: String) -> VerifyResp {
        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None => return VerifyResp { ok: false, matches: false,
                message: "Job not found.".into(), stored_hex: None, computed_hex: None },
        };
        let stored = match &job.math_sha256_hex {
            Some(h) => h.clone(),
            None => return VerifyResp { ok: false, matches: false,
                message: "No math_sha256_hex yet. Run run_confirmation first.".into(),
                stored_hex: None, computed_hex: None },
        };
        let (_, hash) = math_payload(&job.input, &job.input_sha256_hex);
        let computed  = hex::encode(hash);
        let matches   = stored == computed;
        VerifyResp {
            ok: true, matches,
            message: if matches { "math_sha256 verified.".into() }
                     else { "math_sha256 MISMATCH — unexpected.".into() },
            stored_hex: Some(stored), computed_hex: Some(computed),
        }
    }

    /// Recomputes audit JSON and checks against stored output_sha256_hex.
    /// Will NOT match after contract redeployment (code_hash changes).
    pub fn verify_output(&self, run_id: String) -> VerifyResp {
        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None => return VerifyResp { ok: false, matches: false,
                message: "Job not found.".into(), stored_hex: None, computed_hex: None },
        };
        let stored = match &job.output_sha256_hex {
            Some(h) => h.clone(),
            None => return VerifyResp { ok: false, matches: false,
                message: "No output_sha256_hex yet.".into(),
                stored_hex: None, computed_hex: None },
        };
        let code_hash_hex = hex::encode(env::code_hash());
        let block_height  = env::block_height();
        let (bytes, _) = confirmation_audit_json(
            &job.input, &job.input_sha256_hex, &code_hash_hex, block_height
        );
        let computed = sha256_hex(&bytes);
        let matches  = stored == computed;
        VerifyResp {
            ok: true, matches,
            message: if matches { "output_sha256 verified.".into() }
                     else { "Mismatch — expected after redeployment (code_hash/block_height change).".into() },
            stored_hex: Some(stored), computed_hex: Some(computed),
        }
    }

    // ── Kaspa anchoring ───────────────────────────────────────────────────────

    /// Step 1 — compute anchor_commitment and store it.
    /// Requires status == Done.
    /// Returns anchor_commitment_hex to embed as Kaspa tx payload.
    pub fn anchor_job(&mut self, run_id: String, plateau_flag: bool) -> AnchorResp {
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None => return AnchorResp { ok: false, anchor_commitment_hex: None,
                message: "Job not found.".into() },
        };
        if job.status != JobStatus::Done {
            return AnchorResp { ok: false, anchor_commitment_hex: None,
                message: "Job not confirmed yet. Run run_confirmation first.".into() };
        }
        let commit_hex = job.commit_hash_hex.clone().unwrap_or_default();
        let raw  = compute_anchor_commitment(&run_id, &commit_hex, plateau_flag);
        let hex  = hex::encode(raw);

        job.anchor_commitment_hex = Some(hex.clone());
        job.anchor_status         = Some(AnchorStatus::Pending);
        job.anchor_note           = Some("Commitment computed; awaiting Kaspa broadcast.".into());
        self.jobs.insert(run_id, job);

        AnchorResp {
            ok: true,
            anchor_commitment_hex: Some(hex.clone()),
            message: format!(
                "Anchor commitment stored. Embed this hex as Kaspa tx payload, \
then call record_kaspa_txid(run_id, txid). Payload: {}", hex
            ),
        }
    }

    /// Step 2 — record Kaspa txid after external broadcast.
    /// Advances anchor_status: Pending → Recorded.
    pub fn record_kaspa_txid(&mut self, run_id: String, kaspa_txid: String) -> bool {
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None => return false,
        };
        if job.anchor_commitment_hex.is_none() { return false; }
        let txid = kaspa_txid.clone();
        job.kaspa_txid    = Some(kaspa_txid);
        job.anchor_status = Some(AnchorStatus::Recorded);
        job.anchor_note   = Some(format!(
            "Kaspa txid recorded: {}. Call verify_kaspa_anchor_offchain or an oracle.", txid
        ));
        self.jobs.insert(run_id, job);
        true
    }

    /// Step 3 — record result of off-chain (or oracle) Kaspa verification.
    /// NEAR contracts cannot make HTTP calls; verification is done externally
    /// and the result is submitted here by a trusted caller (owner or oracle).
    /// For trustless verification, integrate an oracle contract.
    pub fn submit_kaspa_verification(
        &mut self,
        run_id:       String,
        verified:     bool,
        note:         String,
    ) -> bool {
        // Only contract owner or the submitting account may record verification
        let caller = env::predecessor_account_id();
        assert!(
            caller == self.owner,
            "Only the contract owner can submit Kaspa verification results."
        );

        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None => return false,
        };
        job.anchor_status = Some(if verified { AnchorStatus::Verified } else { AnchorStatus::Failed });
        job.anchor_note   = Some(note);
        self.jobs.insert(run_id, job);
        true
    }

    // ── Cross-chain correlation helper ────────────────────────────────────────

    /// Checks whether this job's math_sha256_hex matches the ICP tier's
    /// math_sha256 that was provided at submit_job time.
    /// Both use the same deterministic formula → should always match
    /// for the same input.
    pub fn verify_icp_cross_chain(&self, run_id: String) -> VerifyResp {
        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None => return VerifyResp { ok: false, matches: false,
                message: "Job not found.".into(), stored_hex: None, computed_hex: None },
        };
        let near_math = match &job.math_sha256_hex {
            Some(h) => h.clone(),
            None => return VerifyResp { ok: false, matches: false,
                message: "No NEAR math_sha256_hex yet.".into(), stored_hex: None, computed_hex: None },
        };
        let icp_math = match &job.icp_math_sha256_hex {
            Some(h) => h.clone(),
            None => return VerifyResp { ok: false, matches: false,
                message: "No ICP math_sha256_hex provided at submit time.".into(),
                stored_hex: Some(near_math), computed_hex: None },
        };
        // NOTE: ICP uses N=16, NEAR uses N=32 — different N means different
        // probe_version and different eigenvalues. math_sha256 will NOT match
        // between tiers unless input is identical AND probe_version is same.
        // This call is informational — it shows the two hashes side-by-side
        // for manual audit. A match would indicate same N (unexpected here).
        let matches = near_math == icp_math;
        VerifyResp {
            ok: true,
            matches,
            message: if matches {
                "Hashes match (same N or same stub formula — investigate).".into()
            } else {
                "Hashes differ as expected: ICP N=16 vs NEAR N=32 produce distinct math_sha256. Manual audit: compare eigenvalue structure.".into()
            },
            stored_hex: Some(near_math),    // NEAR math_sha256
            computed_hex: Some(icp_math),   // ICP math_sha256 (reference)
        }
    }

    // ── Read methods (view) ───────────────────────────────────────────────────

    pub fn get_job(&self, run_id: String) -> Option<JobRecord> {
        self.jobs.get(&run_id).cloned()
    }

    pub fn get_anchor_info(&self, run_id: String) -> Option<AnchorInfoResp> {
        self.jobs.get(&run_id).map(|j| AnchorInfoResp {
            run_id:                j.run_id.clone(),
            anchor_commitment_hex: j.anchor_commitment_hex.clone(),
            kaspa_txid:            j.kaspa_txid.clone(),
            anchor_status:         j.anchor_status.clone(),
            anchor_note:           j.anchor_note.clone(),
        })
    }

    pub fn get_contract_meta(&self) -> ContractMeta {
        ContractMeta {
            contract_version:   "0.1.0".into(),
            probe_version:      "smrk-guegap-near-v1".into(),
            n_matrix:           NMAT as u32,
            code_hash_hex:      hex::encode(env::code_hash()),
            block_height:       env::block_height(),
            block_timestamp_ns: env::block_timestamp(),
        }
    }

    pub fn get_job_count(&self) -> u64 { self.job_count }

    pub fn get_owner(&self) -> AccountId { self.owner.clone() }
}
