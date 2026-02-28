// ══════════════════════════════════════════════════════════════════════════════
// SMRK-GUEGAP Probe v4 — NEAR Confirmation Tier  (v0.3.0)
//
// probe_version : smrk-guegap-near-v3
// N             : 32  (confirmation; spec N=512 → scaled for gas budget)
// Role          : Second-layer confirmation of ICP screening results
//
// ── Changelog from v0.2.0 (smrk-guegap-near-v2) ────────────────────────────
//
//   FIX-4  "Already confirmed" branch returned hardcoded pass/true
//          JobRecord now stores h1_pass, h2_pass, r_mean, delta1, h1_h2_proxy
//          at confirmation time. "Already confirmed" branch returns stored
//          spectral values — never hardcoded constants.
//
//   FIX-5  STORAGE_PRICE_PER_BYTE was hardcoded (cross-network fragility)
//          Now uses env::storage_byte_cost() at runtime — always correct
//          for the network the contract is deployed on.
//
//   FIX-6  submit_job() accepted arbitrary input_sha256_hex without validation
//          Now validates: length == 64, chars ∈ [0-9a-f], lowercase canonical.
//          Also validates: input_len > 0 and input_len <= MAX_INPUT_LEN.
//
//   FIX-7  run_confirmation() accepted unbounded input Vec<u8>
//          Now asserts input.len() <= MAX_INPUT_BYTES (4096).
//          Prevents gas grief / compute grief from oversized inputs.
//
//   FIX-8  NaN in jacobi_eigen() caused panic → revert, no audit trail
//          jacobi_eigen() now returns Result<[f64;N], String>.
//          compute_smrk() returns Result<SmrkResult, String>.
//          run_confirmation() catches Err: sets JobStatus::Failed + error
//          field, saves to state, returns ConfirmResp { ok: false }.
//          Audit trail always present, even for degenerate inputs.
//
//   FIX-9  record_kaspa_txid() had no state-machine guard
//          Now asserts anchor_status == Some(Pending) before writing txid.
//          Enforces: anchor_job() must precede record_kaspa_txid().
//
// ── Carried from v0.2.0 ───────────────────────────────────────────────────────
//
//   FIX-1  verify_output() used live env::block_height()/code_hash() → mismatch
//   FIX-2  Storage griefing: submit_job #[payable], input NOT stored on-chain
//   FIX-3  record_kaspa_txid open write: one-time + caller ACL
//   PERF-1 Single Jacobi per confirmation (SmrkResult passed to both builders)
//   PERF-2 h1_h2_proxy from bool, not string parsing
//   SAFE-1 NaN guard (now returns Result instead of panic — see FIX-8)
//   API-1  compare_icp_reference() (renamed from verify_icp_cross_chain)
//
// ── NEAR gas budget ───────────────────────────────────────────────────────────
//   300 TGas / tx. Jacobi N=32, 300 sweeps ≈ 10M f64 ops ≈ ~30 TGas → safe.
//   Single Jacobi per confirmation (PERF-1).
//
// Spec refs: §2 (operator), §3 (observables), §4 (acceptance), §D.1.3 (anchor)
// ══════════════════════════════════════════════════════════════════════════════

use near_sdk::borsh::{BorshDeserialize, BorshSerialize};
use near_sdk::serde::{Deserialize, Serialize};
use near_sdk::store::UnorderedMap;
use near_sdk::{env, near_bindgen, AccountId, NearToken, PanicOnDefault, Promise};
use sha2::{Digest, Sha256};
use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════════
// SMRK operator constants
// ══════════════════════════════════════════════════════════════════════════════

const NMAT: usize = 32;
const M_SHIFTS: usize = 8;
const PRIMES: [usize; 8] = [2, 3, 5, 7, 11, 13, 17, 19];
const GOLDEN: f64   = 0.6180339887498949; // α = (√5−1)/2  golden gauge
const S_EXP: f64    = 1.0;                // prime-decay exponent s
const ETA: f64      = 0.5;                // Liouville modulation strength η

const H1_BAND: f64    = 0.10;             // |⟨r⟩ − r_GUE| ≤ 0.10
const R_GUE: f64      = 0.5996;           // GUE bulk ratio reference
const H2_MIN_GAP: f64 = 1e-4;            // Δ1 ≥ 1×10⁻⁴

const BULK_I0: usize  = 6;               // ⌊0.2·32⌋
const BULK_I1: usize  = 23;              // ⌊0.8·32⌋ − 2

// FIX-7: input size guard — prevents gas/compute grief
const MAX_INPUT_BYTES: usize = 4_096;

// FIX-6: input_len sanity bounds
const MAX_INPUT_LEN: u32 = 4_096;

// Probe + contract version strings — bump on every release
const PROBE_VERSION:    &str = "smrk-guegap-near-v3";
const CONTRACT_VERSION: &str = "0.3.0";

// ══════════════════════════════════════════════════════════════════════════════
// Complex arithmetic  (stack-allocated, no alloc)
// ══════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy)]
struct C64 { re: f64, im: f64 }

impl C64 {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };
    fn conj(self)                -> Self { Self { re: self.re, im: -self.im } }
    fn norm(self)                -> f64  { (self.re * self.re + self.im * self.im).sqrt() }
    fn scale(self, s: f64)       -> Self { Self { re: self.re * s, im: self.im * s } }
    fn arg(self)                 -> f64  { self.im.atan2(self.re) }
    fn polar(r: f64, theta: f64) -> Self { Self { re: r * theta.cos(), im: r * theta.sin() } }
    fn is_nan(self)              -> bool { self.re.is_nan() || self.im.is_nan() }
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
    // n always ≥ 1 at call sites (n+1 passed). Guard kept for safety.
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

    // Diagonal: D_nn = λ · ln(n+1),  0-indexed → ln(n+2)
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
            let amp   = kappa / (1.0 + pm as f64).ln().powf(S_EXP) * (1.0 + ETA * lv);
            let phi   = 2.0 * PI * GOLDEN * n1 + phi_offset;
            let w     = C64::polar(amp, phi);
            h[n][sigma] = h[n][sigma] + w;
            h[sigma][n] = h[sigma][n] + w.conj();
        }
    }
    h
}

// ══════════════════════════════════════════════════════════════════════════════
// Cyclic complex Jacobi eigenvalue solver  (§6.12)
//
// FIX-8: Returns Result<[f64; NMAT], String> instead of panicking on NaN.
//        Caller (compute_smrk) propagates Err → run_confirmation saves Failed.
// ══════════════════════════════════════════════════════════════════════════════

fn jacobi_eigen(h_in: &Matrix) -> Result<[f64; NMAT], String> {
    let mut a = *h_in;
    let tol   = 1e-11_f64;

    'outer: for _ in 0..300 {
        let mut max_off = 0.0_f64;
        for i in 0..NMAT {
            for j in 0..NMAT {
                if i != j {
                    let v = a[i][j].norm();
                    if v > max_off { max_off = v; }
                }
            }
        }
        if max_off < tol { break 'outer; }

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
                let c   = 1.0 / (1.0 + t * t).sqrt();
                let s   = t * c;
                let ep  = C64::polar(1.0,  phi);
                let enp = C64::polar(1.0, -phi);

                // Right-multiply: A ← A · G
                let mut cp = [C64::ZERO; NMAT];
                let mut cq = [C64::ZERO; NMAT];
                for i in 0..NMAT {
                    cp[i] = a[i][p].scale(c) + (enp * a[i][q]).scale(-s);
                    cq[i] = (ep * a[i][p]).scale(s) + a[i][q].scale(c);
                }
                for i in 0..NMAT { a[i][p] = cp[i]; a[i][q] = cq[i]; }

                // Left-multiply: A ← G† · A
                let rp = a[p];
                let rq = a[q];
                for j in 0..NMAT {
                    a[p][j] = rp[j].scale(c)        + (ep  * rq[j]).scale(-s);
                    a[q][j] = (enp * rp[j]).scale(s) + rq[j].scale(c);
                }

                // Suppress floating-point drift
                a[p][p].im = 0.0;
                a[q][q].im = 0.0;
                a[q][p]    = a[p][q].conj();
            }
        }
    }

    // FIX-8: NaN → Err instead of panic. Saves audit trail.
    for i in 0..NMAT {
        if a[i][i].is_nan() {
            return Err(format!(
                "jacobi_eigen: NaN on diagonal at index {}. \
                 Input produced a degenerate operator. Saved as Failed.",
                i
            ));
        }
    }

    let mut eigs = [0.0_f64; NMAT];
    for i in 0..NMAT { eigs[i] = a[i][i].re; }
    eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    Ok(eigs)
}

// ══════════════════════════════════════════════════════════════════════════════
// SMRK observables  (spec §3)
//
// FIX-8: Returns Result<SmrkResult, String> — Err propagated to caller.
// PERF-1: Called ONCE per confirmation. SmrkResult shared by both JSON builders.
// ══════════════════════════════════════════════════════════════════════════════

struct SmrkResult {
    eigs:    [f64; NMAT],
    delta1:  f64,
    r_mean:  f64,
    r_count: usize,
    h1_pass: bool,
    h2_pass: bool,
}

fn compute_smrk(input: &[u8]) -> Result<SmrkResult, String> {
    let h          = sha256_raw(input);
    let kappa      = 0.05 + (h[0] as f64 / 255.0) * 0.15;
    let lam        = (h[1] as f64 / 255.0) * 0.20;
    let phi_offset = (h[2] as f64 / 255.0) * 0.20;

    let mat  = build_smrk_matrix(kappa, lam, phi_offset);
    let eigs = jacobi_eigen(&mat)?;   // FIX-8: propagate Err

    let delta1 = eigs[1] - eigs[0];

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

    Ok(SmrkResult {
        eigs,
        delta1,
        r_mean,
        r_count,
        h1_pass: (r_mean - R_GUE).abs() <= H1_BAND,
        h2_pass: delta1 >= H2_MIN_GAP,
    })
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
/// spec §D.1.3  — identical formula to ICP tier
fn compute_anchor_commitment(run_id: &str, commit_hash_hex: &str, plateau_flag: bool) -> [u8; 32] {
    let mut data = b"QFC|v4|TWISTOR|".to_vec();
    data.extend_from_slice(run_id.as_bytes());
    data.extend_from_slice(commit_hash_hex.as_bytes());
    data.extend_from_slice(if plateau_flag { b"true" } else { b"false" });
    sha256_raw(&data)
}

// FIX-6: validate that a caller-supplied SHA-256 hex string is well-formed.
// Returns Err(message) so the caller can assert! with a meaningful message.
fn validate_sha256_hex(s: &str) -> Result<(), String> {
    if s.len() != 64 {
        return Err(format!(
            "input_sha256_hex must be exactly 64 hex characters, got {}.", s.len()
        ));
    }
    for c in s.chars() {
        if !matches!(c, '0'..='9' | 'a'..='f') {
            return Err(format!(
                "input_sha256_hex must contain only lowercase hex [0-9a-f], got '{}'.", c
            ));
        }
    }
    Ok(())
}

// ══════════════════════════════════════════════════════════════════════════════
// Canonical JSON builders  (PERF-1: receive &SmrkResult, no second Jacobi)
// Key ordering: lexicographic. Float format: 16-decimal scientific (spec §A.1.1)
// ══════════════════════════════════════════════════════════════════════════════

fn fmt_eigs(eigs: &[f64], k: usize) -> String {
    let parts: Vec<String> = eigs.iter().take(k).map(|e| format!("{:.16e}", e)).collect();
    format!("[{}]", parts.join(","))
}

/// Math payload — NO deployment metadata → SHA-256 stable across upgrades.
fn math_payload_from_result(res: &SmrkResult, input_sha256_hex: &str) -> (Vec<u8>, [u8; 32]) {
    let pass_str  = if res.h1_pass && res.h2_pass { "pass" } else { "fail" };
    let eigs_json = fmt_eigs(&res.eigs, 5);
    let json = format!(
        "{{\"H1_H2_proxy\":\"{}\",\
\"N\":{},\
\"bulk_r\":{{\"count\":{},\"r_mean\":{:.16e}}},\
\"eigs_low\":{},\
\"gap\":{{\"delta1\":{:.16e}}},\
\"input_sha256_hex\":\"{}\",\
\"probe_version\":\"{}\"}}",
        pass_str, NMAT,
        res.r_count, res.r_mean,
        eigs_json,
        res.delta1,
        input_sha256_hex,
        PROBE_VERSION,
    );
    let bytes = json.into_bytes();
    let hash  = sha256_raw(&bytes);
    (bytes, hash)
}

/// Audit output JSON — adds deployment metadata.
///
/// FIX-1: code_hash_hex + block_height are parameters captured at confirmation
///        time and stored in JobRecord. Never read from live env::* here.
fn audit_json_from_result(
    res:              &SmrkResult,
    input_sha256_hex: &str,
    math_sha256_hex:  &str,
    code_hash_hex:    &str,
    block_height:     u64,
) -> (Vec<u8>, [u8; 32]) {
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
\"contract_version\":\"{}\",\
\"input_sha256_hex\":\"{}\",\
\"math_sha256_hex\":\"{}\"\
}},\
\"probe_version\":\"{}\"}}",
        pass_str, NMAT,
        res.r_count, res.r_mean,
        eigs_json,
        res.delta1,
        block_height,
        code_hash_hex,
        CONTRACT_VERSION,
        input_sha256_hex,
        math_sha256_hex,
        PROBE_VERSION,
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

/// State stored per job.
///
/// FIX-4: h1_pass, h2_pass, r_mean, delta1, h1_h2_proxy stored at confirmation.
/// FIX-1: confirmed_block_height, confirmed_code_hash_hex stored at confirmation.
/// FIX-2: raw input bytes NOT stored. Only input_sha256_hex + input_len.
#[derive(BorshSerialize, BorshDeserialize, Serialize, Deserialize, Clone, Debug)]
#[serde(crate = "near_sdk::serde")]
pub struct JobRecord {
    pub run_id:              String,
    pub created_at_ns:       u64,
    pub submitted_by:        AccountId,
    pub status:              JobStatus,

    // Input commitment — raw bytes never stored (FIX-2)
    pub input_sha256_hex:    String,
    pub input_len:           u32,

    // ICP cross-chain reference (optional)
    pub icp_run_id:          Option<String>,
    pub icp_math_sha256_hex: Option<String>,

    // NEAR confirmation output — hashes
    pub output_sha256_hex:   Option<String>,
    pub math_sha256_hex:     Option<String>,
    pub commit_hash_hex:     Option<String>,  // = output_sha256_hex, used in anchor formula
    pub error:               Option<String>,

    // FIX-4: cached spectral result — returned by "Already confirmed" branch
    pub h1_h2_proxy:         Option<String>,  // "pass" | "fail"
    pub h1_pass:             Option<bool>,
    pub h2_pass:             Option<bool>,
    pub r_mean:              Option<f64>,
    pub delta1:              Option<f64>,

    // FIX-1: env snapshot at confirmation time — used by verify_output()
    pub confirmed_block_height:   Option<u64>,
    pub confirmed_code_hash_hex:  Option<String>,

    // Kaspa anchoring
    pub anchor_commitment_hex: Option<String>,
    pub kaspa_txid:            Option<String>,
    pub anchor_status:         Option<AnchorStatus>,
    pub anchor_note:           Option<String>,
}

// ── Response types ────────────────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct SubmitJobResp {
    pub run_id:           String,
    pub input_sha256_hex: String,
    pub input_len:        u32,
    pub storage_used:     u64,
    pub deposit_returned: bool,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct ConfirmResp {
    pub ok:                      bool,
    pub run_id:                  String,
    pub h1_h2_proxy:             String,
    pub h1_pass:                 bool,
    pub h2_pass:                 bool,
    pub r_mean:                  f64,
    pub delta1:                  f64,
    pub math_sha256_hex:         Option<String>,
    pub confirmed_block_height:  u64,
    pub confirmed_code_hash_hex: String,
    pub message:                 String,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct VerifyResp {
    pub ok:           bool,
    pub matches:      bool,
    pub message:      String,
    pub stored_hex:   Option<String>,
    pub computed_hex: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct IcpCompareResp {
    pub ok:                   bool,
    pub near_math_sha256_hex: Option<String>,
    pub icp_math_sha256_hex:  Option<String>,
    pub near_probe_version:   String,
    pub icp_probe_version:    Option<String>,
    pub hashes_match:         bool,
    pub note:                 String,
}

#[derive(Serialize, Deserialize)]
#[serde(crate = "near_sdk::serde")]
pub struct AnchorResp {
    pub ok:                    bool,
    pub anchor_commitment_hex: Option<String>,
    pub message:               String,
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
    pub contract_version:   String,
    pub probe_version:      String,
    pub n_matrix:           u32,
    pub code_hash_hex:      String,
    pub block_height:       u64,
    pub block_timestamp_ns: u64,
}

// ══════════════════════════════════════════════════════════════════════════════
// Contract state
// ══════════════════════════════════════════════════════════════════════════════

#[near_bindgen]
#[derive(BorshSerialize, BorshDeserialize, PanicOnDefault)]
pub struct SmrkGuegapContract {
    owner:     AccountId,
    jobs:      UnorderedMap<String, JobRecord>,
    job_count: u64,
}

// ── Helper: build a failed ConfirmResp without stored data ───────────────────
fn fail_confirm(run_id: String, msg: &str) -> ConfirmResp {
    ConfirmResp {
        ok: false, run_id,
        h1_h2_proxy: "fail".into(),
        h1_pass: false, h2_pass: false,
        r_mean: 0.0, delta1: 0.0,
        math_sha256_hex: None,
        confirmed_block_height: 0,
        confirmed_code_hash_hex: String::new(),
        message: msg.into(),
    }
}

#[near_bindgen]
impl SmrkGuegapContract {

    // ── Constructor ──────────────────────────────────────────────────────────

    #[init]
    pub fn new(owner: AccountId) -> Self {
        Self { owner, jobs: UnorderedMap::new(b"j"), job_count: 0 }
    }

    // ── Job submission ────────────────────────────────────────────────────────
    //
    // FIX-2: #[payable]. Raw input NOT stored. Only hash + len.
    // FIX-5: storage cost via env::storage_byte_cost() — correct on any network.
    // FIX-6: input_sha256_hex validated: length 64, chars [0-9a-f], len > 0.
    // Idempotent: duplicate hash returns existing record, refunds deposit.

    #[payable]
    pub fn submit_job(
        &mut self,
        input_sha256_hex: String,
        input_len:        u32,
        icp_run_id:       Option<String>,
        icp_math_sha256_hex: Option<String>,
    ) -> SubmitJobResp {

        // FIX-6: validate hex format
        if let Err(e) = validate_sha256_hex(&input_sha256_hex) {
            env::panic_str(&e);
        }
        // FIX-6: validate length bounds
        assert!(input_len > 0, "input_len must be > 0.");
        assert!(
            input_len <= MAX_INPUT_LEN,
            "input_len {} exceeds MAX_INPUT_LEN {}.", input_len, MAX_INPUT_LEN
        );

        let run_id = input_sha256_hex.clone();

        // Idempotent: refund deposit and return existing record
        if self.jobs.contains_key(&run_id) {
            let deposit = env::attached_deposit();
            if deposit > NearToken::from_yoctonear(0) {
                Promise::new(env::predecessor_account_id()).transfer(deposit);
            }
            return SubmitJobResp {
                run_id, input_sha256_hex, input_len,
                storage_used: 0, deposit_returned: true,
            };
        }

        let storage_before = env::storage_usage();

        self.jobs.insert(run_id.clone(), JobRecord {
            run_id:               run_id.clone(),
            created_at_ns:        env::block_timestamp(),
            submitted_by:         env::predecessor_account_id(),
            status:               JobStatus::Queued,
            input_sha256_hex:     input_sha256_hex.clone(),
            input_len,
            icp_run_id,
            icp_math_sha256_hex,
            output_sha256_hex:    None,
            math_sha256_hex:      None,
            commit_hash_hex:      None,
            error:                None,
            h1_h2_proxy:          None,   // FIX-4
            h1_pass:              None,   // FIX-4
            h2_pass:              None,   // FIX-4
            r_mean:               None,   // FIX-4
            delta1:               None,   // FIX-4
            confirmed_block_height:  None,  // FIX-1
            confirmed_code_hash_hex: None,  // FIX-1
            anchor_commitment_hex: None,
            kaspa_txid:           None,
            anchor_status:        None,
            anchor_note:          None,
        });
        self.job_count += 1;

        let storage_after = env::storage_usage();
        let storage_used  = storage_after.saturating_sub(storage_before);

        // FIX-5: runtime storage cost — correct on testnet, mainnet, localnet
        let byte_cost    = env::storage_byte_cost().as_yoctonear();
        let required     = storage_used as u128 * byte_cost;
        let attached     = env::attached_deposit().as_yoctonear();

        // Standard NEAR pattern: insert first, then check deposit.
        // On panic the state reverts, so the insert is rolled back.
        assert!(
            attached >= required,
            "Insufficient storage deposit. Required: {} yoctoNEAR ({} bytes × {} per byte). Attached: {}.",
            required, storage_used, byte_cost, attached
        );

        let excess = attached.saturating_sub(required);
        if excess > 0 {
            Promise::new(env::predecessor_account_id())
                .transfer(NearToken::from_yoctonear(excess));
        }

        SubmitJobResp {
            run_id, input_sha256_hex, input_len,
            storage_used, deposit_returned: excess > 0,
        }
    }

    // ── NEAR confirmation ─────────────────────────────────────────────────────
    //
    // PERF-1: single compute_smrk() call.
    // FIX-1:  env snapshot stored in JobRecord.
    // FIX-4:  h1_pass/h2_pass/r_mean/delta1 stored in JobRecord.
    // FIX-7:  input size capped at MAX_INPUT_BYTES.
    // FIX-8:  compute_smrk() returning Err → JobStatus::Failed saved, no panic.

    pub fn run_confirmation(&mut self, run_id: String, input: Vec<u8>) -> ConfirmResp {
        // FIX-7: guard against oversized input before any computation
        assert!(
            input.len() <= MAX_INPUT_BYTES,
            "input length {} exceeds MAX_INPUT_BYTES {}.", input.len(), MAX_INPUT_BYTES
        );

        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None    => return fail_confirm(run_id, "Job not found. Call submit_job first."),
        };

        // FIX-2: verify supplied input matches committed hash
        let supplied_hash = hex::encode(sha256_raw(&input));
        if supplied_hash != job.input_sha256_hex {
            return fail_confirm(run_id,
                "Input mismatch: SHA-256(supplied_input) != stored input_sha256_hex. \
                 Supply the exact bytes used when calling submit_job.");
        }

        // FIX-4: "Already confirmed" returns cached spectral values — not constants
        if job.status == JobStatus::Done {
            return ConfirmResp {
                ok:             true,
                run_id:         run_id.clone(),
                h1_h2_proxy:    job.h1_h2_proxy.clone().unwrap_or_else(|| "fail".into()),
                h1_pass:        job.h1_pass.unwrap_or(false),
                h2_pass:        job.h2_pass.unwrap_or(false),
                r_mean:         job.r_mean.unwrap_or(0.0),
                delta1:         job.delta1.unwrap_or(0.0),
                math_sha256_hex: job.math_sha256_hex.clone(),
                confirmed_block_height:  job.confirmed_block_height.unwrap_or(0),
                confirmed_code_hash_hex: job.confirmed_code_hash_hex.clone().unwrap_or_default(),
                message: "Already confirmed (cached result).".into(),
            };
        }

        if job.status == JobStatus::Failed {
            return fail_confirm(run_id,
                &format!("Job previously failed: {}", job.error.as_deref().unwrap_or("?")));
        }

        job.status = JobStatus::Running;

        // FIX-1: capture env snapshot at this block — stored for verify_output()
        let code_hash_hex = hex::encode(env::code_hash());
        let block_height  = env::block_height();

        // PERF-1: single Jacobi run — FIX-8: handle Err without panic
        let res = match compute_smrk(&input) {
            Ok(r)  => r,
            Err(e) => {
                // FIX-8: save Failed with error message → audit trail preserved
                job.status = JobStatus::Failed;
                job.error  = Some(e.clone());
                self.jobs.insert(run_id.clone(), job);
                return fail_confirm(run_id, &format!("Computation failed: {}", e));
            }
        };

        // Build canonical JSON outputs
        let (_, math_hash)    = math_payload_from_result(&res, &job.input_sha256_hex);
        let math_sha_hex      = hex::encode(math_hash);
        let (_, audit_hash)   = audit_json_from_result(
            &res, &job.input_sha256_hex, &math_sha_hex, &code_hash_hex, block_height,
        );
        let out_sha_hex = hex::encode(audit_hash);

        // PERF-2 + FIX-4: derive all result fields from res booleans
        let h1_pass     = res.h1_pass;
        let h2_pass     = res.h2_pass;
        let h1h2_str    = if h1_pass && h2_pass { "pass" } else { "fail" };
        let r_mean      = res.r_mean;
        let delta1      = res.delta1;

        // FIX-4: persist spectral result
        job.h1_h2_proxy  = Some(h1h2_str.into());
        job.h1_pass      = Some(h1_pass);
        job.h2_pass      = Some(h2_pass);
        job.r_mean       = Some(r_mean);
        job.delta1       = Some(delta1);

        job.status               = JobStatus::Done;
        job.output_sha256_hex    = Some(out_sha_hex.clone());
        job.math_sha256_hex      = Some(math_sha_hex.clone());
        job.commit_hash_hex      = Some(out_sha_hex);
        job.error                = None;
        job.confirmed_block_height  = Some(block_height);   // FIX-1
        job.confirmed_code_hash_hex = Some(code_hash_hex.clone()); // FIX-1

        self.jobs.insert(run_id.clone(), job);

        ConfirmResp {
            ok: true, run_id,
            h1_h2_proxy:             h1h2_str.into(),
            h1_pass, h2_pass, r_mean, delta1,
            math_sha256_hex:         Some(math_sha_hex),
            confirmed_block_height:  block_height,
            confirmed_code_hash_hex: code_hash_hex,
            message: "N=32 confirmation complete.".into(),
        }
    }

    // ── Verification ──────────────────────────────────────────────────────────

    /// Recomputes math_payload and checks against stored math_sha256_hex.
    /// Always matches (stable across redeployments). Caller supplies input.
    pub fn verify_math(&self, run_id: String, input: Vec<u8>) -> VerifyResp {
        assert!(input.len() <= MAX_INPUT_BYTES,
            "input too large: {} > {}", input.len(), MAX_INPUT_BYTES);

        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None    => return VerifyResp { ok: false, matches: false,
                message: "Job not found.".into(), stored_hex: None, computed_hex: None },
        };
        let stored = match &job.math_sha256_hex {
            Some(h) => h.clone(),
            None    => return VerifyResp { ok: false, matches: false,
                message: "No math_sha256_hex yet. Run run_confirmation first.".into(),
                stored_hex: None, computed_hex: None },
        };
        let supplied_hash = hex::encode(sha256_raw(&input));
        if supplied_hash != job.input_sha256_hex {
            return VerifyResp { ok: false, matches: false,
                message: "Input mismatch: SHA-256(supplied) != stored input_sha256_hex.".into(),
                stored_hex: Some(stored), computed_hex: None };
        }
        let res = match compute_smrk(&input) {
            Ok(r)  => r,
            Err(e) => return VerifyResp { ok: false, matches: false,
                message: format!("Recomputation failed: {}", e),
                stored_hex: Some(stored), computed_hex: None },
        };
        let (_, hash) = math_payload_from_result(&res, &job.input_sha256_hex);
        let computed  = hex::encode(hash);
        let matches   = stored == computed;
        VerifyResp {
            ok: true, matches,
            message: if matches { "math_sha256 verified." }
                     else       { "math_sha256 MISMATCH — unexpected." }.into(),
            stored_hex: Some(stored), computed_hex: Some(computed),
        }
    }

    /// Recomputes audit JSON using STORED block_height + code_hash_hex (FIX-1).
    /// Matches on same deployment; intentionally differs after redeployment.
    pub fn verify_output(&self, run_id: String, input: Vec<u8>) -> VerifyResp {
        assert!(input.len() <= MAX_INPUT_BYTES,
            "input too large: {} > {}", input.len(), MAX_INPUT_BYTES);

        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None    => return VerifyResp { ok: false, matches: false,
                message: "Job not found.".into(), stored_hex: None, computed_hex: None },
        };
        let stored = match &job.output_sha256_hex {
            Some(h) => h.clone(),
            None    => return VerifyResp { ok: false, matches: false,
                message: "No output_sha256_hex yet. Run run_confirmation first.".into(),
                stored_hex: None, computed_hex: None },
        };
        // FIX-1: use stored env snapshot, NOT live env::*
        let code_hash_hex = match &job.confirmed_code_hash_hex {
            Some(h) => h.clone(),
            None    => return VerifyResp { ok: false, matches: false,
                message: "No confirmed_code_hash_hex (job confirmed before v0.3.0?).".into(),
                stored_hex: Some(stored), computed_hex: None },
        };
        let block_height = match job.confirmed_block_height {
            Some(h) => h,
            None    => return VerifyResp { ok: false, matches: false,
                message: "No confirmed_block_height (job confirmed before v0.3.0?).".into(),
                stored_hex: Some(stored), computed_hex: None },
        };
        let supplied_hash = hex::encode(sha256_raw(&input));
        if supplied_hash != job.input_sha256_hex {
            return VerifyResp { ok: false, matches: false,
                message: "Input mismatch: SHA-256(supplied) != stored input_sha256_hex.".into(),
                stored_hex: Some(stored), computed_hex: None };
        }
        let res = match compute_smrk(&input) {
            Ok(r)  => r,
            Err(e) => return VerifyResp { ok: false, matches: false,
                message: format!("Recomputation failed: {}", e),
                stored_hex: Some(stored), computed_hex: None },
        };
        let (_, math_hash)   = math_payload_from_result(&res, &job.input_sha256_hex);
        let math_sha_hex     = hex::encode(math_hash);
        let (audit_bytes, _) = audit_json_from_result(
            &res, &job.input_sha256_hex, &math_sha_hex, &code_hash_hex, block_height,
        );
        let computed = sha256_hex(&audit_bytes);
        let matches  = stored == computed;
        VerifyResp {
            ok: true, matches,
            message: if matches {
                "output_sha256 verified (same deployment).".into()
            } else {
                "Mismatch — expected after contract redeployment (code_hash changed).".into()
            },
            stored_hex: Some(stored), computed_hex: Some(computed),
        }
    }

    // ── Kaspa anchoring ───────────────────────────────────────────────────────

    /// Step 1 — compute anchor_commitment. Requires status == Done.
    pub fn anchor_job(&mut self, run_id: String, plateau_flag: bool) -> AnchorResp {
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None    => return AnchorResp { ok: false, anchor_commitment_hex: None,
                message: "Job not found.".into() },
        };
        if job.status != JobStatus::Done {
            return AnchorResp { ok: false, anchor_commitment_hex: None,
                message: "Job not Done. Run run_confirmation first.".into() };
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
                "Anchor commitment stored. Embed as Kaspa tx payload, \
then call record_kaspa_txid(run_id, txid). Payload: {}", hex
            ),
        }
    }

    /// Step 2 — record Kaspa txid after broadcast.
    ///
    /// FIX-3: one-time write + caller ACL.
    /// FIX-9: state-machine guard — anchor_status must be Pending.
    pub fn record_kaspa_txid(&mut self, run_id: String, kaspa_txid: String) -> bool {
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None    => return false,
        };

        // FIX-9: state machine — must be Pending
        assert!(
            job.anchor_status == Some(AnchorStatus::Pending),
            "anchor_status must be Pending. Call anchor_job() first, \
             and do not call record_kaspa_txid() twice."
        );

        // FIX-3(a): one-time write
        assert!(
            job.kaspa_txid.is_none(),
            "kaspa_txid already recorded. Cannot overwrite."
        );

        // FIX-3(b): ACL — submitted_by or owner
        let caller = env::predecessor_account_id();
        assert!(
            caller == job.submitted_by || caller == self.owner,
            "Only the job submitter or contract owner may record kaspa_txid."
        );

        assert!(job.anchor_commitment_hex.is_some(),
            "No anchor_commitment_hex. Call anchor_job() first.");

        let txid = kaspa_txid.clone();
        job.kaspa_txid    = Some(kaspa_txid);
        job.anchor_status = Some(AnchorStatus::Recorded);
        job.anchor_note   = Some(format!("Kaspa txid recorded: {}. Awaiting verification.", txid));
        self.jobs.insert(run_id, job);
        true
    }

    /// Step 3 — record off-chain Kaspa verification result.
    /// Only owner. For trustless verification, integrate an oracle.
    pub fn submit_kaspa_verification(
        &mut self, run_id: String, verified: bool, note: String,
    ) -> bool {
        let caller = env::predecessor_account_id();
        assert!(caller == self.owner,
            "Only the contract owner can submit Kaspa verification results.");
        let mut job = match self.jobs.get(&run_id) {
            Some(j) => j.clone(),
            None    => return false,
        };
        job.anchor_status = Some(if verified { AnchorStatus::Verified } else { AnchorStatus::Failed });
        job.anchor_note   = Some(note);
        self.jobs.insert(run_id, job);
        true
    }

    // ── Cross-chain reference comparison ──────────────────────────────────────
    //
    // API-1: informational only. ICP N=16 vs NEAR N=32 → hashes will differ.

    pub fn compare_icp_reference(&self, run_id: String) -> IcpCompareResp {
        let job = match self.jobs.get(&run_id) {
            Some(j) => j,
            None    => return IcpCompareResp {
                ok: false, near_math_sha256_hex: None, icp_math_sha256_hex: None,
                near_probe_version: PROBE_VERSION.into(), icp_probe_version: None,
                hashes_match: false, note: "Job not found.".into(),
            },
        };
        let near_math = match &job.math_sha256_hex {
            Some(h) => h.clone(),
            None    => return IcpCompareResp {
                ok: false, near_math_sha256_hex: None, icp_math_sha256_hex: None,
                near_probe_version: PROBE_VERSION.into(), icp_probe_version: None,
                hashes_match: false,
                note: "No NEAR math_sha256_hex yet. Run run_confirmation first.".into(),
            },
        };
        let icp_math = match &job.icp_math_sha256_hex {
            Some(h) => h.clone(),
            None    => return IcpCompareResp {
                ok: true,
                near_math_sha256_hex: Some(near_math), icp_math_sha256_hex: None,
                near_probe_version: PROBE_VERSION.into(), icp_probe_version: None,
                hashes_match: false,
                note: "No ICP hash provided at submit time. \
                       Supply icp_math_sha256_hex in submit_job() for cross-chain audit.".into(),
            },
        };
        let hashes_match = near_math == icp_math;
        IcpCompareResp {
            ok: true,
            near_math_sha256_hex: Some(near_math),
            icp_math_sha256_hex:  Some(icp_math),
            near_probe_version:   PROBE_VERSION.into(),
            icp_probe_version:    Some("smrk-guegap-icp-v3".into()),
            hashes_match,
            note: if hashes_match {
                "Hashes match — UNEXPECTED for N=16 vs N=32. Investigate.".into()
            } else {
                "Hashes differ — expected. ICP N=16 vs NEAR N=32 produce distinct \
                 math_sha256. Compare eigs_low / delta1 / r_mean for manual audit.".into()
            },
        }
    }

    // ── Read methods (view, no gas) ───────────────────────────────────────────

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
            contract_version:   CONTRACT_VERSION.into(),
            probe_version:      PROBE_VERSION.into(),
            n_matrix:           NMAT as u32,
            code_hash_hex:      hex::encode(env::code_hash()),
            block_height:       env::block_height(),
            block_timestamp_ns: env::block_timestamp(),
        }
    }

    pub fn get_job_count(&self) -> u64 { self.job_count }
    pub fn get_owner(&self)    -> AccountId { self.owner.clone() }
}
