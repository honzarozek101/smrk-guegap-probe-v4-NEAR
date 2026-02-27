# SMRK–GUEGAP Probe v4 — NEAR Confirmation Tier

**Second-layer confirmation for the NEAR Protocol**,  
as defined in the *SMRK–GUEGAP Probe v4* specification:  
*Deterministic Spectral Gap Probing on ICP and NEAR with Kaspa Time Anchoring.*

---

## Overview

This contract implements the **Quansistor confirmation tier** — a deterministic on-chain spectral experiment that constructs a 32×32 Hermitian matrix with arithmetic prime-shift couplings, diagonalises it using a cyclic complex Jacobi solver, and verifies two spectral hypotheses:

- **H1** — bulk level-spacing ratios ⟨r⟩ consistent with GUE random-matrix statistics
- **H2** — a stable low-energy spectral gap Δ₁

Every result is cryptographically committed via two SHA-256 hashes and time-anchored to the Kaspa proof-of-work blockchain. The computation is fully deterministic and publicly reproducible from the input alone.

---

## Parameters

| Parameter       | Value                                              |
|-----------------|----------------------------------------------------|
| `probe_version` | `smrk-guegap-near-v1`                              |
| Matrix **N**    | **32** (spec N=512 → scaled to NEAR gas budget)    |
| Role            | Confirmation tier — second layer after ICP screening (N=16) |
| H1 acceptance   | \|⟨r⟩ − 0.5996\| ≤ 0.10                           |
| H2 acceptance   | Δ₁ ≥ 1×10⁻⁴                                       |
| Bulk window     | i₀ = 6 … i₁ = 23 → ~18 ratio samples              |
| Gas (confirm)   | ~30–50 TGas (limit: 300 TGas / tx)                 |

---

## Architecture

A **single NEAR contract** unifies both the registry and the compute roles.  
Cross-contract calls on NEAR are asynchronous (promise chains with callbacks) — for this use-case a single-contract design is simpler, fully auditable, and avoids callback debt.

```
  caller
    │  submit_job(input, icp_run_id?, icp_math_sha256_hex?)
    ▼
  ┌────────────────────────────────────────────────┐
  │  SmrkGuegapContract                            │
  │  storage: UnorderedMap<run_id, JobRecord>      │
  │                                                │
  │  run_confirmation()  ← Jacobi N=32 on-chain   │
  │  verify_math()       ← always stable           │
  │  verify_output()     ← changes on redeploy     │
  │  anchor_job()        ← SHA-256 commitment      │
  │  record_kaspa_txid() ← store txid              │
  │  submit_kaspa_verification()  ← owner only     │
  │  verify_icp_cross_chain()     ← audit          │
  └────────────────────────────────────────────────┘
         │
         │  Kaspa (external broadcast + off-chain verify)
         ▼
  api.kaspa.org
```

> **Note on Kaspa HTTP:** NEAR contracts cannot make outbound HTTP calls.  
> Kaspa txid verification is performed off-chain and the result is written back  
> to the contract via `submit_kaspa_verification()` (owner-gated).  
> For trustless production use, integrate an oracle contract.

---

## ICP vs NEAR — Key Differences

| Aspect                      | ICP Screening (N=16)         | NEAR Confirmation (N=32)       |
|-----------------------------|------------------------------|--------------------------------|
| Architecture                | 2 canisters                  | 1 contract                     |
| `probe_version`             | `smrk-guegap-icp-v3`         | `smrk-guegap-near-v1`          |
| Bulk window i₀..i₁          | 3..10 (~8 samples)           | 6..23 (~18 samples)            |
| H1 acceptance band          | ±0.15                        | ±0.10 (tighter)                |
| Kaspa HTTP verification     | In-contract HTTP outcall     | Off-chain / oracle             |
| Audit build metadata        | `git_commit` + `canister_version` | `code_hash_hex` + `block_height` |
| Cross-chain `math_sha256`   | —                            | Differs by design (N=16 ≠ N=32) |

Cross-chain correlation: `verify_icp_cross_chain()` shows both `math_sha256_hex` values side-by-side. Because N=16 ≠ N=32, the hashes intentionally differ — both results remain independently reproducible from the same input.

---

## Operator Family — SMRK–GUEGAP (spec §2)

The 32×32 Hermitian matrix H = D + C is constructed as follows:

| Component          | Definition                                                          |
|--------------------|---------------------------------------------------------------------|
| Diagonal D_nn      | λ · ln(n+1),  n = 1..32                                            |
| Prime shifts p_m   | {2, 3, 5, 7, 11, 13, 17, 19}  (M = 8 shifts)                      |
| Cyclic shift σ_m(n)| (n + p_m) mod 32 — periodic boundary conditions                    |
| Amplitude A_m(n)   | κ / ln(1+p_m) · (1 + 0.5 · λ_L(n))                                |
| Liouville λ_L(n)   | (−1)^Ω(n),  Ω = total prime factors with multiplicity              |
| Phase φ_m(n)       | 2π · α · n + φ_offset,  α = (√5−1)/2 ≈ 0.6180 (golden gauge)     |
| Coupling w_m(n)    | A_m(n) · e^{i·φ_m(n)};  C[n][σ] += w,  C[σ][n] += conj(w)       |
| Parameters         | κ = 0.05 + h[0]/255·0.15,  λ = h[1]/255·0.20,  φ = h[2]/255·0.20 |

Parameters κ, λ, φ_offset are deterministically derived from SHA-256(input) bytes 0–2. Each distinct input explores a unique point in (κ, λ, φ) parameter space.

---

## Build & Deploy

### Prerequisites

```bash
rustup target add wasm32-unknown-unknown
cargo install near-cli-rs        # or: npm install -g near-cli
```

### Build

```bash
cargo build --target wasm32-unknown-unknown --release
# output: target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm
```

### Deploy — Testnet

```bash
near create-account smrk-probe.testnet --useFaucet

near deploy smrk-probe.testnet \
  target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm \
  --initFunction new \
  --initArgs '{"owner": "your-account.testnet"}'
```

### Deploy — Mainnet

```bash
near deploy smrk-probe.near \
  target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm \
  --initFunction new \
  --initArgs '{"owner": "your-account.near"}' \
  --networkId mainnet
```

---

## API Reference

### Initialisation & metadata

| Method                | Type   | Description                                                         |
|-----------------------|--------|---------------------------------------------------------------------|
| `new(owner)`          | init   | Constructor — must be called as `--initFunction` on deploy          |
| `get_owner()`         | view   | Returns the contract owner `AccountId`                              |
| `get_contract_meta()` | view   | `probe_version`, `n_matrix`, `code_hash_hex`, `block_height`, etc.  |
| `get_job_count()`     | view   | Total number of jobs submitted since deploy                         |

### Job lifecycle

| Method                                                      | Type   | Description                                                                 |
|-------------------------------------------------------------|--------|-----------------------------------------------------------------------------|
| `submit_job(input, icp_run_id?, icp_math_sha256_hex?)`      | call   | Derives `run_id = SHA-256(input)`. Idempotent — duplicate returns existing `run_id` without re-insertion. |
| `run_confirmation(run_id)`                                  | call   | Full N=32 Jacobi computation on-chain. Returns `h1_h2_proxy` + `math_sha256_hex`. Gas: ~30–50 TGas. |
| `get_job(run_id)`                                           | view   | Returns the complete `JobRecord`                                            |

### Verification

| Method                           | Type | Description                                                                           |
|----------------------------------|------|---------------------------------------------------------------------------------------|
| `verify_math(run_id)`            | view | Recomputes `math_payload` and compares SHA-256 to stored `math_sha256_hex`. Always matches — stable across redeploys. |
| `verify_output(run_id)`          | view | Recomputes audit JSON (current `code_hash` + `block_height`). Matches only on the same deploy — intentionally changes on redeploy. |
| `verify_icp_cross_chain(run_id)` | view | Shows NEAR `math_sha256` vs ICP `math_sha256` side-by-side for manual audit.         |

### Kaspa anchoring

| Method                                              | Type | Description                                                                 |
|-----------------------------------------------------|------|-----------------------------------------------------------------------------|
| `anchor_job(run_id, plateau_flag)`                  | call | Computes `anchor_commitment` on-chain, stores hex. Returns payload for Kaspa tx. Requires `status = Done`. |
| `record_kaspa_txid(run_id, kaspa_txid)`             | call | Stores `kaspa_txid` after external broadcast. `anchor_status`: Pending → Recorded. |
| `submit_kaspa_verification(run_id, verified, note)` | call | **Owner only.** Records off-chain / oracle verification result. `anchor_status`: → Verified or Failed. |
| `get_anchor_info(run_id)`                           | view | Lightweight read — anchor fields only, without the full output blob.        |

---

## Complete Workflow

```bash
CONTRACT="smrk-probe.testnet"
ACCOUNT="your-account.testnet"

# A — Submit job
near call $CONTRACT submit_job \
  '{"input": [104,101,108,108,111], "icp_run_id": null, "icp_math_sha256_hex": null}' \
  --accountId $ACCOUNT
# → {"run_id": "<hex>", "input_sha256_hex": "<hex>"}

RUN_ID="<run_id from above>"

# B — N=32 confirmation (Jacobi on-chain)
near call $CONTRACT run_confirmation \
  '{"run_id": "'$RUN_ID'"}' \
  --accountId $ACCOUNT --gas 100000000000000

# C — Verify math (always matches, stable across redeploys)
near view $CONTRACT verify_math '{"run_id": "'$RUN_ID'"}'

# D — Verify output (matches only on same deploy)
near view $CONTRACT verify_output '{"run_id": "'$RUN_ID'"}'

# E — Cross-chain audit vs ICP tier
near view $CONTRACT verify_icp_cross_chain '{"run_id": "'$RUN_ID'"}'

# F — Compute Kaspa anchor commitment
near call $CONTRACT anchor_job \
  '{"run_id": "'$RUN_ID'", "plateau_flag": true}' \
  --accountId $ACCOUNT
# → {"anchor_commitment_hex": "<64 hex chars>", ...}

# G — Broadcast to Kaspa (external — use Kaspa SDK)
#     embed anchor_commitment_hex as the transaction payload

# H — Record Kaspa txid
near call $CONTRACT record_kaspa_txid \
  '{"run_id": "'$RUN_ID'", "kaspa_txid": "<kaspa_txid>"}' \
  --accountId $ACCOUNT

# I — Submit verification result (owner only)
near call $CONTRACT submit_kaspa_verification \
  '{"run_id": "'$RUN_ID'", "verified": true, "note": "Verified via api.kaspa.org"}' \
  --accountId $ACCOUNT

# J — Read anchor state
near view $CONTRACT get_anchor_info '{"run_id": "'$RUN_ID'"}'

# — Contract metadata
near view $CONTRACT get_contract_meta '{}'
```

---

## Dual-Hash Audit Model

Each confirmed job produces two independent SHA-256 commitments with different stability guarantees:

| Hash                   | Stable?          | Depends on                                          |
|------------------------|------------------|-----------------------------------------------------|
| `math_sha256_hex`      | ✅ Always         | Input bytes + eigenvalues only. No environment metadata. |
| `output_sha256_hex`    | ❌ Changes on redeploy | `code_hash_hex` (WASM SHA-256) + `block_height` + `math_sha256_hex` |
| `anchor_commitment_hex`| ✅ Always         | `run_id` + `commit_hash_hex` + `plateau_flag`        |

`math_sha256_hex` is the scientific commitment — reproducible forever from the input alone.  
`output_sha256_hex` is the audit trail — it captures the exact runtime environment at the time of computation.

### Kaspa Anchor Commitment Formula (spec §D.1.3)

```
anchor_commitment = SHA256("QFC|v4|TWISTOR|" ‖ run_id ‖ commit_hash_hex ‖ plateau_flag)
```

`plateau_flag` is serialised as the UTF-8 string `"true"` or `"false"` — consistent with the ICP tier, ensuring identical commitments from the same inputs across both chains.

---

## Canonical JSON Formats

### math payload — stable across redeploys

```json
{
  "H1_H2_proxy": "pass",
  "N": 32,
  "bulk_r": { "count": 18, "r_mean": 5.9123456789012345e-01 },
  "eigs_low": [E0, E1, E2, E3, E4],
  "gap":      { "delta1": 3.2100000000000000e-02 },
  "input_sha256_hex": "...",
  "probe_version": "smrk-guegap-near-v1"
}
```

### audit output — changes on redeploy

```json
{
  "H1_H2_proxy": "pass",
  "N": 32,
  "bulk_r": { "count": 18, "r_mean": 5.9123456789012345e-01 },
  "eigs_low": [E0, E1, E2, E3, E4],
  "gap":      { "delta1": 3.2100000000000000e-02 },
  "meta": {
    "block_height":     12345678,
    "code_hash_hex":    "<sha256 of wasm bytecode>",
    "contract_version": "0.1.0",
    "input_sha256_hex": "...",
    "math_sha256_hex":  "..."
  },
  "probe_version": "smrk-guegap-near-v1"
}
```

All keys are lexicographically ordered. All floating-point values use 16-decimal scientific notation (spec §A.1.1) to guarantee bitwise reproducibility.

---

## Gas Estimates

| Operation                          | Estimated gas   |
|------------------------------------|-----------------|
| `submit_job`                       | ~5 TGas         |
| `run_confirmation` (Jacobi N=32)   | ~30–50 TGas     |
| `verify_math`                      | ~30 TGas        |
| `anchor_job`                       | ~5 TGas         |

Use `--gas 100000000000000` (100 TGas) for `run_confirmation`.

---

## Security Notes

- `submit_kaspa_verification` is restricted to the **contract owner**. For production, consider replacing the owner gate with an oracle contract that performs on-chain HTTP verification.
- The full `input` blob is stored in contract state. For large inputs, store only `input_sha256_hex` on-chain and pass the raw data ad-hoc to view functions.
- Jobs are **idempotent**: calling `submit_job` twice with the same input returns the same `run_id` without duplication or state overwrite.
- `anchor_job` requires `status = Done` — it is not possible to anchor an incomplete or failed job.
- `verify_output` depends on the current `block_height` and will return a different hash on every block. Use `verify_math` for reproducible long-term verification.


License
=======

© 2026 101research.group

----------------------------------------------------------------------
Documents and Textual Materials
----------------------------------------------------------------------

All documents, whitepapers, PDFs, and written materials distributed
via this canister are licensed under the Creative Commons
Attribution–NonCommercial–NoDerivatives 4.0 International License
(CC BY–NC–ND 4.0).

You are free to:
- Read, download, and share the materials
- Cite the materials with proper attribution

Under the following conditions:
- Attribution — Appropriate credit must be given to the author(s)
- NonCommercial — The material may not be used for commercial purposes
- NoDerivatives — The material may not be modified, adapted, or remixed

Modification, creation of derivative works, commercial use,
rebranding, or redistribution as part of proprietary systems
requires explicit written permission from the author.

License text:
https://creativecommons.org/licenses/by-nc-nd/4.0/

----------------------------------------------------------------------
Concepts and Frameworks
----------------------------------------------------------------------

The concepts, terminology, and names including (but not limited to):

- Quansistor
- Quansistor Field Computing (QFC)
- SMRK Hamiltonian

are provided for research and educational use with attribution.
Any commercialization, rebranding, or incorporation into
commercial or proprietary systems requires explicit permission
from the author.

----------------------------------------------------------------------
Metadata and AI Navigation Files
----------------------------------------------------------------------

Machine-readable metadata and AI navigation files, including but not
limited to:

- qfc.map.json
- llms.txt
- index.yaml
- related metadata descriptors

are licensed under the Creative Commons Attribution 4.0 International
License (CC BY 4.0), unless explicitly stated otherwise.

----------------------------------------------------------------------
General Note
----------------------------------------------------------------------

This licensing structure reflects the research-oriented and
pre-publication nature of the Quansistor Field Computing project.

For licensing inquiries, contact:
101research.group@seznam.cz
