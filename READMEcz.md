# SMRK–GUEGAP Probe v4 — NEAR Confirmation Tier

Confirmation layer pro **NEAR Protocol** jako definováno ve specifikaci
*SMRK–GUEGAP Probe v4 — Deterministic Spectral Gap Probing on ICP and NEAR with Kaspa Time Anchoring*.

| Parametr | Hodnota |
|---|---|
| `probe_version` | `smrk-guegap-near-v1` |
| N (matice) | **32** (spec N=512 → škálováno pro gas budget) |
| Role | Confirmation tier — druhá vrstva po ICP screening (N=16) |
| H1 band | ±0.10 od r_GUE = 0.5996 |
| H2 min gap | Δ1 ≥ 1×10⁻⁴ |
| Bulk okno | i₀=6, i₁=23 → ~18 ratio vzorků |

---

## Architektura

Jeden NEAR kontrakt (registry + compute unified). Na rozdíl od ICP implementace,
cross-contract calls na NEAR jsou async s promise callbacky — pro tuto use-case
je single-contract design čistší a auditovatelný.

```
  caller
    │  submit_job(input, icp_run_id?, icp_math_sha256_hex?)
    ▼
  ┌─────────────────────────────────────────────┐
  │  SmrkGuegapContract                         │
  │  storage: UnorderedMap<run_id, JobRecord>   │
  │                                             │
  │  run_confirmation()    ← Jacobi N=32        │
  │  verify_math()         ← stabilní hash      │
  │  verify_output()       ← mění se po redeploy│
  │  anchor_job()          ← SHA-256 commitment │
  │  record_kaspa_txid()   ← uloží txid         │
  │  submit_kaspa_verification()  ← owner only  │
  │  verify_icp_cross_chain()     ← audit       │
  └─────────────────────────────────────────────┘
         │
         │  Kaspa (externí broadcast + off-chain verify)
         ▼
  api.kaspa.org  (HTTP — mimo NEAR kontrakt)
```

> **Kaspa HTTP**: NEAR kontrakty nemohou dělat HTTP calls. Verifikace Kaspa txid
> probíhá off-chain a výsledek se submittuje přes `submit_kaspa_verification()`
> vlastníkem kontraktu, nebo přes oracle kontrakt.

---

## Rozdíly ICP vs NEAR

| Aspekt | ICP (N=16) | NEAR (N=32) |
|---|---|---|
| Architektura | 2 canistry | 1 kontrakt |
| probe_version | smrk-guegap-icp-v3 | smrk-guegap-near-v1 |
| Bulk i₀..i₁ | 3..10 (~8 vzorků) | 6..23 (~18 vzorků) |
| H1 band | ±0.15 | ±0.10 (užší) |
| Kaspa HTTP | verify_kaspa_anchor() in-contract | off-chain / oracle |
| Build meta v auditu | git_commit + canister_version | code_hash + block_height |
| math_sha256 cross-chain | ❌ různé N → různé hashes | — |

Cross-chain korelace: `verify_icp_cross_chain()` zobrazí oba math_sha256 hashů
side-by-side. Protože N=16 ≠ N=32, hashes se lišit **záměrně** — oba výsledky
jsou ale deterministicky reprodukovatelné ze stejného vstupního vstupu.

---

## Build & deploy

### Prerekvizity

```bash
rustup target add wasm32-unknown-unknown
cargo install near-cli-rs   # nebo: npm install -g near-cli
```

### Build

```bash
cargo build --target wasm32-unknown-unknown --release
# wasm binary: target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm
```

### Deploy (testnet)

```bash
near create-account smrk-probe.testnet --useFaucet
near deploy smrk-probe.testnet \
  target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm \
  --initFunction new \
  --initArgs '{"owner": "your-account.testnet"}'
```

### Deploy (mainnet)

```bash
near deploy smrk-probe.near \
  target/wasm32-unknown-unknown/release/smrk_guegap_near.wasm \
  --initFunction new \
  --initArgs '{"owner": "your-account.near"}' \
  --networkId mainnet
```

---

## Kompletní workflow

```bash
CONTRACT="smrk-probe.testnet"
ACCOUNT="your-account.testnet"

# 1. Submit (input jako base64 nebo hex přes near-cli)
near call $CONTRACT submit_job \
  '{"input": [104,101,108,108,111], "icp_run_id": null, "icp_math_sha256_hex": null}' \
  --accountId $ACCOUNT
# → {"run_id": "<hex>", "input_sha256_hex": "<hex>"}

RUN_ID="<run_id from above>"

# 2. N=32 confirmation
near call $CONTRACT run_confirmation \
  '{"run_id": "'$RUN_ID'"}' \
  --accountId $ACCOUNT --gas 100000000000000

# 3. Verify math (vždy matches)
near view $CONTRACT verify_math '{"run_id": "'$RUN_ID'"}'

# 4. Anchor — compute commitment
near call $CONTRACT anchor_job \
  '{"run_id": "'$RUN_ID'", "plateau_flag": true}' \
  --accountId $ACCOUNT
# → {"anchor_commitment_hex": "<64 hex chars>", ...}

# 5. Broadcast do Kaspa (externí — použij Kaspa SDK)
#    payload = anchor_commitment_hex

# 6. Record txid
near call $CONTRACT record_kaspa_txid \
  '{"run_id": "'$RUN_ID'", "kaspa_txid": "<kaspa_txid>"}' \
  --accountId $ACCOUNT

# 7. Submit verifikaci (owner only, nebo oracle)
near call $CONTRACT submit_kaspa_verification \
  '{"run_id": "'$RUN_ID'", "verified": true, "note": "Verified via api.kaspa.org"}' \
  --accountId $ACCOUNT

# 8. Read anchor state
near view $CONTRACT get_anchor_info '{"run_id": "'$RUN_ID'"}'

# 9. Cross-chain audit (ICP vs NEAR math hashes)
near view $CONTRACT verify_icp_cross_chain '{"run_id": "'$RUN_ID'"}'

# 10. Contract metadata
near view $CONTRACT get_contract_meta '{}'
```

---

## Canonical JSON

### math payload (stable)
```json
{
  "H1_H2_proxy": "pass|fail",
  "N": 32,
  "bulk_r": {"count": 18, "r_mean": 5.9123456789012345e-01},
  "eigs_low": [E0, E1, E2, E3, E4],
  "gap": {"delta1": 3.2100000000000000e-02},
  "input_sha256_hex": "...",
  "probe_version": "smrk-guegap-near-v1"
}
```

### audit output (mění se po redeploy)
```json
{
  "H1_H2_proxy": "pass|fail",
  "N": 32,
  "bulk_r": {...}, "eigs_low": [...], "gap": {...},
  "meta": {
    "block_height": 12345678,
    "code_hash_hex": "<wasm sha256>",
    "contract_version": "0.1.0",
    "input_sha256_hex": "...",
    "math_sha256_hex": "..."
  },
  "probe_version": "smrk-guegap-near-v1"
}
```

### Hash stabilita

| Hash | Stabilní | Závislý na |
|---|---|---|
| `math_sha256_hex` | ✅ Vždy | Pouze vstup + eigenvalua |
| `output_sha256_hex` | ❌ Po redeploy | code_hash + block_height |
| `anchor_commitment` | ✅ Vždy | run_id + commit_hash_hex + plateau_flag |

---

## Gas odhad

| Operace | Odhadovaný gas |
|---|---|
| `submit_job` | ~5 TGas |
| `run_confirmation` (Jacobi N=32) | ~30–50 TGas |
| `verify_math` | ~30 TGas |
| `anchor_job` | ~5 TGas |

Doporučený `--gas` flag pro `run_confirmation`: `100000000000000` (100 TGas).

---

## Bezpečnostní poznámky

- `submit_kaspa_verification` je omezena na contract owner — pro produkci
  zvažte oracle kontrakt s on-chain HTTP verifikací.
- Vstup `input` je uložen v kontraktním stavu — pro velké vstupy zvažte
  ukládání jen `input_sha256_hex` a předávání vstupu ad-hoc.
- Job je idempotentní: dvojí `submit_job` se stejným vstupem vrátí stejný
  `run_id` bez duplikace.
