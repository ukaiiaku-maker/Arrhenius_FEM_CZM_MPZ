#!/usr/bin/env python3
"""Read-only inventory of the PF final v10.2.30 four-class campaign.

Walks
`PF-fracture-fatigue_v10_2_21_persistent_sites_top1/runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1`
across every rate directory, material-class option, and temperature/seed
leaf case, extracting completion status, first-passage/R-curve summary
statistics, event counts, kernel reference, and content hashes.

This script never writes into the PF repository -- it only reads. Output
goes to PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json and
PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md in this workspace's root, plus a
flat pf_case_inventory.csv under the scripts' own output directory.

Per CLAUDE.md: this campaign is the principal PF physical-correspondence
reference bank going forward (superseding the earlier, now-unavailable
v10.4.1 selective_reuse Peak/1000K case as the primary reference, though
that older case's provenance record remains valid for what it audited at
the time).
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

CAMPAIGN_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
    "/runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1"
)
RATE_DIRS = ["rate0p01x", "rate1x", "rate100x"]

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = WORKSPACE_ROOT / "PF_FINAL_CAMPAIGN_REFERENCE_MANIFEST.json"
SUMMARY_PATH = WORKSPACE_ROOT / "PF_FINAL_CAMPAIGN_REFERENCE_SUMMARY.md"
INVENTORY_CSV_PATH = WORKSPACE_ROOT / "runs" / "pf_final_campaign_v10230" / "pf_case_inventory.csv"

TEMP_SEED_RE = re.compile(r"^T(\d+)K_th(\d+)_seed(\d+)$")


def _sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _safe_json(path: Path) -> tuple[Any, str | None]:
    if not path.is_file():
        return None, "missing"
    try:
        with open(path) as f:
            return json.load(f), None
    except Exception as exc:  # noqa: BLE001 -- deliberately broad for a read-only audit
        return None, str(exc)


def _extract_case(case_dir: Path, rate: str) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "rate": rate,
        "case_dir": str(case_dir),
        "option_dir": case_dir.parent.name,
        "temp_seed_dir": case_dir.name,
        "missing_or_malformed": [],
    }

    m = TEMP_SEED_RE.match(case_dir.name)
    if m:
        entry["temperature_K"] = int(m.group(1))
        entry["theta_deg"] = int(m.group(2))
        entry["seed"] = int(m.group(3))
    else:
        entry["missing_or_malformed"].append("unparsed_temp_seed_dirname")

    complete_path = case_dir / "COMPLETE"
    exit_code_path = case_dir / "exit_code.txt"
    entry["complete_marker_present"] = complete_path.is_file()
    entry["complete_marker_text"] = (
        complete_path.read_text().strip() if complete_path.is_file() else None
    )
    entry["exit_code"] = (
        exit_code_path.read_text().strip() if exit_code_path.is_file() else None
    )
    if not complete_path.is_file():
        entry["missing_or_malformed"].append("COMPLETE_marker_missing")
    if entry["exit_code"] not in ("0", None):
        entry["missing_or_malformed"].append(f"nonzero_exit_code:{entry['exit_code']}")

    status, err = _safe_json(case_dir / "stage3_case_status.json")
    if err:
        entry["missing_or_malformed"].append(f"stage3_case_status.json:{err}")
    if status:
        entry["status_schema"] = status.get("schema")
        entry["status_complete"] = status.get("complete")
        entry["status_string"] = status.get("status")
        entry["first_passage_recorded"] = status.get("first_passage_recorded")
        entry["projected_extension_um"] = status.get("projected_extension_um")
        entry["target_extension_um"] = status.get("target_extension_um")
        entry["target_tolerance_um"] = status.get("target_tolerance_um")
        entry["Kc_first_MPa_sqrt_m"] = status.get("Kc_first_MPa_sqrt_m")
        summ = status.get("summary", {}) or {}
        entry["N_em_init"] = summ.get("N_em_init")
        entry["N_em_final"] = summ.get("N_em_final")
        entry["W_emit_J_per_m"] = summ.get("W_emit_J_per_m")
        entry["sigma_back_init_GPa"] = summ.get("sigma_back_init_GPa")
        entry["r_eff_over_r0_init"] = summ.get("r_eff_over_r0_init")
        entry["n_advances"] = summ.get("n_advances")
        entry["n_advances_primary"] = summ.get("n_advances_primary")
        entry["n_advances_branch"] = summ.get("n_advances_branch")
        entry["n_geometry_events"] = summ.get("n_geometry_events")
        entry["branched"] = summ.get("branched")
        entry["mode"] = summ.get("mode")
        entry["a_final_mm"] = summ.get("a_final_mm")
        entry["hbar_tip_m"] = summ.get("hbar_tip_m")
        entry["n_nodes"] = summ.get("n_nodes")
        shelf = summ.get("shelf", {}) or {}
        entry["material_class_label"] = shelf.get("material_class")
        entry["t_total_s"] = shelf.get("t_total")
        entry["clock_completable"] = shelf.get("clock_completable")
    else:
        entry["missing_or_malformed"].append("stage3_case_status.json_unreadable")

    run_args, err = _safe_json(case_dir / "run_args.json")
    if err:
        entry["missing_or_malformed"].append(f"run_args.json:{err}")
    if run_args:
        entry["dU_m"] = run_args.get("dU")
        entry["dt_s"] = run_args.get("dt")
        entry["n_stagger"] = run_args.get("n_stagger")
        entry["bulk_kinetics_model"] = run_args.get("bulk_kinetics_model")
        entry["bulk_kinetics_model_detail"] = run_args.get("bulk_kinetics_model_detail")
        entry["bulk_mult_frac"] = run_args.get("bulk_mult_frac")
        entry["crack_backend"] = run_args.get("crack_backend")
        entry["target_crack_extension_um_run_args"] = run_args.get("target_crack_extension_um")
        entry["Kdot"] = run_args.get("Kdot")

    entry["kernel_family_path"] = None
    entry["kernel_family_sha256_now"] = None
    entry["kernel_family_sha256_provenance_note"] = None
    cmd_path = case_dir / "command.sh"
    if cmd_path.is_file():
        text = cmd_path.read_text()
        km = re.search(r"--signed-kernel-family\s+(\S+)", text)
        if km:
            kpath = Path(km.group(1))
            entry["kernel_family_path"] = str(kpath)
            if kpath.is_file():
                entry["kernel_family_sha256_now"] = _sha256_of(kpath)
                entry["kernel_family_sha256_provenance_note"] = (
                    "Re-hashed from the LIVE path at inventory time. This "
                    "campaign's case files are dated 2026-07-29; this exact "
                    "kernel_cache path was independently observed to change "
                    "during a later, unrelated session (2026-08-03/04) -- "
                    "this hash is NOT guaranteed to equal what the campaign "
                    "actually consumed at runtime. No runtime-captured "
                    "kernel content hash exists in this case's own audit "
                    "JSON to cross-check against; only shape/policy metadata "
                    "(v10_2_30_hazard_energy_gate_audit.json's "
                    "signed_burgers_shared_physics.kernel block) is recorded "
                    "there, not a byte fingerprint."
                )
            else:
                entry["missing_or_malformed"].append("kernel_family_path_missing_now")
        else:
            entry["missing_or_malformed"].append("signed_kernel_family_flag_not_found_in_command_sh")
    else:
        entry["missing_or_malformed"].append("command.sh_missing")

    steps_candidates = sorted(case_dir.glob("steps_*K.csv"))
    if len(steps_candidates) == 1:
        steps_path = steps_candidates[0]
        entry["steps_csv_path"] = str(steps_path)
        entry["steps_csv_sha256"] = _sha256_of(steps_path)
        entry["steps_csv_size_bytes"] = steps_path.stat().st_size
        try:
            with open(steps_path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            entry["steps_csv_rows"] = len(rows)
            if rows:
                entry["steps_csv_columns"] = list(rows[0].keys())
                fp_idx = None
                for i, r in enumerate(rows):
                    try:
                        if float(r.get("crack_extension_m", 0) or 0) > 0:
                            fp_idx = i
                            break
                    except (TypeError, ValueError):
                        continue
                if fp_idx is not None:
                    fp_row = rows[fp_idx]
                    entry["first_passage_row_index"] = fp_idx
                    entry["first_passage_step"] = fp_row.get("step")
                    entry["first_passage_J_J_per_m2"] = fp_row.get("J_effective_direct_J_per_m2")
                    kj_raw = fp_row.get("KJ_Pa_sqrtm")
                    entry["first_passage_KJ_MPa_sqrt_m"] = (
                        float(kj_raw) / 1.0e6 if kj_raw else None
                    )
                    try:
                        entry["first_passage_time_s"] = sum(
                            float(r.get("dt_cur_s", 0) or 0) for r in rows[: fp_idx + 1]
                        )
                    except Exception:  # noqa: BLE001
                        entry["missing_or_malformed"].append("first_passage_time_computation_failed")
                else:
                    entry["missing_or_malformed"].append("no_first_passage_row_found")

                last_row = rows[-1]
                entry["final_crack_extension_m"] = last_row.get("crack_extension_m")
                kj_final = last_row.get("KJ_Pa_sqrtm")
                entry["final_KJ_MPa_sqrt_m"] = float(kj_final) / 1.0e6 if kj_final else None
                entry["final_J_J_per_m2"] = last_row.get("J_effective_direct_J_per_m2")
                entry["final_B"] = last_row.get("B")
                entry["final_N_em"] = last_row.get("N_em")
                entry["final_W_bulk_plastic_cumulative_J_per_m"] = last_row.get(
                    "W_bulk_plastic_cumulative_J_per_m"
                )
                entry["final_mpz_mobile_count"] = last_row.get("mpz_mobile_count")
                entry["final_mpz_retained_count"] = last_row.get("mpz_retained_count")
                entry["final_step"] = last_row.get("step")

                # Event-length statistics: successive first-passage-like
                # jumps in crack_extension_m across accepted rows (n_fire>0
                # where available, else any positive delta).
                exts = []
                for r in rows:
                    try:
                        exts.append(float(r.get("crack_extension_m", 0) or 0))
                    except (TypeError, ValueError):
                        exts.append(None)
                deltas = []
                prev = None
                for v in exts:
                    if v is None:
                        continue
                    if prev is not None and v > prev + 1e-15:
                        deltas.append(v - prev)
                    prev = v
                if deltas:
                    entry["event_length_count"] = len(deltas)
                    entry["event_length_mean_m"] = sum(deltas) / len(deltas)
                    entry["event_length_min_m"] = min(deltas)
                    entry["event_length_max_m"] = max(deltas)
        except Exception as exc:  # noqa: BLE001
            entry["missing_or_malformed"].append(f"steps_csv_parse_error:{exc}")
    else:
        entry["missing_or_malformed"].append(
            f"steps_csv_ambiguous_or_missing(found={len(steps_candidates)})"
        )

    crack_path_candidates = sorted(case_dir.glob("crack_path_*K.csv"))
    entry["crack_path_csv_path"] = str(crack_path_candidates[0]) if crack_path_candidates else None
    if not crack_path_candidates:
        entry["missing_or_malformed"].append("crack_path_csv_missing")

    avalanche_path = case_dir / "stochastic_avalanche_geometry_events.json"
    if avalanche_path.is_file():
        entry["stochastic_avalanche_geometry_events_path"] = str(avalanche_path)
        entry["stochastic_avalanche_geometry_events_sha256"] = _sha256_of(avalanche_path)
    else:
        entry["missing_or_malformed"].append("stochastic_avalanche_geometry_events_missing")

    for name, key in [
        ("v10_2_27_case_contract.json", "case_contract_path"),
        ("v10_2_22_parameter_selection.json", "parameter_selection_path"),
        ("v10_2_30_hazard_energy_gate_audit.json", "hazard_energy_gate_audit_path"),
        ("v10_2_27_energy_ledger_output_audit.json", "energy_ledger_output_audit_path"),
    ]:
        p = case_dir / name
        entry[key] = str(p) if p.is_file() else None
        if not p.is_file():
            entry["missing_or_malformed"].append(f"{name}_missing")

    return entry


def _material_class_from_option(option_dir: str) -> str:
    lower = option_dir.lower()
    if "peak" in lower:
        return "Peak"
    if "dbtt" in lower:
        return "DBTT"
    if "weakt" in lower or "weak_t" in lower or "weak-t" in lower:
        return "weak-T"
    if "ceramic" in lower:
        return "ceramic"
    return "unknown"


def main() -> None:
    if not CAMPAIGN_ROOT.is_dir():
        raise SystemExit(f"campaign root not found: {CAMPAIGN_ROOT}")

    cases: list[dict[str, Any]] = []
    for rate in RATE_DIRS:
        rate_dir = CAMPAIGN_ROOT / rate
        if not rate_dir.is_dir():
            cases.append({"rate": rate, "missing_or_malformed": ["rate_dir_missing"]})
            continue
        option_dirs = sorted(
            p for p in rate_dir.iterdir() if p.is_dir() and p.name.startswith("v913_")
        )
        for option_dir in option_dirs:
            temp_seed_dirs = sorted(
                p for p in option_dir.iterdir() if p.is_dir() and TEMP_SEED_RE.match(p.name)
            )
            for case_dir in temp_seed_dirs:
                entry = _extract_case(case_dir, rate)
                entry["material_class"] = _material_class_from_option(option_dir.name)
                cases.append(entry)

    manifest = {
        "schema": "pf_final_campaign_v10230_reference_manifest_v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "campaign_root": str(CAMPAIGN_ROOT),
        "rate_dirs": RATE_DIRS,
        "n_cases": len(cases),
        "n_cases_missing_or_malformed": sum(1 for c in cases if c.get("missing_or_malformed")),
        "cases": cases,
    }

    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n")

    INVENTORY_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    flat_keys: list[str] = []
    for c in cases:
        for k in c.keys():
            if k not in flat_keys and k != "missing_or_malformed":
                flat_keys.append(k)
    flat_keys.append("missing_or_malformed")
    with open(INVENTORY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=flat_keys, extrasaction="ignore")
        writer.writeheader()
        for c in cases:
            row = dict(c)
            row["missing_or_malformed"] = ";".join(c.get("missing_or_malformed", []))
            if "steps_csv_columns" in row and isinstance(row["steps_csv_columns"], list):
                row["steps_csv_columns"] = "|".join(row["steps_csv_columns"])
            writer.writerow(row)

    print(f"Wrote manifest: {MANIFEST_PATH} ({len(cases)} cases)")
    print(f"Wrote flat inventory CSV: {INVENTORY_CSV_PATH}")
    n_bad = manifest["n_cases_missing_or_malformed"]
    print(f"Cases with missing/malformed fields: {n_bad} / {len(cases)}")


if __name__ == "__main__":
    main()
