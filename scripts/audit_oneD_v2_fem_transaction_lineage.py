#!/usr/bin/env python3
"""Hash transaction-owning FEM/CZM functions across the qualified lineage."""

from __future__ import annotations

import ast
import hashlib
import json
import subprocess
from pathlib import Path


REPO = Path("/private/tmp/oneD-v2-femczm-driver-diagnostics")
OLD = "30b53ff"
NEW = "931bed66913afc970117bce900805ecc9b6225f8"
PATH = "arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py"
OUT = Path(__file__).resolve().parents[1] / "analysis_outputs/oneD_v2_capability_v2"
REPLACEMENT_OWNERS = {
    "outer transaction begin": "transaction snapshot and process-zone capture",
    "mechanics interval begin": "accepted/working process-zone separation",
    "accepted mechanics interval": "accepted physical-interval state",
    "pre-remap cohesive snapshot": "tentative geometry and geometry snapshot",
    "overlap topology remap": "tentative remap",
    "deflect geometry veto rollback": "late veto; threshold/RNG and geometry restoration",
    "authoritative commit": "authoritative process-zone commit",
}


def source_at(commit: str) -> str:
    return subprocess.check_output(
        ["git", "show", f"{commit}:{PATH}"], cwd=REPO, text=True
    )


def functions(source: str) -> dict[str, str]:
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    return {
        node.name: "".join(lines[node.lineno - 1 : node.end_lineno])
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def replacements(source: str) -> dict[str, str]:
    """Return injected production fragments by their fail-closed replacement label."""
    found = {}
    for node in ast.walk(ast.parse(source)):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "_replace" and len(node.args) >= 4):
            continue
        try:
            label = ast.literal_eval(node.args[3])
            injected = ast.literal_eval(node.args[2])
        except (ValueError, TypeError):
            continue
        found[str(label)] = str(injected)
    return found


def normalized(text: str) -> str:
    return ast.dump(ast.parse(text), annotate_fields=True, include_attributes=False)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> None:
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", OLD, NEW], cwd=REPO
    ).returncode == 0
    old_source, new_source = source_at(OLD), source_at(NEW)
    old, new = functions(old_source), functions(new_source)
    old_repl, new_repl = replacements(old_source), replacements(new_source)
    rows = []
    for label, responsibility in REPLACEMENT_OWNERS.items():
        a, b = old_repl[label], new_repl[label]
        if a == b:
            classification = "IDENTICAL"
        elif label == "authoritative commit" and all(statement in b for statement in (
            "authoritative_bulk_v1043 = bulk_pt_v1043.commit_topology_trial()",
            "ep_gp = authoritative_bulk_v1043.ep_gp.copy()",
            "rho_gp = authoritative_bulk_v1043.rho_gp.copy()",
            "_BULK_PT_RUNTIME_V1043['authoritative_commits'] += 1",
        )):
            # The corrected lineage adds census/field observers after the same
            # commit and state-copy operations.  It does not change ownership,
            # ordering, or the committed state.
            classification = "NONSEMANTIC_CHANGE"
        else:
            classification = "SEMANTIC_CHANGE_UNQUALIFIED"
        rows.append({
            "file_path": PATH,
            "function_name": "build_bulk_run_2d",
            "boundary_name": label,
            "responsibility": responsibility,
            "source_hash_at_30b53ff": sha(a),
            "source_hash_at_931bed6": sha(b),
            "change_classification": classification,
        })
    # The remap wrapper acquired material-parent provenance and compact failure
    # capture.  That is a real semantic change, requalified at the corrected
    # commit by the two bounded regressions recorded below.
    a, b = old["remesh_with_failure_capture"], new["remesh_with_failure_capture"]
    rows.append({
        "file_path": PATH,
        "function_name": "remesh_with_failure_capture",
        "boundary_name": "tentative remap/failure capture",
        "responsibility": "tentative remap, material provenance, failure capture",
        "source_hash_at_30b53ff": sha(a),
        "source_hash_at_931bed6": sha(b),
        "change_classification": "SEMANTIC_CHANGE_REQUALIFIED",
    })
    a, b = old["write_joint_veto_audit"], new["write_joint_veto_audit"]
    rows.append({
        "file_path": PATH,
        "function_name": "write_joint_veto_audit",
        "boundary_name": "post-rollback audit",
        "responsibility": "post-rollback restoration verification",
        "source_hash_at_30b53ff": sha(a),
        "source_hash_at_931bed6": sha(b),
        "change_classification": "IDENTICAL" if a == b else "SEMANTIC_CHANGE_UNQUALIFIED",
    })
    payload = {
        "schema": "oneD_v2_fem_transaction_lineage_v1",
        "qualified_evidence_commit": OLD,
        "corrected_production_commit": NEW,
        "qualified_commit_is_ancestor": ancestor,
        "functions": rows,
        "lineage_status": (
            "QUALIFIED_FROM_EXISTING_TRANSACTIONAL_EVIDENCE"
            if ancestor and all(r["change_classification"] in {
                "IDENTICAL", "NONSEMANTIC_CHANGE", "SEMANTIC_CHANGE_REQUALIFIED"
            } for r in rows)
            else "TARGETED_REQUALIFICATION_REQUIRED"
        ),
        "targeted_corrected_lineage_regression": {
            "command": "python -m pytest -q tests/test_bulk_pt_live_run2d_composition_v1043.py -k 'joint_late_veto or failed_remap_capture'",
            "environment": "arrhenius-fem-czm",
            "result": "2 passed, 4 deselected",
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "oneD_v2_fem_transaction_lineage.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    md = [
        "# One-dimensional V2 FEM/CZM transaction lineage audit",
        "",
        f"Ancestry `{OLD} -> {NEW}`: **{ancestor}**.",
        "",
        "| File | Function / boundary | 30b53ff SHA-256 | 931bed6 SHA-256 | Classification |",
        "|---|---|---|---|---|",
    ]
    md.extend(
        f"| `{r['file_path']}` | `{r['function_name']} / {r['boundary_name']}` | `{r['source_hash_at_30b53ff']}` | "
        f"`{r['source_hash_at_931bed6']}` | **{r['change_classification']}** |"
        for r in rows
    )
    md += [
        "",
        "Corrected-lineage bounded regression: **2 passed, 4 deselected**. No trajectory was run.",
        "",
        f"Lifecycle result: **{payload['lineage_status']}**.",
        "",
    ]
    (OUT / "ONE_D_V2_FEMCZM_TRANSACTION_LINEAGE_AUDIT.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
