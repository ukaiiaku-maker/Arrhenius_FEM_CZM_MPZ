"""Guard against silently breaking the AST/text-patch anchors that ~15
other modules rely on to build specialized variants of
`sharp_front.run_2d` (see the "sharp_front.py architecture" notes in
CLAUDE_PROGRESS.md). Those modules locate their insertion points via
exact literal-substring matching against `inspect.getsource(run_2d)`, so
any edit to `run_2d` that changes whitespace/indentation of a matched
region -- even without changing its logic -- breaks them silently until
something downstream tries to exercise that specific patch chain.

This was discovered the hard way this session: an earlier edit wrapping
the mechanics/dU-selection block in a new `if pf_target_active: ... else:
<reindented original>` structure re-indented text that
mode_i_first_passage_v10_0_5_5_stochastic_vhcf.py's `_CACHE_INIT_ANCHOR`
and `_CACHE_MECHANICS_ANCHOR` depend on matching verbatim, breaking 7
tests that only surfaced when running the FULL suite (not the narrower
focused PF-controller suite this project has mostly been using this
session). This test makes that failure mode visible immediately, on its
own, without needing to run the full suite or notice a downstream test
several import-hops away.

Anchors that are pre-existing failures (confirmed broken already at the
workspace's starting commit ad06e4c, unrelated to any PF-controller work)
are explicitly excluded below with a comment explaining why -- do not add
new exclusions here without the same baseline-commit verification.
"""
from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest

from arrhenius_fracture import sharp_front

PACKAGE_DIR = Path(__file__).resolve().parents[1] / "arrhenius_fracture"

# (module filename, variable name) pairs confirmed broken already at the
# workspace's starting commit (ad06e4c), before any PF-controller work --
# unrelated pre-existing technical debt, not a regression to guard here.
KNOWN_PRE_EXISTING_BROKEN_ANCHORS = {
    ("mode_i_first_passage_v10_0_5_3_fatigue.py", "accepted_anchor"),
    ("mode_i_first_passage_v10_0_5_9_production_j_probe.py", "_AUDIT_INSERT_ANCHOR"),
}


def _discover_anchor_literals() -> list[tuple[str, str, str]]:
    """Find every `<name ending in _anchor/_Anchor> = "..."` top-level or
    function-local string-literal assignment across arrhenius_fracture/.
    """
    found = []
    for py_file in sorted(PACKAGE_DIR.glob("*.py")):
        text = py_file.read_text()
        if "anchor" not in text.lower():
            continue
        try:
            tree = ast.parse(text, filename=str(py_file))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if not (isinstance(target, ast.Name) and "anchor" in target.id.lower()):
                    continue
                try:
                    value = ast.literal_eval(node.value)
                except Exception:
                    continue
                if isinstance(value, str) and value.strip():
                    found.append((py_file.name, target.id, value))
    return found


ANCHOR_LITERALS = _discover_anchor_literals()


def test_at_least_the_known_anchors_were_discovered():
    # Sanity check on the discovery mechanism itself: if this drops to 0,
    # the AST-based scan broke, not the anchors -- don't let that pass
    # silently as "everything's fine".
    assert len(ANCHOR_LITERALS) >= 15


@pytest.mark.parametrize(
    "fname,varname,anchor",
    ANCHOR_LITERALS,
    ids=[f"{f}:{v}" for f, v, _ in ANCHOR_LITERALS],
)
def test_anchor_appears_exactly_once_in_run_2d_source(fname, varname, anchor):
    if (fname, varname) in KNOWN_PRE_EXISTING_BROKEN_ANCHORS:
        pytest.skip(f"{fname}:{varname} is a known pre-existing break, unrelated to this workspace")
    source = inspect.getsource(sharp_front.run_2d)
    count = source.count(anchor)
    assert count == 1, (
        f"{fname}:{varname} expected exactly one match in run_2d source, "
        f"found {count}. An edit to run_2d changed text this anchor "
        "depends on matching verbatim (see this file's module docstring)."
    )
