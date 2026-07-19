from pathlib import Path


def test_legacy_variational_fracture_modules_are_absent():
    package = Path(__file__).resolve().parents[1] / "arrhenius_fracture"
    forbidden = {
        "a" + "t1.py",
        "a" + "t2_overlay.py",
        "phase_" + "field.py",
        "main.py",
        "fatigue_" + "pf2d.py",
        "sn_" + "pf2d.py",
        "sn_" + "pf2d_fullplastic.py",
    }
    assert forbidden.isdisjoint({path.name for path in package.iterdir()})


def test_duplicated_legacy_source_tree_is_absent():
    root = Path(__file__).resolve().parents[1]
    assert not (root / ("legacy_" + "first_passage_v8")).exists()
