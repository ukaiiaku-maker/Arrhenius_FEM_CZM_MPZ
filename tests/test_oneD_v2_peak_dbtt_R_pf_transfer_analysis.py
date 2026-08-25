import pandas as pd

from scripts.analyze_oneD_v2_peak_dbtt_R_pf_transfer import _target_prefix, cases


def test_transfer_plan_has_six_new_and_six_archived_cases():
    planned = cases()
    assert len(planned) == 12
    assert sum(bool(row["new_PF_run"]) for row in planned) == 6
    assert {row["material_class"] for row in planned} == {"Peak", "DBTT"}


def test_target_prefix_stops_at_first_target_state():
    frame = pd.DataFrame({"crack_extension_m": [0.0, 90e-6, 101e-6, 110e-6]})
    result = _target_prefix(frame)
    assert list(result.crack_extension_m) == [0.0, 90e-6, 101e-6]
