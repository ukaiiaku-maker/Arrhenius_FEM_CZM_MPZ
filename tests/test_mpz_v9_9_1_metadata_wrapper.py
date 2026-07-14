import pandas as pd

from continue_mpz_v9_9_1_barrier_scale import idempotent_dataframe_insert


def test_duplicate_metadata_insert_overwrites_existing_column():
    frame = pd.DataFrame(
        {
            "T_K": [300.0, 700.0],
            "target_class": ["old", "old"],
        }
    )
    original = pd.DataFrame.insert
    with idempotent_dataframe_insert():
        frame.insert(0, "continuation_candidate_id", "candidate")
        frame.insert(1, "target_class", "DBTT")
        frame.insert(2, "barrier_scale", 0.6)
    assert pd.DataFrame.insert is original
    assert frame.target_class.tolist() == ["DBTT", "DBTT"]
    assert frame.continuation_candidate_id.tolist() == ["candidate", "candidate"]
    assert frame.barrier_scale.tolist() == [0.6, 0.6]
    assert list(frame.columns).count("target_class") == 1
