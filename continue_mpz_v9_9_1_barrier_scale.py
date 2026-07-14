#!/usr/bin/env python3
"""v9.9.1 compatibility wrapper for barrier continuation metadata.

The v9.9 objective detail table is produced by merging the design-target table
with the calculated temperature response.  The merge already carries the
``target_class`` column, while the output assembly attempted to insert that
column a second time.  pandas 3 raises ``ValueError`` for this duplicate insert.

This wrapper makes metadata insertion idempotent for the isolated continuation
process: an existing metadata column is overwritten, while a new column uses
the normal pandas insertion behavior.  The original pandas method is restored
when the continuation exits.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator

import pandas as pd

import continue_mpz_v9_9_barrier_scale as implementation


@contextmanager
def idempotent_dataframe_insert() -> Iterator[None]:
    """Allow output metadata to overwrite an existing DataFrame column."""
    original = pd.DataFrame.insert

    def safe_insert(self, loc, column, value, *args, **kwargs):
        allow_duplicates = kwargs.get("allow_duplicates", False)
        if column in self.columns and allow_duplicates is not True:
            self[column] = value
            return None
        return original(self, loc, column, value, *args, **kwargs)

    pd.DataFrame.insert = safe_insert
    try:
        yield
    finally:
        pd.DataFrame.insert = original


def main() -> None:
    with idempotent_dataframe_insert():
        implementation.main()


if __name__ == "__main__":
    main()
