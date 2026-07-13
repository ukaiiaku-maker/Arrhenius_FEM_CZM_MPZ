#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest

import numpy
import scipy

from arrhenius_fracture.mixed_mode_first_passage_v8 import MODEL_ID

EXPECTED_MODEL_ID = (
    "FEM_CZM_mixed_mode_first_passage_v8_exact_backend_full_circle_boolean_safe"
)


def main() -> None:
    print("python:", sys.executable)
    print("numpy:", numpy.__version__)
    print("scipy:", scipy.__version__)
    print("model:", MODEL_ID)

    if MODEL_ID != EXPECTED_MODEL_ID:
        raise SystemExit(
            "unexpected v8 model ID:\n"
            f"  installed: {MODEL_ID}\n"
            f"  expected:  {EXPECTED_MODEL_ID}"
        )

    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    for module_name in (
        "tests.test_mixed_mode_first_passage_v8",
        "tests.test_mixed_mode_campaign_v8",
        "tests.test_mixed_mode_probe_boolean_v8",
    ):
        suite.addTests(loader.loadTestsFromName(module_name))

    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)

    print("MIXED_MODE_V8_1 verification OK")


if __name__ == "__main__":
    main()
