#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest

import numpy
import scipy

from arrhenius_fracture.mixed_mode_first_passage_v8 import MODEL_ID


def main():
    print("python:", sys.executable)
    print("numpy:", numpy.__version__)
    print("scipy:", scipy.__version__)
    print("model:", MODEL_ID)
    expected = "FEM_CZM_mixed_mode_first_passage_v8_exact_backend_full_circle"
    if MODEL_ID != expected:
        raise SystemExit(f"unexpected v8 model ID: {MODEL_ID}")
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    suite.addTests(loader.loadTestsFromName("tests.test_mixed_mode_first_passage_v8"))
    suite.addTests(loader.loadTestsFromName("tests.test_mixed_mode_campaign_v8"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    print("MIXED_MODE_V8 verification OK")


if __name__ == "__main__":
    main()
