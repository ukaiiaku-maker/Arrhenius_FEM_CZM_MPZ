#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest

import numpy
import scipy

from arrhenius_fracture.mixed_mode_first_passage_v5 import MODEL_ID


def main():
    print("python:", sys.executable)
    print("numpy:", numpy.__version__)
    print("scipy:", scipy.__version__)
    print("model:", MODEL_ID)
    if MODEL_ID != "FEM_CZM_mixed_mode_first_passage_v5_anisotropic_calibrated_tip":
        raise SystemExit("unexpected v5 model ID")
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    suite.addTests(loader.loadTestsFromName("tests.test_mixed_mode_first_passage_v5"))
    suite.addTests(loader.loadTestsFromName("tests.test_mixed_mode_campaign_v5"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    print("MIXED_MODE_V5 verification OK")


if __name__ == "__main__":
    main()
