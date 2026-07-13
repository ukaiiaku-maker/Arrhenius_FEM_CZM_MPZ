#!/usr/bin/env python3
from __future__ import annotations
import sys, unittest

def main():
    import numpy, scipy, arrhenius_fracture
    from arrhenius_fracture.mixed_mode_first_passage_v3_2 import MODEL_ID
    print('python:',sys.executable);print('numpy:',numpy.__version__);print('scipy:',scipy.__version__);print('model:',MODEL_ID)
    expected='FEM_CZM_mixed_mode_first_passage_v3_2_J_consistent_signed_basis'
    if MODEL_ID!=expected:raise SystemExit(f'wrong v3.2 model: {MODEL_ID}')
    suite=unittest.defaultTestLoader.loadTestsFromName('tests.test_mixed_mode_first_passage_v3_2')
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
    print('MIXED_MODE_V3_2 verification OK')
if __name__=='__main__':main()
