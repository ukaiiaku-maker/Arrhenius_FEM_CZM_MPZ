#!/usr/bin/env python3
from __future__ import annotations
import importlib,sys,unittest

def main():
    import numpy,scipy
    import arrhenius_fracture
    m=importlib.import_module('arrhenius_fracture.mixed_mode_first_passage_v2')
    print('python:',sys.executable);print('numpy:',numpy.__version__);print('scipy:',scipy.__version__);print('model:',m.MODEL_ID)
    suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromName('tests.test_mixed_mode_first_passage_v2'), unittest.defaultTestLoader.loadTestsFromName('tests.test_mixed_mode_control_v2')])
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
    print('MIXED_MODE_V2 verification OK')
if __name__=='__main__':main()
