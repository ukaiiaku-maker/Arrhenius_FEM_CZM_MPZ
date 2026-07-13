#!/usr/bin/env python3
from __future__ import annotations
import subprocess,sys,unittest

def main():
    import numpy,scipy
    from arrhenius_fracture.mixed_mode_first_passage_v3_3 import MODEL_ID
    expected='FEM_CZM_mixed_mode_first_passage_v3_3_J_consistent_circular_phase'
    print('python:',sys.executable);print('numpy:',numpy.__version__);print('scipy:',scipy.__version__)
    print('mechanics model:',MODEL_ID);print('campaign runner: mixed_mode_fem_czm_v3_4_campaign_probe_fix')
    if MODEL_ID!=expected:raise SystemExit(f'wrong v3.3 mechanics model: {MODEL_ID}')
    # Import probe must succeed without invoking the required mixed-mode CLI.
    cp=subprocess.run([sys.executable,'-c','from arrhenius_fracture.mixed_mode_first_passage_v3_3 import MODEL_ID;print(MODEL_ID)'],text=True,capture_output=True)
    if cp.returncode or expected not in cp.stdout:raise SystemExit(cp.stderr or cp.stdout)
    suite=unittest.TestSuite()
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName('tests.test_mixed_mode_first_passage_v3_3'))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromName('tests.test_mixed_mode_campaign_v3_4'))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
    print('MIXED_MODE_V3_4 campaign verification OK')
if __name__=='__main__':main()
