#!/usr/bin/env python3
import sys,unittest
from arrhenius_fracture.mixed_mode_first_passage_v4 import MODEL_ID
print('model:',MODEL_ID)
if MODEL_ID!='FEM_CZM_mixed_mode_first_passage_v4_anisotropic_traction_controlled':raise SystemExit('wrong model')
s=unittest.TestSuite();l=unittest.defaultTestLoader
for n in ('tests.test_mixed_mode_first_passage_v4','tests.test_mixed_mode_campaign_v4'):s.addTests(l.loadTestsFromName(n))
r=unittest.TextTestRunner(verbosity=2).run(s)
if not r.wasSuccessful():raise SystemExit(1)
print('MIXED_MODE_V4 verification OK')
