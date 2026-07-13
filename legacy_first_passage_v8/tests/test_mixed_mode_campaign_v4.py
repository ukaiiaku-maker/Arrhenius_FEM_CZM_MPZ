import unittest
from pathlib import Path
class CommandTests(unittest.TestCase):
 def test_campaign_flags(self):
  s=Path('run_mixed_mode_fem_czm_v4_campaign.py').read_text();
  for x in ('--crystal-aniso','--crystal-compete','--target-traction-phase-deg','--traction-probe-radius-m'):self.assertIn(x,s)
 def test_no_isotropic_partition_for_kinetics(self):
  s=Path('arrhenius_fracture/mixed_mode_first_passage_v4.py').read_text();self.assertNotIn('KJ * math.cos',s);self.assertIn('sigma_cleave_drive_Pa',s)
if __name__=='__main__':unittest.main()
