import unittest
import run_mixed_mode_fem_czm_v3_4_campaign as runner

class CampaignRunnerTests(unittest.TestCase):
    def test_probe_does_not_invoke_required_cli(self):
        cmd = runner.driver_probe_command('/tmp/python')
        self.assertEqual(cmd[0], '/tmp/python')
        self.assertIn('-c', cmd)
        self.assertNotIn('--help', cmd)
        self.assertNotIn('-m', cmd)

    def test_case_arguments_include_target_phase(self):
        args = runner.required_mixed_args(-87.1, -30.0, 1.0)
        i = args.index('--target-mode-phase-deg')
        self.assertEqual(float(args[i+1]), -30.0)
        self.assertIn('--mixity-loading-angle-deg', args)
        self.assertIn('--deterministic-threshold', args)

if __name__ == '__main__':
    unittest.main()
