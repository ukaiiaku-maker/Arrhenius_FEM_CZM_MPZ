from arrhenius_fracture.kinetic_campaign_czm_v1003 import (
    CampaignAwareV1003TipEngineMixin,
    STATE_MODEL,
)


class DummyCampaign(CampaignAwareV1003TipEngineMixin):
    state_model = STATE_MODEL

    def __init__(self):
        self.predict_args = None
        self.step_args = None
        self._mm_prev_Kcleave = None

    def _mm_drives(self, K):
        return 2.0 * K, 3.0 * K, {
            "KJ": K,
            "fc": 2.0,
            "fe": 3.0,
            "candidate_angle_deg": 0.0,
        }

    def predict_clock_increment_drives(self, Kc, Ke, T, dt):
        self.predict_args = (Kc, Ke, T, dt)
        return 7.0

    def step_drives(self, Kc, Ke, T, dt, metadata=None):
        self.step_args = (Kc, Ke, T, dt, metadata)
        return {"fired": False, "n_fire": 0}


def test_campaign_prediction_uses_separated_drive_api():
    eng = DummyCampaign()
    assert eng.predict_clock_increment(4.0, 700.0, 0.25) == 7.0
    assert eng.predict_args == (8.0, 12.0, 700.0, 0.25)


def test_campaign_step_uses_separated_drive_api():
    eng = DummyCampaign()
    out = eng.step(4.0, 700.0, 0.25)
    assert out == {"fired": False, "n_fire": 0}
    Kc, Ke, T, dt, metadata = eng.step_args
    assert (Kc, Ke, T, dt) == (8.0, 12.0, 700.0, 0.25)
    assert metadata["anisotropic_KJ_Pa_sqrt_m"] == 4.0
    assert metadata["anisotropic_Kcleave_Pa_sqrt_m"] == 8.0
    assert metadata["anisotropic_Kemit_Pa_sqrt_m"] == 12.0
