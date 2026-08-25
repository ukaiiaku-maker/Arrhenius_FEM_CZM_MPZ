from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np

from reduced_fracture_v2.native_state_factory import (
    ExactSourceDualLane,
    FEMCZMNativeStateFactory,
    PFNativeStateFactory,
    TypedNativeState,
    canonical_source_value,
    scientific_fingerprint,
)


class _MPZ:
    def __init__(self): self.values=np.asarray([1.0,2.0]); self.time_s=0.0
    def copy(self): return copy.deepcopy(self)


class _PF:
    def __init__(self):
        self.mpz=_MPZ(); self.B=0.0; self.a_adv=0.0; self.n_adv=0
        self.W_emit=0.0; self.t=0.0; self.K_prev=None; self._lambda_c_prev=None
        self.cb=SimpleNamespace(kind="source"); self.eb=SimpleNamespace(kind="source")
        self.f=SimpleNamespace(r0=1e-6); self.mpz_config=SimpleNamespace(n=2)
    def step(self,K,T,dt):
        self.B += dt; self.t += dt; self.mpz.time_s += dt
        fired=self.B >= 1.0
        if fired: self.B-=1.; self.a_adv+=5e-6; self.n_adv+=1
        return {"fired":fired,"n_fire":int(fired),"K":K,"T":T}


class _FEM(_PF):
    def __init__(self):
        super().__init__(); self.mpz_state=self.mpz
        self._hazard_rng=np.random.default_rng(91)
        self.hazard_threshold_action=.7; self.hazard_action_current=0.
    def _capture_state(self):
        return {"mpz_state":self.mpz_state.copy(),"B":self.B,"a_adv":self.a_adv,
                "n_adv":self.n_adv,"t":self.t,
                "hazard_rng_state":copy.deepcopy(self._hazard_rng.bit_generator.state)}
    def step_drives(self,Kc,Ke,T,dt,metadata=None): return self.step(Kc,T,dt)


def test_production_constructor_defines_lossless_independent_initial_state():
    for factory in (PFNativeStateFactory(_PF),FEMCZMNativeStateFactory(_FEM)):
        product=factory.instantiate()
        assert product.independent_lane_objects
        assert product.production_object is not product.v2_object
        assert product.production_state.source_owned_fields == product.v2_state.source_owned_fields
        assert product.future_state_imports == 0


def test_rng_and_threshold_are_part_of_fem_lossless_export():
    state=TypedNativeState.export("FEMCZM",_FEM())
    text=str(state.source_owned_fields)
    assert "hazard_rng_state" in text
    assert state.fingerprint == scientific_fingerprint(state.source_owned_fields)


def test_dual_lane_does_not_import_future_source_state():
    lanes=ExactSourceDualLane(PFNativeStateFactory(_PF).instantiate())
    result=lanes.advance(10.,10.,1000.,.25)
    assert result.exact_state_match
    assert result.field_status == "EXACT_SOURCE_MATCH"
    assert lanes.future_state_imports == 0
    lanes.source.mpz.values[0] = 99.
    assert lanes.v2.mpz.values[0] != 99.


def test_canonical_arrays_preserve_values_dtype_and_shape():
    array=np.asarray([[1,2],[3,4]],dtype=np.int16)
    value=canonical_source_value(array)
    assert value == {"__ndarray__":True,"dtype":"int16","shape":[2,2],"values":[[1,2],[3,4]]}


def test_native_state_fingerprint_is_deterministic():
    a=TypedNativeState.export("PF",_PF()); b=TypedNativeState.export("PF",_PF())
    assert a.fingerprint == b.fingerprint
