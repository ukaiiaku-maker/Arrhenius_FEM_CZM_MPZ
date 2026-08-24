from __future__ import annotations

from dataclasses import dataclass
from .base import MechanicsProvider, MechanicsState


@dataclass
class ReplayProvider(MechanicsProvider):
    pre: tuple[MechanicsState,...]
    post: tuple[MechanicsState,...]
    reload: tuple[MechanicsState,...]
    archived_geometry_fingerprint: str

    def evaluate_pre_event(self,index): return self.pre[index]
    def evaluate_post_event(self,index): return self.post[index]
    def evaluate_reload(self,index): return self.reload[index]

    def assert_geometry(self,fingerprint:str)->None:
        if fingerprint != self.archived_geometry_fingerprint: raise ValueError("replay geometry mutation")
