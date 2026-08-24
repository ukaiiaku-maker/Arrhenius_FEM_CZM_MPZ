from __future__ import annotations
import numpy as np


class GeneratedSequence:
    def __init__(self,seed:int): self._rng=np.random.default_rng(seed)
    def thresholds(self,n:int)->tuple[float,...]: return tuple(float(x) for x in self._rng.exponential(size=n))


class ReplaySequence:
    def __init__(self,values): self.values=tuple(float(x) for x in values)
    def thresholds(self,n:int)->tuple[float,...]:
        if n>len(self.values): raise ValueError("archived threshold sequence exhausted")
        return self.values[:n]
