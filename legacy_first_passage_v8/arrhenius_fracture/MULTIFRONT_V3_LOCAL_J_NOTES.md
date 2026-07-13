# Multifront v3: local J contour and unresolved-branch handoff

This patch addresses the failure mode where branches are born in the multifront
inventory but remain visually and mechanically short.  The underlying issue is
scale separation: a newborn daughter branch is initially O(da_phys), while the
old J-domain was much larger than that.  Reading an independent J-integral on
such a short stub is not meaningful, so the branch starves numerically.

## Main changes

1. **Local J outer radius.**
   Added `--rJ-outer` as the actual desired outer radius of the local J-integral
   contour.  The code passes `ell = rJ_outer/8` to the existing domain-integral
   routine, whose default outer factor is 8.  The legacy `--rJ` is retained, but
   is interpreted as the old ell-like value, so `rJ_outer=8*rJ` for backward
   compatibility.  The default local outer radius is tied to `da_phys`, `L_pz`,
   and tip mesh size.

2. **Unresolved daughter handoff.**
   A newborn daughter front is marked `resolved=False`.  It is drawn as a sharp
   branch segment, but it does not immediately get its own independent J-integral.
   Instead, while it remains inside the parent-tip J/process-zone domain, it
   co-grows with the parent under the parent event budget.  Once its arclength
   and separation exceed `--branch-resolve-length` (default: local J outer radius
   or several `da_phys`), it is marked resolved and begins using its own local
   J-integral and renewal clock.

3. **Local budget conservation, not global starvation.**
   Parent and unresolved daughter split the parent event budget.  Unrelated
   resolved fronts have their own local renewal-clock budgets.  This preserves
   local bifurcation conservation without dividing one global increment among
   all fronts.

4. **Multi-tip remeshing.**
   `mesh.make_tri_mesh(..., tip_center=...)` now accepts multiple tip centers and
   unions radial fine patches around all active tips.  Remeshing passes all
   active front tips, not only the leading x-tip, so separated branches remain
   locally resolved.

5. **Diagnostics.**
   `fronts_<T>K.csv` now includes `resolved` and `branch_len_m` columns so branch
   handoff can be audited.

## Recommended production flags

Start with:

```bash
--tip-h-fine 0.35e-6 \
--rJ-outer 4e-6 \
--branch-resolve-length 4e-6
```

For cheaper pilots use `--tip-h-fine 0.6e-6 --rJ-outer 6e-6`.  If the branch is
still shorter than the contour, reduce `da_phys` and/or refine the mesh; do not
add velocity or advance caps.
