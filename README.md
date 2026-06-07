# Part XXI — The Irreducible Three-Body Congruence Correlation and the Annihilation Lattice

*Volume II of the Arithmetic Geodynamics programme on the 6N skeleton.*

Part XX recast the two-centre correlation as a circular convolution of castellated
indicators on the torus `T¹`. **This paper extends it to three twin centres and
asks whether the three-body correlation factorizes into pairwise terms. It does
not.**

For twins at N, N+j₁, N+j₂ the joint singular series, normalized by the three
marginals, is `R₃(j₁,j₂)`. The **irreducible three-body factor** is

> U(j₁,j₂) = R₃(j₁,j₂) / ( R₂(j₁) · R₂(j₂) · R₂(j₂−j₁) ),  with R₂ the Part XX pairwise correlation.

**U = 1 would mean the three-body correlation is just its three pairwise pieces.**

### Findings

- **Irreducible everywhere.** Across all 870 ordered lag pairs in the window
  1 ≤ j₁,j₂ ≤ 30, `U = 1` occurs **nowhere**. A genuine three-body congruence term
  is generic; finite values range U ≈ 0.37 … 2.24, sharpest at triple-hole
  alignments (e.g. (5,10): 0,5,10 all ≡ 0 mod 5).
- **Annihilation is tiling, not counting.** R = 0 **iff** the shifted dead-residue
  sets cover ℤ/q for some prime q. The naive heuristic "m_A+m_B ≥ q ⟹ collapse" is
  **false**: covering ℤ/q needs the holes to be *complementary*, not merely
  numerous.
- **No Twin–Triplet annihilation.** A triplet contains a twin, so at q=5
  dead_trip = {0,1,4} ⊇ {1,4} = dead_twin; the sets always overlap and never tile.
  Direct check: I₅(j) ∈ {1/5, 2/5}, never 0. The lowest order at which twin-based
  annihilation occurs is therefore **three**-body.
- **The annihilation lattice.** 216 forbidden lag pairs in the window. The defining
  case (j₁,j₂)=(1,2) — three consecutive twin centres — tiles ℤ/5:
  {1,4}∪{0,3}∪{2,4} = ℤ/5, so 6N−1,6N+1,6N+5,6N+7,6N+11,6N+13 always contain a
  multiple of 5. **This recovers the mod-5 inadmissibility of that six-member
  window from pure phase-space tiling.**

The horizontal annihilation here is distinct from the single-centre **killer
prime** (the vertical ω-collapse of Volume I, where a tuple's own member 6N+5 dies
when 5|N).

## Layout

```
.
├── paper/    Chen_6N_Paper21.{tex,pdf} + figure
├── figures/  fig_three_body.{pdf,png}
├── data/     U_grid.csv (j1,j2,U) · annihilation_lattice.csv (216 forbidden pairs)
├── code/
│   ├── three_body_correlation.py   # computes U(j1,j2), the lattice, and the figure
│   └── verify_tiling.py            # checks the tiling theorem & the corollaries
├── CITATION.cff · .zenodo.json · LICENSE (MIT)
```

## Reproducing

```bash
pip install numpy matplotlib
python code/verify_tiling.py            # all checks pass
python code/three_body_correlation.py   # regenerates U, the 216-point lattice, the figure
```

Expected: consecutive twins tile ℤ/5 (annihilation); Twin–Triplet min I₅ = 0.200
(no annihilation); 0 of 870 pairs factorize (U=1 nowhere).

## Scope

Closed-form throughout. At integer lags the phase-space integral equals the
discrete singular-series ratio exactly; no new prime data, no infinitude claim.
Continues Part XX (doi:10.5281/zenodo.20583994).

## License

MIT — see `LICENSE`.
