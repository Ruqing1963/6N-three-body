#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifier for Part XXI: the tiling annihilation criterion.

Theorem (Annihilation = tiling): a joint constellation density vanishes iff the
shifted dead-residue sets cover Z/q for some prime q. This script checks:
  (1) three consecutive twin centres (lags 0,1,2) tile Z/5  => R3 = 0 (annihilation);
  (2) a Twin-Triplet two-body pair NEVER tiles Z/5  => I_5(j) >= 1/5 > 0 (no annihilation);
  (3) the three-body factor U(j1,j2) != 1 for every lag pair in a window
      (no factorization into pairwise terms).
Standard library only.
"""
import math

PRIMES = [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]
def dead(offsets, q):
    inv6 = pow(6, -1, q)
    return {(-a * inv6) % q for a in offsets}
TWIN = (-1, 1)
TRIP = (-1, 1, 5)

def nu_union(deadsets):
    U = set()
    for d in deadsets: U |= d
    return len(U)

def check_consecutive_twins():
    q = 5
    sets = [ {(r - j) % q for r in dead(TWIN, q)} for j in (0,1,2) ]
    cover = nu_union(sets)
    print("(1) three consecutive twins (lags 0,1,2) at q=5:")
    print(f"    shifted dead-sets {[sorted(s) for s in sets]}  union size {cover}/{q}")
    print(f"    tiles Z/5: {cover==q}  ->  R3 = 0 (annihilation)  [expected True]")
    return cover == q

def check_twin_triplet():
    q = 5
    dt = dead(TWIN, q)
    print("(2) Twin(N) + Triplet(N+j) at q=5, all shifts:")
    Imin = 1.0
    for j in range(q):
        st = {(r - j) % q for r in dead(TRIP, q)}
        I = (q - nu_union([dt, st])) / q
        Imin = min(Imin, I)
    print(f"    minimum I_5(j) over j = {Imin:.3f}  ->  never 0  [expected 0.200]")
    return Imin > 0

def R2(j):
    p = 1.0
    for q in PRIMES:
        nu = nu_union([dead(TWIN,q), {(r-j)%q for r in dead(TWIN,q)}])
        p *= (1 - nu/q) / (1 - 2/q)**2
    return p
def R3(j1, j2):
    p = 1.0
    for q in PRIMES:
        nu = nu_union([dead(TWIN,q), {(r-j1)%q for r in dead(TWIN,q)},
                       {(r-j2)%q for r in dead(TWIN,q)}])
        p *= (1 - nu/q) / (1 - 2/q)**3
    return p
def U(j1, j2):
    den = R2(j1)*R2(j2)*R2(abs(j2-j1))
    return R3(j1,j2)/den if den > 0 else 0.0

def check_irreducible(M=30):
    facts = 0; total = 0
    for j1 in range(1, M+1):
        for j2 in range(1, M+1):
            if j1 == j2: continue
            total += 1
            if abs(U(j1,j2) - 1) < 1e-9: facts += 1
    print(f"(3) U(j1,j2) over {total} lag pairs (window {M}): "
          f"exactly factorizing (U=1) = {facts}  [expected 0]")
    return facts == 0

if __name__ == "__main__":
    a = check_consecutive_twins()
    b = check_twin_triplet()
    c = check_irreducible()
    print("\nALL CHECKS PASS:", a and b and c)
