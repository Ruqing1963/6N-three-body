#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three-body phase wind tunnel: the irreducible three-twin correlation on the 6N
lattice. For twin centres at N, N+j1, N+j2, the joint singular series normalised
by the three marginals is R3(j1,j2); the irreducible three-body factor is

    U(j1,j2) = R3(j1,j2) / ( R2(j1) R2(j2) R2(j2-j1) ),

where R2 is the pairwise (Part XX) correlation. U==1 would mean the three-body
correlation factorises into pairwise interactions. We find U != 1 everywhere:
a genuine irreducible three-body term, vanishing (U=0, total annihilation) on a
lattice of lags where the six twin-holes tile Z/q at a small prime.
Requires: numpy, matplotlib.
"""
import math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PRIMES = [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]
def dead(q): a = pow(6,-1,q); return {a,(-a)%q}
def nu(shifts,q):
    D=dead(q); U=set()
    for s in shifts: U |= {(r-s)%q for r in D}
    return len(U)
def R2(j):
    p=1.0
    for q in PRIMES: p *= (1-nu([0,j],q)/q)/((1-2/q)**2)
    return p
def R3(j1,j2):
    p=1.0
    for q in PRIMES: p *= (1-nu([0,j1,j2],q)/q)/((1-2/q)**3)
    return p
def U(j1,j2):
    den = R2(j1)*R2(j2)*R2(abs(j2-j1))
    return R3(j1,j2)/den if den>0 else float('nan')

def main():
    M=30
    G=np.full((M+1,M+1), np.nan)
    for j1 in range(1,M+1):
        for j2 in range(1,M+1):
            if j1==j2: continue
            G[j1,j2]=U(j1,j2)
    # report
    finite=[(G[j1,j2],j1,j2) for j1 in range(1,M+1) for j2 in range(1,M+1)
            if j1!=j2 and not math.isnan(G[j1,j2])]
    nz=[t for t in finite if t[0]>1e-9]
    print(f"pairs scanned: {len(finite)}; exactly factorizing (U=1): "
          f"{sum(1 for u,_,_ in finite if abs(u-1)<1e-9)}")
    print(f"annihilation (U=0): {sum(1 for u,_,_ in finite if u<1e-9)} lattice points")
    nz.sort()
    print("strongest finite irreducible terms:")
    for u,j1,j2 in nz[-3:][::-1]: print(f"   U({j1},{j2})={u:.3f}")
    for u,j1,j2 in nz[:3]:        print(f"   U({j1},{j2})={u:.3f}")

    # ---------- figure ----------
    fig, ax = plt.subplots(1,2,figsize=(13.5,5.6))
    fig.suptitle("Three-body twin correlation $U(j_1,j_2)=R_3/(R_2R_2R_2)$: "
                 "irreducible everywhere, annihilating on a lattice",
                 fontsize=12.5, fontweight="bold")
    # panel A: heatmap of U (log color), annihilation cells black
    a=ax[0]
    Z=np.ma.masked_invalid(G[1:,1:])
    Zlog=np.ma.masked_less_equal(Z, 1e-9)
    cmap=plt.cm.RdBu_r.copy(); cmap.set_bad("0.6")
    im=a.imshow(np.log10(Zlog), origin="lower", extent=[1,M,1,M],
                cmap=cmap, vmin=-0.4, vmax=0.4)
    # mark annihilation (U=0)
    ann_j1=[j1 for j1 in range(1,M+1) for j2 in range(1,M+1)
            if j1!=j2 and not math.isnan(G[j1,j2]) and G[j1,j2]<1e-9]
    ann_j2=[j2 for j1 in range(1,M+1) for j2 in range(1,M+1)
            if j1!=j2 and not math.isnan(G[j1,j2]) and G[j1,j2]<1e-9]
    a.scatter(ann_j2, ann_j1, s=10, c="k", marker="x", label="annihilation $U=0$")
    a.set_xlabel("$j_2$"); a.set_ylabel("$j_1$")
    a.set_title("(A) $\\log_{10}U$;  black $\\times$ = total annihilation ($R_3=0$)")
    a.legend(fontsize=8, loc="upper right")
    fig.colorbar(im, ax=a, fraction=0.046, label="$\\log_{10}U$")
    # panel B: a slice + the triple-alignment mechanism
    b=ax[1]
    js=list(range(1,31))
    for j1 in (5,7,1):
        b.plot(js, [U(j1,j2) if j2!=j1 else np.nan for j2 in js], "o-", ms=3, lw=1,
               label=f"$j_1={j1}$")
    b.axhline(1,color="k",ls="--",lw=.8); b.axhline(0,color="0.6",ls=":",lw=.8)
    b.set_xlabel("$j_2$"); b.set_ylabel("$U(j_1,j_2)$")
    b.set_title("(B) slices: $U=1$ never reached; dips to 0 = annihilation lattice")
    b.legend(fontsize=9); b.grid(alpha=.3)
    fig.tight_layout(rect=[0,0,1,0.94])
    fig.savefig("fig_three_body.png", dpi=200); fig.savefig("fig_three_body.pdf")
    print("wrote fig_three_body.{png,pdf}")

if __name__=="__main__":
    main()
