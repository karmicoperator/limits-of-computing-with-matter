#!/usr/bin/env python3
"""Verification script for "The Limits of Computing with Matter" (Amberg, 2026).

Recomputes every load-bearing derived value in the paper (derived is the
paper's untagged default state) from the postulates P1-P9 and the tagged
inputs, and checks it against the value as printed (tolerances cover the
paper's rounding). Sections below mirror the
paper's sections. Run:  python3 verify.py     (standard library only)

Exit code 0 iff every check passes.

Provenance of inputs:
  [measured]  experimental/industrial fact, verified against live sources
              on 2026-07-18 (see bibliography note in the paper)
  [estimated] engineering judgment with stated range
  [chosen]    free design parameter
  [projected] extrapolation beyond demonstrated practice
[estimated]/[chosen]/[projected] values are INPUTS of this script, not outputs;
their ranges are exercised in the S13/S13.1 sections.
"""
import math
import subprocess
import sys

# ---- physical constants (CODATA) ----
k    = 1.380649e-23     # J/K (exact)
q    = 1.602177e-19     # C
h    = 6.62607e-34      # J s
hbar = 1.054572e-34
me   = 9.10938e-31      # kg
eps0 = 8.8542e-12       # F/m
mu0  = 4e-7 * math.pi
c    = 2.9979e8
NA   = 6.02214e23

FAIL, NCHK = [], [0]
def chk(name, got, printed, tol=0.05):
    ok = abs(got / printed - 1) < tol
    NCHK[0] += 1
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: computed {got:.4g}, printed {printed:.4g}")
    if not ok:
        FAIL.append(name)

def inrange(name, got, lo, hi, slack=0.10):
    """Value printed as a range: computed must fall within [lo,hi] +- slack."""
    ok = lo * (1 - slack) <= got <= hi * (1 + slack)
    NCHK[0] += 1
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: computed {got:.4g}, printed range [{lo:.3g}, {hi:.3g}]")
    if not ok:
        FAIL.append(name)

def order(name, got, printed, dex=0.7):
    """Order-of-magnitude claim: |log10(got/printed)| <= dex (two-sided)."""
    ok = abs(math.log10(got / printed)) <= dex
    NCHK[0] += 1
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: computed {got:.4g}, printed ~{printed:.4g} (order check)")
    if not ok:
        FAIL.append(name)

# Note on check depth: some checks below are one-step arithmetic on tagged
# inputs (e.g. duty x planes, cluster count x 700 W). They exist to catch
# transcription errors between derivation and print, not to prove physics;
# the multi-step chains (energy chain, ladder rungs, de-rating) are the
# substantive recomputations.

# =====================================================================
# S2/S4  Postulates and the energy chain (OP-E: 400 K, 100 MHz, eps=1e-3)
# =====================================================================
T, f, eps = 400.0, 1e8, 1e-3
kT  = k * T
nu0 = kT / h                       # P1 attempt rate [transition-state heuristic]
Lam = math.log(nu0 / (f * eps))
print("S4  Energy chain (400 K, 100 MHz, eps=1e-3)")
chk("nu0 [1/s]", nu0, 8.3e12)
chk("Lambda", Lam, 18.2, 0.01)
Vmin = kT / q * Lam
chk("V_min at n=1 [V]", Vmin, 0.63, 0.01)
Cmin = q * q / kT
chk("C_min [aF]", Cmin * 1e18, 4.65, 0.01)
n_ideal = 1.3                       # [estimated 1.2-1.5], S12 geometry
chk("electron count n*Lambda (upper)", n_ideal * Lam, 24, 0.02)
chk("electron count at n=1 (lower)", Lam, 18, 0.02)
chk("P1 attempt-rate term ln(nu0/f)", math.log(nu0 / f), 11.3, 0.05)
chk("nu0 sensitivity: +-10x -> dLambda", math.log(10), 2.3, 0.02)

# S12 cell: cylindrical annulus, channel outer r=30 nm, hole r=45 nm,
# ONO eps_r ~ 5 [estimated], gate height 35 nm
Ccell = 2 * math.pi * eps0 * 5 * 35e-9 / math.log(45 / 30)
chk("C_cell [aF]", Ccell * 1e18, 24.0, 0.02)
chk("C_cell/C_min", Ccell / Cmin, 5.2, 0.01)
eta = 10                            # charge recovery [projected; demonstrated 2-3]
R   = 1                             # wiring ratio [estimated, systolic]
E_floor = 0.5 * (1 + R) * kT * Lam**2 / eta
E_cell  = E_floor * (Ccell / Cmin) * n_ideal**2
chk("floor 0.5(1+R)kT L^2/eta [eV]", E_floor / q, 1.15, 0.01)
chk("E_op cell [J]", E_cell, 1.60e-18, 0.01)
chk("E_op cell [eV]", E_cell / q, 10.0, 0.01)
# S4.4 bookkeeping deltas that separate the published floors
chk("Poisson/equilibrium floor ratio (= Lambda)", Lam, 18, 0.02)
inrange("no-attempt-rate floor lower by", (Lam / math.log(1 / eps))**2, 5, 7)

# S4.5 waterfall
print("S4.5 Energy waterfall")
chk("Landauer kT ln2 [eV]", kT * math.log(2) / q, 0.024, 0.01)
chk("barrier kT*Lambda [eV]", kT * Lam / q, 0.63, 0.01)
chk("waterfall x26 (barrier/Landauer)", Lam / math.log(2), 26, 0.02)
chk("waterfall x1.8 (floor/barrier)", E_floor / (kT * Lam), 1.8, 0.02)
chk("waterfall x8.7 (cell/floor)", E_cell / E_floor, 8.7, 0.01)

# =====================================================================
# S9.2  Comparator (B300 Blackwell Ultra, 2025) [measured; mapping +-2x]
# =====================================================================
print("S9.2 Comparator (B300 Blackwell Ultra, 2025)")
# The chain runs from the events-per-op mapping (triangulated +-2x) times the
# measured dense-FP8 rate; the 2.25e19 ev/s anchor makes events/op exactly 5e3.
B_P, B_fp8, B_ntr, B_f = 1400.0, 4.5e15, 208e9, 1.8e9    # [measured; f estimated]
HX = 5e3 * B_fp8                        # 'one B300' in switch events/s
chk("anchor events/s", HX, 2.25e19, 0.001)
H_J = B_P / HX
chk("J/event", H_J, 6.2e-17, 0.01)
chk("eV/event", H_J / q, 388, 0.01)
chk("fJ/op", B_P / B_fp8 * 1e15, 311, 0.01)
ev_per_op = (B_P / B_fp8) / H_J
chk("events per FP8 op", ev_per_op, 5e3, 0.001)
chk("implied alpha [%]", HX / (B_ntr * B_f) * 100, 6, 0.02)
# mapping triangulation (S9.2)
chk("alpha from 3e3 ev/op [%]", 3e3 * B_fp8 / (B_ntr * B_f) * 100, 3.6, 0.02)
chk("alpha from 1e4 ev/op [%]", 1e4 * B_fp8 / (B_ntr * B_f) * 100, 12, 0.01)
chk("J/event at alpha x3 [eV]", H_J / q / 3, 130, 0.01)
chk("switched C at alpha /3, 0.75 V rail [fF]", 2 * (3 * H_J) / 0.75**2 * 1e15, 0.66, 0.02)
chk("ONO series eps_eff (5/5/5 nm ox/nit/ox)", 15 / (5/3.9 + 5/7.5 + 5/3.9), 4.6, 0.02)
chk("waterfall x39 (B300/cell)", H_J / E_cell, 39, 0.01)
# vintage cross-check (2022 H100, retained in [20])
chk("H100 fJ/op", 700 / 1.979e15 * 1e15, 354, 0.01)
chk("vintage J/op ratio (H100/B300)", (700 / 1.979e15) / (B_P / B_fp8), 1.14, 0.01)
chk("H100 eV/event at shared mapping", (700 / 1.979e15) / 5e3 / q, 442, 0.01)

# =====================================================================
# S12  Machine counts
# =====================================================================
print("S12 Machine")
a = 27.5e-3
pil_slab = (a / 130e-9)**2
planes = 1800                       # [projected; production ~300 in 2-3 decks]
cells = 183 * pil_slab * planes
chk("pillars/slab", pil_slab, 4.5e10, 0.01)
chk("pillars cube-wide (streams, S7)", 183 * pil_slab, 8.2e12, 0.01)
chk("cells", cells, 1.5e16, 0.02)
chk("cells/slab", pil_slab * planes, 8.1e13, 0.01)
chk("weights/slab [TB]", pil_slab * planes / 8 / 1e12, 10, 0.01)
chk("weights [PB]", cells / 8 / 1e15, 1.8, 0.03)
chk("x B300 HBM (288 GB)", cells / 8 / 288e9, 6300, 0.02)
Lcube = 183 * 100e-6 + 183 * 50e-6      # 183 slabs + gaps (rounds to 2.75 cm)
chk("cube edge [cm]", Lcube * 100, 2.75, 0.01)
chk("cube volume [cm3]", Lcube * a * a * 1e6, 21, 0.02)
Vsolid = 183 * 100e-6 * a * a
chk("solid volume [cm3]", Vsolid * 1e6, 13.8, 0.01)
atoms = 5.0e22 * Vsolid * 1e6           # Si-class atomic density [measured]
chk("atoms", atoms, 6.9e23, 0.01)
chk("atoms [mol]", atoms / NA, 1.15, 0.01)
chk("atoms per switch", atoms / cells, 4.7e7, 0.01)
chk("rail swing n*V_min [V]", n_ideal * Vmin, 0.82, 0.01)
chk("rail current at 85 kW [kA]", 85e3 / 0.9 / 1e3, 94, 0.01)
chk("per-slab current [A]", 85e3 / 0.9 / 183, 520, 0.01)
chk("retention barrier, year-scale at bare nu0 [eV]",
    kT * math.log(nu0 * 3.156e7) / q, 1.65, 0.03)
# S12 cooling
wall = 183 * 2 * a * a                  # both faces of every slab
chk("cooled area [cm2]", wall * 1e4, 2768, 0.01)
chk("water at 85 kW [mL/s]", 85e3 / 2.20e6 * 1e3, 39, 0.01)
flux_ref = 85e3 / wall / 1e4
chk("flux at reference [W/cm2]", flux_ref, 31, 0.01)
inrange("CHF margin at reference", 300 / flux_ref, 10, 33)
inrange("CHF margin at reference (upper)", 1000 / flux_ref, 10, 33)

# =====================================================================
# S3/S9.3  Throughput (illustrative corollary; cell-level accounting)
# =====================================================================
print("S9.3 Design-point performance")
chk("B(1.4 kW, one comparator) [ev/s]", 1400 / E_cell, 8.7e20, 0.01)
chk("x B300 at 1.4 kW", 1400 / E_cell / HX, 39, 0.01)
B85 = 85e3 / E_cell
chk("B(85 kW) [ev/s]", B85, 5.3e22, 0.01)
chk("x B300 at 85 kW", B85 / HX, 2400, 0.02)
chk("alpha at 85 kW [%]", B85 / (cells * f) * 100, 3.6, 0.03)
Bcap = cells * 0.1 * f                  # alpha_max = 0.1 [chosen]
chk("B(cap) [ev/s]", Bcap, 1.5e23, 0.02)
chk("x B300 at cap", Bcap / HX, 6600, 0.01)
chk("P(cap) [kW]", Bcap * E_cell / 1e3, 236, 0.03)
flux_cap = Bcap * E_cell / wall / 1e4
chk("flux at cap [W/cm2]", flux_cap, 85, 0.03)
inrange("CHF margin at cap", 300 / flux_cap, 3.5, 12)
inrange("CHF margin at cap (upper)", 1000 / flux_cap, 3.5, 12)
chk("efficiency x (error-tolerant)", H_J / E_cell, 39, 0.01)
Lam_x = math.log(nu0 / (f * 1e-12))
chk("Lambda(exact, eps=1e-12)", Lam_x, 39, 0.01)
chk("exact-compute energy penalty", (Lam_x / Lam)**2, 4.6, 0.01)
chk("efficiency x (exact)", H_J / E_cell / (Lam_x / Lam)**2, 8.5, 0.01)
chk("idealized floor x", H_J / E_floor, 340, 0.01)
chk("equal-throughput cluster [MW]", B85 / HX * B_P / 1e6, 3.3, 0.01)
chk("cluster power ratio (1/39 electricity)", B85 / HX * B_P / 85e3, 39, 0.01)
chk("-- under 4x pessimistic event count", B85 / HX * B_P / 85e3 / 4, 9.7, 0.01)
chk("slab x B300", B85 / HX / 183, 13, 0.02)
chk("slab power [W]", 85e3 / 183, 465, 0.01)
chk("slab x B300 at cap", Bcap / HX / 183, 36, 0.01)
chk("per-cell power at cap duty [pW]", E_cell * 0.1 * f * 1e12, 16, 0.01)
chk("per-cell power at reference [pW]", E_cell * (B85 / (cells * f)) * f * 1e12, 6, 0.05)
chk("pillar weight column [bytes]", planes / 8, 225, 0.01)

# =====================================================================
# S5  The smallest switch (300 K, 3 GHz, eps=1e-3)
# =====================================================================
print("S5  Smallest switch")
kT3  = k * 300
nu03 = kT3 / h
Lam3 = math.log(nu03 / (3e9 * 1e-3))
chk("Lambda(300K,3GHz)", Lam3, 14.6, 0.01)
Vm3  = kT3 / q * Lam3
kap  = math.sqrt(2 * 0.19 * me * q * Vm3) / hbar   # m*=0.19 me, Si [measured]
chk("L_min [nm]", Lam3 / (2 * kap) * 1e9, 5.3, 0.01)
chk("C_min-sized plate gate, eps_r=10, t=1nm [nm]",
    math.sqrt(Cmin * 1e-9 / (eps0 * 10)) * 1e9, 7, 0.05)

# =====================================================================
# S6  Cooling bound
# =====================================================================
print("S6  Cooling")
chk("internal-flux ceiling, conservative [MW]", 300e4 * wall / 1e6, 0.83, 0.01)
chk("internal-flux ceiling, optimistic [MW]", 1000e4 * wall / 1e6, 2.77, 0.01)

# =====================================================================
# S7  Addressing bound
# =====================================================================
print("S7  Addressing")
Nch = Vsolid / (100e-9)**3
chk("optical regions (100 nm voxels)", Nch, 1.4e16, 0.02)
chk("atoms per optical channel", atoms / Nch, 5e7, 0.01)
E_ph = h * c / 1e-10                     # 1 Angstrom quantum
chk("X-ray quantum [keV]", E_ph / q / 1e3, 12.4, 0.01)
chk("X-ray / cohesive energy (~6 eV)", E_ph / q / 6, 2e3, 0.05)
chk("X-ray deficit at design cell", E_ph / E_cell, 1.2e3, 0.04)
chk("magnet stress at 100 T [GPa]", 100**2 / (2 * mu0) / 1e9, 4.0, 0.01)
chk("magnetic channels/axis", 28.02e9 * 100 / 1e6, 3e6, 0.07)
chk("in-flight per pillar", 0.1 * planes, 180, 0.01)
chk("concurrent events at duty cap", 183 * pil_slab * 180, 1.5e15, 0.02)

# =====================================================================
# S7.5  Assembled bound (Figure 1 rungs; heat at internal-flux ceiling 2.8 MW)
# =====================================================================
print("S7.5 Assembled bound / Figure 1")
Pmax = 2.8e6
chk("naive mole rung [ev/s]", atoms * nu0, 5.8e36, 0.01)
chk("Landauer-capped rung [ev/s]", Pmax / (kT * math.log(2)), 7e26, 0.05)
chk("barrier rung [ev/s]", Pmax / (kT * Lam), 3e25, 0.08)
inrange("charge-floor rung [ev/s]", Pmax / (0.5 * kT * Lam**2), 3e24, 4e24)
chk("channel rung [ev/s]", Nch * f, 1.4e24, 0.02)
chk("Landauer rung x", Pmax / (kT * math.log(2)) / HX, 3e7, 0.1)
chk("barrier rung x", Pmax / (kT * Lam) / HX, 1.2e6, 0.03)
chk("floor rung x", Pmax / (0.5 * kT * Lam**2) / HX, 1.4e5, 0.04)
chk("channel rung x", Nch * f / HX, 6e4, 0.03)
chk("design entry below assembled bound", (Nch * f) / B85, 26, 0.03)
chk("face-feed vs internal flux (~4x)", 12e6 / Pmax, 4, 0.08)
chk("Fig.1 cut: mole -> Landauer", atoms * nu0 / (Pmax / (kT * math.log(2))), 8e9, 0.02)
chk("Fig.1 cut: Landauer -> barrier", Lam / math.log(2), 26, 0.02)
chk("Fig.1 cut: barrier -> floor", Lam / 2, 9, 0.02)
chk("Fig.1 margin: cap to channel rung", Nch * f / Bcap, 9, 0.05)

# =====================================================================
# S9.1  OP-T corner (836 K, 3 GHz, eps=1e-3, R=3, 21 cm3 monolith)
# =====================================================================
print("S9.1 OP-T corner")
kT8  = k * 836
Lam8 = math.log((kT8 / h) / (3e9 * 1e-3))
E_opt = 0.5 * 4 * kT8 * Lam8**2
chk("E_op(OP-T) [kT]", E_opt / kT8, 485, 0.01)
chk("E_op(OP-T) [J]", E_opt, 5.6e-18, 0.01)
N_opt = 21 * 0.3 / (15e-9 * 15e-9 * 30e-9 * 1e6)   # 30% fill, 15x15x30 nm cells
chk("N(OP-T) switches", N_opt, 9.3e17, 0.01)
P_opt = 2.32e6 * 412 * 0.7 * a * a * 25            # Ga single-pass [measured rho*cp]
chk("P(OP-T) [MW]", P_opt / 1e6, 12.6, 0.01)
chk("B(OP-T) [ev/s]", P_opt / E_opt, 2.3e24, 0.03)
chk("alpha(OP-T)", P_opt / E_opt / (N_opt * 3e9), 8e-4, 0.02)
chk("power density [kW/cm3]", P_opt / 1e3 / 21, 610, 0.02)
chk("X-ray deficit at OP-T", E_ph / E_opt, 355, 0.01)

# =====================================================================
# S9.5  Headroom ladder and delivery
# =====================================================================
print("S9.5 Headroom")
chk("x B300 at conservative CHF", 300e4 * wall / E_cell / HX, 23000, 0.01)
B_a1 = cells * 1.0 * f
chk("B(alpha=1) [ev/s]", B_a1, 1.5e24, 0.02)
chk("x B300 at alpha=1", B_a1 / HX, 66000, 0.01)
chk("P(alpha=1) [MW]", B_a1 * E_cell / 1e6, 2.36, 0.03)
chk("x B300 at optimistic CHF", 1000e4 * wall / E_cell / HX, 77000, 0.01)
Lam_1G = math.log(nu0 / (1e9 * eps))
chk("1 GHz energy falls x", (Lam / Lam_1G)**2, 1.3, 0.01)
chk("busbar at 0.83 MW [MA]", 300e4 * wall / 0.9 / 1e6, 0.9, 0.03)
chk("busbar at 2.4 MW [MA]", 2.4e6 / 0.9 / 1e6, 2.6, 0.03)
chk("water at 0.83 MW [L/s]", 300e4 * wall / 2.2e6, 0.4, 0.06)
chk("water at 2.77 MW [L/s]", 1000e4 * wall / 2.2e6, 1.3, 0.04)

# =====================================================================
# S10  Reversible limit / clocking
# =====================================================================
print("S10 Clocking")
Cagg = cells * 2 * Ccell                 # ~48 aF/cell incl. wiring (R=1)
chk("aggregate C [F]", Cagg, 0.7, 0.02)
order("L for 100 MHz resonance, whole body [H]",
      1 / ((2 * math.pi * 1e8)**2 * Cagg), 1e-18)
order("L per 1e-6 domain [H]",
      1 / ((2 * math.pi * 1e8)**2 * Cagg * 1e-6), 1e-12)
chk("body inductance mu0*a [nH]", mu0 * a * 1e9, 35, 0.02)
chk("self-resonance [kHz]",
    1 / (2 * math.pi * math.sqrt(mu0 * a * Cagg)) / 1e3, 1.0, 0.05)
chk("wire Q at 15 nm, 3 GHz", 2 * math.pi * 3e9 * mu0 * (15e-9)**2 / 1e-7, 5e-5, 0.07)

# =====================================================================
# S11  Atomic/cryogenic limit
# =====================================================================
print("S11 Spins")
gam = 1.760859e11                        # electron gyromagnetic ratio [rad/s/T]
r = 20e-9
nu_dd = mu0 * hbar * gam**2 / (4 * math.pi * r**3) / (2 * math.pi)
chk("dipolar coupling at 20 nm [kHz]", nu_dd / 1e3, 6.5, 0.02)
n_spin = 1 / (20e-7)**3                  # cm^-3 at 20 nm spacing
chk("spin density below Mott (3.74e18)", n_spin / 3.74e18, 0.033, 0.02)
spins = Vsolid * 1e6 / (20e-7)**3
chk("spins in solid volume", spins, 1.7e18, 0.02)
chk("copies per optical channel", spins / Nch, 125, 0.02)

# =====================================================================
# S13  Sensitivity ranges (efficiency x, central 43)
# =====================================================================
print("S13 Sensitivity")
eff = H_J / E_cell
chk("alpha x2 (unfavorable): efficiency", eff / 2, 19, 0.03)
chk("alpha /2 (favorable): efficiency", eff * 2, 78, 0.01)
# correlation structure: both headline figures scale as 1/alpha; ratio locked
chk("headline ratio locked at power ratio", 85e3 / B_P, 61, 0.01)
chk("alpha x2 (unfavorable): throughput x", B85 / (2 * HX), 1200, 0.02)
chk("eta_rec = 3 (S13.1 de-rate basis)", eff * 3 / 10, 12, 0.04)
chk("eta_rec = 2 (basis low end)", eff * 2 / 10, 7.8, 0.01)
chk("eta_rec = 5 (Tier-2 gate)", eff * 5 / 10, 19, 0.03)
chk("eta_rec = 16 (step-count ceiling)", eff * 16 / 10, 62, 0.01)
chk("n = 1.2", eff * n_ideal**2 / 1.2**2, 46, 0.015)
chk("n = 1.5", eff * n_ideal**2 / 1.5**2, 29, 0.01)
chk("C ratio 4 (ONO band low)", eff * (Ccell / Cmin) / 4, 50, 0.01)
chk("C ratio 6 (ONO band high)", eff * (Ccell / Cmin) / 6, 33, 0.02)
chk("exact low (alpha x2)", eff / (Lam_x / Lam)**2 / 2, 4.3, 0.02)
chk("exact high (alpha /2)", eff / (Lam_x / Lam)**2 * 2, 17, 0.01)

# =====================================================================
# S13.1  Compounded projections: the demonstrated-practice case
# =====================================================================
print("S13.1 Demonstrated-practice case")
E_der = E_cell * 10 / 3                  # eta_rec 10 -> 3 (demonstrated)
chk("E_op de-rated [eV]", E_der / q, 33.4, 0.01)
B_der = 85e3 / E_der
chk("B de-rated at 85 kW [ev/s]", B_der, 1.6e22, 0.01)
chk("x B300, three projections de-rated", B_der / HX, 710, 0.01)
chk("efficiency de-rated", H_J / E_der, 12, 0.04)
cells_prod = cells * 300 / 1800          # planes 1800 -> 300 (production, 2-3 decks)
chk("cells at 300 planes", cells_prod, 2.5e15, 0.02)
chk("duty cap at 300 planes [ev/s] (does not bind)", cells_prod * 0.1 * f, 2.5e22, 0.02)
assert cells_prod * 0.1 * f > B_der, "duty cap must not bind in the 1,600x row"
chk("weights at 300 planes [PB]", cells_prod / 8 / 1e15, 0.3, 0.03)
cells_cmos = cells_prod / 10             # IGZO -> bonded-CMOS retreat (10x density)
Bcap_cmos = cells_cmos * 0.1 * f
chk("duty cap, bonded-CMOS retreat [ev/s] (binds)", Bcap_cmos, 2.5e21, 0.02)
chk("x B300, all four de-rated", Bcap_cmos / HX, 110, 0.01)
chk("power, all four de-rated [kW]", Bcap_cmos * E_der / 1e3, 13, 0.02)
chk("weights, retreat [TB]", cells_cmos / 8 / 1e12, 31, 0.01)
chk("x B300 HBM, retreat", cells_cmos / 8 / 288e9, 110, 0.04)
chk("fully de-rated floor at 10^2.0", math.log10(Bcap_cmos / HX), 2.0, 0.02)
chk("three-de-rated at 10^2.8", math.log10(B_der / HX), 2.8, 0.02)
chk("projections buy ~20x vs de-rated corner", B85 / HX / (Bcap_cmos / HX), 20, 0.08)
chk("recovery alone buys x3.3 at fixed power", (B85 / HX) / (B_der / HX), 3.3, 0.02)
chk("composite best-demonstrated (430 planes) floor", cells * 430 / 1800 / 10 * 0.1 * f / HX, 160, 0.03)
# graded floor: process-grounded 600-plane proposal (4-5 demonstrated deck
# cycles, ~30 um stack; MSA-CBA bonding route reaches the same count [40])
cells_600 = cells * 600 / 1800
chk("proposal floor: duty cap at 600-plane retreat [ev/s]", cells_600 / 10 * 0.1 * f, 4.9e21, 0.02)
chk("proposal floor: x B300", cells_600 / 10 * 0.1 * f / HX, 220, 0.01)
chk("proposal floor: power [kW]", cells_600 / 10 * 0.1 * f * E_der / 1e3, 26, 0.02)
chk("proposal floor: weights [TB]", cells_600 / 10 / 8 / 1e12, 61, 0.01)
chk("600 planes = 4-5 demonstrated decks", 600 / 143, 4.2, 0.02)
chk("stack height at 600 planes [um]", 600 * 50e-9 * 1e6, 30, 0.01)
chk("vs tallest shipping stack (~2x)", 600 * 50 / (321 * 45), 2.1, 0.03)
import math as _m
chk("proposal floor at 10^2.3", _m.log10(cells_600 / 10 * 0.1 * f / HX), 2.3, 0.02)

# =====================================================================
# S13.3  Proof-of-concept tiers (mini-cube numbers)
# =====================================================================
print("S13.3 PoC tiers")
mini = 16 * (1e-2 / 130e-9)**2 * 180
chk("mini-cube cells", mini, 1.7e13, 0.01)
chk("mini-cube spec [ev/s]", mini * 0.1 * f, 1.7e20, 0.01)
chk("mini-cube spec x B300", mini * 0.1 * f / HX, 7.6, 0.01)
chk("mini-cube spec power [W]", mini * 0.1 * f * E_cell, 270, 0.02)
chk("Tier-3 pass threshold x B300", 2e19 / HX, 0.9, 0.02)
chk("Tier-2 worst-case pass [eV]", (30e-18 / Cmin) * 1.45**2 * kT * Lam**2 / 10 / q, 15.6, 0.02)

# =====================================================================
print()
print("NOT RECOMPUTED (inputs not fully pinned in the paper; tagged there):")
for note in [
    "S6  face-feed ~12 MW and pumping-power <1% (need channel geometry)",
    "S7  word-plane RC: ~ps segment vs ~9 us slab-global (needs sheet R, C')",
    "S9.4 volumetric-density ~1e6x and communication-radius 1e3-1e4x (order-of-magnitude)",
    "S10 loaded-line velocity ~1e2 m/s (needs per-length L', C')",
    "S11 net ~450x spin advantage (chain includes cryo wall-plug range)",
    "S12 GAA natural length ~12 nm and n ~ 1.3 (TCAD-level, PoC Tier 1 measures)",
]:
    print("  -", note)

print()
try:
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True, timeout=5).stdout.strip()
except Exception:
    sha = ""
print(f"checks: {NCHK[0]}, failures: {len(FAIL)}" + (f"  (git {sha})" if sha else ""))
print("RESULT:", "ALL CHECKS PASS" if not FAIL else f"FAILURES: {FAIL}")
sys.exit(0 if not FAIL else 1)
