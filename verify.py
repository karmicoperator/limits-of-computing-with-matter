#!/usr/bin/env python3
"""Verification supplement for The Limits of Computing with Matter, package v4.0.

This program recomputes the derived numerical claims identified in the paper,
checks the formulae through independent identities, exercises the effective-work
bounds with deterministic property tests, and verifies that the LaTeX source
contains the corresponding rounded values.

It verifies mathematical arithmetic and manuscript synchronization. It does not
validate cited source data, fabrication, cooling closure, a physical compute
primitive, trained-model quality, economics, or achieved performance.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

PACKAGE_VERSION = "4.0"
PAPER_TITLE = "The Limits of Computing with Matter"


@dataclass(frozen=True)
class Cube:
    side_m: float = 27.5e-3
    leaves: int = 170
    leaf_thickness_m: float = 100e-6
    gap_m: float = 50e-6
    pitch_m: float = 130e-9
    silicon_density_kg_m3: float = 2329.0
    silicon_molar_mass_kg_mol: float = 28.0855e-3
    avogadro_mol_inv: float = 6.02214076e23

    @property
    def gaps(self) -> int:
        return self.leaves - 1

    @property
    def internal_faces(self) -> int:
        return 2 * self.gaps

    @property
    def cube_volume_m3(self) -> float:
        return self.side_m**3

    @property
    def active_solid_volume_m3(self) -> float:
        return self.side_m**2 * self.leaves * self.leaf_thickness_m

    @property
    def stack_thickness_m(self) -> float:
        return self.leaves * self.leaf_thickness_m + self.gaps * self.gap_m

    @property
    def remaining_margin_m(self) -> float:
        return self.side_m - self.stack_thickness_m

    @property
    def external_area_m2(self) -> float:
        return 6 * self.side_m**2

    @property
    def active_face_area_mm2(self) -> float:
        return self.internal_faces * (self.side_m * 1e3) ** 2

    @property
    def silicon_equivalent_mass_kg(self) -> float:
        return self.active_solid_volume_m3 * self.silicon_density_kg_m3

    @property
    def silicon_equivalent_moles(self) -> float:
        return self.silicon_equivalent_mass_kg / self.silicon_molar_mass_kg_mol

    @property
    def silicon_equivalent_atoms(self) -> float:
        return self.silicon_equivalent_moles * self.avogadro_mol_inv

    @property
    def solid_cube_moles(self) -> float:
        return self.cube_volume_m3 * self.silicon_density_kg_m3 / self.silicon_molar_mass_kg_mol

    @property
    def columns_per_face(self) -> float:
        return (self.side_m / self.pitch_m) ** 2

    @property
    def leaf_shared_columns(self) -> float:
        return self.leaves * self.columns_per_face

    @property
    def face_distinct_columns(self) -> float:
        return self.internal_faces * self.columns_per_face

    def optical_sampling_count(self, wavelength_m: float, numerical_aperture: float) -> float:
        if wavelength_m <= 0 or numerical_aperture <= 0:
            raise ValueError("wavelength and numerical aperture must be positive")
        return 2 * self.external_area_m2 / (wavelength_m / (2 * numerical_aperture)) ** 2

    def site_capacity_bytes(
        self,
        layers_per_face: int,
        usable_fraction: float = 0.10,
        net_bits_per_site: float = 1.0,
    ) -> float:
        if layers_per_face <= 0 or not 0 < usable_fraction <= 1 or net_bits_per_site <= 0:
            raise ValueError("invalid site-capacity parameter")
        return (
            usable_fraction
            * self.columns_per_face
            * self.internal_faces
            * layers_per_face
            * net_bits_per_site
            / 8
        )


@dataclass(frozen=True)
class ThermionicReference:
    temperature_k: float = 400.0
    frequency_hz: float = 100e6
    s: float = 10.0
    ideality: float = 1.0
    wiring_ratio: float = 1.0
    capacitance_ratio: float = 1.0
    recovery_factor: float = 1.0
    boltzmann_j_k: float = 1.380649e-23
    planck_j_s: float = 6.62607015e-34
    electron_charge_c: float = 1.602176634e-19

    @property
    def attempt_frequency_hz(self) -> float:
        return self.boltzmann_j_k * self.temperature_k / self.planck_j_s

    def log_factor_for(self, epsilon: float) -> float:
        if not 0 < epsilon < 1:
            raise ValueError("epsilon must lie in (0,1)")
        return math.log(
            self.attempt_frequency_hz
            / (self.frequency_hz * (-math.log1p(-epsilon)))
        )

    def approximate_log_factor_for(self, epsilon: float) -> float:
        if not 0 < epsilon < 1:
            raise ValueError("epsilon must lie in (0,1)")
        return math.log(self.attempt_frequency_hz / (self.frequency_hz * epsilon))

    @property
    def cq_f(self) -> float:
        return self.electron_charge_c**2 / (self.boltzmann_j_k * self.temperature_k)

    def barrier_j_for(self, epsilon: float) -> float:
        return self.boltzmann_j_k * self.temperature_k * self.log_factor_for(epsilon)

    def gate_swing_v_for(self, epsilon: float, ideality: float | None = None) -> float:
        n = self.ideality if ideality is None else ideality
        return n * self.boltzmann_j_k * self.temperature_k / self.electron_charge_c * self.log_factor_for(epsilon)

    def no_offset_j_for(
        self,
        epsilon: float,
        *,
        ideality: float | None = None,
        wiring_ratio: float | None = None,
        capacitance_ratio: float | None = None,
        recovery_factor: float | None = None,
    ) -> float:
        n = self.ideality if ideality is None else ideality
        r = self.wiring_ratio if wiring_ratio is None else wiring_ratio
        chi = self.capacitance_ratio if capacitance_ratio is None else capacitance_ratio
        eta = self.recovery_factor if recovery_factor is None else recovery_factor
        if n <= 0 or r < 0 or chi <= 0 or eta <= 0:
            raise ValueError("invalid thermionic parameter")
        return (
            0.5
            * (1 + r)
            * n**2
            * self.boltzmann_j_k
            * self.temperature_k
            * self.log_factor_for(epsilon) ** 2
            * chi
            / eta
        )

    def no_offset_j_direct_for(
        self,
        epsilon: float,
        *,
        ideality: float | None = None,
        wiring_ratio: float | None = None,
        capacitance_ratio: float | None = None,
        recovery_factor: float | None = None,
    ) -> float:
        n = self.ideality if ideality is None else ideality
        r = self.wiring_ratio if wiring_ratio is None else wiring_ratio
        chi = self.capacitance_ratio if capacitance_ratio is None else capacitance_ratio
        eta = self.recovery_factor if recovery_factor is None else recovery_factor
        c_cell = chi * self.cq_f
        voltage = self.gate_swing_v_for(epsilon, n)
        return 0.5 * (1 + r) * c_cell * voltage**2 / eta

    def static_threshold_j_for(self, epsilon: float) -> float:
        return (
            0.5
            * (1 + self.wiring_ratio)
            * self.ideality**2
            * self.boltzmann_j_k
            * self.temperature_k
            * math.log(self.s / epsilon) ** 2
            * self.capacitance_ratio
            / self.recovery_factor
        )

    def ev(self, joules: float) -> float:
        return joules / self.electron_charge_c

    def poisson_failure_for_barrier(self, barrier_j: float) -> float:
        gamma = self.attempt_frequency_hz * math.exp(-barrier_j / (self.boltzmann_j_k * self.temperature_k))
        return -math.expm1(-gamma / self.frequency_hz)


@dataclass(frozen=True)
class DensityEnvelope:
    side_mm: float = 27.5
    faces: int = 338
    density_gbit_mm2: float = 28.5
    retention: float = 0.90
    anchor_layers: float = 280.0

    @property
    def area_mm2(self) -> float:
        return self.faces * self.side_mm**2

    def gross_bytes(self, density_gbit_mm2: float) -> float:
        return self.area_mm2 * density_gbit_mm2 * 1e9 / 8

    def net_bytes(self, density_gbit_mm2: float) -> float:
        return self.retention * self.gross_bytes(density_gbit_mm2)

    def linear_density(self, layers: float) -> float:
        return self.density_gbit_mm2 * layers / self.anchor_layers

    def density_for_net_bytes(self, net_bytes: float) -> float:
        return net_bytes * 8 / (self.retention * self.area_mm2 * 1e9)


@dataclass(frozen=True)
class Operator:
    width: int = 8192
    block: int = 128
    rank: int = 128
    top_k: int = 2
    experts: int = 128

    @property
    def dense_macs(self) -> int:
        return self.width**2

    @property
    def structured_base_macs(self) -> int:
        return 2 * self.width * self.block

    @property
    def active_expert_macs(self) -> int:
        return 2 * self.top_k * self.width * self.rank

    @property
    def active_macs(self) -> int:
        return self.structured_base_macs + self.active_expert_macs

    @property
    def work_ratio(self) -> float:
        return self.dense_macs / self.active_macs

    @property
    def base_parameters(self) -> int:
        return 2 * self.width * self.block

    @property
    def expert_parameters(self) -> int:
        return self.experts * 2 * self.width * self.rank

    @property
    def parameters_per_layer(self) -> int:
        return self.base_parameters + self.expert_parameters


@dataclass(frozen=True)
class B300Reference:
    dense_fp8_flops_s: float = 5e15
    dense_nvfp4_flops_s: float = 15e15
    hbm_bytes: float = 288e9
    hbm_bandwidth_b_s: float = 8e12
    maximum_tgp_w: float = 1400.0

    @property
    def dense_fp8_macs_s(self) -> float:
        return self.dense_fp8_flops_s / 2

    @property
    def dense_nvfp4_macs_s(self) -> float:
        return self.dense_nvfp4_flops_s / 2


class Verifier:
    def __init__(self) -> None:
        self.checks = 0
        self.failures: list[str] = []

    def _record(self, name: str, ok: bool, detail: str = "") -> None:
        self.checks += 1
        print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")
        if not ok:
            self.failures.append(name)

    def close(
        self,
        name: str,
        got: float,
        expected: float,
        *,
        rel_tol: float = 1e-11,
        abs_tol: float = 0.0,
    ) -> None:
        ok = math.isfinite(got) and math.isfinite(expected) and math.isclose(
            got, expected, rel_tol=rel_tol, abs_tol=abs_tol
        )
        self._record(name, ok, f"got {got:.15g}; expected {expected:.15g}")

    def equal(self, name: str, got: Any, expected: Any) -> None:
        self._record(name, got == expected, f"got {got!r}; expected {expected!r}")

    def true(self, name: str, condition: bool, detail: str = "") -> None:
        self._record(name, bool(condition), detail)

    def contains(self, name: str, text: str, snippet: str) -> None:
        self._record(name, snippet in text, f"required snippet {snippet!r}")

    def absent(self, name: str, text: str, snippet: str) -> None:
        self._record(name, snippet not in text, f"forbidden snippet {snippet!r}")


def effective_work_arms(
    available_incremental_power_w: float,
    effective_energy_j: float,
    streams: float,
    stream_rate_hz: float,
    effective_units_per_update: float,
    output_bandwidth_b_s: float,
    output_bytes_per_unit: float,
    sites: float,
    site_rate_hz: float,
    site_cycles_per_unit: float,
) -> tuple[float, float, float, float]:
    inputs = (
        available_incremental_power_w,
        effective_energy_j,
        streams,
        stream_rate_hz,
        effective_units_per_update,
        output_bandwidth_b_s,
        output_bytes_per_unit,
        sites,
        site_rate_hz,
        site_cycles_per_unit,
    )
    if min(inputs) <= 0:
        raise ValueError("all effective-work inputs must be positive")
    return (
        available_incremental_power_w / effective_energy_j,
        streams * stream_rate_hz * effective_units_per_update,
        output_bandwidth_b_s / output_bytes_per_unit,
        sites * site_rate_hz / site_cycles_per_unit,
    )


def effective_work_bound(*args: float) -> float:
    return min(effective_work_arms(*args))


def fixed_half_up(value: float, decimals: int) -> str:
    """Format a computed float with conventional decimal half-up rounding."""
    stabilized = round(value, decimals + 8)
    quantum = Decimal(1).scaleb(-decimals)
    return format(Decimal(str(stabilized)).quantize(quantum, rounding=ROUND_HALF_UP), "f")


def sci_tex(value: float, decimals: int) -> str:
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / 10**exponent
    return rf"{mantissa:.{decimals}f}\times10^{{{exponent}}}"


def build_results() -> dict[str, Any]:
    cube = Cube()
    therm = ThermionicReference()
    density = DensityEnvelope()
    op = Operator()
    b300 = B300Reference()

    raw_eps = 1e-3
    unit_delta = 1e-3
    union_eps = unit_delta / op.active_macs
    exact_eps = -math.expm1(math.log1p(-unit_delta) / op.active_macs)
    union_raw_j = therm.no_offset_j_for(union_eps)
    illustrative_transition_contribution_j = op.active_macs * union_raw_j

    sensitivity_eps = [1e-3, 1e-6, 1e-9, union_eps, 1e-12, 1e-15]
    sensitivity = {
        f"{eps:.17g}": {
            "epsilon": eps,
            "lambda": therm.log_factor_for(eps),
            "energy_eV": therm.ev(therm.no_offset_j_for(eps)),
            "rate_at_1p4kW_s_inv": 1400.0 / therm.no_offset_j_for(eps),
        }
        for eps in sensitivity_eps
    }

    optical_nas = [0.5, 1.0, 2.65]
    wavelength = 532e-9
    optical = {str(na): cube.optical_sampling_count(wavelength, na) for na in optical_nas}

    event_powers = [1400.0, 10_000.0, 100_000.0]
    baseline_raw_j = therm.no_offset_j_for(raw_eps)
    raw_power_rows = []
    for power in event_powers:
        rate = power / baseline_raw_j
        raw_power_rows.append(
            {
                "power_w": power,
                "flux_w_cm2": power / (cube.external_area_m2 * 1e4),
                "raw_transitions_s_inv": rate,
                "raw_yield_optical_na1": rate / (optical["1.0"] * therm.frequency_hz),
                "raw_yield_leaf_columns": rate / (cube.leaf_shared_columns * therm.frequency_hz),
            }
        )

    density300 = density.linear_density(300)
    density600 = density.linear_density(600)
    net_pb_threshold_density = density.density_for_net_bytes(1e15)
    threshold_gross_bytes = density.gross_bytes(net_pb_threshold_density)

    intensity_1b = op.active_macs / op.width
    intensity_2b = op.active_macs / (2 * op.width)
    target_macs = 25e15

    result: dict[str, Any] = {
        "metadata": {"package_version": PACKAGE_VERSION, "paper_title": PAPER_TITLE},
        "inputs": {
            "cube": asdict(cube),
            "thermionic": asdict(therm),
            "density": asdict(density),
            "operator": asdict(op),
            "b300": asdict(b300),
        },
        "geometry": {
            "cube_volume_cm3": cube.cube_volume_m3 * 1e6,
            "active_solid_volume_cm3": cube.active_solid_volume_m3 * 1e6,
            "stack_thickness_mm": cube.stack_thickness_m * 1e3,
            "remaining_margin_mm": cube.remaining_margin_m * 1e3,
            "external_area_cm2": cube.external_area_m2 * 1e4,
            "active_face_area_mm2": cube.active_face_area_mm2,
            "silicon_equivalent_mass_g": cube.silicon_equivalent_mass_kg * 1e3,
            "silicon_equivalent_moles": cube.silicon_equivalent_moles,
            "silicon_equivalent_atoms": cube.silicon_equivalent_atoms,
            "solid_cube_moles": cube.solid_cube_moles,
        },
        "thermionic": {
            "attempt_frequency_hz": therm.attempt_frequency_hz,
            "lambda_raw_1e_3": therm.log_factor_for(raw_eps),
            "lambda_approx_raw_1e_3": therm.approximate_log_factor_for(raw_eps),
            "approx_energy_relative_difference": (
                therm.no_offset_j_for(raw_eps)
                - therm.boltzmann_j_k * therm.temperature_k * therm.approximate_log_factor_for(raw_eps) ** 2
            )
            / therm.no_offset_j_for(raw_eps),
            "cq_aF": therm.cq_f * 1e18,
            "barrier_raw_1e_3_eV": therm.ev(therm.barrier_j_for(raw_eps)),
            "static_raw_1e_3_eV": therm.ev(therm.static_threshold_j_for(raw_eps)),
            "no_offset_raw_1e_3_eV": therm.ev(baseline_raw_j),
            "no_offset_static_ratio": baseline_raw_j / therm.static_threshold_j_for(raw_eps),
            "sensitivity": sensitivity,
        },
        "reliability": {
            "unit_transitions": op.active_macs,
            "unit_error_target": unit_delta,
            "union_raw_epsilon": union_eps,
            "exact_independent_raw_epsilon": exact_eps,
            "exact_union_relative_difference": exact_eps / union_eps - 1,
            "union_raw_energy_eV": therm.ev(union_raw_j),
            "exact_raw_energy_eV": therm.ev(therm.no_offset_j_for(exact_eps)),
            "illustrative_transition_contribution_j": illustrative_transition_contribution_j,
            "illustrative_transition_contribution_pj": illustrative_transition_contribution_j * 1e12,
            "illustrative_rate_quotient_at_1p4kw_s_inv": 1400.0 / illustrative_transition_contribution_j,
            "replay_factor_union_epsilon": 1.0 / (1.0 - union_eps),
        },
        "channels": {
            "wavelength_nm": wavelength * 1e9,
            "optical_sampling_counts": optical,
            "columns_per_face": cube.columns_per_face,
            "leaf_shared_columns": cube.leaf_shared_columns,
            "face_distinct_columns": cube.face_distinct_columns,
        },
        "raw_resources": {
            "constituent_cycles_at_100MHz_s_inv": cube.silicon_equivalent_atoms * therm.frequency_hz,
            "power_rows": raw_power_rows,
            "serial_inactive_300": 1 - 1 / 300,
            "serial_inactive_600": 1 - 1 / 600,
        },
        "illustrative_normalization": {
            "energy_arm_s_inv": 1400.0 / illustrative_transition_contribution_j,
            "site_cycle_arm_s_inv": cube.silicon_equivalent_atoms * therm.frequency_hz / op.active_macs,
            "output_arm_per_TBps_1byte_s_inv": 1e12 / op.width,
            "output_arm_per_TBps_2byte_s_inv": 1e12 / (2 * op.width),
        },
        "capacity": {
            "site_300_TB": cube.site_capacity_bytes(300) / 1e12,
            "site_600_TB": cube.site_capacity_bytes(600) / 1e12,
            "site_300_to_b300_hbm": cube.site_capacity_bytes(300) / b300.hbm_bytes,
            "site_600_to_b300_hbm": cube.site_capacity_bytes(600) / b300.hbm_bytes,
            "published_gross_PB": density.gross_bytes(28.5) / 1e15,
            "published_net_PB": density.net_bytes(28.5) / 1e15,
            "linear_density_300_Gb_mm2": density300,
            "linear_density_600_Gb_mm2": density600,
            "linear_net_300_PB": density.net_bytes(density300) / 1e15,
            "linear_net_600_PB": density.net_bytes(density600) / 1e15,
            "net_petabyte_threshold_density_Gb_mm2": net_pb_threshold_density,
            "net_petabyte_threshold_gross_PB": threshold_gross_bytes / 1e15,
            "net_petabyte_threshold_net_PB": density.net_bytes(net_pb_threshold_density) / 1e15,
        },
        "operator": {
            "dense_macs": op.dense_macs,
            "structured_base_macs": op.structured_base_macs,
            "active_expert_macs": op.active_expert_macs,
            "active_macs": op.active_macs,
            "work_ratio": op.work_ratio,
            "base_parameters": op.base_parameters,
            "expert_parameters": op.expert_parameters,
            "parameters_per_layer": op.parameters_per_layer,
            "parameters_100_layers_GB_1byte": op.parameters_per_layer * 100 / 1e9,
            "output_intensity_1byte_MAC_B": intensity_1b,
            "output_intensity_2byte_MAC_B": intensity_2b,
            "required_resident_reuse_at_8TBps_MAC_B": target_macs / b300.hbm_bandwidth_b_s,
            "minimum_output_bandwidth_1byte_TBps": target_macs / intensity_1b / 1e12,
            "minimum_output_bandwidth_2byte_TBps": target_macs / intensity_2b / 1e12,
        },
        "b300": {
            "dense_fp8_macs_s": b300.dense_fp8_macs_s,
            "dense_nvfp4_macs_s": b300.dense_nvfp4_macs_s,
            "fp8_spec_quotient_fJ_MAC": b300.maximum_tgp_w / b300.dense_fp8_macs_s * 1e15,
            "nvfp4_spec_quotient_fJ_MAC": b300.maximum_tgp_w / b300.dense_nvfp4_macs_s * 1e15,
            "tenfold_target_PetaMAC_s": 10 * b300.dense_fp8_macs_s / 1e15,
            "tenfold_energy_ceiling_fJ_MAC": b300.maximum_tgp_w / (10 * b300.dense_fp8_macs_s) * 1e15,
        },
    }
    return result


def verify_formulae(v: Verifier, r: dict[str, Any]) -> None:
    cube = Cube()
    therm = ThermionicReference()
    density = DensityEnvelope()
    op = Operator()

    print("\nINDEPENDENT FORMULA IDENTITIES")
    v.close("cube volume by millimetre conversion", r["geometry"]["cube_volume_cm3"], (27.5**3) / 1000)
    v.close("active volume by leaf sum", r["geometry"]["active_solid_volume_cm3"], 170 * 27.5**2 * 0.1 / 1000)
    v.close("stack plus margin equals side", r["geometry"]["stack_thickness_mm"] + r["geometry"]["remaining_margin_mm"], 27.5)
    v.close("moles recover from atoms", r["geometry"]["silicon_equivalent_atoms"] / cube.avogadro_mol_inv, r["geometry"]["silicon_equivalent_moles"])
    v.close("mass recover from moles", r["geometry"]["silicon_equivalent_moles"] * cube.silicon_molar_mass_kg_mol * 1e3, r["geometry"]["silicon_equivalent_mass_g"])
    v.close("active face area matches cube geometry", r["geometry"]["active_face_area_mm2"], cube.internal_faces * (cube.side_m * 1e3) ** 2)

    for eps_key, row in r["thermionic"]["sensitivity"].items():
        eps = row["epsilon"]
        compact = therm.no_offset_j_for(eps)
        direct = therm.no_offset_j_direct_for(eps)
        v.close(f"thermionic compact/direct energy identity epsilon={eps_key}", compact, direct, rel_tol=2e-14)
        v.close(f"barrier reproduces Poisson error epsilon={eps_key}", therm.poisson_failure_for_barrier(therm.barrier_j_for(eps)), eps, rel_tol=2e-7, abs_tol=1e-18)
        v.close(f"energy-rate inverse identity epsilon={eps_key}", row["rate_at_1p4kW_s_inv"] * compact, 1400.0, rel_tol=2e-14)

    raw_eps = 1e-3
    swing = therm.gate_swing_v_for(raw_eps)
    v.close("thermionic swing creates barrier", therm.electron_charge_c * swing / therm.ideality, therm.barrier_j_for(raw_eps), rel_tol=2e-14)
    static_v = therm.ideality * therm.boltzmann_j_k * therm.temperature_k / therm.electron_charge_c * math.log(therm.s / raw_eps)
    v.close("static threshold current ratio", math.exp(therm.electron_charge_c * static_v / (therm.ideality * therm.boltzmann_j_k * therm.temperature_k)), therm.s / raw_eps, rel_tol=2e-14)

    union_eps = r["reliability"]["union_raw_epsilon"]
    exact_eps = r["reliability"]["exact_independent_raw_epsilon"]
    m = r["reliability"]["unit_transitions"]
    delta = r["reliability"]["unit_error_target"]
    v.true("union bound is sufficient", m * union_eps <= delta * (1 + 1e-15))
    v.close("exact independent composition returns unit target", -math.expm1(m * math.log1p(-exact_eps)), delta, rel_tol=2e-12)
    v.true("union bound is stricter than exact independent requirement", union_eps < exact_eps)
    v.close("illustrative transition contribution is transition count times union energy", r["reliability"]["illustrative_transition_contribution_j"], m * therm.no_offset_j_for(union_eps), rel_tol=2e-14)

    for na_text, count in r["channels"]["optical_sampling_counts"].items():
        na = float(na_text)
        alternate = 8 * cube.external_area_m2 * na**2 / (532e-9) ** 2
        v.close(f"optical count two-form identity NA={na_text}", count, alternate, rel_tol=2e-14)
    v.close("column count from mm/nm units", r["channels"]["columns_per_face"], (27.5e6 / 130) ** 2)

    v.close("site capacity direct area-pitch identity 300", r["capacity"]["site_300_TB"] * 1e12, 0.10 * (cube.side_m / cube.pitch_m) ** 2 * cube.internal_faces * 300 / 8)
    v.close("density area agrees with cube", density.area_mm2, cube.active_face_area_mm2)
    v.close("net-petabyte threshold closes to one PB", r["capacity"]["net_petabyte_threshold_net_PB"], 1.0, rel_tol=2e-14)

    v.equal("active MAC decomposition", op.active_macs, op.structured_base_macs + op.active_expert_macs)
    v.equal("parameter decomposition", op.parameters_per_layer, op.base_parameters + op.expert_parameters)
    v.close("work ratio identity", r["operator"]["work_ratio"] * op.active_macs, op.dense_macs, rel_tol=2e-14)
    v.close("one-byte intensity identity", r["operator"]["output_intensity_1byte_MAC_B"] * op.width, op.active_macs)
    v.close("two-byte bandwidth doubles", r["operator"]["minimum_output_bandwidth_2byte_TBps"], 2 * r["operator"]["minimum_output_bandwidth_1byte_TBps"])


def verify_bounds(v: Verifier) -> None:
    print("\nBOUND AND SCALING PROPERTY TESTS")
    rng = random.Random(20260721)
    for i in range(100):
        args = [10 ** rng.uniform(-3, 6) for _ in range(10)]
        arms = effective_work_arms(*args)
        bound = effective_work_bound(*args)
        if not math.isclose(bound, min(arms), rel_tol=0.0, abs_tol=0.0):
            v.true("effective-work bound returns minimum arm", False, f"case {i}")
            break
    else:
        v.true("effective-work bound returns minimum arm for 100 deterministic cases", True)

    base = [1400.0, 1e-9, 1e9, 1e8, 2.0, 8e12, 16384.0, 1e12, 1e8, 1e3]
    b0 = effective_work_bound(*base)
    for resource_index in (0, 2, 3, 4, 5, 7, 8):
        changed = list(base)
        changed[resource_index] *= 2
        v.true(f"bound monotone in resource input {resource_index}", effective_work_bound(*changed) >= b0)
    for cost_index in (1, 6, 9):
        changed = list(base)
        changed[cost_index] *= 2
        v.true(f"bound antitone in cost input {cost_index}", effective_work_bound(*changed) <= b0)

    def scaling_bound(L: float) -> float:
        return effective_work_bound(
            3 * L**2,
            2.0,
            5 * L**2,
            7.0,
            11.0,
            13 * L**2,
            17.0,
            19 * L**3,
            23.0,
            29.0,
        )

    v.close("surface scaling numerical sanity R(2L)/R(L)", scaling_bound(2.0) / scaling_bound(1.0), 4.0)
    v.close("resident-normalized scaling sanity", (scaling_bound(2.0) / (2.0**3)) / (scaling_bound(1.0) / 1.0**3), 0.5)

    c_t = 3.0
    c_m = 7.0
    tau = 2.0
    dark_10 = max(0.0, 1 - tau * c_t / (c_m * 10))
    dark_100 = max(0.0, 1 - tau * c_t / (c_m * 100))
    v.true("dark-state lower bound increases with size", dark_100 > dark_10)
    v.true("dark-state lower bound lies in [0,1]", 0 <= dark_10 <= 1 and 0 <= dark_100 <= 1)

    # Correct schedule-average output relation. It is intentionally not a pointwise bound on gamma_u.
    actual_updates = 4e6
    interval_s = 2.0
    output_bandwidth = 8e12
    bytes_per_unit = 16384.0
    effective_units = output_bandwidth * interval_s / bytes_per_unit
    average_yield = effective_units / actual_updates
    v.close(
        "schedule-average yield relation at saturated output",
        average_yield,
        output_bandwidth * interval_s / (actual_updates * bytes_per_unit),
        rel_tol=2e-14,
    )


def manuscript_snippets(r: dict[str, Any]) -> list[tuple[str, str]]:
    g = r["geometry"]
    t = r["thermionic"]
    rel = r["reliability"]
    ch = r["channels"]
    rr = r["raw_resources"]
    al = r["illustrative_normalization"]
    cap = r["capacity"]
    op = r["operator"]
    b = r["b300"]
    bi = r["inputs"]["b300"]

    rows = rr["power_rows"]
    sens = list(t["sensitivity"].values())
    return [
        ("cube volume", f"{g['cube_volume_cm3']:.4f}~\\mathrm{{cm^3}}"),
        ("active volume", f"{g['active_solid_volume_cm3']:.4f}~\\mathrm{{cm^3}}"),
        ("stack thickness", f"{g['stack_thickness_mm']:.2f}~mm"),
        ("remaining margin", f"{g['remaining_margin_mm']:.2f}~mm"),
        ("external area", f"{g['external_area_cm2']:.3f}~\\mathrm{{cm^2}}"),
        ("silicon-equivalent mass", f"{g['silicon_equivalent_mass_g']:.3f}~\\mathrm{{g}}"),
        ("silicon-equivalent moles", f"{g['silicon_equivalent_moles']:.4f}~\\mathrm{{mol}}"),
        ("silicon-equivalent atoms", f"{sci_tex(g['silicon_equivalent_atoms'], 3)}\\ \\text{{atoms}}"),
        ("solid cube moles", f"{g['solid_cube_moles']:.3f}~mol"),
        ("attempt frequency", f"{sci_tex(t['attempt_frequency_hz'], 3)}~\\mathrm{{s^{{-1}}}}"),
        ("baseline lambda", f"\\Lambda_\\varepsilon&={t['lambda_raw_1e_3']:.4f}"),
        ("Cq", f"{t['cq_aF']:.3f}~\\mathrm{{aF}}"),
        ("barrier", f"{t['barrier_raw_1e_3_eV']:.5f}~\\mathrm{{eV}}"),
        ("static energy", f"{t['static_raw_1e_3_eV']:.3f}~\\mathrm{{eV}}"),
        ("raw energy", f"{t['no_offset_raw_1e_3_eV']:.4f}~\\mathrm{{eV}}"),
        ("small-error approximation fraction", f"{abs(t['approx_energy_relative_difference'])*1e5:.1f}\\times10^{{-5}}"),
        ("union epsilon", f"{rel['union_raw_epsilon']*1e10:.6f}\\times10^{{-10}}"),
        ("exact epsilon", f"{rel['exact_independent_raw_epsilon']*1e10:.6f}\\times10^{{-10}}"),
        ("union transition energy", fixed_half_up(rel['union_raw_energy_eV'], 4)),
        ("illustrative rate quotient", f"{sci_tex(rel['illustrative_rate_quotient_at_1p4kw_s_inv'], 3)}"),
        ("chi 10 transition rate", sci_tex(rel['illustrative_rate_quotient_at_1p4kw_s_inv'] / 10, 3)),
        ("chi 100 transition rate", sci_tex(rel['illustrative_rate_quotient_at_1p4kw_s_inv'] / 100, 3)),
        ("optical NA 0.5", f"{sci_tex(ch['optical_sampling_counts']['0.5'], 2)}"),
        ("optical NA 1", f"{sci_tex(ch['optical_sampling_counts']['1.0'], 2)}"),
        ("optical NA 2.65", f"{sci_tex(ch['optical_sampling_counts']['2.65'], 2)}"),
        ("columns per face", f"{sci_tex(ch['columns_per_face'], 3)}"),
        ("leaf columns", f"{sci_tex(ch['leaf_shared_columns'], 3)}"),
        ("face columns", f"{sci_tex(ch['face_distinct_columns'], 4)}"),
        ("constituent cycles", f"{sci_tex(rr['constituent_cycles_at_100MHz_s_inv'], 2)}"),
        ("serial 300 inactive", f"{100*rr['serial_inactive_300']:.3f}\\%"),
        ("serial 600 inactive", f"{100*rr['serial_inactive_600']:.3f}\\%"),
        ("raw power 1.4kW rate", f"{sci_tex(rows[0]['raw_transitions_s_inv'], 2)}"),
        ("raw power 10kW rate", f"{sci_tex(rows[1]['raw_transitions_s_inv'], 2)}"),
        ("raw power 100kW rate", f"{sci_tex(rows[2]['raw_transitions_s_inv'], 2)}"),
        ("raw flux 1.4kW", f"{rows[0]['flux_w_cm2']:.1f}"),
        ("raw flux 10kW", f"{rows[1]['flux_w_cm2']:.0f}"),
        ("raw flux 100kW", f"{rows[2]['flux_w_cm2']:.0f}"),
        ("effective output 1byte", sci_tex(al['output_arm_per_TBps_1byte_s_inv'], 3)),
        ("effective output 2byte", sci_tex(al['output_arm_per_TBps_2byte_s_inv'], 3)),
        ("site capacity 300", f"{fixed_half_up(cap['site_300_TB'], 4)}~\\mathrm{{TB}}"),
        ("site capacity 600", f"{fixed_half_up(cap['site_600_TB'], 4)}~\\mathrm{{TB}}"),
        ("HBM ratio 300", f"{cap['site_300_to_b300_hbm']:.2f}"),
        ("HBM ratio 600", f"{cap['site_600_to_b300_hbm']:.2f}"),
        ("face area", f"{int(g['active_face_area_mm2']):,}".replace(",", "{,}") + f"{g['active_face_area_mm2'] % 1:.1f}"[1:] + "~\\mathrm{mm^2}"),
        ("published gross", f"{cap['published_gross_PB']:.4f}~PB"),
        ("published net", f"{cap['published_net_PB']:.4f}~PB"),
        ("linear density 300", f"{cap['linear_density_300_Gb_mm2']:.3f} Gb/mm$^2$"),
        ("linear density 600", f"{cap['linear_density_600_Gb_mm2']:.3f} Gb/mm$^2$"),
        ("linear net 300", f"{cap['linear_net_300_PB']:.4f} PB"),
        ("linear net 600", f"{cap['linear_net_600_PB']:.4f} PB"),
        ("petabyte threshold density", f"{cap['net_petabyte_threshold_density_Gb_mm2']:.3f} Gb/mm$^2$"),
        ("petabyte threshold gross", f"{cap['net_petabyte_threshold_gross_PB']:.4f} PB"),
        ("dense MACs", f"{op['dense_macs']:,}".replace(",", "{,}")),
        ("active MACs", f"{op['active_macs']:,}".replace(",", "{,}")),
        ("work ratio", f"{op['work_ratio']:.4f}"),
        ("router MAC count", f"{8192*128:,}".replace(",", "{,}")),
        ("router-inclusive ratio", f"{op['dense_macs']/(op['active_macs']+8192*128):.4f}"),
        ("expert parameters", f"{op['expert_parameters']:,}".replace(",", "{,}")),
        ("100-layer GB", f"{op['parameters_100_layers_GB_1byte']:.1f}~GB"),
        ("reuse at 8TBps", f"{op['required_resident_reuse_at_8TBps_MAC_B']:.0f} MAC/byte"),
        ("intensity 1byte", f"{op['output_intensity_1byte_MAC_B']:.0f} MAC/byte"),
        ("intensity 2byte", f"{op['output_intensity_2byte_MAC_B']:.0f} MAC/byte"),
        ("output BW 1byte", f"{int(op['minimum_output_bandwidth_1byte_TBps']*10)/10:.1f}"),
        ("output BW 2byte", f"{int(op['minimum_output_bandwidth_2byte_TBps']*10)/10:.1f}"),
        ("B300 FP8", f"{b['dense_fp8_macs_s']/1e15:.1f}~PetaMAC/s"),
        ("B300 NVFP4", f"{b['dense_nvfp4_macs_s']/1e15:.1f}~PetaMAC/s"),
        ("B300 quotient FP8", f"{b['fp8_spec_quotient_fJ_MAC']:.0f}"),
        ("B300 quotient NVFP4", f"{b['nvfp4_spec_quotient_fJ_MAC']:.3f}"),
        ("datasheet FP8 quotient", f"{bi['maximum_tgp_w']/(4.5e15/2)*1e15:.2f}"),
        ("update-product gate", f"{int(b['tenfold_target_PetaMAC_s']*1e15/op['active_macs']/1e8)/10:.1f}" + "\\times10^{9}"),
        ("expert instances to fill net envelope", sci_tex(cap['published_net_PB']*1e15/(op['expert_parameters']/128), 1)),
        ("experts per layer at 1000 layers", sci_tex(cap['published_net_PB']*1e15/(op['expert_parameters']/128)/1000, 1)),
        ("per-expert megabytes", f"{op['expert_parameters']/128/1e6:.2f}~MB"),
        ("egress bits per second", sci_tex(int(op['minimum_output_bandwidth_1byte_TBps']*10)/10*1e12*8, 1)),
        ("egress per sample location", sci_tex(int(op['minimum_output_bandwidth_1byte_TBps']*10)/10*1e12*8/r["channels"]["optical_sampling_counts"]["1.0"], 1)),
        ("datasheet NVFP4 quotient", f"{bi['maximum_tgp_w']/(14e15/2)*1e15:.0f}~fJ/MAC"),
        ("B300 target", f"{b['tenfold_target_PetaMAC_s']:.0f}~PetaMAC/s"),
        ("B300 ceiling", f"{b['tenfold_energy_ceiling_fJ_MAC']:.0f}~\\mathrm{{fJ/MAC}}"),
        ("B300 HBM", f"{bi['hbm_bytes']/1e9:.0f}~GB"),
        ("B300 bandwidth", f"{bi['hbm_bandwidth_b_s']/1e12:.0f}~TB/s"),
        ("B300 TGP", f"{bi['maximum_tgp_w']/1e3:.1f}~kW"),
        ("illustrative transition contribution joules", sci_tex(rel['illustrative_transition_contribution_j'], 5)),
        ("illustrative transition contribution pJ", f"{fixed_half_up(rel['illustrative_transition_contribution_pj'], 3)}~\\mathrm{{pJ}}"),
        ("common site-cycle arm", sci_tex(al['site_cycle_arm_s_inv'], 2)),
        # Reliability sensitivity table rows, generated from computed values.
        ("sensitivity row 1e-3", f"$10^{{-3}}$ & {sens[0]['lambda']:.4f} & {sens[0]['energy_eV']:.3f} & ${sci_tex(sens[0]['rate_at_1p4kW_s_inv'], 2)}$"),
        ("sensitivity row 1e-6", f"$10^{{-6}}$ & {sens[1]['lambda']:.4f} & {sens[1]['energy_eV']:.3f} & ${sci_tex(sens[1]['rate_at_1p4kW_s_inv'], 2)}$"),
        ("sensitivity row 1e-9", f"$10^{{-9}}$ & {sens[2]['lambda']:.4f} & {sens[2]['energy_eV']:.3f} & ${sci_tex(sens[2]['rate_at_1p4kW_s_inv'], 2)}$"),
        ("sensitivity row union", f"${rel['union_raw_epsilon']*1e10:.5f}\\times10^{{-10}}$ & {sens[3]['lambda']:.4f} & {sens[3]['energy_eV']:.3f} & ${sci_tex(sens[3]['rate_at_1p4kW_s_inv'], 2)}$"),
        ("sensitivity row 1e-12", f"$10^{{-12}}$ & {sens[4]['lambda']:.4f} & {sens[4]['energy_eV']:.3f} & ${sci_tex(sens[4]['rate_at_1p4kW_s_inv'], 2)}$"),
        ("sensitivity row 1e-15", f"$10^{{-15}}$ & {sens[5]['lambda']:.4f} & {sens[5]['energy_eV']:.3f} & ${sci_tex(sens[5]['rate_at_1p4kW_s_inv'], 2)}$"),
        # Raw resource table rows, generated from computed values.
        ("raw table 1.4 kW", f"1.4 kW & {rows[0]['flux_w_cm2']:.1f} & ${sci_tex(rows[0]['raw_transitions_s_inv'], 2)}$ & {rows[0]['raw_yield_optical_na1']:.1f} & {rows[0]['raw_yield_leaf_columns']:.2f}"),
        ("raw table 10 kW", f"10 kW & {rows[1]['flux_w_cm2']:.0f} & ${sci_tex(rows[1]['raw_transitions_s_inv'], 2)}$ & {rows[1]['raw_yield_optical_na1']:.0f} & {rows[1]['raw_yield_leaf_columns']:.2f}"),
        ("raw table 100 kW", f"100 kW & {rows[2]['flux_w_cm2']:.0f} & ${sci_tex(rows[2]['raw_transitions_s_inv'], 2)}$ & {rows[2]['raw_yield_optical_na1']:.0f} & {rows[2]['raw_yield_leaf_columns']:.1f}"),
        # Capacity table rows, generated from computed values.
        ("capacity row published", f"28.500 Gb/mm$^2$ & {cap['published_gross_PB']:.4f} PB & {cap['published_net_PB']:.4f} PB"),
        ("capacity row 300", f"{cap['linear_density_300_Gb_mm2']:.3f} Gb/mm$^2$ & {cap['linear_net_300_PB']/0.9:.4f} PB & {cap['linear_net_300_PB']:.4f} PB"),
        ("capacity row threshold", f"{cap['net_petabyte_threshold_density_Gb_mm2']:.3f} Gb/mm$^2$ & {cap['net_petabyte_threshold_gross_PB']:.4f} PB & {cap['net_petabyte_threshold_net_PB']:.4f} PB"),
        ("capacity row 600", f"{cap['linear_density_600_Gb_mm2']:.3f} Gb/mm$^2$ & {cap['linear_net_600_PB']/0.9:.4f} PB & {cap['linear_net_600_PB']:.4f} PB"),
    ]


def verify_manuscript(
    v: Verifier,
    r: dict[str, Any],
    tex_path: Path,
    readme_path: Path,
) -> None:
    print("\nMANUSCRIPT AND PACKAGE SYNCHRONIZATION")
    text = tex_path.read_text(encoding="utf-8")
    for name, snippet in manuscript_snippets(r):
        v.contains(f"manuscript contains {name}", text, snippet)

    v.contains("paper title present", text, PAPER_TITLE)
    v.contains("subtitle present", text, "From thermodynamic bounds to a candidate cubic computer design")
    for retired in ("accepted", "audit"):
        v.absent(f"retired term absent: {retired}", text.lower(), retired)
    v.contains("abstract opens with the mole question", text, "Would a mole of matter, organized as a machine, deliver a mole's worth of computation? Existing answers sit at two extremes")
    v.contains("Feynman attribution follows in the introduction", text, "Feynman's question about manipulating matter at small scales has a direct computational form")
    v.contains("accounting boundary standardized", text, "first externally observable boundary")
    v.contains("same boundary comparison rule present", text, "same completion boundary and acceptance rule")
    v.contains("operational stream definition present", text, "physically distinguishable, separately schedulable update channel")
    v.contains("statistical-independence clarification present", text, "does not mean statistical independence")
    v.contains("common-normalization contribution stated", text, "The contribution is the common normalization")
    v.contains("scaling result demoted to corollary", text, r"\begin{corollary}[Compact surface-limited families]")
    v.contains("dark-state fraction defined before result", text, "Define $D_\\tau(L)$ as the fraction")
    v.contains("serial inactive fraction defined before result", text, "define $D_{\\mathrm{inst}}$ as the instantaneous fraction")
    v.equal("three proposition environments", text.count(r"\begin{proposition}"), 3)
    v.equal("one corollary environment", text.count(r"\begin{corollary}"), 1)
    v.contains("illustrative energy symbol present", text, r"\Eunorm=m_u\Eraw")
    v.contains("illustrative energy caveat present", text, "not asserted as a physical lower bound")
    v.contains("resident-touch intensity stated", text, "one MAC per resident byte")
    v.absent("raw resource ladder removed", text, "fig:ladder")
    v.contains("resource conversion table present", text, "Conversion from inventories to effective work")
    v.contains("cube labeled object of study", text, "object of study and measurement target")
    v.contains("decisive conclusion present", text, "This paper does not demonstrate a molar-scale accelerator")
    v.equal("verdict sentence in abstract and conclusion", text.count("A mole of active matter, organized as a compact irreversible machine, cannot compute like a mole"), 2)
    v.contains("verification limitations stated", text, "does not verify cited source data or physical realization")
    v.contains("verification script referenced generically", text, "accompanying verification script")
    v.contains("code availability statement present", text, "A verification script is published alongside this manuscript")
    for filename in ("verify.py", "results.json", "README.md", "paper.tex", "paper.pdf", "v4.0"):
        v.absent(f"no package filename in manuscript: {filename}", text, filename)
    v.absent("old internal script filename absent", text, "ai_clean")
    v.absent("mole title retired", text, "Computing with a Mole of Matter")
    v.absent("old v2 title absent", text, "Effective-Work Bounds for Compact Physical Computers")
    v.absent("no Unicode em dash", text, chr(0x2014))
    v.absent("no LaTeX em dash token", text, "-" * 3)
    v.contains("AI disclosure present", text, "AI tools were used for drafting and typesetting.")
    v.absent("appendix macro removed", text, "\\appendix")
    v.absent("appendix headings removed", text, "\\section*{Appendix")

    for phrase in (
        "additionally",
        "serves as",
        "stands as",
        "pivotal",
        "underscores",
        "showcases",
        "highlighting",
        "emphasizing",
        "enhancing",
        "delve",
        "tapestry",
        "testament",
        "The design rule is:",
        "The shortest experimental route is:",
    ):
        v.absent(f"AI-style phrase absent: {phrase}", text.lower(), phrase.lower())
    v.absent("no rather-than contrast", text.lower(), "rather than")
    v.absent("no inline bold item headers", text, r"\item \textbf")
    v.absent("no prose bold headers", text, r"\noindent\textbf")
    v.equal("common-normalization sentence appears once", text.count("The contribution is the common normalization"), 1)
    v.true("illustrative wording not overused", text.lower().count("illustrative") <= 12, str(text.lower().count("illustrative")))

    labels = re.findall(r"\\label\{([^}]+)\}", text)
    v.equal("all LaTeX labels unique", len(labels), len(set(labels)))
    refs = re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", text)
    missing = sorted(set(refs) - set(labels))
    v.equal("all internal references resolve in source", missing, [])

    if readme_path.exists():
        readme = readme_path.read_text(encoding="utf-8")
        for snippet in (
            "Script for confirming the math in the paper.",
            "python3 verify.py",
            "paper.tex",
        ):
            v.contains(f"README contains {snippet}", readme, snippet)
    else:
        v.true("README exists", False, str(readme_path))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tex",
        type=Path,
        default=Path(__file__).with_name("paper.tex"),
        help="LaTeX source to check for synchronized displayed values",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path(__file__).with_name("results.json"),
        help="path for the computed result manifest",
    )
    parser.add_argument(
        "--readme",
        type=Path,
        default=Path(__file__).with_name("README.md"),
        help="README file to check for package consistency",
    )
    args = parser.parse_args(argv)

    verifier = Verifier()
    results = build_results()

    print("STRICT NUMERICAL VERIFICATION")
    # Sanity check: a zero result must not pass for a small expected value.
    verifier.true(
        "strict tolerance rejects zero for union epsilon",
        not math.isclose(0.0, results["reliability"]["union_raw_epsilon"], rel_tol=1e-11, abs_tol=0.0),
    )
    verifier.true(
        "strict tolerance rejects zero for illustrative transition contribution",
        not math.isclose(0.0, results["reliability"]["illustrative_transition_contribution_j"], rel_tol=1e-11, abs_tol=0.0),
    )

    verifier.equal("package version metadata", results["metadata"]["package_version"], PACKAGE_VERSION)
    verifier.equal("paper title metadata", results["metadata"]["paper_title"], PAPER_TITLE)
    verify_formulae(verifier, results)
    verify_bounds(verifier)

    if args.tex.exists():
        verify_manuscript(verifier, results, args.tex, args.readme)
    else:
        verifier.true("LaTeX source exists", False, str(args.tex))

    args.json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    verifier.true("JSON result manifest written", args.json.exists() and args.json.stat().st_size > 0, str(args.json))

    print("\nSUMMARY")
    print(f"Checks run: {verifier.checks}")
    if verifier.failures:
        print("Failures:")
        for name in verifier.failures:
            print(f"  - {name}")
        return 1
    print("All numerical, identity, property, manuscript, and package-synchronization checks passed.")
    print("Scope: arithmetic and source synchronization only; no physical validation is implied.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
