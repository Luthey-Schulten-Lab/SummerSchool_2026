"""
Driver for successive btree_chromo / LAMMPS runs (Syn3A replication + looping).

Usage:
    python3 run_btree_chromo.py <seed> <run_name> [start_min] [end_min] [is_restart]
        [v_replication] [n_smc] [tau_basal] [v_translocation]

SMC looping parameters (N, v_translocation, tau_basal, tau_stall, tau_bypass):
    basal_death_prob = 20 / (tau_basal * v_translocation)
    stall_death_prob = 20 / (tau_stall * v_translocation)
    bypass = 20 / (tau_bypass * v_translocation)
    knockoff = 1 - bypass
    numSmc = N

Fork replication and SMC translocation use separate speeds:
    v_replication = 100 bp/s  (transform / map_replication batch amounts)
    v_translocation = 350 bp/s (translocate directive and SMC loop kinetics)

Summer-school defaults: N=50, v_translocation=500 bp/s, tau_basal=100 s, tau_stall=100 s,
    tau_bypass=7 s.

Replication is batched in 60 s of bio time per minute (12 s with BD + 48 s without).
Per-batch amounts scale with batch duration:
    replicate_amount = v_replication/10 * seconds_per_batch  (beads)
    translocate_amount = v_translocation/20 * seconds_per_batch
Minute 30 uses 0.4 s BD batches (30 repeats) for higher-resolution VMD output.
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass

import numpy as np

# --- Paths -------------------------------------------------------------------

base_dir = os.path.dirname(os.path.abspath(__file__))
dna_monomers_dir = os.path.join(base_dir, "../data/coords/")
output_dir = os.path.join(base_dir, "../data/")
template_dir = os.path.join(base_dir, "./")
template_file = os.path.join(template_dir, "template_replicate.inp")

btree_chromo_executable = os.environ.get(
    "BTREE_CHROMO",
    "/ps/btree_chromo_gpu/build/apps/btree_chromo",
)

DEFAULT_LD_LIBRARY_PATH = (
    "/usr/local/lib64:/usr/local/lib:/usr/local/nvidia/lib:"
    "/usr/local/nvidia/lib64:/.singularity.d/libs:/Software/LAMMPS/OMP_GPU_Kokkos/lib"
)

# --- Cell-cycle model parameters ---------------------------------------------
# minutes 0–59: growing sphere; minutes 60–90: twin spheres (division)

R0 = 2000.0
R_FUDGE = 160.0
DR_DT = 10.0  # radius growth per minute

T_DIV = 60
T_END = 90
R_DIV = 2160.0
R_AMP = 600.0
H0 = 30.0
H_SPAN = 2130.0

BD_STEPS = 20_000
BD_OUT_FACTOR = 2
BD_DIV0 = 20_000
BD_DIV1 = 200_000

N_CHROMO = 54338
RIBOS_PER_MIN = 5
LOOP_EQUIL = "translocate:360000,F"  # 1 h bio time at 100 beads/s
LOOP_PARAMS_FILE = os.path.join(template_dir, "loop_params.txt")

# SMC looping parameters (see Minimal_Cell_ChromosomeSegregation submit_jobs_3d_sweep.sh)
N_SMC_DEFAULT = 50
V_REPLICATION_DEFAULT = 100
V_TRANSLOCATION_DEFAULT = 500
TAU_BASAL_DEFAULT = 100
TAU_STALL_DEFAULT = 100
TAU_BYPASS_DEFAULT = 7

# Replication batching: 60 s bio time per minute = 12 s with BD + 48 s without
DEFAULT_SECONDS_PER_BATCH = 2.0
BD_BIO_SECONDS = 12.0
NO_BD_BIO_SECONDS = 48.0
HIGH_RES_MINUTE = 30
HIGH_RES_SECONDS_PER_BATCH = 0.4


@dataclass(frozen=True)
class CellBoundary:
    """Growing or dividing cell boundary geometry."""

    sphere_radius: float
    sphere_height: float
    separation_axis: np.ndarray
    load_boundary: str


@dataclass(frozen=True)
class SmcParams:
    """SMC loop-extrusion parameters and derived death probabilities."""

    n_smc: int
    v_translocation: int
    tau_basal: int
    tau_stall: int
    tau_bypass: int

    @property
    def basal_death_prob(self) -> float:
        return 20.0 / (self.tau_basal * self.v_translocation)

    @property
    def stall_death_prob(self) -> float:
        return 20.0 / (self.tau_stall * self.v_translocation)

    @property
    def bypass(self) -> float:
        return 20.0 / (self.tau_bypass * self.v_translocation)

    @property
    def knockoff(self) -> float:
        return 1.0 - self.bypass

    @property
    def s(self) -> int:
        return self.n_smc * self.v_translocation * self.tau_basal


@dataclass(frozen=True)
class GlobalSimulationArgs:
    seed: int
    run_name: str
    start_time: int
    end_time: int
    is_restart: bool
    v_replication: int
    smc: SmcParams


@dataclass(frozen=True)
class ReplicationBatching:
    """Per-minute replication loop counts and amounts for template_replicate.inp."""

    seconds_per_batch_bd: float
    seconds_per_batch_no_bd: float
    repeat_bd: int
    repeat_no_bd: int
    transform_bd: str
    transform_no_bd: str
    translocate_bd: int
    translocate_no_bd: int


@dataclass(frozen=True)
class PerMinuteDirectives:
    """btree_chromo template placeholders that vary each biological minute."""

    input_state: str
    load_loops: str
    equilibrate_loops: str
    append_string: str
    run_dynamics: str


def parse_args(argv: list[str]) -> GlobalSimulationArgs:
    return GlobalSimulationArgs(
        seed=int(argv[1]),
        run_name=argv[2],
        start_time=int(argv[3]) if len(argv) > 3 else 0,
        end_time=int(argv[4]) if len(argv) > 4 else 90,
        is_restart=(argv[5].lower() == "true") if len(argv) > 5 else False,
        v_replication=int(argv[6]) if len(argv) > 6 else V_REPLICATION_DEFAULT,
        smc=SmcParams(
            n_smc=int(argv[7]) if len(argv) > 7 else N_SMC_DEFAULT,
            v_translocation=int(argv[9]) if len(argv) > 9 else V_TRANSLOCATION_DEFAULT,
            tau_basal=int(argv[8]) if len(argv) > 8 else TAU_BASAL_DEFAULT,
            tau_stall=TAU_STALL_DEFAULT,
            tau_bypass=TAU_BYPASS_DEFAULT,
        ),
    )


def write_loop_params(smc: SmcParams, path: str = LOOP_PARAMS_FILE) -> None:
    """Write loop_params.txt from SMC parameters."""
    content = (
        f"basal_death_prob={smc.basal_death_prob:.6f}\n"
        f"step_prob=1.0\n"
        f"stall_death_prob={smc.stall_death_prob:.6f}\n"
        f"knockoff={smc.knockoff:.6f}\n"
        f"bypass={smc.bypass:.6f}\n"
        f"N={N_CHROMO}\n"
        f"numSmc={smc.n_smc}\n"
        f"smcWidth=1\n"
    )
    with open(path, "w") as f:
        f.write(content)
    print(f"Wrote {path}")


def build_runtime_env() -> dict[str, str]:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = DEFAULT_LD_LIBRARY_PATH
    return env


def replicate_amount(v_replication: float, seconds_per_batch: float) -> int:
    """Beads replicated per batch: (v_replication bp/s) / (10 bp/bead) * batch duration."""
    return int(round(v_replication / 10.0 * seconds_per_batch))


def translocate_amount(v_translocation: float, seconds_per_batch: float) -> int:
    """SMC translocation per batch: (v_translocation bp/s) / (20 bp/bead-step) * batch duration."""
    return int(round(v_translocation / 20.0 * seconds_per_batch))


def transform_directive(v_replication: float, seconds_per_batch: float) -> str:
    """btree_chromo transform string with matched cw/ccw replication amounts."""
    beads = replicate_amount(v_replication, seconds_per_batch)
    return f"transform:m_cw{beads}_ccw{beads}"


def seconds_per_batch_bd(timestep: int) -> float:
    """Batch duration for replication loops that include BD dynamics."""
    if timestep == HIGH_RES_MINUTE:
        return HIGH_RES_SECONDS_PER_BATCH
    return DEFAULT_SECONDS_PER_BATCH


def compute_replication_batching(
    timestep: int, v_replication: int, v_translocation: int
) -> ReplicationBatching:
    """Derive repeat counts and per-batch amounts for one biological minute."""
    spb_bd = seconds_per_batch_bd(timestep)
    spb_no_bd = DEFAULT_SECONDS_PER_BATCH
    repeat_bd = int(round(BD_BIO_SECONDS / spb_bd))
    repeat_no_bd = int(round(NO_BD_BIO_SECONDS / spb_no_bd))
    return ReplicationBatching(
        seconds_per_batch_bd=spb_bd,
        seconds_per_batch_no_bd=spb_no_bd,
        repeat_bd=repeat_bd,
        repeat_no_bd=repeat_no_bd,
        transform_bd=transform_directive(v_replication, spb_bd),
        transform_no_bd=transform_directive(v_replication, spb_no_bd),
        translocate_bd=translocate_amount(v_translocation, spb_bd),
        translocate_no_bd=translocate_amount(v_translocation, spb_no_bd),
    )


def coord_paths(run_name: str, timestep: int) -> tuple[str, str]:
    dna_path = f"{dna_monomers_dir}dna_{run_name}_{timestep}.bin"
    ribo_path = f"{dna_monomers_dir}ribo_{run_name}_{timestep}.bin"
    return dna_path, ribo_path


def _random_point_in_sphere(radius: float) -> np.ndarray:
    while True:
        point = np.random.uniform(-radius, radius, 3)
        if np.linalg.norm(point) <= radius:
            return point


def add_ribosomes_to_bin(
    bin_path: str,
    n_ribosomes: int,
    sphere_radius: float,
    sphere_height: float,
    separation_axis: np.ndarray,
    order: str = "row",
    seed: int | None = None,
) -> None:
    """Append ribosome coordinates inside the current cell boundary."""
    if seed is not None:
        np.random.seed(seed)

    data = np.fromfile(bin_path, dtype=np.float64)
    if data.size % 3 != 0:
        raise ValueError("Binary file size is not divisible by 3")
    n_existing = data.size // 3

    if order == "row":
        coords = data.reshape((n_existing, 3))
    elif order == "col":
        coords = np.column_stack(
            (
                data[0:n_existing],
                data[n_existing : 2 * n_existing],
                data[2 * n_existing : 3 * n_existing],
            )
        )
    else:
        raise ValueError("order must be 'row' or 'col'")

    axis = separation_axis / np.linalg.norm(separation_axis)
    sphere1_center = -sphere_height * axis
    sphere2_center = sphere_height * axis

    new_coords = []
    while len(new_coords) < n_ribosomes:
        center = sphere1_center if np.random.rand() < 0.5 else sphere2_center
        new_coords.append(_random_point_in_sphere(sphere_radius) + center)

    combined = np.vstack([coords, np.array(new_coords)])

    if order == "row":
        combined.astype(np.float64).tofile(bin_path)
    else:
        out = np.concatenate([combined[:, 0], combined[:, 1], combined[:, 2]])
        out.astype(np.float64).tofile(bin_path)

    print(f"Appended {n_ribosomes} new coordinates to '{bin_path}'.")
    print(f"Final total: {combined.shape[0]} coordinates.")


def division_axis_from_dna(dna_monomers_path: str) -> np.ndarray:
    """Unit vector between sister-chromosome centers of mass (division axis)."""
    with open(dna_monomers_path, "rb") as f:
        dna_bin = np.fromfile(f, dtype=np.float64, count=-1)
    dna_coords = dna_bin.reshape((3, dna_bin.shape[0] // 3), order="F").T
    left = dna_coords[:N_CHROMO]
    right = dna_coords[N_CHROMO:]
    axis = np.average(left, axis=0) - np.average(right, axis=0)
    return axis / np.linalg.norm(axis)


def _spherical_boundary(sphere_radius: float) -> str:
    return f"spherical_bdry:{sphere_radius:.1f}, 0.0, 0.0, 0.0"


def _overlapping_spheres_boundary(
    sphere_height: float, sphere_radius: float, separation_axis: np.ndarray
) -> str:
    axis_string = ", ".join(f"{x:.3f}" for x in separation_axis)
    return (
        f"overlapping_spheres_bdry:{sphere_height:.1f}, {sphere_radius:.1f}, "
        f"0.0, 0.0, 0.0, {axis_string}"
    )


def _soft_harmonic_dynamics(bd_steps: int) -> str:
    output_steps = bd_steps * BD_OUT_FACTOR
    return f"simulator_run_soft_harmonic:{bd_steps},1000,{output_steps},append,nofirst"


def compute_run_dynamics(timestep: int) -> str:
    """BD step count for the post-replication equilibration (scales up during division)."""
    div_span = T_END - T_DIV
    if timestep >= T_DIV:
        dt = timestep - T_DIV
        bd_steps = int(BD_DIV0 + dt * (BD_DIV1 - BD_DIV0) / div_span)
    else:
        bd_steps = BD_STEPS
    return _soft_harmonic_dynamics(bd_steps)


def compute_cell_boundary(timestep: int, dna_monomers_path: str) -> CellBoundary:
    """Return boundary geometry for one biological minute."""
    div_span = T_END - T_DIV
    sphere_radius = R0 + R_FUDGE + DR_DT * timestep
    sphere_height = 0.0
    separation_axis = np.array([1.0, 0.0, 0.0])

    if timestep >= T_DIV:
        dt = timestep - T_DIV
        sphere_radius = R_DIV + R_AMP * (timestep - T_END) ** 2 / div_span**2
        sphere_height = H0 + np.sqrt(dt) * H_SPAN / np.sqrt(div_span)
        separation_axis = division_axis_from_dna(dna_monomers_path)

    load_boundary = (
        _overlapping_spheres_boundary(sphere_height, sphere_radius, separation_axis)
        if timestep >= T_DIV
        else _spherical_boundary(sphere_radius)
    )

    return CellBoundary(
        sphere_radius=sphere_radius,
        sphere_height=sphere_height,
        separation_axis=separation_axis,
        load_boundary=load_boundary,
    )


def directives_for_minute(run_name: str, timestep: int) -> PerMinuteDirectives:
    """Replication, looping, and dynamics settings for one biological minute."""
    run_dynamics = compute_run_dynamics(timestep)

    if timestep == 0:
        return PerMinuteDirectives(
            input_state=f"input_state:{template_dir}rep_state_initial.txt",
            load_loops="",
            equilibrate_loops=LOOP_EQUIL,
            append_string="noappend,first",
            run_dynamics=run_dynamics,
        )

    return PerMinuteDirectives(
        input_state=f"input_state:{output_dir}rep_states/rep_state_{run_name}_{timestep}.txt",
        load_loops=f"load_loops:{output_dir}loops/loops_{run_name}_{timestep}.txt",
        equilibrate_loops="",
        append_string="append,nofirst",
        run_dynamics=run_dynamics,
    )


def create_directives(
    run_name: str, seed: int, timestep: int, smc: SmcParams, v_replication: int
) -> str:
    """Build one biological-minute btree_chromo directive file from template_replicate.inp."""
    dna_monomers_path, ribos_path = coord_paths(run_name, timestep)
    boundary = compute_cell_boundary(timestep, dna_monomers_path)
    minute = directives_for_minute(run_name, timestep)
    batching = compute_replication_batching(
        timestep, v_replication, smc.v_translocation
    )

    add_ribosomes_to_bin(
        ribos_path,
        RIBOS_PER_MIN,
        boundary.sphere_radius,
        boundary.sphere_height,
        boundary.separation_axis,
    )

    with open(template_file) as file:
        directives_content = file.read()

    directives_content = directives_content.format(
        sim_prng_seed=seed,
        base_dir=base_dir,
        run_name=run_name,
        input_state=minute.input_state,
        dna_monomers_path=dna_monomers_path,
        ribos_path=ribos_path,
        output_dir=output_dir,
        timestep=timestep,
        next_timestep=timestep + 1,
        dna_monomers_dir=dna_monomers_dir,
        load_boundary=boundary.load_boundary,
        load_loops=minute.load_loops,
        equilibrate_loops=minute.equilibrate_loops,
        append_string=minute.append_string,
        run_dynamics=minute.run_dynamics,
        repeat_bd=batching.repeat_bd,
        repeat_no_bd=batching.repeat_no_bd,
        transform_bd=batching.transform_bd,
        transform_no_bd=batching.transform_no_bd,
        translocate_bd=batching.translocate_bd,
        translocate_no_bd=batching.translocate_no_bd,
    )

    directives_filename = f"{output_dir}/btree_chromo_directives_{timestep}_{run_name}.inp"
    with open(directives_filename, "w") as f:
        f.write(directives_content)

    return directives_filename


def run_btree_chromo(directives_file: str, env: dict[str, str]) -> None:
    subprocess.run([btree_chromo_executable, directives_file], check=True, env=env)
    print(f"Successfully ran btree_chromo with directives: {directives_file}")


def run_simulation_minute(
    run_name: str,
    seed: int,
    timestep: int,
    env: dict[str, str],
    smc: SmcParams,
    v_replication: int,
) -> None:
    """Create directives, run btree_chromo for one minute, and remove the temp file."""
    directives_file = create_directives(
        run_name, seed, timestep, smc=smc, v_replication=v_replication
    )
    try:
        run_btree_chromo(directives_file, env)
    finally:
        os.remove(directives_file)


def main() -> None:
    args = parse_args(sys.argv)

    if not os.path.isfile(btree_chromo_executable):
        raise FileNotFoundError(
            f"{btree_chromo_executable} not found. "
            "Run build.sh on a login node before submitting the simulation."
        )

    env = build_runtime_env()
    smc = args.smc

    write_loop_params(smc)

    print(f"Running simulation: {args.run_name}")
    print(f"Start time: {args.start_time}, End time: {args.end_time}, Restart: {args.is_restart}")
    print(f"Replication fork speed: {args.v_replication} bp/s")
    print(
        f"SMC parameters: N={smc.n_smc}, v_translocation={smc.v_translocation} bp/s, "
        f"tau_basal={smc.tau_basal} s, tau_stall={smc.tau_stall} s, "
        f"tau_bypass={smc.tau_bypass} s, S={smc.s}"
    )
    print(
        f"  basal_death_prob={smc.basal_death_prob:.6f}, "
        f"stall_death_prob={smc.stall_death_prob:.6f}, "
        f"bypass={smc.bypass:.6f}, knockoff={smc.knockoff:.6f}"
    )
    print(
        f"Replication batching (default): {int(BD_BIO_SECONDS / DEFAULT_SECONDS_PER_BATCH)} BD "
        f"+ {int(NO_BD_BIO_SECONDS / DEFAULT_SECONDS_PER_BATCH)} no-BD batches @ "
        f"{DEFAULT_SECONDS_PER_BATCH:g} s"
    )
    print(
        f"  minute {HIGH_RES_MINUTE} BD section: "
        f"{int(BD_BIO_SECONDS / HIGH_RES_SECONDS_PER_BATCH)} batches @ "
        f"{HIGH_RES_SECONDS_PER_BATCH:g} s"
    )
    default_beads = replicate_amount(args.v_replication, DEFAULT_SECONDS_PER_BATCH)
    print(
        f"  per-batch amounts @ {DEFAULT_SECONDS_PER_BATCH:g} s: "
        f"{default_beads} beads ({transform_directive(args.v_replication, DEFAULT_SECONDS_PER_BATCH)}), "
        f"translocate {translocate_amount(smc.v_translocation, DEFAULT_SECONDS_PER_BATCH)}"
    )

    for timestep in range(args.start_time, args.end_time + 1):
        if timestep == HIGH_RES_MINUTE:
            batching = compute_replication_batching(
                timestep, args.v_replication, smc.v_translocation
            )
            print(
                f"Minute {timestep}: high-res BD — {batching.repeat_bd} batches @ "
                f"{batching.seconds_per_batch_bd:g} s "
                f"({batching.transform_bd}, translocate {batching.translocate_bd})"
            )
        run_simulation_minute(
            args.run_name, args.seed, timestep, env, smc, args.v_replication
        )


if __name__ == "__main__":
    main()
