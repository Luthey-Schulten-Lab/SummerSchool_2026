#!/usr/bin/env python3
"""
Compute and plot daughter-chromosome partitioning vs time.

Usage (from full_cell_simulation/):
    python3 plot_partitioning.py
    python3 plot_partitioning.py --from-bins
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

N_PARENT = 54338
MIDPOINT_1BASED = 27169

TRJ = Path("DNA_SummerSchool_2026/data/summerschool.lammpstrj")
COORDS_DIR = Path("DNA_SummerSchool_2026/data/coords")
RUN_NAME = "summerschool"
OUT_TXT = Path("partitioning_summerschool.txt")
OUT_PNG = Path("partitioning_summerschool.png")
END_MINUTE = 90.0
CUTOFF = 110.0
TITLE = "Daughter chromosome partitioning vs time"
YLIM = (0.0, 1.0)


def sanitize_partitioning(values: np.ndarray) -> np.ndarray:
    """Treat NaN as 0 (e.g. minute 0 before replication has no daughter beads)."""
    return np.nan_to_num(values, nan=0.0)


def daughter_partitioning(coords: np.ndarray, cutoff: float = CUTOFF) -> float:
    """Return partitioning metric in [0, 1] for one frame; 0 before replication."""
    n = len(coords)
    if n <= N_PARENT:
        return 0.0

    n_repl = n - N_PARENT
    right_daughter = np.arange(N_PARENT, n)
    left_start = (MIDPOINT_1BASED - 1) - (n_repl - 1) // 2
    left_end = (MIDPOINT_1BASED - 1) + n_repl // 2 + 1
    left_daughter = np.arange(left_start, left_end)

    label = np.zeros(n, dtype=np.int8)
    label[left_daughter] = 1
    label[right_daughter] = 2

    marked: set[int] = set()
    tree = cKDTree(coords)
    for i in left_daughter:
        for j in tree.query_ball_point(coords[i], cutoff):
            if i != j and label[j] == 2:
                marked.add(int(i))
                break
    for i in right_daughter:
        for j in tree.query_ball_point(coords[i], cutoff):
            if i != j and label[j] == 1:
                marked.add(int(i))
                break

    return 1.0 - len(marked) / (2 * n_repl)


def read_bin(path: Path) -> np.ndarray:
    data = np.fromfile(path, dtype=np.float64)
    if data.size % 3 != 0:
        raise ValueError(f"{path} size not divisible by 3")
    return data.reshape((data.size // 3, 3))


def discover_bin_files(coords_dir: Path, run_name: str) -> list[tuple[int, Path]]:
    pattern = re.compile(rf"dna_{re.escape(run_name)}_(\d+)\.bin$")
    files: list[tuple[int, Path]] = []
    for path in coords_dir.glob(f"dna_{run_name}_*.bin"):
        match = pattern.match(path.name)
        if match:
            files.append((int(match.group(1)), path))
    return sorted(files, key=lambda item: item[0])


def iter_lammpstrj_frames(path: Path):
    frame_index = 0
    with open(path, "r") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if "ITEM: TIMESTEP" not in line:
                continue
            handle.readline()

            if "ITEM: NUMBER OF ATOMS" not in handle.readline():
                break
            n_atoms = int(handle.readline().strip())

            if "ITEM: BOX" not in handle.readline():
                break
            for _ in range(3):
                handle.readline()

            header_line = handle.readline()
            if "ITEM: ATOMS" not in header_line:
                break
            header = header_line.strip().split()[2:]
            idx = {name: header.index(name) for name in ("x", "y", "z", "c_id_track", "c_type_track")}

            rows: list[tuple[int, float, float, float]] = []
            for _ in range(n_atoms):
                atom_line = handle.readline()
                if not atom_line:
                    break
                parts = atom_line.split()
                if int(parts[idx["c_type_track"]]) < 3:
                    continue
                rows.append(
                    (
                        int(parts[idx["c_id_track"]]),
                        float(parts[idx["x"]]),
                        float(parts[idx["y"]]),
                        float(parts[idx["z"]]),
                    )
                )

            if rows:
                rows.sort(key=lambda row: row[0])
                max_id = rows[-1][0]
                coords = np.zeros((max_id, 3), dtype=np.float64)
                for c_id, x, y, z in rows:
                    coords[c_id - 1] = (x, y, z)
            else:
                coords = np.empty((0, 3), dtype=np.float64)

            yield frame_index, coords
            frame_index += 1


def compute_partitioning_from_bins(coords_dir: Path, run_name: str) -> tuple[np.ndarray, np.ndarray]:
    bin_files = discover_bin_files(coords_dir, run_name)
    if not bin_files:
        raise ValueError(
            f"No dna_{run_name}_*.bin files found in {coords_dir} "
            f"(simulation may still be running)"
        )

    minutes: list[int] = []
    partitioning: list[float] = []
    for i, (minute, path) in enumerate(bin_files):
        partitioning.append(daughter_partitioning(read_bin(path)))
        minutes.append(minute)
        if (i + 1) % 25 == 0:
            print(f"  processed minute {minute} ({i + 1}/{len(bin_files)})")

    return np.array(minutes, dtype=int), np.array(partitioning, dtype=np.float64)


def compute_partitioning_from_lammpstrj(trj_path: Path) -> tuple[np.ndarray, np.ndarray]:
    frames: list[int] = []
    partitioning: list[float] = []

    for i, (frame_idx, coords) in enumerate(iter_lammpstrj_frames(trj_path)):
        frames.append(frame_idx)
        partitioning.append(daughter_partitioning(coords))
        if (i + 1) % 25 == 0:
            print(f"  processed {i + 1} frames")

    if not frames:
        raise ValueError(f"No frames read from {trj_path}")

    return np.array(frames, dtype=int), np.array(partitioning, dtype=np.float64)


def save_partitioning_txt(
    time_values: np.ndarray,
    partitioning: np.ndarray,
    source: str,
    source_path: Path,
    time_column: str,
) -> None:
    with open(OUT_TXT, "w") as handle:
        handle.write(f"# partitioning time series from {source}\n")
        handle.write(f"# cutoff={CUTOFF}\n")
        handle.write(f"# source={source_path.resolve()}\n")
        handle.write(f"# time_column={time_column}\n")
        handle.write(f"#\n# {time_column}\tpartitioning\n")
        for time_val, value in zip(time_values, partitioning):
            handle.write(f"{int(time_val)}\t{value}\n")


def load_partitioning_txt() -> tuple[np.ndarray, np.ndarray, str]:
    time_values, values = [], []
    time_column = "frame"
    with open(OUT_TXT) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("# time_column="):
                time_column = line.split("=", 1)[1]
                continue
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                time_values.append(int(parts[0]))
                values.append(float(parts[1]))
    return np.array(time_values, dtype=int), np.array(values, dtype=np.float64), time_column


def plot_partitioning(
    time_values: np.ndarray,
    partitioning: np.ndarray,
    time_is_minutes: bool,
) -> None:
    order = np.argsort(time_values)
    time_values = time_values[order]
    partitioning = partitioning[order]

    if time_is_minutes:
        time_min = time_values.astype(float)
    else:
        max_frame = int(time_values.max()) if len(time_values) else 0
        time_min = (
            time_values.astype(float) * END_MINUTE / max_frame
            if max_frame > 0
            else time_values.astype(float)
        )

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(time_min, partitioning, color="C0", linewidth=1.5, label=RUN_NAME)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Partitioning")
    ax.set_title(TITLE)
    ax.set_xlim(0, END_MINUTE if END_MINUTE > 0 else (time_min.max() if len(time_min) else 1))
    ax.set_ylim(*YLIM)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {OUT_PNG}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute and plot daughter-chromosome partitioning vs time."
    )
    parser.add_argument(
        "--from-bins",
        action="store_true",
        help="Use per-minute dna_<run_name>_<minute>.bin files instead of .lammpstrj",
    )
    use_bins = parser.parse_args().from_bins

    if use_bins:
        if not COORDS_DIR.is_dir():
            print(f"ERROR: coords directory not found: {COORDS_DIR}")
            return 1
        source_path = COORDS_DIR
        source_label = "bins"
        time_column = "minute"
    else:
        if not TRJ.exists():
            print(f"ERROR: trajectory not found: {TRJ}")
            return 1
        source_path = TRJ
        source_label = "lammpstrj"
        time_column = "frame"

    time_values: np.ndarray | None = None
    partitioning: np.ndarray | None = None
    time_is_minutes = use_bins

    if OUT_TXT.exists():
        time_values, partitioning, cached_column = load_partitioning_txt()
        time_is_minutes = cached_column == "minute"
        if time_is_minutes != use_bins:
            print(
                f"WARNING: cached data is {cached_column}-based but "
                f"{'--from-bins' if use_bins else 'lammpstrj'} was requested; recomputing."
            )
            time_values = None

    if time_values is None:
        if use_bins:
            print(f"Computing partitioning from {COORDS_DIR}/dna_{RUN_NAME}_*.bin ...")
            time_values, partitioning = compute_partitioning_from_bins(COORDS_DIR, RUN_NAME)
        else:
            print(f"Computing partitioning from {TRJ} ...")
            time_values, partitioning = compute_partitioning_from_lammpstrj(TRJ)
        partitioning = sanitize_partitioning(partitioning)
        save_partitioning_txt(time_values, partitioning, source_label, source_path, time_column)
        time_is_minutes = use_bins
        print(f"Saved partitioning trace to {OUT_TXT}")
    else:
        print(f"Loading existing partitioning data from {OUT_TXT}")
        partitioning = sanitize_partitioning(partitioning)

    plot_partitioning(time_values, partitioning, time_is_minutes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
