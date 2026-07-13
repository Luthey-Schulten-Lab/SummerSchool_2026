#!/usr/bin/env python3
"""
Find small (recently born) SMC loops near a given cell-cycle minute.

Minute identity is inferred from replication-fork positions in the concatenated
loops file (loops_summerschool.txt). Only loops that land on an actual
lammpstrj dump frame are reported.

Usage (from full_cell_simulation/):
    python3 find_small_loops.py
    python3 find_small_loops.py --max-length 50 --minute 30
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

# Syn3A bead constants (same as plot_partitioning.py / run_btree_chromo.py)
N_PARENT = 54338
MIDPOINT = 27169

V_REPLICATION_BPS = 100  # bp/s
BP_PER_BEAD = 10
SECONDS_PER_MINUTE = 60
BEADS_PER_MIN_PER_FORK = V_REPLICATION_BPS // BP_PER_BEAD * SECONDS_PER_MINUTE  # 600

# Matches run_btree_chromo.py high-res + dump schedule:
#   BD dump interval = BD_STEPS * BD_OUT_FACTOR = 20000 * 2 = 40000 steps
#   → one lammpstrj frame every 2 BD run_dynamics calls
BD_OUT_FACTOR = 2
BD_DUMP_EVERY = BD_OUT_FACTOR  # BD batches per lammpstrj frame

DEFAULT_SECONDS_PER_BATCH = 2.0
HIGH_RES_SECONDS_PER_BATCH = 0.4
BD_BIO_SECONDS_EARLY = 20.0  # minutes 0–59 BD budget

HIGH_RES_MINUTE = 30
DEFAULT_MAX_LENGTH = 100
DEFAULT_LOOPS = Path("DNA_SummerSchool_2026/data/loops/loops_summerschool.txt")


@dataclass(frozen=True)
class LoopSnapshot:
    index: int
    left_fork: int
    right_fork: int
    loops: list[tuple[int, int]]


@dataclass(frozen=True)
class SmallLoopHit:
    start: int
    end: int
    length: int
    frame: int


def beads_replicated_per_batch(seconds_per_batch: float) -> int:
    return int(round(V_REPLICATION_BPS / BP_PER_BEAD * seconds_per_batch))


def expected_left_fork(minute: int) -> int:
    return MIDPOINT - BEADS_PER_MIN_PER_FORK * minute


def minute_fork_window(minute: int) -> tuple[int, int]:
    """Inclusive left-fork range for [start of minute, end of minute]."""
    left_start = expected_left_fork(minute)
    left_end = expected_left_fork(minute + 1)
    return left_end, left_start


def loop_length(start: int, end: int) -> int:
    return abs(end - start)


def parse_loops_file(path: Path) -> list[LoopSnapshot]:
    snapshots: list[LoopSnapshot] = []
    current_loops: list[tuple[int, int]] = []
    left_fork: int | None = None
    right_fork: int | None = None

    def flush() -> None:
        nonlocal current_loops, left_fork, right_fork
        if left_fork is None or right_fork is None:
            current_loops = []
            left_fork = right_fork = None
            return
        snapshots.append(
            LoopSnapshot(len(snapshots), left_fork, right_fork, list(current_loops))
        )
        current_loops = []
        left_fork = right_fork = None

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("Number of loops:"):
                flush()
            elif line.startswith("Replication forks:"):
                vals = line.split(":", 1)[1].strip().split(",")
                left_fork = int(vals[0].strip())
                right_fork = int(vals[1].strip())
            elif line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        current_loops.append((int(parts[0]), int(parts[1])))
                    except ValueError:
                        continue
    flush()
    return snapshots


def find_minute_snapshots(
    snapshots: list[LoopSnapshot], minute: int
) -> list[LoopSnapshot]:
    """Last contiguous block whose left fork sits in the minute's window."""
    low, high = minute_fork_window(minute)
    in_window = [
        i for i, snap in enumerate(snapshots) if low <= snap.left_fork <= high
    ]
    if not in_window:
        raise ValueError(
            f"No snapshots with left fork in [{low}, {high}] "
            f"(expected start-of-minute-{minute} fork ≈ {expected_left_fork(minute)})"
        )

    end = in_window[-1]
    start = end
    while start > 0 and low <= snapshots[start - 1].left_fork <= high:
        start -= 1

    block = snapshots[start : end + 1]
    trimmed: list[LoopSnapshot] = [block[0]]
    for snap in block[1:]:
        if snap.left_fork <= trimmed[-1].left_fork:
            trimmed.append(snap)
        else:
            trimmed = [snap]
    return trimmed


def classify_bd_section(
    minute_snaps: list[LoopSnapshot], minute: int
) -> list[LoopSnapshot]:
    """High-res BD snapshots (fork step ≈ 4 beads), skipping the minute-start write."""
    high_res_step = beads_replicated_per_batch(HIGH_RES_SECONDS_PER_BATCH)
    normal_step = beads_replicated_per_batch(DEFAULT_SECONDS_PER_BATCH)
    left_start = expected_left_fork(minute)

    bd: list[LoopSnapshot] = []
    for i, snap in enumerate(minute_snaps):
        if i == 0:
            continue
        step = minute_snaps[i - 1].left_fork - snap.left_fork
        if abs(step - high_res_step) <= 1 and snap.left_fork < left_start:
            bd.append(snap)
        elif bd and abs(step - normal_step) <= 1:
            break
        elif bd and step > high_res_step + 1:
            break

    if not bd:
        raise ValueError(
            f"Could not locate high-res BD section for minute {minute} "
            f"(expected fork steps of {high_res_step} beads)."
        )
    return bd


def first_high_res_lammpstrj_frame(minute: int = HIGH_RES_MINUTE) -> int:
    """0-based lammpstrj index of the first dump in the high-res minute (30 → 151)."""
    repeat_bd_normal = int(round(BD_BIO_SECONDS_EARLY / DEFAULT_SECONDS_PER_BATCH))
    dumps_per_normal_minute = repeat_bd_normal // BD_DUMP_EVERY
    return 1 + dumps_per_normal_minute * minute


def frame_range(n_bd_snaps: int, high_res_frame0: int) -> tuple[int, int]:
    n_dumps = (n_bd_snaps + BD_DUMP_EVERY - 1) // BD_DUMP_EVERY
    return high_res_frame0, high_res_frame0 + max(n_dumps - 1, 0)


def lammpstrj_frame(bd_batch: int, high_res_frame0: int) -> int | None:
    """Return the dump frame for this BD batch, or None if it is not a dump step."""
    completed = bd_batch + 1
    if completed % BD_DUMP_EVERY != 0:
        return None
    return high_res_frame0 + completed // BD_DUMP_EVERY - 1


def find_small_loops(
    bd_snaps: list[LoopSnapshot],
    max_length: int,
    high_res_frame0: int,
) -> list[SmallLoopHit]:
    hits: list[SmallLoopHit] = []
    for bd_batch, snap in enumerate(bd_snaps):
        frame = lammpstrj_frame(bd_batch, high_res_frame0)
        if frame is None:
            continue
        for start, end in snap.loops:
            length = loop_length(start, end)
            if 0 < length < max_length:
                hits.append(SmallLoopHit(start, end, length, frame))
    hits.sort(key=lambda h: (h.frame, h.length, h.start, h.end))
    return hits


def _match_cost(a: SmallLoopHit, b: SmallLoopHit) -> int:
    """Endpoint displacement between two observations of a candidate same loop."""
    return abs(a.start - b.start) + abs(a.end - b.end)


def track_loops(
    hits: list[SmallLoopHit], *, max_match_cost: int = 80
) -> list[list[SmallLoopHit]]:
    """
    Group hits into per-loop trajectories across consecutive dump frames.

    A hit is appended to an open track when its endpoints are close to that
    track's latest observation (cost = |Δstart| + |Δend| ≤ max_match_cost)
    and no closer unmatched hit claims it.
    """
    if not hits:
        return []

    by_frame: dict[int, list[SmallLoopHit]] = {}
    for hit in hits:
        by_frame.setdefault(hit.frame, []).append(hit)

    frames = sorted(by_frame)
    tracks: list[list[SmallLoopHit]] = [[h] for h in by_frame[frames[0]]]

    for frame in frames[1:]:
        open_tracks = {i: tracks[i][-1] for i in range(len(tracks))}
        candidates: list[tuple[int, int, int]] = []  # cost, track_i, hit_j
        frame_hits = by_frame[frame]
        for ti, last in open_tracks.items():
            for hj, hit in enumerate(frame_hits):
                cost = _match_cost(last, hit)
                if cost <= max_match_cost:
                    candidates.append((cost, ti, hj))
        candidates.sort()

        used_tracks: set[int] = set()
        used_hits: set[int] = set()
        for cost, ti, hj in candidates:
            if ti in used_tracks or hj in used_hits:
                continue
            tracks[ti].append(frame_hits[hj])
            used_tracks.add(ti)
            used_hits.add(hj)

        for hj, hit in enumerate(frame_hits):
            if hj not in used_hits:
                tracks.append([hit])

    # Newest / smallest births first (best for visualization).
    tracks.sort(key=lambda t: (t[0].frame, t[0].length, t[0].start, t[0].end))
    return tracks


def format_report(
    high_res_frame0: int,
    last_frame: int,
    tracks: list[list[SmallLoopHit]],
) -> str:
    lines = [f"Frame range: {high_res_frame0}-{last_frame}", ""]
    if not tracks:
        lines.append("(none)")
        return "\n".join(lines)

    for i, track in enumerate(tracks, start=1):
        lines.append(f"Loop {i}:")
        lines.append(f"{'start':>7}  {'end':>7}  {'len':>4}  {'frame':>5}")
        lines.append("-" * 30)
        for hit in track:
            lines.append(
                f"{hit.start:7d}  {hit.end:7d}  {hit.length:4d}  {hit.frame:5d}"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find small SMC loops during a high-res cell-cycle minute."
    )
    parser.add_argument(
        "--loops-file",
        type=Path,
        default=DEFAULT_LOOPS,
        help="Concatenated loops dump (default: %(default)s)",
    )
    parser.add_argument(
        "--minute",
        type=int,
        default=HIGH_RES_MINUTE,
        help="Biological minute to search (default: %(default)s)",
    )
    parser.add_argument(
        "--max-length",
        "-N",
        type=int,
        default=DEFAULT_MAX_LENGTH,
        help="Keep loops with abs(end-start) < N (default: %(default)s)",
    )
    parser.add_argument(
        "--lammpstrj-frame0",
        type=int,
        default=None,
        help="Override 0-based index of first high-res lammpstrj frame "
        f"(default: derived as {first_high_res_lammpstrj_frame()})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.loops_file.is_file():
        raise FileNotFoundError(args.loops_file)

    snapshots = parse_loops_file(args.loops_file)
    minute_snaps = find_minute_snapshots(snapshots, args.minute)
    bd_snaps = classify_bd_section(minute_snaps, args.minute)
    high_res_frame0 = (
        args.lammpstrj_frame0
        if args.lammpstrj_frame0 is not None
        else first_high_res_lammpstrj_frame(args.minute)
    )
    first_frame, last_frame = frame_range(len(bd_snaps), high_res_frame0)
    hits = find_small_loops(bd_snaps, args.max_length, high_res_frame0)
    tracks = track_loops(hits)
    print(format_report(first_frame, last_frame, tracks))


if __name__ == "__main__":
    main()
