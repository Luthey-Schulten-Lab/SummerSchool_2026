#!/usr/bin/env python3
"""Generate Minecraft animation schematics from a 4DWCM MinCell.lm trajectory."""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
from mcschematic import Version
from mcschematic_plus import MCSchematicPlus

COLORS = {
    1: "air",  # cytoplasm
    2: "air",  # outer_cytoplasm
    3: "red_stained_glass",  # ribosomes
    4: "red_concrete",  # ribo_centers
    5: "yellow_concrete",  # DNA
    6: "lime_stained_glass",  # membrane
}


def resolve_lm_path(*, group: str | None, run_name: str, user: str) -> str:
    if group == "shared":
        return "/projects/bgvl/SummerSchool_2026/4DWCM/trajectory/Mar31_1/MinCell.lm"
    if group:
        return (
            f"/projects/bgvl/Groups_4DWCM/{group}/4dwcm_run/"
            f"Optimize_4DWCM_Minimal_Cell/Data/{run_name}/MinCell.lm"
        )
    return (
        f"/projects/bgvl/{user}/4dwcm_run/"
        f"Optimize_4DWCM_Minimal_Cell/Data/{run_name}/MinCell.lm"
    )


def generate_schematics(lm_path: str, output_path: str, frame_num: int) -> None:
    if not os.path.isfile(lm_path):
        raise FileNotFoundError(f"MinCell.lm not found: {lm_path}")

    os.makedirs(output_path, exist_ok=True)

    print(f"Reading: {lm_path}")
    with h5py.File(lm_path, "r") as traj:
        sim = traj["Simulations"]["0000001"]
        max_time = int(sim["LatticeTimes"][-1])
        jump = max(1, int(max_time / frame_num))
        print(f"max_time={max_time}, frames={frame_num}, jump={jump}")

        for k in range(frame_num + 1):
            frame_key = f"000000{str(k * jump).zfill(4)}"
            voxels = np.array(sim["Sites"][frame_key])
            mask_voxels = voxels[..., None] == np.arange(np.max(voxels) + 1)
            mask_voxels = mask_voxels[..., 1:]  # drop extracellular (air)

            schem = MCSchematicPlus()
            for i, component in enumerate(COLORS.keys()):
                vol = mask_voxels[..., i]
                for j in range(i):
                    vol = vol & ~mask_voxels[..., j]
                schem.placeVolume(vol, COLORS[component])

            out_file = os.path.join(output_path, f"lm_sim{k}.nbt")
            schem.saveNBT(out_file, Version.JE_1_20_1, shifted=False)
            print(f"  wrote {out_file}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build Minecraft animation schematics from MinCell.lm"
    )
    parser.add_argument("--lm-path", help="Full path to MinCell.lm")
    parser.add_argument(
        "--group",
        choices=["Group1", "Group2", "Group3", "shared"],
        help="Read from Groups_4DWCM/<group> or shared Mar31_1 demo",
    )
    parser.add_argument(
        "--run-name",
        default="4dwcm_7200",
        help="Data subfolder name (default: 4dwcm_7200)",
    )
    parser.add_argument(
        "--user",
        default=os.environ.get("USER", ""),
        help="NCSA username for personal 4dwcm_run (default: $USER)",
    )
    parser.add_argument(
        "--output-path",
        required=True,
        help="Directory for lm_sim*.nbt files",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=30,
        help="Number of animation frames (default: 30)",
    )
    args = parser.parse_args()

    if args.lm_path:
        lm_path = args.lm_path
    elif args.group:
        lm_path = resolve_lm_path(group=args.group, run_name=args.run_name, user=args.user)
    elif args.user:
        lm_path = resolve_lm_path(group=None, run_name=args.run_name, user=args.user)
    else:
        parser.error("Set --lm-path, --group, or --user")

    try:
        generate_schematics(lm_path, args.output_path, args.frames)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Done. Schematics in: {args.output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
