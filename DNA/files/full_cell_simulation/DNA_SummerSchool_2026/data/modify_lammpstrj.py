#!/usr/bin/env python3
"""
Relabel DNA monomer types in a LAMMPS trajectory for VMD daughter coloring.

Per frame, uses the total DNA monomer count (max c_id_track among type-3 beads)
to assign left/right daughter beads — same geometry as plot_partitioning.py:
  - type 3  : mother chromosome (unreplicated parent beads)
  - type 13 : left daughter
  - type 14 : right daughter (DNA_new)

This replaces the old symmetric window around the terminus (width_per_frame * frame),
which did not match the train-track replication geometry and ignored the extra
trajectory frames dumped during minute 30 (HIGH_RES_MINUTE in run_btree_chromo.py).
"""

from __future__ import annotations

from pathlib import Path

N_PARENT = 54338
MIDPOINT_1BASED = 27169
MONO_TYPE = 3
LEFT_TYPE = 13
RIGHT_TYPE = 14

INPUT_FILE = Path("summerschool.lammpstrj")
OUTPUT_FILE = Path("modified.lammpstrj")


def daughter_bead_ids_1based(n_beads: int) -> tuple[set[int], set[int]]:
    """Return (left, right) daughter genomic ids (1-based c_id_track values)."""
    if n_beads <= N_PARENT:
        return set(), set()

    n_repl = n_beads - N_PARENT
    left_start = (MIDPOINT_1BASED - 1) - (n_repl - 1) // 2
    left_end = (MIDPOINT_1BASED - 1) + n_repl // 2 + 1
    left = {bead + 1 for bead in range(left_start, left_end)}
    right = {bead + 1 for bead in range(N_PARENT, n_beads)}
    return left, right


def modify_lammpstrj(input_file: Path, output_file: Path) -> None:
    with open(input_file) as fin, open(output_file, "w") as fout:
        line = fin.readline()
        while line:
            if not line.startswith("ITEM: TIMESTEP"):
                fout.write(line)
                line = fin.readline()
                continue

            fout.write(line)
            fout.write(fin.readline())  # timestep value
            fout.write(fin.readline())  # ITEM: NUMBER OF ATOMS
            num_atoms_line = fin.readline()
            fout.write(num_atoms_line)
            num_atoms = int(num_atoms_line)
            fout.write(fin.readline())  # ITEM: BOX BOUNDS
            for _ in range(3):
                fout.write(fin.readline())

            header_line = fin.readline()
            fout.write(header_line)
            header_columns = header_line.strip().split()[2:]
            id_index = header_columns.index("c_id_track")
            type_index = header_columns.index("c_type_track")

            atom_lines = [fin.readline() for _ in range(num_atoms)]
            mono_max = 0
            for atom_line in atom_lines:
                fields = atom_line.split()
                if int(fields[type_index]) == MONO_TYPE:
                    mono_max = max(mono_max, int(fields[id_index]))

            left_ids, right_ids = daughter_bead_ids_1based(mono_max)

            for atom_line in atom_lines:
                fields = atom_line.split()
                c_id = int(fields[id_index])
                c_type = int(fields[type_index])
                if c_type == MONO_TYPE:
                    if c_id in left_ids:
                        fields[type_index] = str(LEFT_TYPE)
                    elif c_id in right_ids:
                        fields[type_index] = str(RIGHT_TYPE)
                fout.write(" ".join(fields) + "\n")

            line = fin.readline()


if __name__ == "__main__":
    modify_lammpstrj(INPUT_FILE, OUTPUT_FILE)
    print(f"Wrote {OUTPUT_FILE}")
