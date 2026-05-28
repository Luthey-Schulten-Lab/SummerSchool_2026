"""Build three Jupyter notebooks for the Bentopy / Martini summer-school
tutorial and drop each one into its own subdirectory of the SummerSchool
repo's Martini module:

    Martini/tutorial_1/tutorial_1_basic_packing.ipynb
    Martini/tutorial_2/tutorial_2_membrane_packing.ipynb
    Martini/tutorial_3/tutorial_3_multi_compartment.ipynb

The notebooks are tutorial-relative: every cell uses the working directory
(`tutorial_<N>/`) where `structures/`, `topology/`, `mdp_files/` are made
available via symlinks to `../data/`.

Run:    python scripts/build_notebooks.py
"""
from __future__ import annotations
import json
import pathlib
import uuid

HERE = pathlib.Path(__file__).resolve().parent      # Martini/scripts/
MARTINI = HERE.parent                                # Martini/

KSPEC = {
    "kernelspec": {
        "display_name": "Python 3 (ipykernel)",
        "language": "python",
        "name": "python3",
    },
    "language_info": {"name": "python", "version": "3.12"},
}


def _new_id() -> str:
    return uuid.uuid4().hex[:12]


def md(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": _new_id(),
        "metadata": {},
        "source": text.rstrip() + "\n",
    }


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "id": _new_id(),
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.rstrip() + "\n",
    }


def write_notebook(subdir: str, name: str, cells: list[dict]) -> pathlib.Path:
    nb = {
        "cells": cells,
        "metadata": KSPEC,
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out_dir = MARTINI / subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / name
    out.write_text(json.dumps(nb, indent=1))
    return out


# ---------------------------------------------------------------------------
# Common preamble used in each notebook
# ---------------------------------------------------------------------------

PREAMBLE_CODE = """\
import os, pathlib, shutil, sys

# Make the bentopy CLI tools available no matter which Python kernel the
# user picks. The shared venv lives in the workshop's allocation.
VENV_BIN = "/projects/bgvl/alfiaparvez/bentopy_tutorial/.venv/bin"
if VENV_BIN not in os.environ["PATH"]:
    os.environ["PATH"] = VENV_BIN + ":" + os.environ["PATH"]

# This notebook is meant to be opened from inside its tutorial folder,
# e.g.  /projects/bgvl/$USER/SummerSchool_2026/Martini/{0}/ .  Outputs of
# every cell below (.bent, .json, .gro, .top, ...) land in this folder.
NB_DIR = pathlib.Path.cwd()

missing = [d for d in ("structures", "topology", "mdp_files")
           if not (NB_DIR / d).exists()]
if missing:
    raise FileNotFoundError(
        f"Missing input directories in {{NB_DIR}}: {{missing}}.\\n"
        f"Run  bash Martini/copy.sh  first, then open the notebook from "
        f"/projects/bgvl/$USER/SummerSchool_2026/Martini/{0}/."
    )

print("Working dir :", NB_DIR)
print("bentopy-pack:", shutil.which("bentopy-pack"))
"""

VIS_BLOCK = """\
## 5. Visualize

Open the solvated system with VMD (run this in a Delta OnDemand desktop
session, not in the notebook itself):

```bash
export PATH=/projects/bgvl/alfiaparvez/software/vmd/bin:$PATH
cd /projects/bgvl/$USER/SummerSchool_2026/Martini/{subdir}
vmd -e vis_t{n}.tcl solvated_system.gro
```

The QuickSurf rendering script `vis_t{n}.tcl` is in this folder already,
so VMD will pick it up directly.
"""


# ===========================================================================
# Tutorial 1
# ===========================================================================

t1_cells: list[dict] = []

t1_cells.append(md("""\
# Bentopy Tutorial 1 — Basic Packing in Empty Space

This notebook follows the [Bentopy tutorial](https://cgmartini.nl/docs/tutorials/Martini3/Bentopy/)
section *Tutorial 1: Basic Packing in Empty Space*.

Workflow:

1. Write a `.bent` recipe describing the system.
2. Run **`bentopy-pack`** to pack 650 lysozymes into a 40 × 40 × 40 nm³ box.
3. Run **`bentopy-render`** to convert the placement list into a `.gro` file
   plus a GROMACS topology.
4. Run **`bentopy-solvate`** to add water and 0.15 M NaCl.
5. Inspect the final solvated system with VMD.
"""))

t1_cells.append(md("""\
## 0. Verify the workspace

Before opening this notebook you should have already run

```bash
bash /projects/bgvl/alfiaparvez/bentopy_tutorial/notebooks/setup_martini.sh
```

That script created `/projects/bgvl/$USER/Martini/tutorial_1/` with a
copy of this notebook, a copy of `structures/`, `topology/`, `mdp_files/`,
and the visualization script `vis_t1.tcl`.  When you opened *this* copy
of the notebook, JupyterLab set the kernel's working directory to that
folder, which is where every output below will land.

The next cell checks the layout and ensures `bentopy-pack` is on `$PATH`.
"""))

t1_cells.append(code(PREAMBLE_CODE.format("tutorial_1")))

t1_cells.append(md("""\
## 1. Create the packing configuration (`simple_packing.bent`)

The `.bent` file is a small recipe describing the box, the input topologies,
the compartments and which segments to pack into them.
"""))

t1_cells.append(code("""\
bent = '''[ general ]
title "Proteins in a box"
seed 0

[ space ]
# All dimensions in bentopy are given in nanometers.
dimensions 40, 40, 40
resolution 0.5

[ includes ]
"topology/martini_v3.0.0.itp"
"topology/martini_v3.0.0_ions_v1.itp"
"topology/martini_v3.0.0_solvents_v1.itp"
"topology/lysozyme.itp"

[ compartments ]
system is all

[ segments ]
LYZ 650 from "structures/lysozyme.pdb" in system
'''
pathlib.Path("simple_packing.bent").write_text(bent)
print(bent)
"""))

t1_cells.append(md("""\
## 2. Pack the system

`bentopy-pack` reads the recipe and writes a *placement list* — a
lightweight JSON describing where each molecule was placed. It does **not**
yet write atomic coordinates.
"""))

t1_cells.append(code("!bentopy-pack simple_packing.bent placements.json"))

t1_cells.append(md("""\
## 3. Render the placements to a coordinate file

`bentopy-render` expands the placement list into a real `.gro` coordinate
file together with a `.top` topology suitable for GROMACS.
"""))

t1_cells.append(code(
    "!bentopy-render placements.json system.gro -t topol.top\n"
    "print('--- topol.top ---')\n"
    "print(pathlib.Path('topol.top').read_text())"
))

t1_cells.append(md("""\
## 4. Solvation

We add water plus 150 mM NaCl, and let `bentopy-solvate` compute the extra
counter-ions needed to make the system electrically neutral.
"""))

t1_cells.append(code(
    "!bentopy-solvate -i system.gro -o solvated_system.gro \\\n"
    "    -s NA:0.15M -s CL:0.15M --charge neutral -t topol.top\n"
    "print('--- final topol.top ---')\n"
    "print(pathlib.Path('topol.top').read_text())"
))

t1_cells.append(md(VIS_BLOCK.format(subdir="tutorial_1", n=1)))

t1_cells.append(md("""\
### Expected result (Figure 4 of the tutorial)

* 650 lysozymes uniformly distributed throughout the box.
* No overlaps between proteins.
* Cytoplasmic densities; the `bentopy-pack` summary should report **100.0 %**
  of the requested copies placed.
"""))

write_notebook("tutorial_1", "tutorial_1_basic_packing.ipynb", t1_cells)


# ===========================================================================
# Tutorial 2
# ===========================================================================

t2_cells: list[dict] = []

t2_cells.append(md("""\
# Bentopy Tutorial 2 — Packing Around an Existing Structure

This notebook follows the [Bentopy tutorial](https://cgmartini.nl/docs/tutorials/Martini3/Bentopy/)
section *Tutorial 2: Packing Around Existing Structures*.

We add a single POPC bilayer as **excluded volume**, place 300 lysozymes
in the bulk solvent, and 100 ubiquitins **within 5 nm of the membrane
surface** using a proximity rule.

New tools introduced:

* **`bentopy-mask`** — convert a structure into a voxel mask
* **`bentopy-merge`** — concatenate `.gro` files
* The `combines` and `around` rules in `.bent` files
"""))

t2_cells.append(md("""\
## 0. Verify the workspace

This notebook is meant to be opened from
`/projects/bgvl/$USER/Martini/tutorial_2/` (created earlier by
`setup_martini.sh`).  All outputs of the cells below land in that folder.
"""))

t2_cells.append(code(PREAMBLE_CODE.format("tutorial_2")))

t2_cells.append(md("""\
## 1. Inspect the compartments identified in the membrane

`bentopy-mask --visualize-labels` writes a `.gro` file with one bead per
voxel, named after its compartment. Negative labels (`-1`, `-2`, …) are
the open / solvent regions, positive labels (`1`, `2`, …) the excluded
solid structures.
"""))

t2_cells.append(code(
    "!bentopy-mask structures/membrane.gro -b labels.gro\n"
    "print('  labels.gro size:', pathlib.Path('labels.gro').stat().st_size, 'bytes')"
))

t2_cells.append(md("""\
You can open `labels.gro` in VMD and use the selections

* `name "-1"` → the outside solvent region (available for packing)
* `name 1`    → the membrane region (excluded)
"""))

t2_cells.append(md("""\
## 2. Create the membrane mask

We turn label `1` (the membrane voxels) into a `.npz` mask that the next
step will load.
"""))

t2_cells.append(code("!bentopy-mask structures/membrane.gro -l 1:membrane_mask.npz"))

t2_cells.append(md("""\
## 3. Write the new recipe `membrane_packing.bent`

Three new ideas:

* `membrane from "membrane_mask.npz"` — load a compartment from a mask file.
* `solvent combines not membrane` — boolean combination of compartments.
* `close-to-membrane around 5 of membrane` — proximity rule (5 nm shell).

The `:lyz` / `:ubq` suffixes are *tags* that override the residue names in
the output `.gro`, which makes selecting them later in VMD trivial.
"""))

t2_cells.append(code("""\
bent = '''[ general ]
title "Proteins around a membrane"
seed 0

[ space ]
dimensions 40, 40, 40
resolution 0.5

[ includes ]
"topology/martini_v3.0.0.itp"
"topology/martini_v3.0.0_ions_v1.itp"
"topology/martini_v3.0.0_solvents_v1.itp"
"topology/martini_v3.0.0_phospholipids_v1.itp"
"topology/lysozyme.itp"
"topology/ubiquitin.itp"

[ compartments ]
membrane           from "membrane_mask.npz"
solvent            combines not membrane
close-to-membrane  around 5 of membrane

[ segments ]
LYZ:lyz 300 from "structures/lysozyme.pdb"  in solvent
UBQ:ubq 100 from "structures/ubiquitin.pdb" in close-to-membrane
'''
pathlib.Path("membrane_packing.bent").write_text(bent)
print(bent)
"""))

t2_cells.append(md("""\
## 4. Pack the system

Bentopy will pack the larger structure first (lysozyme is bigger than
ubiquitin under the default *moment* heuristic).
"""))

t2_cells.append(code("!bentopy-pack membrane_packing.bent placements.json"))

t2_cells.append(md("""\
## 5. Render and merge with the original membrane

`bentopy-render` writes only the packed proteins. `bentopy-merge`
concatenates the protein `.gro` with the original `membrane.gro` to
produce the full system. The lipid count is then appended to `topol.top`.
"""))

t2_cells.append(code(
    "!bentopy-render placements.json packed_proteins.gro -t topol.top\n"
    "!bentopy-merge packed_proteins.gro structures/membrane.gro -o system.gro\n"
    "with open('topol.top', 'a') as fh:\n"
    "    fh.write('POPC    5408\\n')\n"
    "print(pathlib.Path('topol.top').read_text())"
))

t2_cells.append(md("## 6. Solvation"))

t2_cells.append(code(
    "!bentopy-solvate -i system.gro -o solvated_system.gro -t topol.top \\\n"
    "    -s NA:0.15M -s CL:0.15M --charge neutral\n"
    "print('--- final topol.top ---')\n"
    "print(pathlib.Path('topol.top').read_text())"
))

t2_cells.append(md(VIS_BLOCK.format(subdir="tutorial_2", n=2)))

t2_cells.append(md("""\
### Expected result (Figure 5 of the tutorial)

* The membrane appears as an excluded slab of POPC.
* Lysozymes (`resname lyz`) are spread through the aqueous regions.
* Ubiquitins (`resname ubq`) cluster within 5 nm of the membrane surface.
* No protein overlaps with the membrane.
"""))

write_notebook("tutorial_2", "tutorial_2_membrane_packing.ipynb", t2_cells)


# ===========================================================================
# Tutorial 3
# ===========================================================================

t3_cells: list[dict] = []

t3_cells.append(md("""\
# Bentopy Tutorial 3 — Multi-Compartment Systems

This notebook follows the [Bentopy tutorial](https://cgmartini.nl/docs/tutorials/Martini3/Bentopy/)
section *Tutorial 3: Multi-Compartment Systems with Placement Rules*.

We use a **double POPC membrane**, which creates two distinct open
compartments **A** and **B** (separated by the periodic membrane bilayers).
We place

* 200 lysozymes only in compartment **A**, and
* 100 ubiquitins only in compartment **B** *and* close to the membrane.

This combines several earlier ideas: multiple masks, named compartments,
boolean combinations, and proximity rules — all driven from a single
`.bent` recipe.
"""))

t3_cells.append(md("""\
## 0. Verify the workspace

This notebook is meant to be opened from
`/projects/bgvl/$USER/Martini/tutorial_3/` (created earlier by
`setup_martini.sh`).  All outputs of the cells below land in that folder.
"""))

t3_cells.append(code(PREAMBLE_CODE.format("tutorial_3")))

t3_cells.append(md("""\
## 1. Inspect the compartments of the double membrane

The containment algorithm should now find **four** compartments: two
inter-membrane spaces (`-1`, `-2`) and the two bilayers (`1`, `2`).
"""))

t3_cells.append(code(
    "!bentopy-mask structures/double_membrane.gro -b compartment_labels.gro"
))

t3_cells.append(md("""\
Useful VMD selections for `compartment_labels.gro`:

* Compartment **A**: `name "-1"`
* Compartment **B**: `name "-2"`
* Membranes (excluded): `name 1 2`
"""))

t3_cells.append(md("## 2. Create one mask per compartment"))

t3_cells.append(code(
    "!bentopy-mask structures/double_membrane.gro \\\n"
    "    -l  -1:A_mask.npz                   \\\n"
    "    -l  -2:B_mask.npz                   \\\n"
    "    -l 1,2:membrane_mask.npz"
))

t3_cells.append(md("""\
## 3. Write the recipe `compartment_packing.bent`

Note the boolean combination

```
B-close-to-membrane combines membrane-neighborhood and B
```

which restricts ubiquitins to the *intersection* of compartment B and the
4-nm shell around the membranes.
"""))

t3_cells.append(code("""\
bent = '''[ general ]
title "Proteins in different compartments"
seed 0

[ space ]
dimensions 40, 40, 40
resolution 0.5

[ includes ]
"topology/martini_v3.0.0.itp"
"topology/martini_v3.0.0_ions_v1.itp"
"topology/martini_v3.0.0_solvents_v1.itp"
"topology/martini_v3.0.0_phospholipids_v1.itp"
"topology/lysozyme.itp"
"topology/ubiquitin.itp"

[ compartments ]
membrane                from "membrane_mask.npz"
A                       from "A_mask.npz"
B                       from "B_mask.npz"
membrane-neighborhood   around 4 of membrane
B-close-to-membrane     combines membrane-neighborhood and B

[ segments ]
LYZ:lyz 200 from "structures/lysozyme.pdb"  in A
UBQ:ubq 100 from "structures/ubiquitin.pdb" in B-close-to-membrane
'''
pathlib.Path("compartment_packing.bent").write_text(bent)
print(bent)
"""))

t3_cells.append(md("## 4. Pack, render, merge"))

t3_cells.append(code(
    "!bentopy-pack compartment_packing.bent placements.json\n"
    "!bentopy-render placements.json packed_proteins.gro -t topol.top\n"
    "!bentopy-merge packed_proteins.gro structures/double_membrane.gro -o system.gro\n"
    "with open('topol.top', 'a') as fh:\n"
    "    fh.write('POPC    10816\\n')"
))

t3_cells.append(md("## 5. Solvation"))

t3_cells.append(code(
    "!bentopy-solvate -i system.gro -o solvated_system.gro -t topol.top \\\n"
    "    -s NA:0.15M -s CL:0.15M --charge neutral\n"
    "print('--- final topol.top ---')\n"
    "print(pathlib.Path('topol.top').read_text())"
))

t3_cells.append(md(VIS_BLOCK.format(subdir="tutorial_3", n=3)))

t3_cells.append(md("""\
### Expected result (Figure 6 of the tutorial)

* Two POPC bilayers slice the box into compartments **A** and **B**.
* Lysozymes (`resname lyz`) live exclusively in compartment **A**.
* Ubiquitins (`resname ubq`) live exclusively in compartment **B**, hugging
  the membrane surfaces (within 4 nm).
* No protein overlaps with either membrane.
"""))

t3_cells.append(md("""\
## 6. (Optional) Run a GROMACS simulation

The Bentopy tutorial finishes with energy minimisation, an equilibration
and a 1 µs production run on the system you just built.  On Delta this
should be submitted through SLURM, **not** run inside the notebook.

`setup_martini.sh` already placed `run_md.slurm` next to this notebook.
Submit it from a terminal:

```bash
cd /projects/bgvl/$USER/Martini/tutorial_3
sbatch run_md.slurm
squeue -u $USER
```

The template targets the `gpuH200x8` partition under the
`bgvl-delta-gpu` account.  Edit the `#SBATCH` headers if you'd rather
use `gpuA100x4` or `gpuA40x4` (usually shorter queues).
"""))

t3_cells.append(code(
    "# Quick sanity check: confirm the SLURM template is in place.\n"
    "slurm = NB_DIR / 'run_md.slurm'\n"
    "if slurm.exists():\n"
    "    print(f'OK: {slurm} ({slurm.stat().st_size} bytes)')\n"
    "else:\n"
    "    print(f'MISSING: {slurm}')\n"
    "    print('Re-run setup_martini.sh to copy it in.')"
))

write_notebook("tutorial_3", "tutorial_3_multi_compartment.ipynb", t3_cells)


print("Wrote:")
for nb in sorted(HERE.glob("tutorial_*.ipynb")):
    print(f"  {nb}  ({nb.stat().st_size} bytes)")
