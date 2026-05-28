# Martini 3 Coarse-Grained Models with Bentopy

## Description:

In the ***Martini 3 Coarse-Grained Models with Bentopy*** tutorial you will
learn how to build large coarse-grained biomolecular systems with the
[Bentopy](https://github.com/marrink-lab/bentopy) packing tool[^bentopy] and
the [Martini 3 force field](http://cgmartini.nl/)[^martini3]. The three
tutorials follow the official Bentopy walkthrough on
[cgmartini.nl][^cgmartini-tutorial] and progressively introduce the workflow:

1. Pack a small protein in a water box.
2. Embed a transmembrane protein in a lipid bilayer.
3. Build a multi-compartment system (a POPC vesicle).

Every tutorial ends with a solvated `.gro` file that can be visualized in
[VMD](https://www.ks.uiuc.edu/Research/vmd/) and, in the third tutorial,
energy-minimized / equilibrated / run with [GROMACS](https://www.gromacs.org/).

*This tutorial was prepared for the third edition of the STC-QCB Summer
School (2026).*

## Outline:

1. Set up the tutorial on Delta
2. Tutorial 1 - Basic packing of lysozyme in a water box
3. Tutorial 2 - Membrane packing of Aquaporin Z in a POPE/POPC bilayer
4. Tutorial 3 - Multi-compartment POPC vesicle + GROMACS MD on a GPU node

## 1. Set up the tutorial on Delta

You will SSH into [NCSA Delta](https://docs.ncsa.illinois.edu/systems/delta/en/latest/quick_start.html)
to run the tutorials.

### Log in to Delta

```bash
ssh USERNAME@login.delta.ncsa.illinois.edu
```
> [!WARNING]
> ***Replace*** `USERNAME` with your Delta username.

### Copy the tutorial materials into your own directory

The workshop master copy lives at
`/projects/bgvl/alfiaparvez/SummerSchool_2026/Martini/`. Copy it into your own
user space with the helper script:

```bash
bash /projects/bgvl/alfiaparvez/SummerSchool_2026/Martini/copy.sh
```

This creates `/projects/bgvl/$USER/SummerSchool_2026/Martini/` containing
three self-contained tutorial folders plus a shared `data/` directory:

```
Martini/
├── data/                     # shared inputs (structures, topology, mdp_files)
├── tutorial_1/               # tutorial_1_basic_packing.ipynb + vis_t1.tcl
├── tutorial_2/               # tutorial_2_membrane_packing.ipynb + vis_t2.tcl
├── tutorial_3/               # tutorial_3_multi_compartment.ipynb + vis_t3.tcl + run_md.slurm
├── copy.sh
├── scripts/build_notebooks.py
└── README.md
```

Inside each `tutorial_<N>/` folder, `structures/`, `topology/` and `mdp_files/`
are symlinks back to `data/`, so you only carry one copy of the inputs.

### Bentopy environment

A shared Python virtual environment with `bentopy` (and all its CLI tools:
`bentopy-pack`, `bentopy-render`, `bentopy-solvate`, `bentopy-mask`,
`bentopy-merge`) is preinstalled at:

```
/projects/bgvl/alfiaparvez/bentopy_tutorial/.venv/
```

The first cell of every notebook prepends
`/projects/bgvl/alfiaparvez/bentopy_tutorial/.venv/bin` to `$PATH`, so the
notebooks work out of the box on Delta with the default `python3` kernel -
no custom kernel registration required.

If you want to recreate the venv yourself (e.g. on a different machine):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install bentopy
```

### Launch Jupyter Notebook on Delta

Submit a Jupyter session to a Delta GPU (or CPU) node. Tutorials 1 and 2 are
fast and only need a CPU; tutorial 3 calls `bentopy-render` and benefits from
a GPU node, and the optional MD step requires a GPU.

> [!NOTE]
> The simplest path is to use the **Delta Open OnDemand** Jupyter app at
> <https://openondemand.delta.ncsa.illinois.edu/>:
>
> 1. *Interactive Apps -> Jupyter Notebook (Project)*
> 2. Account: `bgvl-delta-cpu` (Tut. 1, 2) or `bgvl-delta-gpu` (Tut. 3)
> 3. Working directory: `/projects/bgvl/$USER/SummerSchool_2026/Martini`
> 4. Open `tutorial_1/tutorial_1_basic_packing.ipynb` (etc.) and run cells.

If you prefer the command-line route (an `srun` + SSH tunnel), follow the
recipe described in the [CME README](../CME/README.md#launch-jupyter-notebook-on-delta)
and replace `/projects/beyi` with `/projects/bgvl`.

## 2. Tutorial 1 - Basic packing

**Go to [tutorial_1/](tutorial_1/)** and open
[`tutorial_1_basic_packing.ipynb`](tutorial_1/tutorial_1_basic_packing.ipynb).

You will

- write the YAML *recipe* that tells `bentopy-pack` to drop 50 lysozyme
  monomers into a 30 nm cubic box,
- render the placements into a Martini `system.gro`,
- solvate with Martini water (regular and antifreeze beads),
- visualize the result in VMD with `vis_t1.tcl`.

## 3. Tutorial 2 - Membrane packing

**Go to [tutorial_2/](tutorial_2/)** and open
[`tutorial_2_membrane_packing.ipynb`](tutorial_2/tutorial_2_membrane_packing.ipynb).

You will

- pack a POPE / POPC bilayer with a known *membrane mask*,
- insert four copies of Aquaporin Z (AQPZ) into the bilayer,
- solvate the system above and below the membrane,
- visualize the lipid headgroups + transmembrane proteins with
  `vis_t2.tcl`.

## 4. Tutorial 3 - Multi-compartment + GROMACS MD

**Go to [tutorial_3/](tutorial_3/)** and open
[`tutorial_3_multi_compartment.ipynb`](tutorial_3/tutorial_3_multi_compartment.ipynb).

You will

- build a spherical POPC vesicle with `bentopy-pack`'s mask feature,
- merge an inner and an outer compartment with `bentopy-merge`,
- solvate, visualize with `vis_t3.tcl`,
- run an energy minimization, NPT equilibration and 1 ns production MD on
  a Delta GPU node by submitting [`run_md.slurm`](tutorial_3/run_md.slurm).

## Visualization with VMD

A VMD installation is available on Delta at
`/projects/bgvl/alfiaparvez/software/vmd/`. From a Delta OnDemand desktop
session:

```bash
export PATH=/projects/bgvl/alfiaparvez/software/vmd/bin:$PATH
cd /projects/bgvl/$USER/SummerSchool_2026/Martini/tutorial_1
vmd -e vis_t1.tcl solvated_system.gro
```

(replace `tutorial_1`/`vis_t1.tcl` with the corresponding files for the
other tutorials.)

## References:

[^bentopy]: Mol. M. P., Tsanai, M., Marrink, S. J., & Wassenaar, T. A. (2024).
Bentopy: a versatile coarse-grained Martini system builder.
*ChemRxiv* (preprint). <https://github.com/marrink-lab/bentopy>

[^martini3]: Souza, P. C. T., Alessandri, R., Barnoud, J., Thallmair, S.,
Faustino, I., Gruenewald, F., Patmanidis, I., Abdizadeh, H., Bruininks, B.
M. H., Wassenaar, T. A., Kroon, P. C., Melcr, J., Nieto, V., Corradi, V.,
Khan, H. M., Domanski, J., Javanainen, M., Martinez-Seara, H., Reuter, N.,
... Marrink, S. J. (2021). Martini 3: a general purpose force field for
coarse-grained molecular dynamics. *Nature Methods*, **18**, 382-388.
<https://doi.org/10.1038/s41592-021-01098-3>

[^cgmartini-tutorial]: Marrink Lab, *Bentopy / Martini 3 tutorial*.
<https://cgmartini.nl/docs/tutorials/Martini3/Bentopy/>
