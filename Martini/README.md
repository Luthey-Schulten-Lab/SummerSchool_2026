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

1. Run the tutorials on the QCB Delta Gateway
2. Tutorial 1 - Basic packing of lysozyme in a water box
3. Tutorial 2 - Membrane packing of Aquaporin Z in a POPE/POPC bilayer
4. Tutorial 3 - Multi-compartment POPC vesicle + GROMACS MD on a GPU node
5. Visualize the results in Open OnDemand

## 1. Run the tutorials on the QCB Delta Gateway

You run these tutorials in **JupyterLab on the QCB Delta Gateway** (the
browser environment you already use for the other Summer School tutorials) and
then visualize the results in **Delta Open OnDemand**. No SSH is needed to run
the packing tutorials.

### Open the tutorial in the Gateway

1. Open the **QCB Delta Gateway** and start a JupyterLab session.
2. Make sure the `SummerSchool_2026` repository is in your Gateway workspace
   (clone it, or `git pull` if you already have it).
3. Navigate to `SummerSchool_2026/Martini/tutorial_1/` and open
   [`tutorial_1_basic_packing.ipynb`](tutorial_1/tutorial_1_basic_packing.ipynb).

Each `tutorial_<N>/` folder ships the inputs it needs:

```
Martini/
├── data/                     # shared inputs (structures, topology, mdp_files)
├── tutorial_1/               # tutorial_1_basic_packing.ipynb + vis_t1.tcl
├── tutorial_2/               # tutorial_2_membrane_packing.ipynb + vis_t2.tcl
├── tutorial_3/               # tutorial_3_multi_compartment.ipynb + vis_t3.tcl + run_md.slurm
└── README.md
```

Inside each `tutorial_<N>/` folder, `structures/`, `topology/` and `mdp_files/`
point back to `data/`, so you only carry one copy of the inputs.

### Bentopy environment

The notebooks need the `bentopy` CLI tools (`bentopy-pack`, `bentopy-render`,
`bentopy-solvate`, `bentopy-mask`, `bentopy-merge`). The first cell of every
notebook puts them on `$PATH`: it uses a `bentopy` virtual environment if one
is already active, and otherwise falls back to the shared workshop venv at
`/projects/bgvl/alfiaparvez/bentopy_tutorial/.venv/`.

To create your own venv instead:

```bash
python3 -m venv bentopy-venv
source bentopy-venv/bin/activate
pip install bentopy
```

### Where your results are saved

You can't write into `/projects/bgvl/<your-NCSA-name>` from the Gateway, so the
first cell of every notebook writes all outputs into your **Gateway project
home** instead:

```
/projects/bgvl/$USER/Martini/tutorial_<N>/
```

On the Gateway, `$USER` is your Gateway account name (which differs from your
NCSA login name). The setup cell creates that folder, copies the inputs and the
VMD script next to the outputs, switches into it, and **prints the exact path**
so you know where everything lands. The folder is created group-readable
(`umask 002`), which is what lets your real Delta account open the results in
Open OnDemand for visualization (see Section 5).

Then work through the notebook cell by cell.

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

## 5. Visualize the results in Open OnDemand

The Gateway can't run an interactive 3-D viewer, so visualization is done in
**Delta Open OnDemand**, where the session runs as your own Delta account and
can read the files the Gateway wrote into your Gateway project home.

1. Open <https://openondemand.delta.ncsa.illinois.edu/> and start an
   **Interactive Apps → Desktop** session.
2. Open a terminal in that desktop and run the command **printed by the setup
   cell of the notebook** (it points at
   `/projects/bgvl/<your-gateway-home>/Martini/tutorial_<N>`):

```bash
export PATH=/projects/bgvl/alfiaparvez/software/vmd/bin:$PATH
cd /projects/bgvl/<your-gateway-home>/Martini/tutorial_1
vmd -e vis_t1.tcl solvated_system.gro
```

(replace `tutorial_1`/`vis_t1.tcl` with the corresponding files for the other
tutorials). A VMD installation with GPU acceleration is available on Delta at
`/projects/bgvl/alfiaparvez/software/vmd/`.

> [!NOTE]
> `<your-gateway-home>` is your Gateway account name, **not** your NCSA login
> name. Copy the exact path from the setup cell's output.

### Optional: GROMACS MD (Tutorial 3)

Tutorial 3 finishes with an optional energy minimization / equilibration /
production MD run. That step uses SLURM on a GPU node, which the Gateway can't
submit, so run it from a **Delta login shell**:

```bash
ssh USERNAME@login.delta.ncsa.illinois.edu
cd /projects/bgvl/<your-gateway-home>/Martini/tutorial_3
sbatch run_md.slurm
squeue -u $USER
```

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
