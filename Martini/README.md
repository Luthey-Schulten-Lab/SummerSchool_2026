# Modeling Cells with the Martini Force Field

<table>
<tr>
<td width="60%" valign="top">

Coarse-grained molecular dynamics simulations make it possible to study
biological systems across a wide range of length and time scales, with
whole-cell simulations only recently coming within reach.

In this module, each tutorial focuses on a class of cellular components and
how to build and simulate them in Martini: lipid bilayers, proteins, vesicles
with embedded membrane proteins, polymers and nucleic acids, and crowded
molecular systems. The final tutorial integrates these elements into a toy
model of a bacterial cell.

**Tools covered:** `martinize2` · `TS2CG` · `Polyply` · `Bentopy` · `GROMACS`

</td>
<td width="40%" align="center" valign="top">

<img src="./figures/cell.png" width="320">

<sub><i>Snapshot from a molecular dynamics simulation of the JCVI-syn3A minimal bacterial cell (Stevens et al., 2026).</i></sub>

</td>
</tr>
</table>

---

### Prerequisites

A basic understanding of molecular dynamics and prior familiarity with GROMACS
is assumed. For details on GROMACS usage, force fields, and run parameters, see
the GROMACS [user guide](https://manual.gromacs.org/current/user-guide/index.html)
and [website](www.gromacs.org). A collection of well-written GROMACS tutorials is
also available [here](https://tutorials.gromacs.org/).

---

### Tutorials

|  #  | Tutorial | Topic |
| :-: | :--- | :--- |
| I   | [Bilayer Self-Assembly](01_bilayer_self_assembly/tutorial.md)        | Self-assemble a lipid bilayer from a random initial configuration. |
| II  | [Protein Basics](02_protein_basics/tutorial.md)                      | Coarse-grain and simulate a protein with `martinize2`. |
| III | [Membranes and Vesicles](03_membranes_and_vesicles/tutorial.md)      | Build a vesicle with embedded membrane proteins using `TS2CG`. |
| IV  | [Polymers and DNA](04_polymers_and_DNA/tutorial.md)                  | Generate polymers and single-stranded DNA with `Polyply`. |
| V   | [Packing Biomolecular Systems](05_packing_biomolecular_systems/tutorial.md) | Assemble crowded molecular systems with `Bentopy`. |
| VI  | [Martini Cell](06_martini_cell/tutorial.md)                          | Integrate all components into a toy model of a bacterial cell. |

---

### Setup

**1. Launch a *Desktop* session on Delta** through
[Open OnDemand](https://openondemand.delta.ncsa.illinois.edu/):

- Log in with your NCSA credentials.
- From the menu, open **Desktop**.
- Fill in the session form with the settings below (defaults are fine for
  anything not listed; *Name of the reservation* you can leave empty):

  | Field | Value |
  | :--- | :--- |
  | Container image | `RHEL 9 w/ CUDA 12.8` |
  | Account         | `bgvl-delta-gpu` |
  | Partition       | `gpuA40x4` |
  | Duration of job | `4:00:00` |
  | Number of CPUs  | `8` |
  | Amount of RAM   | `64G` |
  | Number of GPUs  | `1` |

- Click **Launch** and wait for the session to be allocated (usually under
  a minute).
- Once the session card shows **Running**, click **Launch Desktop** to open
  the remote desktop in your browser. Open a terminal inside the desktop
  to continue.

**2. Clone the workshop repository:**

```sh
git clone https://github.com/Luthey-Schulten-Lab/SummerSchool_2026.git
cd SummerSchool_2026/Martini
```

Each tutorial lives in its own folder, with all input files prepared for use.

**3. Set up the workshop environment** with micromamba. First install
micromamba (accept all default settings when prompted):

```sh
"${SHELL}" <(curl -L micro.mamba.pm)
```

Then create the workshop environment from the provided file and activate it:

```sh
micromamba create -n workshop -f ./files/environment.yml
micromamba activate workshop
```

---

These tutorials were written by *Jan Stevens and Marieke Westendorp*, with
several tutorials adapted from the
[*2025 Martini online workshop*](https://cgmartini.nl/docs/tutorials/Martini3/workshop.html).
