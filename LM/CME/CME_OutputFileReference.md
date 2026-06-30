# CME Output File Reference #

Lattice Microbes CME simulations produce a single **`.lm`** file per run. Despite the `.lm` extension the file is a standard **HDF5** file (format version 4) and can be opened with any HDF5-compatible tool (h5py, HDFView, etc.).

## Opening the File

```python
import h5py
import numpy as np

f = h5py.File('MySimulation.lm', 'r')
```

---

## Top-Level Structure

```
MySimulation.lm
├── Model/
│   └── Reaction/
│       ├── DependencyMatrix
│       ├── InitialSpeciesCounts
│       ├── ReactionRateConstants
│       ├── ReactionTypes
│       └── StoichiometricMatrix
├── Parameters/
│   └── SpeciesNames
└── Simulations/
    ├── 0000001/
    │   ├── FirstPassageTimes/
    │   │   └── 00/
    │   │       ├── Counts
    │   │       └── Times
    │   ├── SpeciesCountTimes
    │   └── SpeciesCounts
    ├── 0000002/
    │   └── ...
    └── 000000N/   (one group per replicate)
```

---

## `Model/Reaction/` Datasets

### `SpeciesNames`
| Property | Value |
|----------|-------|
| Path | `Parameters/SpeciesNames` |
| Shape | `(n_species, 1)` |
| dtype | `object` (byte strings) |

Ordered list of species names. The index of each name corresponds to the column index used in `SpeciesCounts`.

```python
names = f['Parameters/SpeciesNames'][:, 0]
# e.g. [b'gene', b'mRNA', b'protein']
```

---

### `InitialSpeciesCounts`
| Property | Value |
|----------|-------|
| Path | `Model/Reaction/InitialSpeciesCounts` |
| Shape | `(n_species,)` |
| dtype | `uint32` |

Particle counts at $t = 0$ for each species, in the same order as `SpeciesNames`.

```python
x0 = f['Model/Reaction/InitialSpeciesCounts'][...]
```

---

### `StoichiometricMatrix`
| Property | Value |
|----------|-------|
| Path | `Model/Reaction/StoichiometricMatrix` |
| Shape | `(n_species, n_reactions)` |
| dtype | `int32` |

Each column describes the net change in particle count for every species when that reaction fires. Positive entries indicate production; negative entries indicate consumption; zero indicates no change.

```
       rxn0  rxn1  rxn2  rxn3
gene  [  0    0    0    0  ]
mRNA  [  1   -1    0    0  ]
ptn   [  0    0    1   -1  ]
```

---

### `ReactionTypes`
| Property | Value |
|----------|-------|
| Path | `Model/Reaction/ReactionTypes` |
| Shape | `(n_reactions,)` |
| dtype | `uint32` |

Integer code identifying the kinetic order of each reaction:

| Code | Meaning |
|------|---------|
| `1` | First-order (unimolecular) |
| `2` | Second-order (bimolecular) |

---

### `ReactionRateConstants`
| Property | Value |
|----------|-------|
| Path | `Model/Reaction/ReactionRateConstants` |
| Shape | `(n_reactions, 10)` |
| dtype | `float64` |

Each row holds the rate constant(s) for one reaction. For first-order reactions only `[0]` is used; all other entries are `NaN`. Units depend on reaction order: s⁻¹ for first-order, M⁻¹ s⁻¹ (particle-based) for second-order.

```python
k = f['Model/Reaction/ReactionRateConstants'][:, 0]
```

---

### `DependencyMatrix`
| Property | Value |
|----------|-------|
| Path | `Model/Reaction/DependencyMatrix` |
| Shape | `(n_species, n_reactions)` |
| dtype | `uint32` |

Boolean matrix (0 or 1). Entry `[i, j] = 1` means species `i` appears in the propensity function of reaction `j`. Used internally by the solver to determine which propensities must be recomputed after a reaction fires.

---

## `Simulations/` — Per-Replicate Data

Replicates are stored as zero-padded 7-digit keys (`0000001`, `0000002`, …). Each group contains identical datasets; the number of replicates equals the number of groups.

```python
reps = sorted(f['Simulations'].keys())   # ['0000001', '0000002', ...]
n_reps = len(reps)
```

---

### `SpeciesCountTimes`
| Property | Value |
|----------|-------|
| Path | `Simulations/XXXXXXX/SpeciesCountTimes` |
| Shape | `(n_timepoints,)` |
| dtype | `float64` |

Simulation time (seconds) corresponding to each row of `SpeciesCounts`. Timepoints are written at the interval set by `writeInterval` during simulation setup.

```python
ts = f['Simulations/0000001/SpeciesCountTimes'][:]
```

---

### `SpeciesCounts`
| Property | Value |
|----------|-------|
| Path | `Simulations/XXXXXXX/SpeciesCounts` |
| Shape | `(n_timepoints, n_species)` |
| dtype | `int32` |

Integer particle counts for every species at every recorded timepoint. Columns are ordered identically to `SpeciesNames`.

```python
counts = f['Simulations/0000001/SpeciesCounts'][:]
mRNA_traj = counts[:, 1]   # index matches SpeciesNames order
```

To read all replicates into a single array:

```python
names  = [n.decode() for n in f['Parameters/SpeciesNames'][:, 0]]
ts     = f['Simulations/0000001/SpeciesCountTimes'][:]
all_counts = np.stack(
    [f[f'Simulations/{r}/SpeciesCounts'][:] for r in reps],
    axis=0
)   # shape: (n_reps, n_timepoints, n_species)
```

---

### `FirstPassageTimes/`
| Property | Value |
|----------|-------|
| Path | `Simulations/XXXXXXX/FirstPassageTimes/00/` |
| Subkeys | `Counts` `(1,) uint32`, `Times` `(1,) float64` |

Records the first time each species reaches a threshold count. In basic GIP simulations a single threshold (`00`) is defined. If no threshold is crossed during the simulation, `Times` contains `0.0`.

```python
fpt_time = f['Simulations/0000001/FirstPassageTimes/00/Times'][0]
```

---

## Quick-Access Cheatsheet

```python
import h5py, numpy as np

f     = h5py.File('MySimulation.lm', 'r')
reps  = sorted(f['Simulations'].keys())

names = [n.decode() for n in f['Parameters/SpeciesNames'][:, 0]]
ts    = f['Simulations/0000001/SpeciesCountTimes'][:]

# All replicates: shape (n_reps, n_timepoints, n_species)
data  = np.stack([f[f'Simulations/{r}/SpeciesCounts'][:] for r in reps])

# Population mean trajectory for species 'mRNA'
idx   = names.index('mRNA')
mean_mRNA = data[:, :, idx].mean(axis=0)

f.close()
```
