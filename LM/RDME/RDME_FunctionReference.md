# RDME Function Reference — Lattice Microbes / jLM #

This reference covers every jLM function used across the RDME tutorial notebooks. Functions are organized by workflow stage so you can follow the same order when building your own simulations.

---

## Workflow Overview ##

```
1. Create simulation object  →  RDMESim(...)
2. Build geometry            →  RegionBuilder, build.ellipsoid(), build.dilate(), build.compose()
3. Define regions            →  sim.region()
4. Define species            →  sim.species()
5. Place particles           →  sim.placeNumber(), sim.distributeNumber()
6. Define reactions          →  sim.rateConst(), region.addReaction()  [or sim.reaction()]
7. Set diffusion             →  sim.diffusionConst(), sim.transitionRate()
8. Finalize and run          →  sim.finalize(), sim.run()
9. Analyze (whole-system)    →  RDMEFile(), traj.getNumberTrajectory()
10. Analyze (region-specific) →  traj._speciesIdParse(), traj._regionIdParse(), LMLattice.latticeStatsAll_h5fmt()
```

---

## Key Constraints — Read This First ##

| Constraint | Detail |
|---|---|
| **Lattice dimensions divisible by 32** | Each of `nx`, `ny`, `nz` must be a multiple of 32 (GPU block-tiling requirement). |
| **Timestep upper bound** | `dt ≤ latticeSpacing² / (6 × D_max)`. Violating this makes hop probabilities > 1. Use `sim.maxDiffusionRate()` to check. |
| **Particles per site limit** | Set at LM compile time (`sim.pps`; typically 4 for Int lattice). Exceeding it causes an error at `finalize()`. |
| **Diffusion table completeness** | Every `(species, fromRegion, toRegion)` combination must be assigned a rate before `sim.finalize()`. The standard idiom: zero everything first with `sim.transitionRate(None, None, None, sim.diffusionZero)`, then selectively enable. |
| **Reaction units are macroscopic** | Supply rates in M/s, s⁻¹, or M⁻¹ s⁻¹. jLM converts to stochastic units internally. |
| **Solver must match lattice type** | `latticeType="Int"` → `IntMpdRdmeSolver()`. `latticeType="Byte"` → `MpdRdmeSolver()`. Mismatching raises a runtime error. |
| **`finalize()` must be called before `run()`** | `finalize()` writes the `.lm` file; `run()` calls the solver on that file. |
| **Diffusion transitions are directional** | Setting `A → B` diffusion does **not** automatically set `B → A`. Each direction must be specified explicitly. |
| **Reactions are per-region** | A reaction added with `region.addReaction()` only fires in that region. Call again on each region where it should occur. |

---

## 1. Simulation Setup ##

### `RDMESim` ###

Creates the central simulation object that holds all model data — geometry, species, reactions, diffusion constants, and particle placements.

```python
sim = RDMESim(name, filename, dimensions, latticeSpacing, regionName,
              latticeType=None, dt=None)
```
+ `name` (*String*): human-readable simulation label;
+ `filename` (*String*): path to the output `.lm` (HDF5) file;
+ `dimensions` (*list of 3 Ints*): number of voxels along x, y, z — **each must be divisible by 32**;
+ `latticeSpacing` (*Float*): side length of each cubic voxel in **metres**;
+ `regionName` (*String*): name of the default region applied to all voxels at initialisation;
+ `latticeType` (*String, optional*): `"Int"` for up to ~4 million species (requires `IntMpdRdmeSolver`); `"Byte"` for up to 255 species (default, requires `MpdRdmeSolver`);
+ `dt` (*Float, optional*): simulation timestep in **seconds** — must satisfy `dt ≤ latticeSpacing² / (6 × D_max)`;

**After construction, set these attributes before finalizing:**

| Attribute | Description | Units |
|---|---|---|
| `sim.simulationTime` | Total biological time to simulate | seconds |
| `sim.latticeWriteInterval` | Timesteps between full lattice snapshots | timesteps |
| `sim.speciesWriteInterval` | Timesteps between species count writes | timesteps |

**Convenient namespaces created automatically:**
- `sim.sp.<name>` — access a species by name
- `sim.reg.<name>` — access a region by name
- `sim.rc.<name>` — access a rate constant by name
- `sim.dc.<name>` — access a diffusion constant by name

**Usage:**
```python
sim = RDMESim(name         = "RDME_GIP",
              filename      = "./Results/output.lm",
              dimensions    = [64, 64, 64],
              latticeSpacing = 8e-9,          # 8 nm voxels
              regionName    = "extracellular",
              latticeType   = "Int",
              dt            = 50e-6)          # 50 µs timestep

sim.simulationTime       = 6300    # seconds
sim.latticeWriteInterval = 20000   # write every 20000 × 50 µs = 1 s
sim.speciesWriteInterval = 20000
```

---

## 2. Building the Geometry ##

Geometry is built using the `RegionBuilder` class, which provides methods for creating boolean 3-D voxel masks. Masks are combined with `compose()` to assign region types to the site lattice.

### `RegionBuilder` ###

Creates a helper object for designing region masks.

```python
build = RegionBuilder(sim)
```
+ `sim` (*RDMESim object*): the simulation object from which lattice dimensions are read;

**Returns:** a `RegionBuilder` instance with the following useful structuring elements:

| Property | Shape | Description |
|---|---|---|
| `build.se6` | `(3,3,3) bool` | 6-connectivity (face-adjacent voxels only) |
| `build.se26` | `(3,3,3) bool` | 26-connectivity (all 26 neighbors including edges and corners) |

---

### `build.ellipsoid` ###

Constructs a boolean voxel mask for an ellipsoidal region (pass a scalar radius for a sphere).

```python
mask = build.ellipsoid(radius, center=None, angles=None)
```
+ `radius` (*Float or list of 3 Floats*): semiaxis length(s) in **voxels**; a single value gives a sphere; a 3-element list `[rx, ry, rz]` gives an ellipsoid; convert from metres with `radius_voxels = int(np.ceil(radius_m / lattice_spacing))`;
+ `center` (*list of 3 Floats, optional*): voxel coordinates of the centroid; defaults to the lattice midpoint;
+ `angles` (*list of 3 Floats, optional*): ZXZ Euler rotation angles `[alpha, beta, gamma]` in **radians**; defaults to no rotation;

**Returns:** `numpy.ndarray(shape=(nx,ny,nz), dtype=bool)` — `True` inside the ellipsoid.

---

### `build.dilate` ###

Morphological binary dilation — expands a boolean mask outward by one structuring-element step.

```python
dilated_mask = build.dilate(binaryMask, se=None, radius=None)
```
+ `binaryMask` (*numpy.ndarray, bool*): the mask to dilate;
+ `se` (*numpy.ndarray, bool*): structuring element (e.g. `build.se26` or `build.se6`); use `se26` to guarantee membrane connectivity;
+ `radius` (*Int, optional*): number of 6-connected dilation iterations; alternative to `se`;

**Returns:** `numpy.ndarray(shape=(nx,ny,nz), dtype=bool)` — the dilated mask.

**Usage (building shell and membrane regions):**
```python
cyto_dilation = build.dilate(cytosol_mask, se=build.se26)
shell_mask    = cyto_dilation & ~cytosol_mask         # 1-voxel-thick inner layer

cyto_dilation = build.dilate(cyto_dilation, se=build.se26)
membrane_mask = cyto_dilation & ~shell_mask & ~cytosol_mask
```

---

### `build.compose` ###

Writes a set of `(Region, mask)` pairs into the simulation site lattice. **Later entries overwrite earlier ones** where masks overlap.

```python
build.compose(
    (sim.region('extracellular'), extracellular_mask),
    (sim.region('cytosol'),       cytosol_mask),
    (sim.region('DNA'),           dna_mask),
    ...
)
```
+ Each positional argument is a 2-tuple of `(Region object, boolean mask)`;
+ Regions must already exist (created via `sim.region()`) before being passed here;
+ Order matters — if any two masks overlap, the later one wins;

**Returns:** None (modifies `sim.siteLattice` in place).

---

### Other `RegionBuilder` geometry methods ###

| Method | Signature | What it creates |
|---|---|---|
| `build.cylinder` | `(radius, length, center=None, angles=None)` | Cylindrical boolean mask |
| `build.capsule` | `(length, width, center=None, angles=None)` | Spherocylinder (cylinder + hemispherical endcaps) |
| `build.box` | `(lx, ly, lz, center=None, angles=None)` | Rectangular cuboid |
| `build.erode` | `(binaryMask, radius=None, se=None)` | Morphological erosion |
| `build.closing` | `(binaryMask, radius=None, se=None)` | Dilate then erode |
| `build.opening` | `(binaryMask, radius=None, se=None)` | Erode then dilate |

**NOTE:** for a description of each of these methods, see the [RegionBuilderCheatSheet](./RegionBuilderCheatSheet.md).

---

## 3. Defining Regions ##

### `sim.region` ###

Looks up an existing region by name, or creates and registers a new one.

```python
reg_obj = sim.region(name, texRepr=None, annotation=None)
```
+ `name` (*String*): unique region name;
+ `texRepr` (*String, optional*): LaTeX representation for display;
+ `annotation` (*String, optional*): free-text description;

**Returns:** a `Region` object (also accessible as `sim.reg.<name>`).

**Note:** the default region (named in the `RDMESim` constructor) is created automatically. You only need to call `sim.region()` for additional regions.

---

## 4. Defining Species ##

### `sim.species` ###

Registers a chemical species in the simulation.

```python
sp_obj = sim.species(name, texRepr=None, annotation=None)
```
+ `name` (*String*): unique species name;
+ `texRepr` (*String, optional*): LaTeX math representation (e.g. `r'D_{mRNA}'`);
+ `annotation` (*String, optional*): free-text description;

**Returns:** a `Species` object (also accessible as `sim.sp.<name>`).

**Usage:**
```python
sim.species("gene", texRepr="gene", annotation="gene in GIP process")
sim.species("mRNA", texRepr=r"mRNA")
sim.species("P",    texRepr=r"P_{protein}")
```

---

## 5. Placing Particles ##

### `sim.placeNumber` ###

Places exactly `n` particles of a species at a specific voxel. Used for deterministic placement (e.g., a single gene at a known DNA locus).

```python
sim.placeNumber(sp, x, y, z, n)
```
+ `sp` (*Species object*): species to place;
+ `x`, `y`, `z` (*Int*): voxel indices of the target site;
+ `n` (*Int*): number of particles to place at that voxel;

**Returns:** None. Placement is queued and resolved during `sim.finalize()`.

**Usage:**
```python
gene_pos = DNA_pos[0]
sim.placeNumber(sp=sp.gene, x=gene_pos[0], y=gene_pos[1], z=gene_pos[2], n=1)
```

---

### `sim.distributeNumber` ###

Distributes `count` particles of a species uniformly and randomly across all voxels of a region.

```python
sim.distributeNumber(sp, reg, count)
```
+ `sp` (*Species object*): species to distribute;
+ `reg` (*Region object*): target region;
+ `count` (*Int*): total number of particles to distribute;

**Returns:** None. Placement is resolved during `sim.finalize()`.

**Usage:**
```python
sim.distributeNumber(sp=sp.mRNA, reg=reg.cytosol, count=1)
```

---

### `sim.distributeConcentration` ###

Converts a molar concentration to a particle count and distributes particles across a region.

```python
sim.distributeConcentration(sp, reg, conc)
```
+ `sp` (*Species object*): species to distribute;
+ `reg` (*Region object*): target region;
+ `conc` (*Float*): initial concentration in **mol/L (M)**;

**Returns:** None. Count computed as `int(N_A × V_voxel × conc × n_voxels_in_region)`.

---

## 6. Defining Reactions ##

### Rate constant units ###

The RDME simulation object expects rate constants in **macroscopic concentration-based units**:

| Reaction order | Example | Units |
|---|---|---|
| 0th order | $\emptyset \rightarrow A$ | M / s |
| 1st order | $C \rightarrow A + B$ | s⁻¹ |
| 2nd order | $A + B \rightarrow C$ | M⁻¹ s⁻¹ |

jLM converts these automatically using:
$$k_{\text{stochastic}} = k_{\text{det}} \cdot \left( N_A \cdot V_{\text{voxel}} \right)^{1 - \text{order}}$$

---

### `sim.rateConst` ###

Creates (or retrieves) a named rate constant.

```python
rc_obj = sim.rateConst(rate, value, order, texRepr=None, annotation=None)
```
+ `rate` (*String*): unique name for this rate constant; retrieved later as `sim.rc.<rate>`;
+ `value` (*Float*): numerical value in macroscopic concentration units (see table above);
+ `order` (*Int*): reaction order: `0`, `1`, or `2`;
+ `texRepr` (*String, optional*): LaTeX representation;
+ `annotation` (*String, optional*): free-text description;

**Returns:** a `RateConst` object (also stored as `sim.rc.<rate>`).

**Usage:**
```python
kf = sim.rateConst('kf',    1.07e5, order=2)    # M⁻¹ s⁻¹
kr = sim.rateConst('kr',    0.351,  order=1)    # s⁻¹
k_trans  = sim.rateConst('trans',   6.14e-4, order=1)
k_transl = sim.rateConst('transl',  7.20e-2, order=1)
```

---

### `region.addReaction` ###

Adds a chemical reaction that occurs **only in the region on which the method is called**.

```python
region.addReaction(reactants, products, rate, value=None)
```
+ `reactants` (*list of Species objects or* `[]`): species consumed; pass `[]` for a 0th-order creation reaction;
+ `products` (*list of Species objects or* `[]`): species produced; pass `[]` for pure degradation;
+ `rate` (*RateConst object or String*): rate constant — use the object returned by `sim.rateConst()`, or provide a string name with `value`;
+ `value` (*Float, optional*): if `rate` is a string, creates a new `RateConst` on the fly;

**Returns:** the `Region` object (for chaining).

**Usage (Tutorial 1 style):**
```python
cyt = sim.region('cytosol')
cyt.addReaction([sp.A, sp.B], [sp.C],      kf)
cyt.addReaction([sp.C],       [sp.A, sp.B], kr)
```

---

### `sim.reaction` ###

Alternative to `region.addReaction()` — creates a reaction and optionally assigns it to one or more regions in a single call (Tutorial 2 style).

```python
rxn = sim.reaction(reactants, products, rate, value=None, regions=None, annotation=None)
```
+ `reactants` (*Species, String, or list*): reactant(s);
+ `products` (*Species, String, or list*): product(s);
+ `rate` (*RateConst or String*): rate constant;
+ `value` (*Float, optional*): value if `rate` is a string name;
+ `regions` (*Region, String, or list, optional*): region(s) where reaction is active;
+ `annotation` (*String, optional*): free-text description;

**Returns:** a `BuilderReaction` object.

**Usage (Tutorial 2 style):**
```python
sim.reaction([sp.gene], [sp.gene, sp.mRNA], rc.trans,   regions=[reg.DNA],       annotation="transcription")
sim.reaction([sp.mRNA], [],                 rc.degrad_m, regions=[reg.shell],     annotation="mRNA degradation")
sim.reaction([sp.mRNA], [sp.mRNA, sp.P],   rc.transl,   regions=[reg.ribosomes], annotation="translation")
sim.reaction([sp.P],    [],                 rc.degrad_p, regions=[reg.shell],     annotation="protein degradation")
```

---

## 7. Setting Diffusion Coefficients ##

Diffusion rules determine how likely a particle is to hop between neighbouring voxels each timestep. The hop probability for species α moving from site type $s$ to $s'$ is:

$$q = \frac{D_\alpha(s \to s') \cdot \tau}{\lambda^2}$$

Every `(species, fromRegion, toRegion)` combination must be assigned a rate. The standard workflow is to zero all transitions first, then selectively enable the ones you want.

---

### `sim.diffusionConst` ###

Creates (or retrieves) a named diffusion constant.

```python
dc_obj = sim.diffusionConst(rate, value, texRepr=None, annotation=None)
```
+ `rate` (*String*): unique name; retrieved later as `sim.dc.<rate>`;
+ `value` (*Float*): diffusion coefficient in **m² / s**;
+ `texRepr` (*String, optional*): LaTeX representation;
+ `annotation` (*String, optional*): free-text description;

**Returns:** a `DiffusionConst` object (also stored as `sim.dc.<rate>`).

**Usage:**
```python
d_mRNA    = sim.diffusionConst('mrna',    4.13e-14, annotation="mRNA diffusion in Syn3A")
d_protein = sim.diffusionConst('protein', 1.0e-13)
```

---

### `sim.transitionRate` ###

Assigns a diffusion rate for particles transitioning between two region types. `None` is a wildcard that applies to all species or all regions.

```python
sim.transitionRate(sp, rFrom, rTo, rate, value=None)
```
+ `sp` (*Species object or* `None`): species to configure; `None` applies to **all** species;
+ `rFrom` (*Region object or* `None`): source region; `None` applies to **all** source regions;
+ `rTo` (*Region object or* `None`): destination region; `None` applies to **all** destination regions;
+ `rate` (*DiffusionConst object*): diffusion constant to assign (e.g. `sim.diffusionZero` or a constant from `sim.diffusionConst()`);
+ `value` (*Float, optional*): if `rate` is a string name, creates a new `DiffusionConst` on the fly;

**Returns:** None.

**Standard pattern — zero all first, then enable selectively:**
```python
# Step 1: zero everything
sim.transitionRate(None, None, None, sim.diffusionZero)

# Step 2: enable transitions you want
sim.transitionRate(sp.mRNA, reg.cytosol,   reg.cytosol,   dc.mrna)
sim.transitionRate(sp.mRNA, reg.cytosol,   reg.ribosomes, dc.mrna)
sim.transitionRate(sp.mRNA, reg.ribosomes, reg.cytosol,   dc.mrna)
sim.transitionRate(sp.mRNA, reg.ribosomes, reg.ribosomes, dc.mrna)
```

---

### Built-in diffusion constants ###

| Property | Value | Description |
|---|---|---|
| `sim.diffusionZero` | 0.0 m²/s | Particle cannot hop — effectively immobile in that transition |
| `sim.diffusionFast` | `latticeSpacing² / (6 × dt)` | Maximum allowed diffusion constant (timestep upper bound) |

---

### `sim.maxDiffusionRate` ###

Computes the maximum physically valid diffusion coefficient for the current lattice spacing and timestep.

```python
d_max = sim.maxDiffusionRate(latticeSpacing=None, dt=None)
```
+ `latticeSpacing` (*Float, optional*): voxel edge length in metres; defaults to `sim.latticeSpacing`;
+ `dt` (*Float, optional*): timestep in seconds; defaults to `sim.timestep`;

**Returns:** `float` in m²/s. If any diffusion constant exceeds this value, the simulation will produce unphysical results.

---

## 8. Finalizing and Running ##

### `sim.finalize` ###

Serializes the complete model into the `.lm` HDF5 file. **Must be called before `sim.run()`.**

```python
sim.finalize()
```

**Returns:** None.

**What it does:**
- Deletes any existing file at `sim.filename` and creates a fresh one
- Resolves all queued `placeNumber` / `distributeNumber` calls
- Validates that every `(species, fromRegion, toRegion)` diffusion entry is set — raises `RuntimeError` if any is missing
- Prints `nodiffregion` — a list of region IDs that have no non-zero diffusion rates (for sanity checking)
- Writes all metadata needed for `RDMEFile` to reload the model

---

### `sim.run` ###

Launches the RDME solver and returns a file handle to the results.

```python
result = sim.run(solver=None, replicate=1, seed=None, cudaDevices=None,
                 checkpointInterval=0, sample_frame=False, max_frames=100)
```
+ `solver` (*RDMESolver object, optional*): use `IntMpdRdmeSolver()` for `"Int"` lattice; `MpdRdmeSolver()` for `"Byte"` lattice; defaults to `MpdRdmeSolver()`;
+ `replicate` (*Int, optional*): replicate index written into the output file; increment to add multiple independent trajectories to the same `.lm` file without overwriting (default: `1`);
+ `seed` (*Int, optional*): fixed random seed for reproducibility; `None` uses time-based seeding (default: `None`);
+ `cudaDevices` (*list of Int, optional*): CUDA GPU device indices (e.g. `[0]` for the first GPU); defaults to `[0]`;
+ `checkpointInterval` (*Int, optional*): wall-clock seconds between checkpoint saves; `0` disables (default: `0`);
+ `sample_frame` (*Bool, optional*): if `True`, lattice frames are sub-sampled when the result file is loaded back (default: `False`);
+ `max_frames` (*Int, optional*): target frame count when `sample_frame=True` (default: `100`);

**Returns:** a `jLM.RDME.File` handle to the output file.

**Usage:**
```python
sim.finalize()
sim.run(solver=IntMpdRdmeSolver(), cudaDevices=[0])
```

---

## 9. Visualization / Introspection ##

### `sim.displayGeometry` ###

Renders an interactive 3-D view of the site lattice in the Jupyter notebook.

```python
sim.displayGeometry(filterFunctions=None, mode="widget")
```
+ `filterFunctions` (*dict, optional*): per-region clipping functions `{region_name: f(x,y,z) -> bool}` to hide parts of the lattice;
+ `mode` (*String, optional*): `"widget"` (inline), `"download_x3d"` (download X3D file), or `"download_html"` (download standalone HTML);

---

### `sim.showAllSpecies` ###

Displays a summary of all species — counts by region, diffusion tables, and associated reactions.

```python
sim.showAllSpecies()
```

---

## 10. Analysis ##

### `RDMEFile` ###

Opens an existing `.lm` result file for analysis.

```python
traj = RDMEFile(fname, replicate=1, latticeType=None,
                sample_frame=False, max_frames=100)
```
+ `fname` (*String*): path to the `.lm` HDF5 output file;
+ `replicate` (*Int, optional*): which replicate to load initially (1-based, default: `1`);
+ `latticeType` (*String, optional*): `"Int"` or `"Byte"`; read from file if not specified;
+ `sample_frame` (*Bool, optional*): if `True`, sub-sample the lattice trajectory frames;
+ `max_frames` (*Int, optional*): target frame count when `sample_frame=True`;

**Returns:** a `jLM.RDME.File` instance with:
- `traj.h5` — raw `h5py.File` handle for direct HDF5 access
- `traj.speciesList` — list of all species in the model
- `traj.regionList` — list of all regions in the model

---

### `traj.getNumberTrajectory` ###

Returns the total particle count time series for one or more species across the entire simulation volume. Reads from `SpeciesCounts` — **fast**.

```python
ts, counts = traj.getNumberTrajectory(species=None, regex=None,
                                       replicate=None,
                                       frameStart=None, frameEnd=None,
                                       timeStart=None, timeEnd=None)
```
+ `species` (*Species object, String, or list, optional*): species to retrieve; lists are summed into a single trajectory; mutually exclusive with `regex`;
+ `regex` (*String, optional*): regular expression matched against species names; all matching species are summed; mutually exclusive with `species`;
+ `replicate` (*Int, optional*): replicate to read (1-based); defaults to the current replicate;
+ `frameStart` (*Int, optional*): index of the first frame (inclusive); defaults to start;
+ `frameEnd` (*Int, optional*): index of the last frame (exclusive); defaults to end;
+ `timeStart` (*Float, optional*): simulation time in seconds of the first frame;
+ `timeEnd` (*Float, optional*): simulation time in seconds of the last frame;

**Returns:** two `numpy.ndarray` of shape `(nt,)` — `ts` (evaluation times in seconds) and `counts` (particle counts).

**Usage:**
```python
ts, mRNAs    = traj.getNumberTrajectory(species="mRNA")
ts, proteins = traj.getNumberTrajectory(species="P")
```

---

### `traj._speciesIdParse` ###

Resolves a species name to its numeric ID as stored in the HDF5 lattice arrays.

```python
ids = traj._speciesIdParse(regex, name, startIndex=1)
```
+ `regex` (*String or* `None`): regular expression to match species names; pass `None` when using `name`;
+ `name` (*String or* `None`): exact species name to look up; pass `None` when using `regex`;
+ `startIndex` (*Int, optional*): offset added to the 0-based index; use `startIndex=1` because particle lattice arrays store species IDs starting at 1, not 0 (default: `1`);

**Returns:** list of integer IDs. Take `[0]` for a single species.

**Usage:**
```python
mRNA_id = traj._speciesIdParse(None, "mRNA", startIndex=1)[0]
```

---

### `traj._regionIdParse` ###

Resolves a region name to its numeric ID as stored in the HDF5 site-type lattice array.

```python
ids = traj._regionIdParse(regex, name)
```
+ `regex` (*String or* `None`): regular expression to match region names; pass `None` when using `name`;
+ `name` (*String or* `None`): exact region name to look up;

**Returns:** list of integer IDs. Take `[0]` for a single region.

**Usage:**
```python
cyto_id  = traj._regionIdParse(None, "cytosol")[0]
shell_id = traj._regionIdParse(None, "shell")[0]
```

---

### `LMLattice.latticeStatsAll_h5fmt` ###

Low-level Cython function that counts particles by region and species ID across the entire lattice in a single pass.

> **Background:** `traj.getNumberTrajectoryFromRegion()` contains a bug in jLM 2.5.0 — its internal array is allocated with size equal to the number of *requested* regions (1), but then indexed by the actual region numeric ID (which can be > 0), causing an `IndexError`. This function is the correct workaround: it returns the full count matrix indexed by actual IDs with no truncation.

```python
from jLM import Lattice as LMLattice

pCount, sCount = LMLattice.latticeStatsAll_h5fmt(plattice, slattice)
```
+ `plattice` (*numpy.ndarray, shape (nx, ny, nz, pps), dtype uint32*): particle lattice for a single frame — each entry holds the species ID of the particle occupying that slot (0 = empty);
+ `slattice` (*numpy.ndarray, shape (nx, ny, nz), dtype uint8*): site-type lattice — each entry holds the region ID of that voxel;

**Returns:**
- `pCount` (*numpy.ndarray, shape (16384, 16384)*): `pCount[region_id, species_id]` gives the number of particles of that species in that region;
- `sCount` (*numpy.ndarray, shape (16384,)*): `sCount[region_id]` gives the number of voxels assigned to that region;

**Complete region-specific trajectory pattern:**

```python
from jLM import Lattice as LMLattice
import warnings

# Resolve numeric IDs
mRNA_id  = traj._speciesIdParse(None, "mRNA",    startIndex=1)[0]
cyto_id  = traj._regionIdParse(None,  "cytosol")[0]
shell_id = traj._regionIdParse(None,  "shell")[0]

# Load static site-type lattice (same for every frame)
site_lattice = traj.h5['/Model/Diffusion/LatticeSites'][...]

# Iterate over saved lattice frames
rep_key  = list(traj.h5['Simulations'].keys())[0]
frames   = sorted(traj.h5[f'Simulations/{rep_key}/Lattice'].keys())
ts_r     = traj.h5[f'Simulations/{rep_key}/LatticeTimes'][:]

mRNA_cyto_counts  = np.zeros(len(frames), dtype=np.int64)
mRNA_shell_counts = np.zeros(len(frames), dtype=np.int64)

particle_lattice = None
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i, frame in enumerate(frames):
        if particle_lattice is not None:
            traj.h5[f'Simulations/{rep_key}/Lattice/{frame}'].read_direct(particle_lattice)
        else:
            particle_lattice = traj.h5[f'Simulations/{rep_key}/Lattice/{frame}'][...]
        pCount, _ = LMLattice.latticeStatsAll_h5fmt(particle_lattice, site_lattice)
        mRNA_cyto_counts[i]  = pCount[cyto_id,  mRNA_id]
        mRNA_shell_counts[i] = pCount[shell_id, mRNA_id]
```

**Notes:**
- `read_direct` reuses the pre-allocated array rather than allocating a new one each frame — important for memory efficiency with large lattices;
- The `warnings` context suppresses HDF5-level deprecation warnings that appear in jLM 2.5.0 when reading lattice frames;
- Prefer `traj._speciesIdParse(..., startIndex=1)` over hard-coding species IDs — the index depends on the order species were registered;

---

## Quick-Reference Table ##

| Function | Stage | What it does |
|---|---|---|
| `RDMESim(name, filename, dimensions, latticeSpacing, regionName, latticeType, dt)` | Setup | Create the simulation object |
| `sim.simulationTime = t` | Setup | Set total simulation time (s) |
| `sim.latticeWriteInterval = n` | Setup | Set lattice write interval (timesteps) |
| `sim.speciesWriteInterval = n` | Setup | Set species write interval (timesteps) |
| `RegionBuilder(sim)` | Geometry | Create geometry helper |
| `build.ellipsoid(radius, center, angles)` | Geometry | Sphere / ellipsoid boolean mask |
| `build.dilate(mask, se)` | Geometry | Expand mask outward by one step |
| `build.compose((reg1, mask1), ...)` | Geometry | Write masks into site lattice (later overwrites earlier) |
| `build.se6` / `build.se26` | Geometry | 6- or 26-connected structuring element |
| `sim.region(name)` | Regions | Create or retrieve a region |
| `sim.species(name)` | Species | Create or retrieve a species |
| `sim.rateConst(name, value, order)` | Reactions | Create a named rate constant |
| `region.addReaction(reactants, products, rate)` | Reactions | Add reaction restricted to one region |
| `sim.reaction(reactants, products, rate, regions=[...])` | Reactions | Add reaction assigned to one or more regions |
| `sim.diffusionConst(name, value)` | Diffusion | Create a named diffusion constant (m²/s) |
| `sim.transitionRate(sp, rFrom, rTo, rate)` | Diffusion | Assign diffusion rate for a (sp, from, to) triple |
| `sim.diffusionZero` | Diffusion | Built-in zero diffusion constant |
| `sim.diffusionFast` | Diffusion | Built-in maximum-allowed diffusion constant |
| `sim.maxDiffusionRate()` | Diffusion | Compute the timestep upper bound for D |
| `sim.placeNumber(sp, x, y, z, n)` | Particles | Place n particles at a specific voxel |
| `sim.distributeNumber(sp, reg, count)` | Particles | Distribute particles uniformly across a region |
| `sim.distributeConcentration(sp, reg, conc)` | Particles | Place particles at a molar concentration |
| `sim.finalize()` | Run | Write model to `.lm` file — required before `run()` |
| `sim.run(solver, replicate, seed, cudaDevices, ...)` | Run | Launch the RDME solver |
| `sim.displayGeometry()` | Inspect | 3-D interactive lattice viewer |
| `sim.showAllSpecies()` | Inspect | Summary table of species, reactions, diffusion |
| `RDMEFile(fname, replicate)` | Analysis | Open a `.lm` result file |
| `traj.getNumberTrajectory(species, ...)` | Analysis | Whole-system particle count time series (fast) |
| `traj._speciesIdParse(None, name, startIndex=1)` | Analysis | Resolve species name → numeric lattice ID |
| `traj._regionIdParse(None, name)` | Analysis | Resolve region name → numeric lattice ID |
| `LMLattice.latticeStatsAll_h5fmt(plattice, slattice)` | Analysis | Count particles by region and species across a single lattice frame |
