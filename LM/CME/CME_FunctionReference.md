# CME Function Reference — Lattice Microbes / jLM #

This reference covers every function used across the CME tutorial notebooks. Functions are organized by workflow stage so you can follow the same order when building your own simulations.

---

## Workflow Overview ##

```
1. Create simulation object  →  CMESimulation()
2. Define species            →  defineSpecies(), addParticles()
3. Define reactions          →  addReaction()
4. Set parameters            →  setSimulationTime(), setWriteInterval()
5. Save to disk              →  sim.save()
6. Run                       →  sim.run()
7. Analyze                   →  openLMFile(), getTimesteps(), plotTrace(), ...
```

---

## 1. Setting Up the Simulation ##

### `CMESimulation` ###

Constructs the central simulation object that holds all species, reactions, and parameters.

```python
sim = CME.CMESimulation(volume=None, name="unnamed")
```
+ `volume` (*Float, optional*): reaction vessel volume in **litres**; when provided, jLM automatically converts macroscopic rate constants to stochastic units — you do not need to divide by $N_A$ or $V$ yourself; pass `None` only if you have already done this conversion manually;
+ `name` (*String, optional*): human-readable label stored on the object as `sim.name`;

**Returns:** a `CMESimulation` instance. Key attributes created:
- `sim.particleMap` — `dict` mapping species name → 1-based integer index
- `sim.species_id` — ordered list of species names
- `sim.initial_counts` — `dict` of initial particle counts (all start at 0)

**Usage:**
```python
sim = CME.CMESimulation(name='Bimolecular Reaction')
```

---

## 2. Defining Species ##

### `sim.defineSpecies` ###

Registers one or more chemical species in the simulation and initialises their counts to zero.

```python
sim.defineSpecies(species)
```
+ `species` (*list of Strings*): species names to register; call once with a list, or multiple times for individual names; each name must be unique;

**Usage:**
```python
sim.defineSpecies(['A', 'B', 'C'])
sim.defineSpecies(['gene', 'mRNA', 'protein'])
```

---

### `sim.addParticles` ###

Sets the initial particle count for a previously defined species.

```python
sim.addParticles(species='unknown', count=1)
```
+ `species` (*String*): name of an already-defined species;
+ `count` (*Int*): number of molecules at time zero; calling this function multiple times on the same species accumulates the total count;

**Usage:**
```python
sim.addParticles(species='A', count=50)
sim.addParticles(species='gene', count=1)
```

---

## 3. Defining Reactions ##

### `sim.addReaction` ###

Adds a chemical reaction to the model. Rate constants must be in **stochastic (particle-based) units** — i.e., s⁻¹ for all orders — unless `volume` was passed to the constructor (in which case jLM converts macroscopic units automatically).

```python
sim.addReaction(reactant, product, rate)
```
+ `reactant` (*String, 2-tuple of Strings, or* `''`): reactant species; use `''` for a 0th-order production reaction; use a 2-tuple `('A', 'B')` for a bimolecular reaction; use the same name twice `('A', 'A')` for a self-reaction;
+ `product` (*String, tuple of Strings, or* `''`): product species; use `''` for a pure degradation reaction with no products; use a tuple `('gene', 'mRNA')` for multiple products;
+ `rate` (*Float*): stochastic rate constant in **s⁻¹**;

| Reaction type | `reactant` form | `product` form | Example |
|---|---|---|---|
| 0th order (production) | `''` | `'A'` | `addReaction('', 'A', k)` |
| 1st order (conversion) | `'C'` | `('A', 'B')` | `addReaction('C', ('A','B'), kr)` |
| 1st order (degradation) | `'mRNA'` | `''` | `addReaction('mRNA', '', k_deg)` |
| 2nd order | `('A', 'B')` | `'C'` | `addReaction(('A','B'), 'C', kf)` |
| 2nd order (self) | `('A', 'A')` | `'B'` | `addReaction(('A','A'), 'B', k)` |

**Usage:**
```python
# Bimolecular association
sim.addReaction(reactant=('A', 'B'), product='C',        rate=kf)
# Dissociation
sim.addReaction(reactant='C',        product=('A', 'B'), rate=kr)
# Transcription (gene acts as catalyst)
sim.addReaction(reactant='gene', product=('gene', 'mRNA'), rate=k_trans)
# mRNA degradation
sim.addReaction(reactant='mRNA', product='',               rate=k_deg_mRNA)
```

---

## 4. Setting Simulation Parameters ##

### `sim.setSimulationTime` ###

Sets the total biological simulation time.

```python
sim.setSimulationTime(time)
```
+ `time` (*Float*): duration of simulation in **seconds**;

**Usage:**
```python
sim.setSimulationTime(30.0)    # Bimolecular tutorial: 30 s
sim.setSimulationTime(6300)    # GIP tutorial: one full cell cycle (~6300 s)
```

---

### `sim.setWriteInterval` ###

Sets how frequently the system state (particle counts) is sampled and written to the output file.

```python
sim.setWriteInterval(time)
```
+ `time` (*Float*): interval between writes in **seconds**; smaller values give finer temporal resolution but produce larger output files;

**Usage:**
```python
sim.setWriteInterval(30e-6)                 # every 30 µs
sim.setWriteInterval(units.microsecond(30)) # equivalent using the units helper
sim.setWriteInterval(1)                     # every 1 s
```

---

## 5. Saving the Simulation ##

### `sim.save` ###

Serialises the complete model (species, reactions, parameters) to an HDF5-based `.lm` file.

```python
sim.save(filename)
```
+ `filename` (*String*): path to the output `.lm` file; the file is created fresh, so delete any pre-existing file first;

**Note:** always remove any pre-existing file before saving to avoid appending to stale data.

**Usage:**
```python
import os
os.system("rm -rf %s" % filename)
sim.save(filename)
```

---

## 6. Running the Simulation ##

### `sim.run` ###

Executes the stochastic simulation and writes each replicate into the `.lm` file.

```python
sim.run(filename, method, replicates=1, seed=None, cudaDevices=None, checkpointInterval=0)
```
+ `filename` (*String*): path to the `.lm` file created by `sim.save()`;
+ `method` (*String*): solver to use; for the Gillespie SSA use `"lm::cme::GillespieDSolver"`;
+ `replicates` (*Int, optional*): number of independent stochastic replicates to run sequentially; each is seeded differently (default: `1`);
+ `seed` (*Int, optional*): fixed random seed for reproducibility; `None` uses time-based seeding (default: `None`);
+ `cudaDevices` (*list of Int, optional*): CUDA device indices to run on; defaults to `[0]`;
+ `checkpointInterval` (*Int, optional*): interval at which checkpoints are written; `0` disables checkpointing (default: `0`);

**Usage:**
```python
sim.run(filename=filename, method="lm::cme::GillespieDSolver", replicates=100)
```

---

## 7. Analysis ##

All post-processing functions are accessed via `jLM.CMEPostProcessing` (commonly imported as `PostProcessing`).

---

### `PostProcessing.openLMFile` ###

Opens an LM output file for reading. The returned handle must be passed to all subsequent analysis functions.

```python
fileHandle = PostProcessing.openLMFile(filename)
```
+ `filename` (*String*): path to the `.lm` HDF5 output file;

**Returns:** `h5py.File` handle (read-only). Always call `closeLMFile` when done.

---

### `PostProcessing.closeLMFile` ###

Closes an open file handle. Always call this when you are finished reading.

```python
PostProcessing.closeLMFile(fileHandle)
```
+ `fileHandle` (*h5py.File*): handle returned by `openLMFile`;

---

### `PostProcessing.getTimesteps` ###

Extracts the array of time points at which data was recorded.

```python
timestep = PostProcessing.getTimesteps(fileHandle)
```
+ `fileHandle` (*h5py.File*): open file handle;

**Returns:** `numpy.ndarray` of shape `(n_timepoints,)` in **seconds**.

---

### `PostProcessing.getNumTrajectories` ###

Returns the number of completed replicates stored in the file.

```python
n = PostProcessing.getNumTrajectories(fileHandle)
```
+ `fileHandle` (*h5py.File*): open file handle;

**Returns:** `int`.

---

### `PostProcessing.getSpecieTrace` ###

Extracts the particle-count time trace for a single species.

```python
trace = PostProcessing.getSpecieTrace(fileHandle, specie, replicate=None)
```
+ `fileHandle` (*h5py.File*): open file handle;
+ `specie` (*String*): species name;
+ `replicate` (*Int, optional*): 1-indexed replicate number; if `None`, returns a `dict` keyed by replicate number;

**Returns:** `numpy.ndarray` of shape `(n_timepoints,)` for a single replicate, or `dict[int, ndarray]` for all replicates.

---

### `PostProcessing.getAvgVarTrace` ###

Computes the mean and variance of a species' particle count across all replicates at every time point.

```python
avg, var, times = PostProcessing.getAvgVarTrace(fileHandle, specie)
```
+ `fileHandle` (*h5py.File*): open file handle;
+ `specie` (*String*): species name;

**Returns:** three `numpy.ndarray` objects of shape `(n_timepoints,)` — mean counts, variance, and time points.

---

### `PostProcessing.plotTrace` ###

Plots raw particle-count trajectories for specified species from a single replicate.

```python
PostProcessing.plotTrace(fileHandle, species=None, replicate=1, filename=None, **kwargs)
```
+ `fileHandle` (*h5py.File*): open file handle;
+ `species` (*list of Strings, optional*): species names to plot;
+ `replicate` (*Int, optional*): 1-indexed replicate number (default: `1`);
+ `filename` (*String, optional*): if provided, saves the figure to this path; `None` displays inline in the notebook;
+ `**kwargs`: passed directly to `matplotlib.pyplot.plot`;

**Returns:** `matplotlib.figure.Figure`.

**Usage:**
```python
PostProcessing.plotTrace(fileHandle, species=['A', 'C'], replicate=5)
```

---

### `PostProcessing.plotAvgVar` ###

Plots the mean and variance of each species over time across all replicates, in a two-panel figure.

```python
PostProcessing.plotAvgVar(fileHandle, species=None, filename=None, **kwargs)
```
+ `fileHandle` (*h5py.File*): open file handle;
+ `species` (*list of Strings, optional*): species names to include;
+ `filename` (*String, optional*): output path for saving; `None` displays inline;
+ `**kwargs`: passed to `matplotlib.pyplot.plot`;

**Returns:** `matplotlib.figure.Figure`.

**Usage:**
```python
PostProcessing.plotAvgVar(fileHandle, species=['A', 'C'])
```

---

## 8. Unit Conversion Helpers ##

The `jLM.units` module provides convenience functions that convert common physical units to SI base units (metres or seconds). Each function accepts one or more values and returns the converted result.

```python
import jLM.units as units
value_in_seconds = units.microsecond(30)   # → 3e-5
```

| Function | Input unit | SI output |
|---|---|---|
| `units.angstrom(*qty)` | Å | metres |
| `units.nm(*qty)` | nanometres | metres |
| `units.micron(*qty)` | micrometres | metres |
| `units.mm(*qty)` | millimetres | metres |
| `units.cm(*qty)` | centimetres | metres |
| `units.ns(*qty)` | nanoseconds | seconds |
| `units.microsecond(*qty)` | microseconds | seconds |
| `units.ms(*qty)` | milliseconds | seconds |
| `units.second(*qty)` | seconds | seconds |
| `units.minute(*qty)` | minutes | seconds |
| `units.hr(*qty)` | hours | seconds |
| `units.day(*qty)` | days | seconds |

**Note:** the function name is `microsecond` (singular), not `microseconds`.

---

## Quick-Reference Table ##

| Function | Stage | What it does |
|---|---|---|
| `CME.CMESimulation(volume, name)` | Setup | Create the simulation object |
| `sim.defineSpecies(species)` | Species | Register species names |
| `sim.addParticles(species, count)` | Species | Set initial particle counts |
| `sim.addReaction(reactant, product, rate)` | Reactions | Add a reaction to the model |
| `sim.setSimulationTime(time)` | Parameters | Set total simulation duration (s) |
| `sim.setWriteInterval(time)` | Parameters | Set data write frequency (s) |
| `sim.save(filename)` | Save | Write model to `.lm` file |
| `sim.run(filename, method, replicates, ...)` | Run | Execute the simulation |
| `PostProcessing.openLMFile(filename)` | Analysis | Open output file for reading |
| `PostProcessing.closeLMFile(f)` | Analysis | Close the output file |
| `PostProcessing.getTimesteps(f)` | Analysis | Get array of recorded time points |
| `PostProcessing.getNumTrajectories(f)` | Analysis | Get number of completed replicates |
| `PostProcessing.getSpecieTrace(f, specie, replicate)` | Analysis | Get count trace for one species |
| `PostProcessing.getAvgVarTrace(f, specie)` | Analysis | Get mean and variance across replicates |
| `PostProcessing.plotTrace(f, species, replicate, ...)` | Analysis | Plot raw trajectory for one replicate |
| `PostProcessing.plotAvgVar(f, species, ...)` | Analysis | Plot mean ± variance across replicates |
| `units.microsecond(qty)` | Units | Convert µs → s |
