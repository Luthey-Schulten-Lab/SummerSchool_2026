# A Well-Mixed Bimolecular Reaction System Solved Both Deterministically and Stochastically #

## Introduction ##
In tutorial 1, we will use deterministic and stochastic methods to solve for the trajectories of each chemical species in a bimolecular reaction. We will use the following reaction scheme as our system:

```math
A + B \leftrightarrow C
```

This system includes two reactions, one forward reaction where A and B combine to form C, and one reverse reaction where C dissociates into A and B. 

In Tutorial 1.1, we will follow a standard deterministic treatment of this system and will use the following system of equations:

```math
\frac{d[A]}{dt} = k_{r}[C] - k_{f}[A][B]
```

```math
\frac{d[B]}{dt} = k_{r}[C] - k_{f}[A][B]
```

```math
\frac{d[C]}{dt} =  k_{f}[A][B] - k_{r}[C]
```

Because we are primarily interested in modeling biological systems, we will start with small initial molecular abundances; A: 100, B: 100, and C: 0. To model this deterministically, we will convert convert these counts to concentrations and use standard concentration-based reaction coefficients. Here, we will set $k_f$=1.07 $\times10^6$ M<sup>-1</sup>s<sup>-1</sup> and $k_r$=0.351 s<sup>-1</sup>. Continuing with the theme of modeling biological processes, we will use a reaction volume of 1 fL (1e-15 L), the size of an average *E. coli* cell.

In contrast, in Tutorial 1.2, we will treat the same exact chemical system but use stochastic methods. When modeling a system with stochastic methods, we usually work with discrete molecular counts instead of concentrations. To make the forward rate constant compatible with molecule-based stochastic simulation, we convert it by dividing by Avogadro’s number and the system volume, yielding units of reactions per molecule per second.

## How to Run the Jupyter Notebooks ##

You will use the notebook `Tut1.1-ODEBimol.ipynb` to simulate the bimolecular reaction using ordinary differential equations (ODEs) solved with SciPy, and `Tut1.2-CMEBimol.ipynb` to simulate the same reaction using the Chemical Master Equation (CME) with jLM in Lattice Microbes.

To run a Jupyter Notebook, you can use short-cut **Control+Enter** to execute the selected cell. Or, you can click **Cell-Run All** to run the entire script.

## Using Lattice Microbes to Solve a Well-Mixed Chemical System Stochastically - Outline of Tutorial 1.2 ##

The following section can be used as a reference guide for setting up and solving well-mixed chemical systems stochastically in Lattice Microbes.

As introduced earlier, the jLM Problem Solving Environment is now widely used to construct CME systems. To get started, we first import the necessary `jLM` modules.

Next, we define the system parameters which usually consist of the system volume and Avogadro's number. 

Next, we define our empty simulation object with the line `sim = CME.CMESimulation()`. The next lines define the chemical species. In `jLM.CME`, species are represented by Python strings and must be registered with the simulation using the `defineSpecies` command. The following lines with `addParticles` define the initial counts of each species accumulatively. 

Rate constants are then defined and converted into particle-based values. Reactions are then added using the `addReaction` function. In Lattice Microbes, both the forward and reverse reactions must be specified separately. The first and second arguments of `addReaction` can be either a tuple of reactants or a single string if only one reactant is involved. `jLM.CME` currently supports 0th-, 1st-, and 2nd-order reactions, and reaction rate constants must be given in stochastic units. For 0th-order reactions, use the empty string `""` as the reactant. Similarly, annihilation reactions can be specified by passing `""` as the product.

Next, we define the simulation parameters. In this tutorial, we write simulation output every 30 microseconds and run the simulation for a total of 30 seconds. We will perform 10 replicates using the same initial conditions. An important detail: the simulation must be saved to a file before it can be executed.

Finally, we call the `run(...)` command on the simulation object, specifying the filename, the simulation method, and the number of independent trajectories (replicates) to run. In this tutorial, we use the **Direct Gillespie Algorithm** to sample the stochastic dynamics of the bimolecular reaction system.

When the simulation runs, you will see a standard output listing the number of completed replicates. In general, CME simulations finish relatively quickly.

We use the built-in module `jLM.CMEPostProcessing` for basic post-processing tasks. These include plotting time-dependent trajectories for individual replicates and calculating population averages and variances across replicates.

Additionally, we also demonstrate a more flexible analysis approach: the simulation output is first serialized into a 3D NumPy array (with dimensions of species, time, and replicates) and then visualized using custom plotting functions. You will apply this method in CME tutorials 2 and 3 for more advanced analyses.

## Further Investigation ##
### 1. How do deterministic and stochastic treatments of the same chemical system differ from one another? ###
Once you have completed both tutorials, you'll notice that the counts of each species in the deterministic treatment change smoothly over time. However, in reality, molecule counts must be integers due to the discrete nature of matter. Since reactions occur through collisions between individual molecules, the changes in molecular counts also happen in integer steps. For systems with low molecule counts, this discreteness becomes significant, and ODEs no longer provide accurate descriptions. Stochastic modeling was developed to capture this behavior.

You may notice that the behavior in the stochastic (CME) results is qualitatively similar to the deterministic ODE results. However, the CME trajectories exhibit considerable fluctuations, even after the system reaches equilibrium. These fluctuations arise from the intrinsic randomness of the process, where reactions can transiently deviate from the equilibrium state.

Additionally, you might observe that the particle counts change in discrete steps between time points. This feature becomes even more apparent when molecule counts are lower.

 - At what initial copy numbers do the behaviors of the deterministic and stochastic treatment begin to look fairly identical?
 - Is the final behavior of the deterministic treatment significantly different from the stochastic treatment? How could you test this?
 - If it is different, at what molecule counts does this difference disappear?

### 2. How does sample size affect the average behavior of stochastically treated species trajectories? ###
In Tutorial 1.2, try increasing the number of replicates from 10 to 100—or even more. To do this, change the variable `reps` and restart the Jupyter Notebook kernel by clicking **Kernel-Restart & Clear Output** to begin a new CME simulation. 

 - What measurement quantities are affected by the change in sample size? Which quantities are not affected?


### 3. How do the magnitudes of reaction rates affect the average behavior of stochastically treated species trajectories? ###

Multiply both the forward and backward rate constants by 10 or 100 by modifying the variable `fold`. Then, restart and rerun both the ODE and CME Jupyter Notebooks.

- In the ODE simulation, does the system reach equilibrium more quickly? Does the equilibrium concentration change?
- In the CME simulation, do the fluctuations in a single replicate occur more rapidly? How does the ensemble average respond?
