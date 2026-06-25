# Chemical Kinetics in Spatially Heterogeneous Systems #
Our analysis of chemical kinetics now continues with systems that are treated as spatially heterogeneous. In contrast to well-mixed systems and the assumption that the probability of chemical reactions taking place does not depend on the locations of molecules within the system, we will now consider spatial heterogeneity. The inclusion of spatial heterogeneity in our treatment of a chemical system can take on many forms (e.g., continuous vs voxelized, see [Earnest et al., 2018](https://doi.org/10.1088/1361-6633/aaae2c)). Here, we will subdivide our system volume into $N$ number of subvolumes, thereby "voxelizing" the chemical system. Then, we will treat each subvolume as though it is well mixed, thereby enabling us to use techniques similar to the stochastic techniques we used previously. The differences here are as follows: 

1. The system is split into an arbitrary number of subsystems (or "voxels").
2. Each voxel, rather than the whole system, is treated as being well-mixed.

As stated above, there are other ways to treat spatially heterogeneous systems stochastically such as particle-based methods, however, Lattice Microbes does not implement these methods. These methods can be useful for systems where more spatial accuracy is necessary, but are more difficult to simulate for longer periods of time. This becomes important when creating whole-cell models because we often want to simulate the entirety of the cell cycle.

As with spatially-homogeneous systems, the assumption that a system is spatially-heterogeneous does not tell us whether we will treat the chemical events as being deterministic or stochastic. Both types of treatments are possible and are briefly described below.

## 1. Deterministic Chemical Kinetics in Spatially Heterogeneous Systems ##
Deterministic chemical kinetics in spatially heterogeneous systems are often modeled using systems of partial differential equations (PDEs), where each chemical species in the system receives its own individual equation. This type of treatment is similar to using ODEs to solve a deterministic spatially homogeneous system, but now includes information about the concentration of molecules across space and how these concentrations change over time. We will not be covering this type of treatment in our present tutorials, but if you are interested in learning more about these methods, the following [Wikipedia page](https://en.wikipedia.org/wiki/Method_of_lines) and [Sympy tutorial](https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/60-chemical-kinetics-reaction-diffusion.html) may be helpful.

## 2. Stochastic Chemical Kinetics in Spatially Heterogeneous Systems ##
Alternatively, stochastic chemical kinetics in spatially heterogeneous systems can be treated using a specific form of the generalized master equation called the reaction-diffusion master equation (RDME): 

```math
\begin{aligned}
\frac{dP(\mathbf{x}, t)}{dt} 
&= \mathbf{R} P(\mathbf{x}, t) + \mathbf{D} P(\mathbf{x}, t) \\
&= \sum_{v}^{V} \sum_{r}^{R} 
\Big[ -a_r(\mathbf{x}_v) P(\mathbf{x}_v, t) 
+ a_r(\mathbf{x}_v - \mathbf{S}_r) P(\mathbf{x}_v - \mathbf{S}_r, t) \Big] \\
&\quad + \sum_{v}^{V} \sum_{\xi}^{i,j,k} \sum_{\alpha}^{N} 
\Big[ -d^\alpha x^\alpha_v P(\mathbf{x}, t) 
+ d^\alpha \big( x^\alpha_{v+\xi} + 1 \big) 
P\big( \mathbf{x} + 1^\alpha_{v+\xi} - 1^\alpha_v, t \big) \Big]
\end{aligned}
```

For a full treatment of this equation see either [Earnest et al., 2018](https://doi.org/10.1088/1361-6633/aaae2c) or [Isaacson & Isaacson, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3405976/). In essence, this equation extends the previously-used CME by adding terms for diffusion across voxels in the x, y, and z directions.

As with deterministic treatments, analytical solutions to the RDME for anything but very simple systems become intractable. Here, we will also resort to using numeric solvers. Lattice Microbes uses what is called the multiparticle diffusion (MPD) operator to enable numerical solutions to the RDME, as discussed in [Roberts et al., 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3762454/). In essence, this operator treats diffusion in each spatial direction as independent, and for each particle in each voxel, it computes the probability of diffusing forward or backward in the given direction. After each spatial direction diffusion kernel is performed, the lattice is updated with new particle positions and the kernel for the next spatial direction is computed. After the diffusion kernels for all spatial directions have been finished, the algorithm performs the Gillespie algorithm on each voxel, independently. As a result, we are able to model both diffusion and reaction processes using the same algorithm.

Similar to the trajectory solutions we obtained from the spatially homogeneous treatment of chemical systems, the solution to the spatially heterogeneous system is a set of species trajectories (abundances over time). However, because we have also included spatial considerations in these systems, we can decompose these trajectories both by species type and by any arbitrary region of the system. For example, we may want to know how many molecules of species A there are across time *in the right half* of the system. Thinking biologically, we may want to partition our system into various compartments/organelles such as the cytosol, plasma membrane, nucleus, etc. Then, we can generate trajectories of each chemical species, *in each compartment*. However, to do this, we must first specify the architecture or geometry of our system

## 3. Including Architectures or Geometries in Chemical Systems ##
In Lattice Microbes, we begin designing our spatial architecture with a uniform cube-like structure, subdivided into evenly-spaced voxels in each direction. The number of voxels in each direction and the length of each voxel side are set by the user (see constraints on lattice spacing, diffusion coefficients, and time steps in Section 5). This system cube is stored as a 3D array, with each entry in the array representing a single subvolume or voxel. Next, each discretized subvolume can be assigned a predefined region type. The most simple architecture one can create is a system in which every subvoxel is of the same region type. Regions can be specified by the user and can be named arbitrarily. Subvoxels are assigned a default region, but the user can create additional 3D arrays for each added region, with the array specifying whether a given subvoxel belongs to the respective region type. Eventually, each region type receives its own 3D array indicating whether a voxel is considered part of that region type. Regions are then sequentially added to the RDME simulation object with overlapping regions being set to the most recently added region type.

Region types are used to restrict chemical reactions and specify diffusion coefficients for individual chemical species. Each reaction that is added to the RDME simulation object must specify the regions in which this reaction is allowed to take place. Some reactions may be permitted in all regions and other in only one or two predetermined regions. A biological example of this would be a reaction being allowed to take place in the nucleus, but not the cytosol (e.g., transcription). Region types also enable the user to define diffusion coefficients for each species with a defined startng and ending region type. For example, if species S can diffuse freely in region A, but is not able to enter region B, one can set the diffusion coefficient of this species in for A-A region diffusion to be some positive number, and for A-B region diffusion to be zero. A biological example of this may be a protein that can freely diffuse in the cytosol, but is not allowed to enter the membrane. More information on diffusion rules will be given below (see Section 5). Finally, region types allow the user to more easily visualize the chemical system and its unique compartments using visualization software such as VMD.

System architectures can be created arbitrarily using idealized geometric shapes (see [RegionBuilderCheatSheet.md](./RegionBuilderCheatSheet.md)). However, if one has experimental data such as tomograms from a cryo-ET experiment, Lattice Microbes allows the incorporation of this data into system geometries. All the user must do is convert these experimental data into boolean 3D arrays of a size corresponding to the chemical system, and then supply these arrays to the simulation object for a given region type.

## 4. Defining Chemical Species, Initial Abundances, and Species Locations ##
After building the system geometry, and as with the spatially homogeneous treatment of chemical systems, we must define all chemical species in the system. This is done in a similar way to spatially homogeneous systems, however, when defining the initial abundances of these species in spatially heterogeneous systems, we must also supply information about their locations. In Lattice Microbes, the default option for placement of chemical species is in a uniform random distribution across the region type of choice. However, we are also able to manually select the initial abundance of each species in each subvoxel of our system if this information is known. 

Lattice Microbes solvers are preprogramed with the number of particles that are allowed to occupy a single lattice site. In the MPD solver for spatially heterogeneous systems, the standard number of particles allowed is 16. If, due to diffusion or reaction events, the particle number exceeds this for any given lattice, Lattice Microbes will flag an "overflow" event, and will then distribute the particles randomly to neighboring voxels of the same region type.

## 5. Defining Diffusion Rules for Each Chemical Species ##
After adding chemical species along with their initial abundances and locations, we need to define each species' diffusion rules. In contrast to molecular dynamics simulations, where the movement of particles is governed by force fields and energy potentials, Lattice Microbes treats diffusion as a Brownian process where a general diffusion coefficient is given and probabilities for diffusion events are determined using the following equation:


```math
q_- = \frac{D_\alpha(s_\text{current} \to s_\text{minus}) \cdot \tau}{\lambda^2}
```

```math
q_+ = \frac{D_\alpha(s_\text{current} \to s_\text{plus}) \cdot \tau}{\lambda^2}
```

```math
p_\text{stay} = 1 - q_- - q_+
```

with the following variable definitions:

| Symbol | Meaning |
|--------|---------|
| $\lambda$ | Lattice spacing (voxel edge length) |
| $\tau$ | Simulation timestep |
| $s, s_{minus}, s_{plus}$ | Source and destination site types |
| $\alpha$ | Chemical species index |
| $D_\alpha(s \to s')$ | Diffusion coefficient for species $\alpha$ from site type $s$ to $s'$ |
| $q_\pm$ | Hop probability in the ± direction along the current axis |
| $p_\text{stay}$ | Probability of not hopping: $1 - q_- - q_+$ |



Lattice Microbes performs these computations for each particle, in each voxel, in each spatial direction, at every timestep.

Importantly, each chemical species will need a generalized diffusion coefficient to be supplied by the user for all possible region type transitions. For a system containing a single region type (region A), the only diffusion coefficient required is from a region A voxel to another region A voxel. For a system with two region types (region A and region B), a diffusion coefficient is required for A to A, A to B, B to A, and B to B transitions. Users can prohibit particles of specific chemical species from diffusing into specific region types by setting all diffusion coefficients into this region type to zero.


## 6. Defining the Chemical Reaction Network ##
As with spatially homogeneous systems, all possible chemical reactions must be specified by the user. However, in spatially heterogeneous systems, we are also allowed to specify which region types a given chemical reaction is allowed to take place. Biologically, this may mean that we can allow a reaction to take place in the nucleus, but not in the cytosol or plasma membrane. Inherent in the spatially heterogeneous treatment of a chemical system is the fact that for a reaction to take place, the required reactant particles must exist in the same voxel.

An important note about spatially heterogeneous simulations in Lattice Microbes is that only a single chemical reaction can take place within each voxel across the duration of a single timestep. Lattice Microbes computes the overall reaction propensity for each voxel to determine if a reaction will take place at any given timestep. Then, if a reaction does take place, a single reaction is randomly picked weighted by the given reaction propensities within the voxel.

## 7. Defining the Simulation Parameters ##
Finally, we must specify various system parameters for a spatially heterogeneous system. First, as described above, we must specify the number of voxels to subdivide our system into along each axis and the physical length of the voxel edges. The length scales applied to voxel edges are arbitrary, thereby allowing the user to create a system as physically small or large as desired. The overall size of the system therefore does not dictate how many subvoxels into which the system is divided, however, each length dimension must be divided into a multiple of 32 voxels due to GPU constraints. 

We must also specify the duration of the timestep taken between each set of diffusion-reaction kernel computations. The choice in timestep is not arbitrary, but must satisfy the following equation:

```math
\tau \leq \frac{\lambda^2}{2\, D_\text{max}}
```

where $\tau$ is the timestep duration, $\lambda$ is the edge length of each subvoxel, and $D_{max}$ is the largest diffusion coefficient within the simulation system. Lattice Microbes will not allow the user to specify a timestep larger than dictated by the above equation. The reasoning behind this constraint is that because we only allow particles to diffuse a single voxel in any spatial direction during the length of a timestep, we do not want to allow timesteps to become large enough to physically allow a particle to diffuse across multiple subvoxels. Using a large timestep and limiting the particle diffusion to a single voxel in each direction becomes physically unrealistic.

In addition to the timestep, we must also specify the total amount of biological time we want to simulate, the lattice write interval and the species write interval. The lattice write interval is a parameter that defines the number of simulation timesteps that take place before writing the state of the system lattice to the results file. The lattice state gives a description of the number of every species at every subvoxel in the system. This type of data is often used either for region-specific abundances or for visualization of the system trajectories with software such as VMD. The species write interval is a parameter that defines the number of simulation timesteps that take place before writing the overall abundance of each species to the results file. This is a more coarse-grained representation of the system state, and takes up much less memory when it is saved to the results file. Therefore, it is common to have the species interval be substantially smaller than the lattice interval, especially for large or complex systems.

## 8. Using Lattice Microbes to Model Stochastic Chemical Kinetics in Spatially Heterogeneous Systems ##
As with the spatially homogeneous system treatment, **Figure 1** gives a simple overview of how to setup a spatially heterogeneous system in Lattice Microbes. In our current tutorial, we will not discuss custom events. This topic will be covered in tutorial 3 where we give a brief overview of how to use Lattice Microbes to run whole-cell model simulations. The specific syntax to complete each of these steps will be described within the tutorials and a comprehensive overview of all RDME-related functions used in the tutorials can be found in the [RDME_FunctionReference](./RDME_FunctionReference.md) document.

<p align="center">
  <img src="./SupplementaryFigures/jLM_Flowchart.png" width="450" alt="Schematic diagram of the LM architecture"> <br>
  <b>Figure 1. Lattice Microbes Flowchart</b>
</p>