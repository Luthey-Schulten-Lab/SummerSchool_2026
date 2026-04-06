# Simulating DNA with LAMMPS

## Description:
<img align="right" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/spotlight.png">

We will walk you through how to set up and run a LAMMPS simulation using GPUs on the Delta HPC cluster. We will simulate the DNA dynamics of the Minimal Cell JCVI-syn3A, including DNA replication, disentanglement of daughter chromosomes, and partitioning of daughter chromosomes into their respective daughter volumes. The coarse-grained model of the DNA, ribosomes and cell membrane will be discussed, as well as the use of LAMMPS to perform energy minimizations and Brownian dynamics. We will also go into greater detail about how we model biological mechanisms such as topoisomerase-induced strand crossing and SMC looping. After the runs complete, you will get a chance to visualize your trajectory in VMD.

*This tutorial was prepared for the second edition of the STC QCB Summer School, held July 21-25, 2025.*

## Outline of tutorial:

1. Introduction to DNA simulation with LAMMPS
2. **Setting up and submitting your job to Delta (do this on Monday)**
3. Generating an initial structure
4. Modeling DNA replication
5. Modeling chromosome dynamics
6. A closer look at SMC dynamics
7. Understanding btree_chromo commands
8. Visualization with VMD

Most of the content of this tutorial, including the implementation of energy terms for the DNA polymer, DNA disentanglement, and general procedure for simulating Brownian dynamics and energy minimization with LAMMPS on a GPU, is also explained in our recent manuscript[^thornburg2025] which you can check out on bioRxiv. The content on SMC blocking/bypassing and daughter chromosome partitioning without the need for an additional fictitious force is a work in progress. 

## 1. Introduction to DNA simulation with LAMMPS

Here, we simulate DNA replication and dynamics using [LAMMPS](https://www.lammps.org/#gsc.tab=0) (Large-scale Atomic/Molecular Massively Parallel Simulator), a molecular dynamics program from Sandia National Laboratories. We will not have to worry about writing our own LAMMPS input scripts. Instead, we will be running the C++ program `btree_chromo`, available online at https://github.com/Luthey-Schulten-Lab/btree_chromo_gpu/tree/btree_chromo_gpu_SummerSchool2025. This program was created mainly for the purposes of simulating the minimal cell chromosome, but it can be used to simulate any circular chromosome. The main purpose of the program is to model replication states of the chromosome, as well as perform simulation of chromosome dynamics by calling LAMMPS. 

<img align="center" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/DNA_model_0.png">

The DNA that btree_chromo simulates is coarse-grained at a 10 bp resolution. This means that a single, 3.4 nm diameter bead is used to represent 10 base pairs. We also use beads to represent ribosomes (10 nm), and the cell membrane. The program emulates the effects of SMC (structural maintenence of chromosomes) proteins that extrude loops of DNA to effect chromosome organization, as well as type II topoisomerases which allow for DNA strand crossings when they become tangled.

Today you will run a simulation using a variant of LAMMPS which utilizes the GPUs on the Delta HPC cluster. We will simulate the cell cycle of the minimal cell including the effects of SMC proteins, topoisomerase, and Brownian dynamics. We will start by generating an initial configuration for the DNA and ribosomes of the minimal cell in a spherical cell membrane. The DNA will replicate, disentangle, and partition, and the cell membrane will grow and divide. At the end of the simulation we should have two cells that each look roughly like the cell we started with.

At the core of the simulation we use LAMMPS for simulating the DNA dynamics, but their are several layers code wrapped around LAMMPS that are designed to make our lives easier. These layers are as follows: we will submit a slurm script that launches an Apptainer container. This container runs the fortran program to generate the initial DNA/ribosome coordinates, and then runs a python script which writes and executes 90 `btree_chromo` input scripts, each of which correspond to 1 biological minute of the Syn3A cell cycle. Each `btree_chromo` script involves writing and executing 6 LAMMPS input scripts (each of which correspond to 2 biological seconds of DNA replication and SMC looping).

At one level even deeper, LAMMPS uses Kokkos, a library that lets the same LAMMPS code run efficiently on different types of hardware, like AMD and NVIDIA GPUs. In our case, Kokkos lets us perform the force calculations for the energy minimizations and Brownian dynamics on the GPU. Running on the GPU is around an order of magnitude faster than running on the CPU, and some GPUs can be much faster than others - for example, the A100 GPUs on Delta are around 2.5 times as fast as the RTX A5000 GPUs on my office desktop computer.

## 2. Setting up submitting your job to Delta
In this section, we will log on to Delta and launch a container which has btree_chromo and LAMMPS already installed. Then, we will start running a simulation of the minimal cell chromosome. 

> [!NOTE]
The reason we are doing this first, is so that the simulation will be left running throughout the rest of day and tomorrow morning. The simulations take around 14 hours to run, even using the A100 GPUs on Delta which are very fast. Tomorrow morning/afternoon we can visualize the results of our simulations.
> 

**Step 1: Log in to Delta**

Connect to the Delta HPC using the following command. You will need to type your password and do 2FA.

```bash
ssh $USERNAME@login.delta.ncsa.illinois.edu
```

> [!WARNING]
You will need to replace $USERNAME with your own.
> 

**Step 2: Copy the workspace folder and the launch scripts**

```bash
bash /projects/beyi/prelaunch_btree_chromo.sh

```

This bash script copies all the necessary files into your personal directory, `/projects/beyi/$USERNAME`.

**Step 3: Run the simulation**

We will run our btree_chromo simulation via a bash script that submits a job to Delta.

You will obtain output in the folder `/projects/beyi/${USER}/DNA_SummerSchool_2025`.

```bash
/projects/beyi/$USERNAME/
├── launch_simulation.sh
└── DNA_SummerSchool_2025/
    ├── data/
    │   ├── coords/
    │   ├── loops/
    │   ├── rep_states/
    │   └── run_name.lammpstrj  # etc.
    ├── scripts/
    │   ├── run_btree_chromo.py
    │   ├── template.inp
    │   ├── run_sc_chain_generation.sh 
    │   ├── Syn3A_chromosome_init.inp
    │   └── BD_lengths.txt
```

Run the command
```bash
cd /projects/beyi/${USER}/
sbatch launch_simulation.sh
```

If you would like to monitor the progress of your job, you can do the command

```bash
squeue -u ${USER}
```
which will show you your jobid, how long the job has been running, as well as what partition and GPU node it is running on.


## 3. Generating an initial structure 

At the start of every simulation, we need initial configuration, i.e. coordinates for the DNA, ribosomes, and cell membrane. Initial configurations for the chromosome are generated using a midpoint-displacement algorithm that creates three-dimensional, closed curves formed from overlapping spherocylinder segments. We assume a spherical cell with a known ribosome distribution (nearly randomly distributed, according to Cryo-ET), and "grow in" a self-avoiding chain of these spherocylinders. This process involves adding segments iteratively while avoiding overlaps with ribosomes and preventing knots. Spherical monomers are then interpolated along the spherocylinders. This method accurately models the "fractal globule" chromosome configuration present in Syn3A cells.

See the figure below for  a schematic of algorithm used to generate initial conditions of the chromosome. The beads coordinates generated by this algorithm are somewhat jagged, but an energy minimization will relax the structure.

<img align="center" width="600" src="./figures/3. Modeling the minimal cell/sc_growth_composite_0.png">

The very first thing that the job we submitted to Delta does is to generate initial configurations (i.e. coordinates) for the DNA and ribosomes using the program `sc_chain_generation` which was written by a previous graduate student Ben Gilbert. The code is available at [github.com/brg4/sc_chain_generation](https://github.com/brg4/sc_chain_generation). If one would like to generate coordinates for DNA, as well as ribosomes, for some other purpose, one should download, compile and run `sc_chain_generation` from the github link above. One can specify the number of ribosomes easily: just vary the number of obstacles `N_o` in the input script (`.inp` file) for `sc_chain_generation`. There are also parameters for sphere diameter, chromosome length, etc. The program `sc_chain_generation` can output the coordinates in either a `.bin`, `.dat`, or `.xyz` file format, the first of which is meant to be read by btree_chromo, and the last of which is human readable and easily read by VMD. To run an input file, you do  `/path/to/sc_chain_generation/src/gen_sc_chain --i_f=${input_fname} --o_d=${outputDirectory} --o_l=Syn3A_chromosome_init --s=10 --l=${log_fname} --n_t=8 --bin --xyz`.

## 4. Modeling DNA replication

The JCVI-syn3A minimal cell has a 543379 bp (543 kbp) genome comprised of 493 genes. This means an unreplicated chromosome is represented as a circular polymer of 54338 beads. Replication begins at a location on the genome called the origin (_Ori_), proceeds along the DNA in the clockwise and counterclockwise directions with Y-shaped structures (Fork), and ends at the terminal site, also called the terminus (_Ter_).  It turns out the replication states of the minimal cell aren't that interesting: it undergoes one replication initiation event per cell cycle, which means it starts with one unreplicated circular chromosome, and replication proceeds from _Ori_ to _Ter_ until we have two complete circular chromosomes.

<img align="center" width="700" src="./figures/4. Modeling chromosome dynamics/rep_state.png">

**Figure 1: Representing replication states.**  _Ori_, _Ter_, and Forks given in red, orange, and violet respectively.

As replication proceeds (starting from the _Ori_, along the Forks), the mother chromosome splits into two chromosomes, which we call the left and right daughter chromosomes. The structure is now no longer circular; it is now called a "theta structure" due to its resemblance to the Greek letter $\theta$. Both the left and right daughters have their own _Ori_'s, so in principle, they could begin to replicate too. However, DNA sequencing of the minimal cell indicates that we have only one replication initiation event per cell cycle.

For our simulations, we implement the "train-track" model of bacterial DNA replication[^gogou2021], where replisomes independently move along the opposite arms of the mother chromosome at each replication fork, replicating the DNA. There is another model called the "replication factory" model, but since Syn3A has so few regulatory mechanisms, this second one unlikely. (Plus, the train track model is also more consistent with our understanding of replication initiation[^thornburg2022].) In our implementation, new monomers are added to the left and right daughter chromosomes during replication by creating pairs of monomers centered around the corresponding position of the mother chromosome's monomers. 

<img align="center" width="800" src="./figures/4. Modeling chromosome dynamics/traintrack_updated.png">

**Figure 2: Replication with train-track model.**  Starting with an unreplicated Syn3A chromosome (543,379 bp) inside a 200 nm radius cell with 500 ribosomes (not shown), the 90% of the chromosome was replicated using the train-track model (refer to the schematic on the right). Spotlights are used to highlight the origins of replication (Oris), and a third spotlight shows the termination site (Ter) and replication forks, which have nearly reached the Ter. The magnifying glass magnifies one of the replication forks in the spotlight. Lime and magenta arrows indicate chromosome partitioning.

In our simulations, we update the replication state every 2 seconds of biological time. We assume the maximum replication rate of 100 bp/s for each of the replication forks, which means the replication forks each move 20 beads every time we update the replication state until they hit the Ter.

## 5. Modeling chromosome dynamics

### Energy of the system

The total potential energy for the chromosome/ribosome system is

$$U= \sum_{i=1}^{N_{\mathrm{DNA}}}\left[U_i^b+U_i^s\right] +\sum_{i=1}^{N_{\mathrm{DNA}}-1} \sum_{j=i+1}^{N_{\mathrm{DNA}}} U_{i j}^{\mathrm{DNA}-\mathrm{DNA}}+\sum_{i=1}^{N_{\mathrm{DNA}}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\mathrm{DNA}-\text { ribo }} +\sum_{i=1}^{N_{\text {ribo }}-1} \sum_{j=i+1}^{N_{\text {ribo }}} U_{i j}^{\text {ribo-ribo }} +\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\mathrm{DNA}}} U_{i j}^{\text {bdry-DNA }}+\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\text {bdry-ribo }}.$$

The energies for the bending, stretching and excluded volume interactions are shown below.

|**Bending:** | **Stretching** | **Excluded Volume** |
|:--:|:--:|:--:|
| <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_bending_0.png/" width="300"/>  |  <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_stretching_0.png/" width="200"/> |  <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_LJ_0.png/" width="120"/> |
|$U_i^b=\kappa_b\left[1-\cos \left(\pi-\theta_i\right)\right]$|$U_i^s= -\frac{\kappa_s L_0^2}{2} \log \left[1-\left(l_i / L_0\right)^2\right]$ $+4 \epsilon_s\left[\left(\frac{\sigma_s^s}{l_i}\right)^{12}-\left(\frac{\sigma_s}{l_i}\right)^6\right]$ $\times \Theta\left(2^{\frac{1}{6}} \sigma_s-l_i\right)$  |$U_{i j}^{e . v .}=  4 \epsilon_{e . v}\left[\left(\frac{\sigma_{e . v}}{r_{i j}}\right)^{12}-\left(\frac{\sigma_{e . v}}{r_{i j}}\right)^6\right]$ $\times  \Theta\left(2^{\frac{1}{6}} \sigma_{{e.v. }}-r_{i j}\right)$|

The excluded volume interaction between the DNA/ribosomes and "boundary" beads (cell membrane) ensure the DNA and ribosomes stay inside the spherical/overlapping sphere shaped volume.

<img align="center" width=250 src="./figures/4. Modeling chromosome dynamics/bdry_alt.svg">

### Brownian dynamics

In btree_chromo, we simulate dynamics of the DNA and ribosomes using a GPU-accelerated version of the Brownian dynamics integrator, which performs time-integration to update the coordinates of each of the beads. Let's quickly go through how the Brownian dynamics integrator works. First, recall how the acceleration of each bead is related to the net force on that bead and it's mass, which is given by Newton's second law, $F_{\text{net}}=ma$. 

For bead $i$, which has mass $m_i$, there are three forces acting on the bead: the "system" force associated with the interactions between beads, $F_{\text{system}} = -\nabla U_i $, as well as the drag force $F_{\text{drag}}$ and random force $F_{\text{rand}}$. So, Newton's second law reads

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Newton_eq.png">

The drag force is proportional to the bead's velocity, $\mathbf{v} = \text{d}\mathbf{x}/\text{d}t$, as well as its translational damping constant $\gamma_i$, given by the Einstein-Stokes equation: 

<img align="center" height=25 src="./figures/4. Modeling chromosome dynamics/Stokes-Einstein_eq.png">

Plugging in $F_{\text{drag}} = -\gamma_i v_i$ and dividing by $\gamma_i$, we obtain the Langevin equation of motion:

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Langevin_eq.png">

In the large friction limit, i.e. where $\gamma_i$ is very large, we can neglect the left hand side of the Langevin equation, and we obtain the equation for Brownian motion:

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Brownian_eq.png">

LAMMPS uses this equation to update positions, according to a simple first order integration scheme known as the Euler-Maruyama method (basically, the Euler method).

Both Langevin and Brownian dynamics can be used to correctly sample the NVT ensemble, but Brownian dynamics is preferred in our case since it allows us to take comparatively large time steps. Brownian dynamics is also sometimes called overdamped Langevin dynamics. This approximation is valid for timesteps that satisfy $\Delta t \gg m_i/\gamma_i$.

### SMC looping and topoisomerases

During the genome reduction process of Syn3A, guided by transposon mutagenesis studies on the original JCVI-syn1.0 genome and its intermediate reduced versions, it was found that structural maintenence of chromosomes (SMC) proteins were essential. Ganji et al. directly visualized the process by which condensin (aka, an SMC dimer) complexes extrude DNA into loops[^ganji2018]. They demonstrated that a single condensin can pull in DNA from one side at a force-dependent rate, supporting the loop extrusion model as a mechanism for chromosome organization. This finding provides strong evidence that SMC protein complexes like condensin actively shape the spatial arrangement of the genome.Magnetic tweezer experiments have been done to determine loop extrusion step size of ~200 bp[^ryu2022], and simulations indicate an extrusion frequency of ~2.5 steps/s[^nomidis2022].

<img align="center" height=300 src="./figures/4. Modeling chromosome dynamics/ganji.png">

**Figure 3: DNA loop extrusion by SMC complex.**  A series of snapshots shows DNA loop extrusion intermediates caused by SMC on a SxO-stained double-tethered DNA strand. A constant flow at a large angle to DNA axis stretches extruded loop and maintains DNA in imaging plane. Adapted from Ganji et al[^ganji2018].


https://github.com/user-attachments/assets/e46c7361-aabe-438c-96b5-40ca38abdc43

**Movie 1: Real-time imaging of DNA loop extrusion by SMC complex.**  Movie corresponding to Figure 6 showing DNA loop extrusion by SMC on a SxO-stained double-tethered DNA strand. A constant flow at a large angle to DNA axis stretches extruded loop and maintains DNA in imaging plane. From Ganji et al[^ganji2018].

(Most of the experiments on SMCs have been with eukaryotic SMCs like condensin. Syn3A actually has the bacterial SMC, SMC-ScpA/B, consisting of an SMC homodimer plus an ScpA kleisin protein and two accesory ScpB proteins. Recently [a preprint went up on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.11.653314v1) which claims to have repeated the Ganji et al. experiments with SMC-ScpA/B from Ureaplasma parvum, which is a not-too-distant relative of Mycoplasma mycoides which is the organism from which Syn3A is derived.)

The simulation methodology we use for SMC looping is that of Bonato and Michieletto, in which DNA loops are created by adding harmonic bonds bewteen two DNA monomers [^bonato2021] which reprent DNA bound by SMC. Many studies indicate that SMC loop extrusion in vivo is bidirectional, in other words both sides of the SMC translocate. Therefore in our model, both DNA monomers are updated. 

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/extrusion.png">

We take an extrusion rate of 500 bp/s, split into 250 bp/s for each side. This seems *very* fast: for example, RNAP transcribes at a maximum rate of 85 bp/s, and for DNA replication the forks have a maximum rate of 100 bp/s. However, the 0.5-1 kbp value has been corroborated by several studies.

There is also the question of the step size of the SMC's. It turns out the exact step size does not matter much in our simulations with respect to partitioning of the daughter chromosomes, but here update the SMC positions in our LAMMPS simulations every 2 biological seconds, which corresponds to extruding 500 bp on each side.

| Parameter | Description |
| --- | --- |
| Total number of loops | Number of active anchor+hinge pairs that are extruding loops. We know that there are ~100 SMC dimers[^gilbert2023]. In our simulations we estimate around half of them are bound at one time. |
| Loop extrusion frequency (s^-1) | How often does loop extrusion occur? Our best estimate is around every 0.4 s[^nomidis2022]. In our simulations we update the loops every 2 s.|
| Unbind/Rebind frequency (s^-1) | How fast are SMC unbind and rebinding? In our simulations, dwell times are on the order of minutes. We assume the unbinding and rebinding frequencies are equal, so that half of the SMCs are bound at one time|
| Extrusion step size (bp) | ~200 bp[^ryu2022]. In our simulations we update loops every 5 extrusion steps, corresponding to 1kbp or 500 bp on each side for bidirectional extrusion.|

In order to get the daughter chromosomes to partition, it turns out it is necessary to model another type of SMC behavior, namely blocking/bypassing. When SMCs encounter each other in our simulations, they block each other from translocating any further, and there is some rate for bypassing each other. Similarly, there is some rate for SMCs to bypass replication forks, but for these simulations we set that to zero.

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/extrude_block_bypass.png">

Also found to be essential were topoisomerases. There is evidence for coordination between topoisomerases and SMC complexes[^zawadzki2015]. For our simulations, topoisomerase is modeled by periodically running a set of minimizations and Brownian dynamics steps with DNA-DNA pair interactions replaced by soft potentials, which permits strand-crossings.

We don't have a great way of keeping track of strand crossings, but they usually happen when the SMC loops update, pulling strands of DNA taught against one another.

<img align="center" width=250 src="./figures/4. Modeling chromosome dynamics/topo.png">

<img align="center" width=250 src="./figures/4. Modeling chromosome dynamics/topo2.png">

## 6. A closer look at SMC dynamics

Consider the following toy example. Suppose we have five SMC’s that load uniformly on a segment of DNA. The SMCs will bind to the DNA and start to form loops, and as they do they will bridge progressively distant genomic sites. 

Below, we represent the looping state of the DNA in three ways: the physical structure, an arc diagram, and a contact map. 
- The physical structure of the DNA shows how the SMC contacts naturally leads to loop formation
- The arc diagram is a 1D line representing genomic locations where arc between i and j shows a contact between genomic locations i and j
- The contact map is a matrix where each axis represents genomic locations and a point corresponds to a contact between two locations. Although the matrix is symmetric, usually the elements both above and below the main diagonal are shown, as we do here. An SMC that bridges genomic locations i and j will be represented on the map by the points (i,j) and (j,i). For the uniform loading case, each of the 5 SMC’s spawn on the main diagonal, and then move diagonally away as i decreases and j increases at the same rate. We represent the growing SMCs with a green dot.

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/uniform_loading_1.png">

Eventually, there comes a point where the SMC’s encounter each other. When this happens, the middle 3 SMC’s are blocked on both sides, while the outer 2 SMC’s are blocked on one side. On the contact map, we represent the 3 stationary SMC’s with a red square, and we represent the two outer SMCs with a yellow triangle which points in the direction in which it is still free to move.

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/uniform_loading_2.png">

We can also draw contact maps for theta structures. Here, we gray square represents the mother chromosome, which spans the genomic positions between forks, and with the terminus in the center. The lime and magenta squares represent the left and right daughter chromosomes, each with their own origin.

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/rep_loop_state.png">

For such theta structures, we can have blocking at the forks. This will appear on contact maps as points at the edges of the squares. If the other side of the SMC is not blocked, the points will move along the edges of the square. In the figure below, can you match each of the loops in the arc diagram to the loops in the physical structure, as well as identify the corresponding points in the contact map?

<img align="center" width=600 src="./figures/4. Modeling chromosome dynamics/rep_loop_state_example.png">

Below, I have attached a movie that shows how the replication/loop state of chromosome changes throughout an ~90 minute cell cycle. Each axis represents the genomic position, which spans 543379 bp × 2 = 1086758 bp since eventually we will have two chromosomes. Starting with the mother chromosome (gray square), the left and right daughter chromosomes grow until we have two complete chromosomes (lime and magenta squares). For these simulations, we start with 50 loops, and add loops proportional to the total amount of DNA until we have 100 loops. Notice how for this set of parameters, most of the loops are fully blocked. You can also see the formation of "+" signs made from several loops, which corresponds to many loops that are blocked in a "traffic jam".

https://github.com/user-attachments/assets/1c6cfac8-d1d2-4148-9c79-f8d8c935bb40

It seems pretty plausible that SMC looping can help with chromosome partitioning, since it scrunches up DNA of the same type (left, right, mother). If you have something that is mixed up and you want to seperate them into distinct volumes, you can de-mix them by introducing bonds between components of the same type (think oil and water): this is exactly what SMC proteins do, and it close to the idea of enthalpy-driven phase separation. You might wonder if there is also an entropic component to the seperation. This is less clear. In organisms with clear nucleoid regions, the phase seperation is aided by crowding and depletion forces which are entropic in nature, but Syn3A does not have a defined nucleoid region. Entropy has been shown to aid with chromosome partitioning in cylindrically shaped bacteria, but Syn3A is spherical during replication.

In our simulations, SMC's start loading onto daughter chromosomes immediately after replication initiation. As the replication forks proceed along the mother chromosome, the daughter chromosomes lay right on top of each other. However, eventually the SMC's will set up "traffic jams" or "condensation centers" on each chromosome, which look like a cluster of SMC's with loops emanating from the clusters -- this is sometimes called a "bottle brush" structure. Then, looping is mostly done by the unblocked SMC's on the edges of the cluster, which reel in newly replicated DNA into the bottle brush, as well as DNA that strays into the other bottle brushes. You might imagine that if the loops are too long, then entropy actually would hinder chromosome separation, because loops would naturally diffuse into the wrong subvolumes. However, this seems to not be the case. If we assume there are 50 SMC's bound on the DNA, the average loop length will be around 543379bp / 50 = 10kbp. If you calculate the radius of gyration of each of those loops, the size of each loop is roughly equal to the radius of the cell. It turns out that the SMC cluster idea is consistent with the idea of a diffuse nucleoid region. 

<img align="center" width=800 src="./figures/4. Modeling chromosome dynamics/cluster_new.png">

In the movie below, I took a chromosome state which is about 2/3 of the way replicated, and ran dynamics for 3 seconds. The system is equilibrated, so it is free to mix if it wanted to, but despite this we only see repositioning within the cell of each of the chromosome domains (L lime, R magenta, and mother violet) and we do not see mixing of those domains. At the center of each of of these domains lies a cluster of SMC's.

https://github.com/user-attachments/assets/8387c708-43c8-486a-8082-b665156d4bbf

## 7. Understanding btree_chromo Commands

On Delta, you can `vim /projects/beyi/${USER}/DNA_SummerSchool_2025/scripts/template.inp` to see the template for the input scripts we run with `btree_chromo`. The most important lines are:

```bash
# run for 12 seconds of bio time (6 batches of 2 seconds) with dynamics
repeat:6
sync_simulator_and_system
set_initial_state
transform:m_cw20_ccw20
set_final_state
output_state:{output_dir}/rep_states/rep_state_{run_name}.txt
map_replication
write_loops:{output_dir}/loops/loops_{run_name}.txt
sys_write_sim_read_LAMMPS_data:{output_dir}data_{run_name}.lammps
translocate:50,T
sys_write_sim_read_LAMMPS_data:{output_dir}data_{run_name}.lammps
simulator_form_loops:F
simulator_minimize_topoDNA_harmonic:1000
simulator_set_delta_t:2.5E+7
{run_dynamics}
end_repeat

sync_simulator_and_system

# run for 48 seconds of bio time (24 batches of 2 seconds) without dynamics
repeat:24
set_initial_state
transform:m_cw20_ccw20
set_final_state
output_state:{output_dir}/rep_states/rep_state_{run_name}.txt
map_replication
write_loops:{output_dir}/loops/loops_{run_name}.txt
translocate:50,T
end_repeat
```

As indicated by the comments, we only run dynamics with LAMMPS for 20% of the biological time, but the other 80% of the time we still perform the replication and SMC blocking/bypassing dynamics.

The `repeat` commands are exactly like for loops. Each iteration represents 2 seconds of biological time, during which we replicate on both sides by 20 beads (`transform:m_cw20_ccw20`) and translocate SMCs by 50 beads (`translocate:50,T`. What this really means is we attempt 50 moves; whether a loop updates or not depends on whether its blocked by another loop, and the parameters we set for the blocking, bypassing, and basal birth/death rates for the SMCs).

The command `simulator_form_loops:F` reads in the loop state from `btree_chromo` into the LAMMPS simulation object. The command `simulator_minimize_topoDNA_harmonic:1000` runs a minimization with strand crossing permitted, `simulator_set_delta_t:2.5E+7` sets the timestep to 25 ns, and `{run_dynamics}` runs Brownian dynamics with strand crossings forbidden.

## 8. Visualization with VMD

Here we follow the [VMD Guide by Tianyu](../RDME/vmd_guide.md) copied below for convenience with some small modifications.

### 1. Initialize the OOD Interactive Session
1. Navigate to the [Open OnDemand dashboard](https://openondemand.delta.ncsa.illinois.edu/pun/sys/dashboard).

2. Log in through CILogon with your NCSA username, password, and Duo MFA.

3. Open the Interactive Apps menu and click Desktop.

4. Configure the job settings and click Launch:
   - Container image: keep default
   - Account: `beyi-delta-gpu`
   - Partition: `cpu-interactive`
   - Duration: `00-00:30:00`
   - Reservation: leave empty if none
   - CPUs: `16`
   - RAM: `64GB`
   - GPUs: `1`

5. Wait for the job status to change from "starting" to "running" in My Interactive Sessions. 

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-starting.png" alt="starting" width="300">

   Click "Connect to Desktop" to access the Linux graphical interface.

   <img src="https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-connect.png" alt="running" width="300">

### 2. Preprocess Trajectory and Load VMD Module
 Open a terminal and run:

```bash
cd /projects/beyi/$USER/DNA_SummerSchool_2025/data/
```

Preprocess the trajectory by doing
```bash
python3 modify_lammpstrj.py
```

This allows the DNA for the left and right daughters and mother to be colored differently. It should take ~3 minutes to run.

Update the load_btree_chromo.tcl file by copying it from my directory (I made some changes to it from yesterday because it turns out Delta doesn't allow GPU accelerated rendering :-()
```bash
cp /projects/beyi/amaytin/DNA_SummerSchool_2025/data/load_btree_chromo.tcl .
```

Next, load and open vmd by doing:


```bash
module load vmd
vmd
```

### 3. Load the LAMMPS trajectory file
   In the VMD "Main" window, click on "Extensions" and then "TkConsole". In the "TkConsole" window, do 
```bash
source load_btree_chromo.tcl
```
This script will take ~2 minutes to run.
You should see a representation of the trajectory for the Minimal Cell growth and division!

**Important considerations for Windows Users:**
For those using a Windows machine, you will need to make sure your environmental variables for LAMMPS are set correctly by setting them via command line outside of VMD. Before starting VMD, in Windows Command Shell, please do the following:
```bash
setx LAMMPSDUMMYPOS "$xd,$yd,$zd"
setx LAMMPSMAXATOMS "200000"
setx LAMMPSREMAPFIELDS "vx=c_id_track,vy=c_type_track"
```
Each entry should produce “SUCCESS: Specified value was saved.”

This workaround is only needed on Windows VMD (i.e. not on Linux and Mac VMD). This issue will be addressed in upcoming VMD releases. As of writing this, the latest VMD is Version 1.9.4.

| Monomer type | Color | Bead Size |
| --- | --- | --- |
| DNA | gray(M), lime(L), magenta(R) | 13.0 |
| Ori | red | 39.0 |
| Ter | orange | 39.0 |
| Fork | violet | 39.0 |
| Ribosome | mauve | 70.0 |
| Boundary | silver | 32.5 |
| SMC1 | black | 19.5 |
| SMC2 | white | 19.5 |

For the presentation on the last day, it would be nice to have a movie of the trajectory on one of your slides. For those who are interested in making a movie, I will work with whoever is interested in our extra time to make it. I have some commands to render each of the frames, and then to compile it into a `.mp4`.

You can use these commands to generate a higher quality movie than allowed by VMD Movie Maker:

Set directory to save frames: set this to something reasonable
```
set outdir "/tmp"
```

```bash
file mkdir $outdir
```

Set the number of frames in the trajectory
```
set nframes [molinfo top get numframes]
```


```bash
# Loop over each frame
for {set i 0} {$i < $nframes} {incr i} {
    # Set the frame
    animate goto $i

    # Format the output filename with zero-padded index
    set fname [format "%s/frame%04d.tga" $outdir $i]

    # Render with TachyonInternal at high resolution and antialiasing
    render TachyonInternal $fname -res 1920 1080 -aa 12

    puts "Rendered frame $i to $fname"
}
```

```bash
ffmpeg -framerate 30 -i frame%04d.tga -c:v libx264 -pix_fmt yuv420p -crf 18 high_quality_movie.mp4
```

## References
[^gilbert2023]: Gilbert, Benjamin R., Zane R. Thornburg, Troy A. Brier, Jan A. Stevens, Fabian Grünewald, John E. Stone, Siewert J. Marrink, and Zaida Luthey-Schulten. “Dynamics of Chromosome Organization in a Minimal Bacterial Cell.” Frontiers in Cell and Developmental Biology 11 (August 9, 2023). https://doi.org/10.3389/fcell.2023.1214962.
[^gogou2021]: Gogou, Christos, Aleksandre Japaridze, and Cees Dekker. “Mechanisms for Chromosome Segregation in Bacteria.” Frontiers in Microbiology 12 (June 2021). https://doi.org/10.3389/fmicb.2021.685687.
[^thornburg2022]: Thornburg, Zane R., David M. Bianchi, Troy A. Brier, Benjamin R. Gilbert, Tyler M. Earnest, Marcelo C. R. Melo, Nataliya Safronova, et al. “Fundamental Behaviors Emerge from Simulations of a Living Minimal Cell.” Cell 185, no. 2 (January 20, 2022): 345-360.e28. https://doi.org/10.1016/j.cell.2021.12.025.
[^thornburg2025]: Thornburg, Zane R., Andrew Maytin, Jiwoong Kwon, Troy A. Brier, Benjamin R. Gilbert, Enguang Fu, Yang-Le Gao, Jordan Quenneville, Tianyu Wu, Henry Li, Talia Long, Weria Pezeshkian, Lijie Sun, John I. Glass, Angad Mehta, Taekjip Ha, and Zaida Luthey-Schulten. “Bringing the Genetically Minimal Cell to Life on a Computer in 4D.” bioRxiv, June 10, 2025. https://doi.org/10.1101/2025.06.10.658899.
[^ryu2022]: Ryu, Je-Kyung, Sang-Hyun Rah, Richard Janissen, Jacob W J Kerssemakers, Andrea Bonato, Davide Michieletto, and Cees Dekker. “Condensin Extrudes DNA Loops in Steps up to Hundreds of Base Pairs That Are Generated by ATP Binding Events.” Nucleic Acids Research 50, no. 2 (January 25, 2022): 820–32. https://doi.org/10.1093/nar/gkab1268.
[^nomidis2022]: Nomidis, Stefanos K, Enrico Carlon, Stephan Gruber, and John F Marko. “DNA Tension-Modulated Translocation and Loop Extrusion by SMC Complexes Revealed by Molecular Dynamics Simulations.” Nucleic Acids Research 50, no. 9 (May 20, 2022): 4974–87. https://doi.org/10.1093/nar/gkac268.
[^ganji2018]: Ganji, Mahipal, Indra A. Shaltiel, Shveta Bisht, Eugene Kim, Ana Kalichava, Christian H. Haering, and Cees Dekker. “Real-Time Imaging of DNA Loop Extrusion by Condensin.” Science 360, no. 6384 (April 2018): 102–5. https://doi.org/10.1126/science.aar7831.
[^bonato2021]: Bonato, Andrea, and Davide Michieletto. “Three-Dimensional Loop Extrusion.” Biophysical Journal 120, no. 24 (December 2021): 5544–52. https://doi.org/10.1016/j.bpj.2021.11.015.
[^zawadzki2015]: Zawadzki, Pawel, Mathew Stracy, Katarzyna Ginda, Katarzyna Zawadzka, Christian Lesterlin, Achillefs N. Kapanidis, and David J. Sherratt. “The Localization and Action of Topoisomerase IV in Escherichia Coli Chromosome Segregation Is Coordinated by the SMC Complex, MukBEF.” Cell Reports 13, no. 11 (December 22, 2015): 2587–96. https://doi.org/10.1016/j.celrep.2015.11.034.
