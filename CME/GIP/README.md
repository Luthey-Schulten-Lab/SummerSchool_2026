# Stochastic Genetic Information Processes in CME

Open [ChatGPT](https://chatgpt.com/) and ask the following question:

**Give an example of large differences in rate constants that lead to significant fluctuations.**

Please also review this previous response: [Stochastic Genetic Information Processes](https://chatgpt.com/c/91d87e61-bae5-4b89-b078-d1d4cfa44274), which corresponds to this tutorial.

## Classic Genetic Information Process Reactions

This classic and simplified genetic information process (GIP) model consists of 3 species and 4 **first-order reactions**. It begins with the transcription of a gene into mRNA. The mRNA can be translated into protein or degraded into its monomers. The resulting protein can also undergo degradation. The reaction scheme and associated rate constants are shown below.

<p align="center">
  <img src="../figs/figs_GIP/GIP_withCMEs.png" width="600" alt="Simple GIP model">  <br>
  <b>Figure 1. Genetic information processing model and its corresponding chemical master equation, where <i>m</i> and <i>n</i> are the counts of mRNAs and proteins.</b>
</p>

The rate constants are for the protein DnaA from the minimal cell. To calculate the first three rate constants, we used the concentrations of nucleotides and aminoacyl-tRNAs, length of gene/mRNA/protein, and the active ratio of RNAP, ribosome, and degradasome as reported in the *Cell* paper[^thornburg_cell]. 

We assume the gene copy number is fixed at 1 and the initial mRNA count is 1, and the initial protein count is 148 from the proteomics study. The protein degradation rate is estimated based on a half-life of 25 hours[^thornburg_kinetic].

**Table 1. Four reactions with their rate constants**

| **Names**              | **Reaction**                          | **Rate Constant (s<sup>-1</sup>)**                              | **Propensity (s<sup>-1</sup>)**                              |
|------------------------|----------------------------------------|------------------------------------------------------|---------------------------------------------------|
| Transcription          | Gene → mRNA                            | *k*<sub>transcription</sub> = 6.41×10<sup>-4</sup>   | *k*<sub>transcription</sub>                       |
| Degradation of mRNA    | mRNA → ∅                               | *k*<sub>deg,m</sub> = 2.59×10<sup>-3</sup>           | *k*<sub>deg,m</sub> · *N*<sub>mRNA</sub>          |
| Translation            | mRNA → mRNA + Protein                  | *k*<sub>translation</sub> = 7.20×10<sup>-2</sup>     | *k*<sub>translation</sub> · *N*<sub>mRNA</sub>    |
| Degradation of Protein | Protein → ∅                            | *k*<sub>deg,p</sub> = 7.70×10<sup>-6</sup>           | *k*<sub>deg,p</sub> · *N*<sub>ptn</sub>           |

## Run the Jupyter Notebook

Open the Jupyter Notebook interface and navigate to the directory `/CME/GIP/`. Run the notebook `Tut2.1-GeneticInformationProcess.ipynb` to simulate this toy model of genetic information processing (GIP). 

By default:
- The total simulation time `simtime` is set to 6300 seconds, representing the full cell cycle of the minimal cell.
- We simulate 10 independent cells (`reps = 10`).
- Trajectories are recorded at intervals of 1 second (`writeInterval = 1`).

## Stochastic Protein Synthesis

We begin by examining the average and variation in mRNA and protein abundances across the population of 10 simulated cells in **Figure 2 Left**. The population-averaged mRNA abundance fluctuates below 1 throughout the cell cycle, while protein levels increase steadily due to translation. Protein degradation is minor under this set of kinetic parameters.

Protein synthesis is a stochastic process that occurs in individual cells. By plotting the time traces of mRNA and protein for several replicates, you can observe the characteristic stair-step pattern of mRNA production and burst-like synthesis of proteins. You will see that protein levels increase or "burst" when mRNAs are present, and may plateau or decrease when mRNAs are absent. We encourage you to compare the patterns shown here with your own simulation results.

<p align="center">
  <img src="../figs/plots_GIP/GIP_mRNA_Protein_10Replicates.png" width="300" alt="mRNA Protein 10 reps"> 
  <img src="../figs/plots_GIP/GIP_mRNA_Protein_Cell1.png" width="300" alt="CME replicate 1"> <br>
  <b>Figure 2. Left: Population average (solid line) and full range (shaded area) of mRNA (red) and protein (blue) abundances across 10 cell replicates. <br> 
  Right: Stair-step mRNA trace and burst-like protein synthesis in a single cell replicate.</b>
</p>

## Discussion

### 1. Steady-State

Do mRNA and protein levels reach a steady state during the 6300-second simulation? How can you tell from the plots? If the fluctuations are large, try increasing the number of replicates `reps` from 10 to 100.

### 2. Doubling the Initial Abundance of Protein for Cell Division

The initial count of protein P\_0001 (DnaA) from experimental proteomics data is 148. In the histogram below, the average DnaA count at the end of the cell cycle is approximately 270. Compare the mean protein count at the end of the simulation to this experimental value. Does the simulation produce roughly 148 new proteins over the full cell cycle? Why is this important? Consider the implications of cell division into two daughter cells.

<p align="center">
  <img src="../figs/plots_GIP/GIP_Proteins_CycleEnd_100replicates.png" width="450" alt="Protein end-cycle distribution"> <br>
  <b>Figure 3. Distribution of protein abundances among 100 cell replicates at the end of the cell cycle.</b>
</p>

## References:
[^thornburg_cell]: Thornburg, Z. R., Bianchi, D. M., Brier, T. A., Gilbert, B. R., Earnest, T. M., Melo, M. C., Safronova, N., Sáenz, J. P., Cook, A. T., Wise, K. S., Hutchison, C. A., Smith, H. O., Glass, J. I., & Luthey-Schulten, Z. (2022). Fundamental behaviors emerge from simulations of a living minimal cell. Cell, 185(2), 345-360.e28. https://doi.org/10.1016/j.cell.2021.12.025

[^thornburg_kinetic]: Thornburg, Z. R., Melo, M. C. R., Bianchi, D., Brier, T. A., Crotty, C., Breuer, M., Smith, H. O., Hutchison, C. A., Glass, J. I., & Luthey-Schulten, Z. (2019). Kinetic modeling of the genetic information processes in a minimal cell. Frontiers in Molecular Biosciences, 6. https://doi.org/10.3389/fmolb.2019.00130
