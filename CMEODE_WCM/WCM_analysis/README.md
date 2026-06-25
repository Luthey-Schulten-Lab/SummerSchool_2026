### Analysis Package for Whole-Cell Model (WCM)
#### Workflow of analysis
1. Serialize output CSV files into a single Python pickel file. The pickel file stores all infos from the simulation, including traces of metabolites, GIP species, and fluxes, also species names, number of replicates.
2. Extract traces or plot from pickel files using functions in script *WCM_analysis.py* with the help of functions in other *WCM_.py* scripts
#### Architecture of Package
The scripts could be divided into functional (with prefix *WCM_*) and analyzing. Core script *WCM_analysis.py* has functions to serialize CSV files into pkl file (with *WCM_traj.py*), extract traces, and plot figures from pkl file. Other functional scripts (with suffix *gene, math, metabolites, mRNA, ptn, ribo, cplx, diagnosis*) contain custom methods supporting the biological analysis. The analyzing scripts are divided by subsystems, like *gene, flux, growth, intrametabolites, mRNA, protein, ribo, tRNAcharging* etc.. These codes **should be changed intensively** for the interests of analyzing.
