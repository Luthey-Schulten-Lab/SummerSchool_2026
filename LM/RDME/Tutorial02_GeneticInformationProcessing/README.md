# A Spatially Heterogeneous Genetic Information Processing System Solved Stochastically #

## Introduction ##
In tutorial 2, we will apply the spatial techniques used in tutorial 1 to a simplified genetic information processing (GIP) system. We will use the same GIP system we used in the spatially homogeneous tutorials, but here we will impose this reaction system onto a lattice which we design to resemble the standard sized (200nm radius) minimal cell (JCVI-Syn3A). In addition to modeling basic GIP reactions, we will also make our system architecture more realistic by including the following region types: the extracellular space, the plasma membrane, the inner membrane region, the cytosol, ribosomes, and DNA. Each of these site types will be given a mask during our `RegionBuilder` steps, some of which (the cytosol, the plasma membrane, and the inner membrane region) we will create using idealized geometries and others (the ribosomes and DNA) we will load into our system using predefined masks. 

Treating the ribosomes, DNA, and the inner membrane space as their own regions enable us to require specific GIP reactions occur at predefined region types. For example, using the ribosome region type allows us to place "effective" ribosomes throughout the cell and only allow translation reactions to take place when an mRNA molecule has diffused to a ribosomal region type. Furthermore, it allows us to require transcription reactions to take place within only DNA region type lattices. Finally, because the mRNA and protein degradation complexes in the minimal cell are bound to the plasma membrane, we create an inner plasma membrane region just inside the true plasma membrane region to mimic the location of these complexes. We then only allow mRNA and protein degradation reactions to take place in this inner plasma membrane region site type.

The creation of a ribosome region type mask will be done using a previously written [python script](./Utils/T2_RibosomeAndDNAMasks.py). This script randomly chooses 500 cytosolic lattice sites and designates them as ribosome sites. Then, because we will use a voxel size of 8nm, we will extend these regions to also include their 6-nearest neighbor lattices to mimic the actual size of the minimal cell's ribosomes. The DNA mask will be created using the coordinates from a single replicate of the minimal cell DNA simulations (see [Gilbert et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8339304/) and [Gilbert et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10445541/)). 

Importantly, we will define diffusion rules for each species to try to mimic how these species would behave in a biological context. For example, we will place the *DnaA* gene in a specific DNA lattice to mimic the gene being located at a unique place on the chromosome and we will not allow this particle to diffuse by setting its diffusion coefficients to zero for all region type transitions. The RNA species in our system will be allowed to diffuse from the DNA to the cytosol, but not the other way around. This will allow it to be synthesized on the DNA lattice where the *DnaA* gene is located, diffuse into the cytosol and ribosomes, but never be able to return to the DNA lattice sites. We will provide a similar rule for the protein species, but instead of being synthesized on the DNA lattice, it will be synthesized on the ribosome lattice. It will then be able to freely diffuse into the cytosol, but not back onto the ribosome lattice. 

## Further Investigation ##
### 1. How are the trajectories of each chemical species different in the spatially homogeneous and the spatially heterogeneous system treatments? ###
 - Do mRNA molecules seem to exist for more or less time than in the spatially homogeneous simulation of the GIP system? Why?
 - Do protein molecules seem to exist for more or less time than in the spatially homogeneous simulation of the GIP system? Why?

### 2. What are the time-dependent abundances of the RNA species in 1) the cytosol and 2) the inner membrane region? ##
 - Using the precomputed trajectory file that lasts 6300s, paste the following code into a new chunk at the end of the tutorial notebook to extract out cytosol-specific and shell-specific abundances of the RNA species. What can you learn from this information?

```python
# Region-specific mRNA trajectories: cytosol vs shell

# Resolve species and region IDs the same way jLM does internally
mRNA_id  = traj._speciesIdParse(None, "mRNA",    startIndex=1)[0]
cyto_id  = traj._regionIdParse(None,  "cytosol")[0]
shell_id = traj._regionIdParse(None,  "shell")[0]

# Load the static site lattice (maps each voxel to its region type ID)
site_lattice = traj.h5['/Model/Diffusion/LatticeSites'][...]

# Iterate over every saved lattice frame and accumulate region-specific mRNA counts
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

# ── Plot ──────────────────────────────────────────────────────────────────────
# Step plots correctly represent discrete molecule counts: horizontal segments
# show occupancy and vertical jumps mark birth/death events.
# Top panel: full 6300 s overview with zoom region highlighted.
# Bottom panel: tight zoom around the birth (~1366 s) and death (~1376 s) event.

BUFFER    = 10
ZOOM_LO   = 1366 - BUFFER
ZOOM_HI   = 1376 + BUFFER

sns.set(style="ticks")
palette   = sns.color_palette("husl", 4)
cyto_col  = palette[1]
shell_col = palette[3]

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(12, 7),
    gridspec_kw={'height_ratios': [1, 2]},
)
fig.patch.set_facecolor('white')

# ── Top panel: full 6300 s overview ──
stride = max(1, len(ts_r) // 1000)
ax_top.step(ts_r[::stride], mRNA_cyto_counts[::stride],
            where='post', label='cytosol', color=cyto_col,  linewidth=1.2)
ax_top.step(ts_r[::stride], mRNA_shell_counts[::stride],
            where='post', label='shell',   color=shell_col, linewidth=1.2, linestyle='--')
ax_top.axvspan(ZOOM_LO, ZOOM_HI, color='gray', alpha=0.12, label='zoomed region')
ax_top.set_xlim(ts_r[0], ts_r[-1])
ax_top.set_ylabel('mRNA count')
ax_top.set_title('Full trajectory (0 – 6300 s)', fontsize=11)
ax_top.legend(frameon=False, loc='upper right', fontsize=9)
ax_top.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
sns.despine(ax=ax_top)

# ── Bottom panel: tight zoom around single mRNA birth/death event ──
zoom_mask = (ts_r >= ZOOM_LO) & (ts_r <= ZOOM_HI)
ax_bot.step(ts_r[zoom_mask], mRNA_cyto_counts[zoom_mask],
            where='post', label='mRNA (cytosol)', color=cyto_col,  linewidth=1.2)
ax_bot.step(ts_r[zoom_mask], mRNA_shell_counts[zoom_mask],
            where='post', label='mRNA (shell)',   color=shell_col, linewidth=1.2, linestyle='--')
ax_bot.set_xlim(ZOOM_LO, ZOOM_HI)
ax_bot.set_xlabel('Time (s)')
ax_bot.set_ylabel('mRNA count')
ax_bot.set_title(f'Zoomed view: t = {ZOOM_LO}–{ZOOM_HI} s  (birth ≈ 1366 s, death ≈ 1376 s)', fontsize=11)
ax_bot.legend(frameon=False, loc='upper right', fontsize=9)
ax_bot.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
sns.despine(ax=ax_bot)

plt.suptitle('Region-Specific mRNA Abundance: Cytosol vs Shell', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('./Plots/TutR2_mRNA_regions.png', facecolor='white', bbox_inches='tight')
plt.show()
```

## 3. How do geometry choices affect the system behavior? ##
 - How would decreasing the size of the cell affect the species' trajectories?
 - If species were only allowed to be placed and react only in the plasma membrane, how might this change the trajectories?
 - If `build.se26` was used as the SE for the dilation of the membrane, would this have any affect on the membrane-only reacting system? If so, how?