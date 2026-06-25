# A Spatially Heterogeneous Bimolecular Reaction System Solved Stochastically #

## Introduction ##
In tutorial 1, we will learn how to setup and simulate a spatially stochastic chemical system using Lattice Microbes. We will use the same chemical system we used in previous tutorials (the bimolecular reaction), but here we will impose this reaction system onto a lattice which we design to resemble a large version (400nm radius) of the minimal cell (JCVI-Syn3A). 

One of the biggest differences in spatially heterogeneous simulations is that we now need to define the architecture of the system. A detailed description of how to use Lattice Microbes to define spatial architectures using the `RegionBuilder` class is given in the [RegionBuilderCheatSheet.md](../RegionBuilderCheatSheet.md) and in the present tutorial. System architectures can be visualized using jLM built-in functions or other stand-alone software such as `VMD`. 

## Further Investigation ##
### 1. How are the trajectories of each chemical species different in the spatially homogeneous and the spatially heterogeneous system treatments? ###
 - Do the final/steady-state abundances differ in these treatments?
 - Does the rate at which steady-state is achieved differ?

### 2. What would the trajectories look like if there was an initial spatial bias in species' locations in the simulation? ##
In cells, chemical species are not always distributed uniformly throughout the cytosol. How might this affect the behavior of the system? Replace the chunk of code in tutorial 1 where the particle abundances and locations were set with the following chunk of code that initially biases species A and B to either side of the cytosol:

**Note:** Before running this code, make sure to change the output file name in the simulation object and the plotting section to ensure you do not write over your previous results!

```python
# Define the center of the cytosol #
cyto_cent = [32,32,32]
# Define cytosol left (x-axis) region #
cyto_left_cent = [32-int(radius_voxel/2),32,32]
# Define cytosol right (x-axis) region #
cyto_right_cent = [32+int(radius_voxel/2),32,32]

# Create chemical species objects #
# Species A #
spA = sim.species('A')
# Species B #
spB = sim.species('B')
# Species C #
spC = sim.species('C')

# Define initial abundances #
total_A = 1000
total_B = 1000
total_C = 0

# Define a function to place initially place particles with a bias #
from itertools import product
from collections import deque

def placeNumbersInCyto(sp, x, y, z, n):
    pps = sim.pps  # max particles allowed per voxel (compile-time constant, usually 16)

    # BFS starting from voxel (x, y, z), expanding outward through all 26 face/edge/corner neighbors
    q = deque([(x, y, z)])   # queue of voxels yet to be filled
    enqueued = {(x, y, z)}   # set of voxels already added to the queue (O(1) lookup)
    visited = []              # ordered list of voxels that were actually filled
    remaining = n             # number of particles still to be placed

    while q and remaining > 0:
        cx, cy, cz = q.popleft()          # take the next voxel in BFS order
        visited.append((cx, cy, cz))

        # Fill this voxel up to its capacity or until we run out of particles
        to_place = min(pps, remaining)
        sim.placeNumber(sp=sp, x=cx, y=cy, z=cz, n=to_place)
        remaining -= to_place

        # Enqueue all 26 neighbors that are inside the lattice and not yet seen
        for dx, dy, dz in product((-1, 0, 1), repeat=3):
            if dx == dy == dz == 0:       # skip the center voxel itself
                continue
            nb = (cx+dx, cy+dy, cz+dz)
            if nb not in enqueued and all(0 <= nb[i] < Ldim[i] for i in range(3)):
                enqueued.add(nb)
                q.append(nb)

    print(f"Placed {n} {sp.name} in the cytoplasm.")
    return visited

# Place species A with a bias #
ini_A_reg = placeNumbersInCyto(spA, cyto_left_cent[0], cyto_left_cent[1], cyto_left_cent[2], total_A)
# Place species B with a bias #
ini_B_reg = placeNumbersInCyto(spB, cyto_right_cent[0], cyto_right_cent[1], cyto_right_cent[2], total_B)
```
 - After running this code, what do you notice about the output plot?
 - How would changing the diffusion coefficient affect the new system behavior?

## 3. How do geometry choices affect the system behavior? ##
 - How would decreasing the size of the cell affect the species' trajectories?
 - If species were only allowed to be placed and react only in the plasma membrane, how might this change the trajectories?
 - If `build.se26` was used as the SE for the dilation of the membrane, would this have any affect on the membrane-only reacting system? If so, how?