"""
Generate Minecraft .schem and .nbt structures from data.
"""

from .mcschematic_plus import MCSchematicPlus
from .data_loaders import (
    to_mc_volume,
    read_tiff,
    read_npy,
    read_mesh,
    voxelize_mesh,
)
from .block_colormap import BlockColormap, get_block_colormap
