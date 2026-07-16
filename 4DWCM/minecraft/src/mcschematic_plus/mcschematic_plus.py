from __future__ import annotations

import numpy as np
import os
from collections import defaultdict
from pathlib import Path
import pyvista as pv
from typing import BinaryIO, Sequence
from mcschematic import MCSchematic, MCStructure, Version
from nbtlib.tag import *
from nbtlib import File, parse_nbt
from .block_colormap import BlockColormap, get_block_colormap

class MCSchematicPlus(MCSchematic):
    def __init__(self, schematicToLoadPath_or_mcStructure: str | os.PathLike | MCStructure = None, version: 'Version' = None):
        if isinstance(schematicToLoadPath_or_mcStructure, (str, os.PathLike)):
            schematicToLoadPath = str(schematicToLoadPath_or_mcStructure)
            if not os.path.isfile(schematicToLoadPath):
                raise FileNotFoundError(f"Schematic file not found: {schematicToLoadPath}")
            if schematicToLoadPath.lower().endswith(".schem"):
                self._initFromFile(schematicToLoadPath)
            else:
                try:
                    import amulet
                except ImportError:
                    raise ImportError("To load non .schem files, please install the 'amulet-core' package")
                level = amulet.load_level(os.path.abspath(schematicToLoadPath))
                dim = "main"
                sel = level.bounds(dim)
                self._defaultInit()
                version = ("java", version.value) if version is not None else ("java", self.getLatestVersion().value)
                for x in range(sel.min_x, sel.max_x):
                    for y in range(sel.min_y, sel.max_y):
                        for z in range(sel.min_z, sel.max_z):
                            block, blockEntity = level.get_version_block(x, y, z, dim, version=version)
                            if block.namespaced_name != "minecraft:air":
                                self.setBlock((x, y, z), block.full_blockstate)
                level.close()
        else:
            super().__init__(schematicToLoadPath_or_mcStructure)

    def placeVolume(self, volumeMask: np.ndarray, blockData_or_color: str | np.ndarray | None, blockColormap: str | BlockColormap | None = None, placePosition: tuple[int, int, int] = (0, 0, 0)):
        """Add blocks for every True voxel in volume_mask. The volume axes should be (x, y, z) = (East, Up, South)."""
        volumeMask = np.asarray(volumeMask, dtype=bool)
        positions = np.argwhere(volumeMask)  # shape (N, 3) with columns [x, y, z]
        if positions.size == 0:
            return
        positions += np.array(placePosition)  # offset coordinates
        if blockData_or_color is not None:
            _bdoc_arr = np.asarray(blockData_or_color)
            if _bdoc_arr.ndim >= 3:
                if _bdoc_arr.shape[:3] != volumeMask.shape:
                    raise ValueError("If blockData_or_color is an array, its first 3 dimensions must match the volume_mask shape")
                blockData_or_color = _bdoc_arr[volumeMask]  # filter to only True voxels
        self.setBlocks(positions, blockData_or_color, blockColormap)
    
    def setBlocks(self, positions: np.ndarray, blockData_or_color: str | Sequence[str] | np.ndarray | None, blockColormap: str | BlockColormap | None = None):
        """Add blocks for every (x,y,z) in points."""
        positions = np.asarray(positions)
        if positions.size == 0:
            return
        positions = np.asarray(positions)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("positions must be an array of shape (N, 3)")
        if blockData_or_color is None:
            for x, y, z in positions:
                pos = (int(x), int(y), int(z))
                if pos in self._structure._blockStates:
                    del self._structure._blockStates[pos]
                if pos in self._structure._blockEntities:
                    del self._structure._blockEntities[pos]
            return
        if blockColormap is not None:
            cmap = get_block_colormap(blockColormap)
            blockData = cmap.get_block_state(blockData_or_color)
        else:
            if isinstance(blockData_or_color, str) or isinstance(blockData_or_color[0], str): # single block or list of blocks
                blockData = blockData_or_color
            else:
                raise ValueError("If blockColormap is not provided, blockData_or_color must be a blockData string or a sequence of blockData strings.")
        if isinstance(blockData, str):
            for x, y, z in positions:
                self.setBlock((int(x), int(y), int(z)), blockData) # TODO slightly inefficient but simple
        else:
            for (x, y, z), bd in zip(positions, blockData):
                self.setBlock((int(x), int(y), int(z)), str(bd))

    def getBlocks(self):
        """Return a dictionary of positions to block names."""
        block_dict : dict[tuple[int, int, int], str] = {}
        for pos in self._structure.getBlockStates().keys():
            block_dict[pos] = self.getBlockDataAt(pos)
        return block_dict
    
    def getBlockName(self, blockData: str | None) -> str | None:
        """Get the block name from blockData string. If blockData is None, return None."""
        if blockData is None:
            return None
        return blockData.split("[")[0]  # Remove properties if present

    def replaceBlocks(self, oldBlock: str, newBlock: str | None):
        """Replace all instances of oldBlock with newBlock. If newBlock is None, remove the block."""
        block_states_copy = dict(self._structure._blockStates)  # avoid modifying dict during iteration
        if newBlock is None:
            for pos, blockPaletteId in block_states_copy.items():
                if self._structure._blockPalette[blockPaletteId] == oldBlock:
                    if pos in self._structure._blockEntities:
                        del self._structure._blockEntities[pos]
                    del self._structure._blockStates[pos]
            # remove from palette
            blockPaletteId = self._structure._blockPalette[oldBlock]
            del self._structure._blockPalette[oldBlock]
            del self._structure._blockPalette[blockPaletteId]  
        
        else:
            palette_copy = dict(self._structure._blockPalette)  # avoid modifying dict during iteration
            # change the palette entry
            blockPaletteId = self._structure._blockPalette[oldBlock]
            self._structure._blockPalette[blockPaletteId] = newBlock
            self._structure._blockPalette[newBlock] = self._structure._blockPalette.pop(oldBlock)

            for pos, blockPaletteId in block_states_copy.items():
                if palette_copy[blockPaletteId] == oldBlock:
                    if pos in self._structure._blockEntities:
                        del self._structure._blockEntities[pos]
                    self.setBlock(pos, newBlock)

    def splitByBlock(self):
        """Return a dictionary of block names and corresponding MCSchematicPlus objects."""
        block_dict : dict[str, MCSchematicPlus] = {}
        for pos, blockPaletteId in self._structure.getBlockStates().items():
            block_name = self._structure._blockPalette[blockPaletteId]
            if block_name not in block_dict:
                block_dict[block_name] = MCSchematicPlus()
            block_dict[block_name].setBlock(pos, block_name) # TODO slightly inefficient but simple
        return block_dict
    
    def getLatestVersion(self) -> 'Version':
        """Get the latest supported Minecraft version."""
        return max(Version, key=lambda v: v.value)

    def save(self, filepath: str | os.PathLike | BinaryIO, version: 'Version' = None, fastSaving: bool = False):
        """Save the structure as a Sponge schematic file at `filepath`."""
        # if filepath is a string, ensure directory exists
        if isinstance(filepath, (str, os.PathLike)):
            filepath = str(filepath)
            directory = os.path.dirname(filepath)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            # add .schem extension if not present
            if not filepath.lower().endswith('.schem'):
                filepath += '.schem'
        if version is None:
            version = self.getLatestVersion()
        #################### Start copied from MCSchematic #####################
        ## Setup
        schemBounds = self._structure.getBounds()
        schemDims = self._structure.getStructureDimensions(schemBounds)
        # The vector amount by which minBounds is offsetted from 0 0 0
        schemOffset = schemBounds[0]

        ## BLOCK PALETTE
        ## We're doing the block palette early because it's gonna
        ## be useful in knowing which algorithm to use in when saving
        ## the blocks to the schematic
        # Get a cleaner version of the block palette, without the duplicates used
        # for back and forth O(1) access
        self._structure.getBlockPalette()
        cleanBlockPalette = self._structure.getBlockPalette()


        ## BLOCK DATA
        encodedBlockStates = self._getEncodedBlockStates(len(cleanBlockPalette),
                                                         schemDims,
                                                         schemOffset,
                                                         fastSaving)


        ## BLOCK ENTITIES
        blockEntitiesCompounds = \
            [
                self._blockEntityStringToSchemCompound(
                    (
                        position[0] - schemOffset[0],
                        position[1] - schemOffset[1],
                        position[2] - schemOffset[2]
                    ), blockEntityStr
                )
                for position, blockEntityStr in self._structure.getBlockEntities().items()
            ]


        ## Generate the schematic file from the byte array
        ## From Fearless's code which was taken from someone else lOl
        schematic = File({

            'Version': Int(2),
            'DataVersion': Int(version.value),
            'Metadata': Compound({
                'WEOffsetX': Int(schemOffset[0]),
                'WEOffsetY': Int(schemOffset[1]),
                'WEOffsetZ': Int(schemOffset[2]),
                'MCSchematicMetadata': Compound({
                    'Generated': String("Generated with love using Sloimay's MCSchematic Python Library, "
                                        "itself dependant on Valentin Berlier's nbtlib library.")
                })
            }),

            'Height': Short(schemDims[1]),
            'Length': Short(schemDims[2]),
            'Width': Short(schemDims[0]),

            'PaletteMax': Int(len(cleanBlockPalette)),
            'Palette': Compound(cleanBlockPalette),
            'BlockData': ByteArray(encodedBlockStates),

            'BlockEntities': List(blockEntitiesCompounds),

        }, gzipped=True, root_name='Schematic')
        #################### End copied from MCSchematic #####################

        schematic.save(filepath)
        

    def saveNBT(self, filepath: str | os.PathLike, version : 'Version' = None, maxSize: int | tuple[int, int, int] | None = None, filenameMode: str = "auto", removeAir = False, shifted=True, rewriteExisting=True, ):
        """
        Save the structure as one or more Minecraft schematic .nbt files in <directory>.
        If the structure exceeds maxSize in any dimension, it will be split into multiple files.

        Parameters
        ----------
        filepath : str or os.PathLike
            Directory to save the .nbt files. If multiple files are created, they will be named
            <filepath>_x_y_z.nbt where x,y,z are the tile indices.
        version : 'Version', optional
            The Version to include in the NBT file. If None, the latest version will be used.
        maxSize : int | tuple | None, optional
            Maximum size in each dimension for a single .nbt file. If None, a single file will be created. Default is None.
        filenameMode : str, optional
            "auto" (default): use base_name.nbt if only one file is needed, otherwise use indexed names.
            "indexed": always use indexed names.
        shifted : bool, optional
            If True (default), block positions are written relative to the structure's
            minimum bound, so the saved schematic starts at (0, 0, 0). This is the usual
            schematic behaviour. If False, positions keep their original world coordinates
            (the box is sized from 0 so they still fit); requires non-negative bounds.
        rewriteExisting : bool, optional
            If True (default), existing files will be overwritten. If False, existing files will not be overwritten.
        """
        directory = os.path.dirname(filepath)
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        if not os.path.exists(directory):
            os.makedirs(directory)
        if version is None:
            version = self.getLatestVersion()
        

        x_min, y_min, z_min = self._structure.getBounds()[0]
        x_max, y_max, z_max = self._structure.getBounds()[1]
        # getBounds() returns the *inclusive* max corner. Everything below treats the
        # max as an exclusive [x_min, x_max) upper bound, so bump it by 1; otherwise the
        # entire max-edge layer of blocks is dropped and `size` is short by 1 per axis.
        x_max, y_max, z_max = x_max + 1, y_max + 1, z_max + 1

        # `shifted` only affects the origin that block positions are written relative
        # to (applied per tile below); it does NOT change how blocks are binned, which
        # always happens in world coordinates. When unshifted we pin the tiling origin
        # to 0 so the saved box spans (0,0,0)..(x_max,y_max,z_max) and still contains
        # the blocks sitting at their world positions.
        if not shifted:
            # The vanilla structure format cannot express negative positions, so
            # preserving world coordinates only works for non-negative bounds.
            if x_min < 0 or y_min < 0 or z_min < 0:
                raise ValueError(
                    f"Negative minimum bounds are not supported when shifted=False. "
                    f"Found minimum bound of ({x_min}, {y_min}, {z_min}). Consider "
                    f"adjusting the structure to have non-negative bounds. Or fix the code."
                )
            x_min, y_min, z_min = 0, 0, 0

        size_x, size_y, size_z = x_max - x_min, y_max - y_min, z_max - z_min
        if maxSize is None:
            maxSize = (size_x, size_y, size_z)
        elif isinstance(maxSize, int):
            maxSize = (maxSize, maxSize, maxSize)
        nx = (size_x + maxSize[0] - 1) // maxSize[0]
        ny = (size_y + maxSize[1] - 1) // maxSize[1]
        nz = (size_z + maxSize[2] - 1) // maxSize[2]

        # Bin every block into its target tile in a single O(N) pass. Previously
        # getBlocks() (a full structure rebuild) was called once per tile, and each
        # tile re-scanned all blocks just to filter by its bounds, making the whole
        # save O(tiles * N) -- catastrophic for a small maxSize. Now it is O(N + tiles).
        # Binning is done entirely in world coordinates; the `shifted` origin is applied
        # per tile below, so no separate shifting of (x, y, z) is needed here.
        blocks_by_tile = defaultdict(list)
        for (x, y, z), block in self.getBlocks().items():
            # x_max/y_max/z_max are exclusive (bumped above), so the max-edge blocks pass.
            if not (x_min <= x < x_max and y_min <= y < y_max and z_min <= z < z_max):
                continue
            ix = (x - x_min) // maxSize[0]
            iy = (y - y_min) // maxSize[1]
            iz = (z - z_min) // maxSize[2]
            blocks_by_tile[(ix, iy, iz)].append((x, y, z, block))

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):

                    # Code to find if file exists already:
                    if not rewriteExisting:
                        if any(Path(directory).glob(f"*{ix}_{iy}_{iz}.nbt")):
                            print(f"Skipping existing tile {ix}_{iy}_{iz}")
                            continue


                    x0, y0, z0 = ix*maxSize[0] + x_min, iy*maxSize[1] + y_min, iz*maxSize[2] + z_min
                    x1, y1, z1 = min(x0+maxSize[0], x_max), min(y0+maxSize[1], y_max), min(z0+maxSize[2], z_max)
                    # Origin block positions are written relative to:
                    #   shifted=True  -> this tile's corner, so the file starts at (0,0,0)
                    #   shifted=False -> world origin, so positions keep their world coords
                    # Sizing from `ox` (rather than `x0`) guarantees every position lands
                    # inside [0, size), keeping the NBT valid in both modes.
                    ox, oy, oz = (x0, y0, z0) if shifted else (0, 0, 0)
                    tile_size = (x1-ox, y1-oy, z1-oz)

                    root = Compound()
                    root["DataVersion"] = Int(version.value)
                    root["size"] = List[Int]([Int(tile_size[0]), Int(tile_size[1]), Int(tile_size[2])])

                    # Palette
                    palette = List[Compound]()
                    palette_index = {}
                    # Always include air
                    palette.append(Compound({"Name": String("minecraft:air")}))
                    palette_index["minecraft:air"] = 0
                    next_index = 1

                    # Blocks
                    blocks = List[Compound]()
                    has_non_air = False

                    for (x, y, z, block) in blocks_by_tile.get((ix, iy, iz), ()):
                        rel = (x - ox, y - oy, z - oz)

                        if block not in palette_index:
                            block_name = block.split("[")[0]
                            entry = Compound({"Name": String(block_name)})
                            if block.find("[") != -1:
                                props_str = block[block.find("[")+1:block.find("]")]
                                props = Compound()
                                for prop in props_str.split(","):
                                    if "=" in prop: # Safety check for malformed strings
                                        key, value = prop.split("=")
                                        props[key] = String(value)
                                entry["Properties"] = props
                            palette.append(entry)
                            palette_index[block] = next_index
                            next_index += 1

                        if block != "minecraft:air":
                            has_non_air = True

                        state = palette_index[block]
                        btag = Compound()
                        btag["state"] = Int(state)
                        btag["pos"] = List[Int]([Int(rel[0]), Int(rel[1]), Int(rel[2])])

                        # Block entities will be messed up if shifted is false
                        if (x,y,z) in self._structure._blockEntities:
                            blockEntityString = self._structure._blockEntities[(x,y,z)]
                            if "{" in blockEntityString:
                                nbtPortion = blockEntityString[blockEntityString.find("{"):]
                                btag["nbt"] = parse_nbt(nbtPortion)

                        blocks.append(btag)

                    root["palette"] = palette

                    # Added by Joshua: skip this tile if it contains only air (or nothing).
                    if removeAir and not has_non_air:
                        if nx>1 or ny>1 or nz>1 or filenameMode == "indexed":
                            fname = f"{base_name}_{ix}_{iy}_{iz}.nbt"
                        else:
                            fname = f"{base_name}.nbt"
                        print(f"Skipped saving {fname} because it contains only air blocks.")
                        continue


                    root["blocks"] = blocks
                    root["entities"] = List[Compound]() # TODO: Add support for entities

                    # Save
                    if nx>1 or ny>1 or nz>1 or filenameMode == "indexed":
                        fname = f"{base_name}_{ix}_{iy}_{iz}.nbt"
                    elif filenameMode == "auto":
                        fname = f"{base_name}.nbt"
                    path = os.path.join(directory, fname)
                    
                    nbt_file = File(root, gzipped=True, root_name="Schematic")
                    nbt_file.save(path)
                    print(f"Saved {path} ({tile_size})")

    def toMesh(self, blockColormap: str | BlockColormap | None = "all") -> pv.UnstructuredGrid:
        """
        Convert to a pyvista UnstructuredGrid mesh for visualization.

        Parameters
        ----------
        blockColormap : str, BlockColormap, or None, optional
            Colormap to use for converting block states to colors. Valid strings are 'standard', 'all', 'smooth' 
            or a path to a custom colormap CSV file. Default is "all".
        
        Returns
        -------
        pyvista.UnstructuredGrid
        """
        bounds = np.array(self.getStructure().getBounds())
        shape = bounds[1] - bounds[0] + 1 # +1 because bounds are inclusive
        blocks = self.getBlocks()
        cmap = get_block_colormap(blockColormap) if blockColormap is not None else None
        mask_3d = np.zeros(shape, dtype=bool)
        if cmap is not None:
            rgba_3d = np.zeros(tuple(shape) + (4,), dtype=np.uint8)
        for pos, block in blocks.items():
            if block == "minecraft:air":
                continue
            pos = tuple(np.array(pos) - bounds[0])
            mask_3d[pos] = True
            if cmap is not None:
                rgba_3d[pos] = cmap.get_rgba(block)

        grid = pv.ImageData(dimensions=np.array(shape)+1, spacing=(1, 1, 1))

        grid.cell_data["mask"] = mask_3d.flatten(order='F')
        if cmap is not None:
            grid.cell_data["colors"] = rgba_3d.reshape(-1, 4, order='F')
        mesh: pv.UnstructuredGrid = grid.threshold(0.5, scalars="mask", invert=False)
        mesh.cell_data.pop("mask")
        mesh.translate(bounds[0], inplace=True) # translate back to world coordinates
        return mesh
    
    def show(self, blockColormap: str | BlockColormap | None = "all", plotter: pv.Plotter | None = None, **kwargs):
        """
        Visualize using pyvista.

        Parameters
        ----------
        blockColormap : str, BlockColormap, or None, optional
            Colormap to use for converting block states to colors. Valid strings are 'standard', 'all', 'smooth' 
            or a path to a custom colormap CSV file. Default is "all".
        plotter : pyvista.Plotter, optional
            An existing pyvista Plotter to use for visualization. If None, a new Plotter will be created. Default is None.
        **kwargs
            Additional keyword arguments to pass to pyvista.Plotter.show()
        
        Returns
        -------
        pyvista.Plotter
            The Plotter used for visualization.
        """
        mesh = self.toMesh(blockColormap=blockColormap)
        p = plotter if plotter is not None else pv.Plotter()
        if blockColormap is None:
            p.add_mesh(mesh)
        else:
            p.add_mesh(mesh, scalars="colors", rgb=True)
        p.camera_position = 'xy'
        p.camera.elevation = 45
        p.show_axes()
        p.show(**kwargs)
        return p