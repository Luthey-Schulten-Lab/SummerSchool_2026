from __future__ import annotations

import csv
from typing import Sequence
import numpy as np
import os
from scipy.spatial import KDTree
from importlib.resources import files
from PIL import ImageColor

_loaded_internal_colormaps = {}
_INTERNAL_COLORMAPS = ["standard", "all", "smooth"]
_HEADER_ITEMS = ["block_state", "r", "g", "b", "a"]

class BlockColormap:
    def __init__(self, csv_path: str | os.PathLike):
        """
        Initializes the BlockColormap by loading block states and their corresponding colors from a CSV file.
        The CSV must have the following columns: 'block_state', 'r', 'g', 'b', 'a'. The RGBA values should be in the range [0, 255].
        """
        bs : list[str] = []
        cs : list[tuple[int, int, int]] = []
        alphas : list[int] = []
        self._dict : dict[str, tuple[int, int, int, int]] = {} # reverse lookup for get_rgb, get_alpha
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            header_indexes = {}
            for item in _HEADER_ITEMS:
                if item not in header:
                    raise ValueError(f"Colormap CSV is missing required column '{item}'")
                header_indexes[item] = header.index(item)
            for row in reader:
                block_state = row[header_indexes["block_state"]]
                r = int(float(row[header_indexes["r"]]))
                g = int(float(row[header_indexes["g"]]))
                b = int(float(row[header_indexes["b"]]))
                a = int(float(row[header_indexes["a"]]))
                bs.append(block_state)
                cs.append((r, g, b))
                alphas.append(a)
                self._dict[block_state] = (r, g, b, a)
        self._block_states = np.array(bs)
        self._colors = np.array(cs)
        self._tree = KDTree(self._colors)
        self._alphas = np.array(alphas) # for now we don't query based on alpha
        self._map_cache : dict[tuple[int, int, int], str] = {} # cache for get_block to speed up repeated lookups

    def get_block_state(self, rgb : Sequence[int] | np.ndarray | str) -> str:
        """
        Returns the closest matching block_state for the given (R, G, B) color or array of colors.
        Caching is used for single color lookups to speed up repeated queries.
        """
        rgb = np.asarray(rgb)
        if np.issubdtype(rgb.dtype, np.number):
            rgb = rgb.astype(int)
        else:
            vfunc = np.vectorize(ImageColor.getrgb) # convert color names to RGB tuples
            rgb = np.asarray(vfunc(rgb))
        rgb = rgb[..., :3]
        if rgb.ndim == 1:
            try:
                return self._map_cache[tuple(rgb)]
            except KeyError:
                pass
        # ignore alpha channel if present
        _dist, index = self._tree.query(rgb, k=1)
        if rgb.ndim == 1:
            self._map_cache[tuple(rgb)] = self._block_states[index]
        return self._block_states[index]
    
    def get_rgba(self, block_state: str) -> np.ndarray:
        """
        Returns the (R, G, B, A) color for the given block_state.
        """
        if block_state not in self._dict:
            raise ValueError(f"Block state {block_state} not found in colormap.")
        return np.array(self._dict[block_state])
    
    def get_rgb(self, block_state: str) -> np.ndarray:
        """
        Returns the (R, G, B) color for the given block_state.
        """
        return self.get_rgba(block_state)[:3]
    
    def get_alpha(self, block_state: str) -> int:
        """
        Returns the alpha (transparency) value for the given block_state.
        """
        return self.get_rgba(block_state)[3]
    
def get_block_colormap(name: str | os.PathLike | BlockColormap) -> BlockColormap:
    """
    Returns a BlockColormap object from a provided colormap or a custom CSV.

    Parameters
    ----------
    name : str, os.PathLike, or BlockColormap
        Name of the colormap: 'standard', 'all', 'smooth' or a path to a custom 
        colormap CSV file. If a BlockColormap object is provided, it is returned directly.

    Returns
    -------
    colormap : BlockColormap
        The requested BlockColormap object.
    """
    if isinstance(name, BlockColormap):
        return name
    if ".csv" in str(name):
        return BlockColormap(name)
    if name not in _INTERNAL_COLORMAPS:
        raise ValueError(f"Colormap '{name}' not found. Available colormaps: {_INTERNAL_COLORMAPS}")
    if name not in _loaded_internal_colormaps:
        _loaded_internal_colormaps[name] = BlockColormap(files("mcschematic_plus").joinpath(f"data/block_colormaps/{name}.csv"))
    return _loaded_internal_colormaps[name]