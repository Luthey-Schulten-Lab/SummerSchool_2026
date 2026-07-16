from __future__ import annotations

import numpy as np
import pyvista as pv
import tifffile
from scipy.spatial import KDTree

def to_mc_volume(arr: np.ndarray) -> np.ndarray:
    """
    Convert raw array in numpy/ImageJ format (z, y, x) = (Up, South, East) to Minecraft coordinates (x, y, z) = (East, Up, South).
    """
    return np.moveaxis(arr, 2, 0)

def read_tiff(path: str) -> np.ndarray:
    """
    Read a multi-page TIFF file (z, y, x) = (Up, South, East) and return a volume in Minecraft coordinates (x, y, z) = (East, Up, South).
    """
    return to_mc_volume(tifffile.imread(path))

def read_npy(path: str) -> np.ndarray:
    """
    Read a NumPy .npy or .npz file (z, y, x) = (Up, South, East) and return a volume in Minecraft coordinates (x, y, z) = (East, Up, South).
    """
    return to_mc_volume(np.load(path))

def read_mesh(filename: str,
              spacing: float = 1.0,
              minimum: tuple = None,
              maximum: tuple = None,
              origin: tuple = None,
              ignore_clip=False,
              fill: bool = True,
              flip_z: bool = True,
              edge_mode: str = "center",
              method: str = "implicit_distance",
              compute_scalars: bool = True,
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Load a mesh (x, y, z) = (East, North, Up) and convert to a boolean voxel volume in Minecraft coordinates (x, y, z) = (East, Up, South).

    Parameters
    ----------
    filename : str
        Path to the mesh file.
    spacing : float or (3,) tuple, optional
        Absolute voxel spacing in mesh coordinates.
    minimum : (3,) tuple, optional
        Minimum bounds of the voxel grid in mesh coordinates. If None, uses mesh minimum bounds.
    maximum : (3,) tuple, optional
        Maximum bounds of the voxel grid in mesh coordinates. If None, uses mesh maximum bounds.
    origin : (3,) tuple, optional
        Origin in mesh coordinates. This is used to calculate the `position` output. If None, uses `minimum`. Default is None.
    ignore_clip : bool, optional
        If True, does not check if the mesh fits completely within the voxel grid. Default is False.
    fill : bool, optional
        If True, fill interior voxels. If False, only surface voxels. Default is True.
    flip_z : bool, optional
        If True, flip the z-axis in the output volume to match Minecraft conventions such that (x, y, z) = (East, Up, South). Default is True.
    edge_mode : {'center', 'inner', 'outer'}, optional
        'center' selects voxels whose centers are within the mesh surface,
        'inner' selects voxels whose centers are inside or touching the surface from inside,
        'outer' selects voxels whose centers are outside or touching the surface from outside.
    method : {'mesh_to_sdf', 'implicit_distance'}, optional
        Method to compute signed distance field. 'mesh_to_sdf' is more robust for complex meshes. 
        'implicit_distance' uses PyVista's built-in method. Default is "implicit_distance".
    compute_scalars : bool, optional
        If True, compute the nearest mesh scalar value at each voxel. For example, RGB colors. Default is True.
    
    Returns
    -------
    voxel_mask : (x, y, z) ndarray
        3D boolean mask indicating voxel occupancy in Minecraft coordinates.
    position : (3,) ndarray
        Position of the voxel grid origin (zero index) in Minecraft coordinates relative to `origin`.
    scalars : (x, y, z, N) ndarray or None
        Scalar value at each voxel, if `compute_scalars` is True. Otherwise, None.
    """
    vol, position, scalars = voxelize_mesh(
        filename,
        spacing=spacing,
        minimum=minimum,
        maximum=maximum,
        origin=origin,
        ignore_clip=ignore_clip,
        fill=fill,
        flip_y=flip_z,
        edge_mode=edge_mode,
        method=method,
        compute_scalars=compute_scalars,
    )
    voxel_mask = to_mc_volume(vol)
    position = position[[2, 0, 1]]  # (Up, South, East) to (East, Up, South)
    scalars = to_mc_volume(scalars) if scalars is not None else None
    return voxel_mask, position, scalars

def voxelize_mesh(mesh: pv.PolyData | str,
                   spacing: float | tuple = 1.0,
                   minimum: tuple = None,
                   maximum: tuple = None,
                   origin: tuple = None,
                   ignore_clip=False,
                   fill: bool = True,
                   flip_y: bool = True,
                   edge_mode: str = "center",
                   method: str = "implicit_distance",
                   compute_scalars: bool = True,
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Load a mesh (x, y, z) = (East, North, Up) and convert to a boolean voxel volume (z, y, x) = (Up, South, East).
    
    Parameters
    ----------
    mesh : pyvista.PolyData or str
        The surface mesh or path to the mesh file.
    spacing : float or (3,) tuple, optional
        Absolute voxel spacing in mesh coordinates.
    minimum : (3,) tuple, optional
        Minimum bounds of the voxel grid in mesh coordinates. If None, uses mesh minimum bounds.
    maximum : (3,) tuple, optional
        Maximum bounds of the voxel grid in mesh coordinates. If None, uses mesh maximum bounds.
    origin : (3,) tuple, optional
        Origin in mesh coordinates. This is used to calculate the `position` output. If None, uses `minimum`. Default is None.
    ignore_clip : bool, optional
        If True, does not check if the mesh fits completely within the voxel grid. Default is False.
    fill : bool, optional
        If True, fill interior voxels. If False, only surface voxels. Default is True.
    flip_y : bool, optional
        If True, flip the y-axis in the output volume to match numpy conventions such that (z, y, x) = (Up, South, East). Default is True.
    edge_mode : {'center', 'inner', 'outer'}, optional
        'center' selects voxels whose centers are within the mesh surface,
        'inner' selects voxels whose centers are inside or touching the surface from inside,
        'outer' selects voxels whose centers are outside or touching the surface from outside.
    method : {'mesh_to_sdf', 'implicit_distance'}, optional
        Method to compute signed distance field. 'mesh_to_sdf' is more robust for complex meshes. 
        'implicit_distance' uses PyVista's built-in method. Default is "implicit_distance".
    compute_scalars : bool, optional
        If True, compute the nearest mesh scalar value at each voxel. For example, RGB colors. Default is True.

    Returns
    -------
    voxel_mask : (z, y, x) ndarray
        3D boolean mask indicating voxel occupancy.
    position : (3,) ndarray
        Position of the voxel grid origin (zero index) in voxel coordinates relative to `origin`.
    scalars : (z, y, x, N) ndarray or None
        Scalar value at each voxel, if `compute_scalars` is True. Otherwise, None.
    """
    # Ensure triangular mesh TODO: is this necessary?
    mesh = pv.read(mesh) if isinstance(mesh, str) else mesh
    if isinstance(mesh, pv.MultiBlock):
        mesh = mesh.extract_surface(algorithm=None)
    elif not isinstance(mesh, pv.PolyData):
        mesh = pv.PolyData(mesh) # last resort
    mesh = mesh.triangulate()
    mins, maxs = np.array(mesh.bounds[::2], dtype=float), np.array(mesh.bounds[1::2], dtype=float)
    min_mesh = mins.copy()
    if minimum is not None:
        for i in range(3):
            if minimum[i] is not None:
                min_mesh[i] = minimum[i]
    max_mesh = maxs.copy()
    if maximum is not None:
        for i in range(3):
            if maximum[i] is not None:
                max_mesh[i] = maximum[i]
        
    if not ignore_clip:
        if np.any(mins < min_mesh):
            raise ValueError(f"minimum {min_mesh} must be <= the minimum bounds of the mesh {mins}. Use ignore_clip=True to override.")
        if np.any(maxs > max_mesh):
            raise ValueError(f"maximum {max_mesh} must be >= the maximum bounds of the mesh {maxs}. Use ignore_clip=True to override.")
        
    spacing = np.array((spacing, spacing, spacing), float) if np.isscalar(spacing) else np.array(spacing, float)
    nxyz = np.ceil((max_mesh - min_mesh) / spacing).astype(int) + 1  # include endpoint
    grid = pv.ImageData(
        dimensions=nxyz,
        spacing=spacing,
        origin=min_mesh,
    )
    if method == "mesh_to_sdf":
        try:
            import trimesh
            import mesh_to_sdf
        except ImportError:
            raise ImportError("trimesh and mesh_to_sdf are required for 'mesh_to_sdf' method. Please install them via 'pip install trimesh mesh_to_sdf'.")
        trimesh_mesh = trimesh.Trimesh(vertices=mesh.points, faces=mesh.faces.reshape((-1, 4))[:, 1:4])
        points = grid.points
        sdf_values = mesh_to_sdf.mesh_to_sdf(trimesh_mesh, points)
    elif method == "implicit_distance":
        distancegrid = grid.compute_implicit_distance(mesh)
        sdf_values = distancegrid.point_data["implicit_distance"]
    sdf = sdf_values.reshape(grid.dimensions[::-1])  # (z,y,x)
    voxel_diagonal = np.linalg.norm(spacing)

    epsilon = 0.1
    if edge_mode == "center":
        mask = np.abs(sdf) <= voxel_diagonal * (0.5)
    elif edge_mode == "inner":
        mask = (sdf <= 0.0) & (sdf >= -voxel_diagonal * (0.5 + epsilon))
    elif edge_mode == "outer":
        mask = (sdf >= 0.0) & (sdf <= voxel_diagonal * (0.5 + epsilon))
    else:
        raise ValueError("edge_mode must be one of {'center','inner','outer'}")
    if fill:
        inside = sdf < 0.0
        mask = inside | mask
    scalars = None
    if compute_scalars:
        if mesh.active_scalars is None:
            raise ValueError("Mesh must have active scalars to use compute_scalars=True.")
        # get nearest point on the mesh for each voxel and get its scalar values
        mask_grid_points = grid.points[mask.flatten()]
        _, idx = KDTree(mesh.points).query(mask_grid_points)
        scalars = np.zeros(mask.shape + mesh.active_scalars.shape[1:], dtype=mesh.active_scalars.dtype)
        scalars[mask] = mesh.active_scalars[idx]
    position = np.zeros(3, dtype=int)
    if origin is not None:
        vox_min = np.ceil(min_mesh / spacing).astype(int)
        vox_origin = np.ceil(np.array(origin, dtype=float) / spacing).astype(int)
        position = vox_min - vox_origin
        position = position[::-1] # to (z, y, x)
    if flip_y:
        mask = mask[:, ::-1, :]
        position[1] = - ((mask.shape[1] - 1) + position[1])
        if scalars is not None:
            scalars = scalars[:, ::-1, :]
    voxel_mask = mask.astype(bool)

    return voxel_mask, position, scalars