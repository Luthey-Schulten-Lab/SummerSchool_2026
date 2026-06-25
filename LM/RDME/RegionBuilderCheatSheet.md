# `RegionBuilder` Cheat Sheet

`RegionBuilder` is a helper class for designing what is termed the "site type lattice" — the 3D array that labels every voxel with a region type (cytoplasm, membrane, nucleus, etc.). Within the class are methods termed "geometry methods", which return a `numpy` boolean mask of shape `(nx, ny, nz)`. After masks have been created, the user may combine the masks into the actual site lattice at the end with the `compose` function.

---

## Instantiation - Define the RegionBuilder Object ##
When using the `RegionBuilder` package, the user must define a `RegionBuilder` object within their python environment. When instantiating this object, the user can either pass an RDME simulation object, or explicit dimensions to the function. The `__init__()` function for the `RegionBuilder` class is given below: 

**`__init__(net=None, dims=None)`**  
Takes either an RDME simulation object (extracts `nx, ny, nz` from the simulation object) or the user can provide an explicit `dims=[nx, ny, nz]` flag. Exactly one must be provided. Initialization creates a voxel index grid `self.xs` (shape `[3, nx, ny, nz]`), `self.center`, and `self.origin` variables for use by the geometry methods.

---

## Geometry Methods — Create `(nx, ny, nz)` Boolean Masks ##
After defining a `RegionBuilder` object, the user can now use any number of geometry methods to create a boolean mask of the exact dimensions of the RDME simulation object. Every geometry method works by testing every single voxel in the lattice against a mathematical condition and returning `True` for voxels that are inside the shape. Basic features of each gemoetry method are described below:

**`ellipsoid(radius, center=None, angles=None, xs=None)`**  
Constructs an ellipsoid mask. `radius` can be a scalar (sphere) or a 3-element list of semiaxes `[rx, ry, rz]`. Centered at `self.center` by default. Supports arbitrary orientation via ZXZ Euler angles (see [Wiki page](https://en.wikipedia.org/wiki/Euler_angles) for more information). 

**`cylinder(radius, length, center=None, angles=None, xs=None)`**  
Constructs a cylinder aligned along the Z axis (before rotation). `radius` is the circular cross-section radius, `length` is the total height. The cylinder extends from `-length/2` to `+length/2` along the axis. Supports translation and ZXZ Euler rotation.

**`capsule(length, width, center=None, angles=None, xs=None)`**  
Constructs a spherocylinder (pill shape) — a cylinder with hemispherical endcaps. `length` is the total tip-to-tip length, `width` is the diameter of the cylindrical section. Implemented as the union of a cylinder and two spheres placed at each end.

**`box(lx, ly, lz, center=None, angles=None, xs=None)`**  
Constructs a rectangular cuboid of dimensions `lx × ly × lz`. Supports ZXZ Euler rotation. The box is oriented with its corner at the origin in the rotated frame (note the `>= 0` bounds — it is not centered).

**`spoke(x0, length, spoke_radius, r, phi, theta)`**  
Constructs a cylinder placed on the surface of a sphere of radius `r` centered at `x0`, at spherical coordinates `(phi, theta)`. The cylinder is automatically oriented radially (normal to the sphere surface) and protrudes `length/2` inward and outward. Useful for modeling flagella anchors or similar surface-attached structures.

---

## Morphological Operations — Modify Existing Masks

All morphological operations use `scipy.ndimage` binary morphology internally to restructure an already-existing boolean mask. The operations work by sliding a structuring element (SE) — a small 3D boolean array that acts as a "brush shape" — over every voxel in the mask. Depending on the method (`dialate`, `erode`, etc.), the center voxel and/or all voxels within the SE are tested for being either `True` or `False`. Then, according to the method's specific criteria, the center voxel is changed if a certain condition has been met. For example, with `dialate`, if the voxel being visited in the boolean array is `True`, then every other voxel within the SE will be changed to also be `True`. This is performed on every voxel in the lattice, and a new boolean mask is then returned.

Every morphological function requires an SE configuration and accepts either a `radius` integer or an explicit `se` array — but not both. If `radius` is given, an octahedral SE of that radius is auto-generated (see Structuring Elements section below). If `se` is given, that array is used directly. The following five functions are the available morphological operation methods:

**`dilate(binaryMask, radius=None, se=None)`**  
Expands the `True` region outward by adding a shell of voxels around the boundary of the mask. At each voxel position, the SE is centred there; if the SE overlaps with *any* `True` voxel in the input, the output voxel is set `True`. The size and shape of the SE controls how many voxels are added and in what direction.

- `radius` — grow by `radius` voxels using the default octahedral SE
- `se` — grow using a custom SE array (e.g. `rb.se6`, `rb.se26`, `RegionBuilder.sphereStructElem(r)`)

*Use case:* growing a cytoplasm mask outward to create a membrane layer on top of it.

**`erode(binaryMask, radius=None, se=None)`**  
Shrinks the `True` region inward by removing voxels from the boundary. A voxel is kept `True` in the output only if the SE, centred there, lies *entirely* within the `True` region of the input. Erosion is the inverse of dilation.

- `radius` — shrink by `radius` voxels using the default octahedral SE
- `se` — shrink using a custom SE array

*Use case:* peeling the outer shell off a region to expose its interior, or thinning a membrane mask.

**`closing(binaryMask, radius=None, se=None, radius1=None, se1=None)`**  
Dilation followed by erosion using the same SE. Fills small holes and gaps in the interior of a mask and bridges narrow gaps between nearby regions, while approximately preserving the overall size and shape.

- `radius` / `se` — SE used for the dilation step
- `radius1` / `se1` — if provided, a *different* SE is used for the erosion step; otherwise the same SE is used for both

*Use case:* smoothing a rough or porous region boundary imported from image data.

**`opening(binaryMask, radius=None, se=None, radius1=None, se1=None)`**  
Erosion followed by dilation using the same SE. Removes small protrusions, thin bridges, and isolated specks from a mask while approximately preserving the overall size and shape. It is the complement of closing.

- `radius` / `se` — SE used for the erosion step
- `radius1` / `se1` — if provided, a *different* SE is used for the dilation step; otherwise the same SE is used for both

*Use case:* cleaning up noise or small artifacts around a geometry boundary.

**`convexHull(binaryMask)`**  
Computes the 3D convex hull of all `True` voxels using `scipy.spatial.ConvexHull` and `Delaunay` triangulation. Returns a mask where every voxel inside the convex hull is `True`. Takes no SE — it is a purely geometric operation with no `radius` parameter.

*Use case:* filling concave regions, for example wrapping a mesh-imported organelle shape so that its interior is fully filled regardless of surface concavities.

---

## Structuring Elements — Shape of Morphological Operations

A structuring element (SE) is the "brush" that determines the neighborhood used by each morphological operation. All four options below are 3D boolean arrays. You can pass any of them as the `se` argument to `dilate`, `erode`, `closing`, or `opening`.

**`se6`** *(property)*  
A 3×3×3 array with only the 6 face-adjacent neighbors set to `True` (±x, ±y, ±z faces). No diagonal neighbors are included. One application grows or shrinks the mask by exactly 1 voxel along each axis — the most conservative option.

```
z=-1:    z=0:    z=1:
0 0 0   0 1 0   0 0 0
0 1 0   1 1 1   0 1 0
0 0 0   0 1 0   0 0 0
```

**`se26`** *(property)*  
A 3×3×3 array with all 26 neighbors (faces + edges + corners) set to `True`. One application grows or shrinks by 1 voxel in all directions including diagonals — the most aggressive single-step option.

```
z=-1:    z=0:    z=1:
1 1 1   1 1 1   1 1 1
1 1 1   1 1 1   1 1 1
1 1 1   1 1 1   1 1 1
```

**`octoStructElem(r)`** *(static — the default)*  
Returns an octahedral SE (3D diamond shape) built by iterating the 6-connected `se6` with itself `r` times. All voxels within Manhattan (taxicab) distance `r` of the center are `True`. This is the SE automatically used when you pass `radius=r` to any morphological function.

The octahedral shape is not a perfect sphere — it has a "pointy" diamond appearance at the corners. For most geometry work this is acceptable and it is the fastest option.

| `r` | Voxels affected per side | Shape |
|---|---|---|
| 1 | 1 voxel (same as `se6`) | Small cross |
| 2 | 2 voxels | Larger diamond |
| 3 | 3 voxels | Even larger diamond |

**`sphereStructElem(r)`** *(static)*  
Returns a true Euclidean sphere SE of radius `r`: a voxel is included if its Euclidean distance from the center is strictly less than `r`. More isotropic than the octahedral SE — produces a rounder boundary after morphological operations — but slightly more expensive to compute.

*When to prefer it:* when the octahedral SE produces visible "diamond-like" flat facets on curved surfaces such as cell membranes or organelle boundaries.

### SE Comparison ###

| SE | Shape | Reach | Diagonal growth | Best for |
|---|---|---|---|---|
| `se6` | Axis-aligned cross | 1 voxel, faces only | No | Minimal single-step ops |
| `se26` | Full cube | 1 voxel, all directions | Yes | Fast single-step with diagonals |
| `octoStructElem(r)` | Octahedron (diamond) | `r` voxels, Manhattan | No | General purpose — the default |
| `sphereStructElem(r)` | Sphere | `r` voxels, Euclidean | Partial | Curved/biological geometries |

---

## Grid Utilities ##

**`transformGrid(xs, x0, alpha, beta, gamma)`** *(static)*  
Rotates and translates an index grid `xs` around center `x0` using ZXZ Euler angles `(alpha, beta, gamma)`. Called internally by all geometry methods. Can be used directly to apply the same transform to multiple primitives before combining them.

**`emptyLatticeMask()`**  
Returns a zero-initialized boolean array of shape `(nx, ny, nz)`. Convenience method for building a mask from scratch.

---

## Lattice Assembly ##
After all region type masks are created, adding them to the RDME simulation object can be done using the `compose` function. Before using this function, make sure that the regions types are defined within the RDME simulation object.

**`compose(*siteSpec, net=None)`**  
The final step — writes multiple `(region, mask)` pairs into the simulation's site lattice in order. Later entries overwrite earlier ones where masks overlap, so ordering matters:

```python
rb.compose(
    (cytoplasm, cytoMask),
    (membrane,  membMask),   # overwrites cytoplasm at the surface
    (nucleus,   nucMask),
)
```

---

## Visualization ##
After composing a lattice structure, it can be helpful to visualize it in either 2D or 3D. The following functions aid the user in visualizing the pregenerated lattice:

**`showBinaryLattices(binLattices, manualColor=None, filterFunctions=None, mode="widget")`**
Renders one or more boolean masks as a 3D interactive scene inside a Jupyter notebook. `filterFunctions` can clip the display to a half-volume or plane for inspecting interior geometry. `mode` can be `"widget"`, `"download_x3d"`, or `"download_html"`.

**`showStack(binLattices, plane='xz', scl=None, maxWidth=600, maxHeight=600)`** 
Displays an interactive 2D slice viewer for one or more masks in a Jupyter notebook. `plane` can be `'xy'`, `'yz'`, or `'xz'`. Useful for quickly checking geometry cross-sections without a full 3D render.

---

## Quick Reference ##

| Method | Category | Returns |
|---|---|---|
| `ellipsoid` | Geometry | mask |
| `cylinder` | Geometry | mask |
| `capsule` | Geometry | mask |
| `box` | Geometry | mask |
| `spoke` | Geometry | mask |
| `dilate` | Morphology | mask |
| `erode` | Morphology | mask |
| `closing` | Morphology | mask |
| `opening` | Morphology | mask |
| `convexHull` | Morphology | mask |
| `compose` | Assembly | writes to site lattice |
| `transformGrid` | Utility | rotated index grid |
| `emptyLatticeMask` | Utility | empty mask |
| `octoStructElem` | Structuring element | SE array |
| `sphereStructElem` | Structuring element | SE array |
| `se6` | Structuring element | SE array |
| `se26` | Structuring element | SE array |
| `showBinaryLattices` | Visualization | Jupyter display |
| `showStack` | Visualization | Jupyter display |
