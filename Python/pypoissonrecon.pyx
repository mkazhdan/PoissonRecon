#cython: language_level=3str

## 2021.07.21: code adapted from https://github.com/mmolero/pypoisson
## https://github.com/mmolero/pypoisson/blob/4f6a651d20ff2bab5ec0d0f8d0a14cc6cef3ef51/src/pypoisson.pyx

cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdlib cimport malloc, free
from libcpp cimport bool

np.import_array()

cdef extern from "../Src/PoissonReconLib.h":
    cdef int PoissonReconLib(int argc, char* argv[])
    cdef vector[double] double_data
    cdef vector[int] int_data
    cdef vector[double] mem_data

def poisson_reconstruction(points, normals,
                           depth=8, full_depth=5, scale=1.1,
                           samples_per_node=1.0,
                           enable_polygon_mesh=False, enable_density=False,
                           nthreads=8, verbose=False):

    """
    Python Binding of Poisson Surface Reconstruction
    Usage:

        faces, vertices = poisson_reconstruction(points, normals, depth=10)


    Parameters
    ----------

    points: array-like

        list of oriented vertices with the x-, y-, and z-coordinates of the positions

    normals: array-like

        list of the x-, y-, and z-coordinates of the normals

    depth: Integer

        This integer is the maximum depth of the tree that will be used for surface reconstruction.
        Running at depth d corresponds to solving on a voxel grid whose resolution is no larger than 2^d x 2^d x 2^d.
        Note that since the reconstructor adapts the octree to the sampling density, the specified reconstruction depth is only an upper bound.
        The default value for this parameter is 8.

    full_depth: Integer

        This integer specifies the depth beyond depth the octree will be adapted.
        At coarser depths, the octree will be complete, containing all 2^d x 2^d x 2^d nodes.
        The default value for this parameter is 5.

    scale: float

        This floating point value specifies the ratio between the diameter of the cube used for reconstruction and the diameter of the samples' bounding cube.
        The default value is 1.1.

    samples_per_node: float

        This floating point value specifies the minimum number of sample points that should fall within an octree node as the octree construction is adapted to sampling density.
        For noise-free samples, small values in the range [1.0 - 5.0] can be used.
        For more noisy samples, larger values in the range [15.0 - 20.0] may be needed to provide a smoother, noise-reduced, reconstruction.
        The default value is 1.0.

    enable_polygon_mesh: Bool
        Enabling this flag tells the reconstructor to output a polygon mesh (rather than triangulating the results of Marching Cubes).
        The default value for this parameter is False.

    enable_density: Bool
        Enabling this flag tells the reconstructor to output the estimated depth values of the iso-surface vertices
        The default value for this parameter is False.

    nthreads: int
        This parameter specifies the number of threads across which the solver should be parallelized.
        The default value for this parameter is 8.

    verbose: Bool
        Enable verbose mode.



    Returns
    -------

    faces: array-like
        faces of the reconstructed mesh

    vertices: array-like

        vertices of the reconstructed mesh




    """

    return _poisson_reconstruction(np.ascontiguousarray(np.float64(points)),
                                   np.ascontiguousarray(np.float64(normals)),
                                   depth, full_depth, scale, samples_per_node,
                                   enable_polygon_mesh, enable_density,
                                   nthreads, verbose)



cdef _poisson_reconstruction(np.float64_t[:, ::1] points,
                             np.float64_t[:, ::1] normals,
                             int depth=8,
                             int full_depth=5,
                             double scale=1.10,
                             double samples_per_node=1.0,
                             bool enable_polygon_mesh=False,
                             bool enable_density=False,
                             int nthreads=0,
                             bool verbose=False):

    cdef:
        char **c_argv
        string arg_depth = str(depth).encode()
        string arg_full_depth = str(full_depth).encode()
        string arg_scale = str(scale).encode()
        string arg_samples_per_node = str(samples_per_node).encode()
        string arg_nthreads = str(nthreads).encode()

    int_data.clear()
    double_data.clear()
    mem_data.clear()

    point_nrows, point_ncols = np.shape(points)
    normal_nrows, normal_ncols = np.shape(normals)

    mem_data.resize(point_ncols * point_nrows + normal_ncols * normal_nrows)

    for i in range(point_nrows):
        for j in range(point_ncols):
            mem_data[j +  i*(point_ncols + normal_ncols)] = points[i,j]
            mem_data[j + point_ncols + i *(point_ncols + normal_ncols)] = normals[i,j]


    args = [b"PoissonRecon", b"--in", b"none", b"--out", b"none",
                             b"--depth", arg_depth.c_str(),
                             b"--fullDepth", arg_full_depth.c_str(),
                             b"--scale",   arg_scale.c_str(),
                             b"--samplesPerNode",  arg_samples_per_node.c_str(),
                             b"--threads", arg_nthreads.c_str()
                             ]

    if verbose == True:
        args += [b"--verbose"]
    if enable_polygon_mesh:
        args += [b"--polygonMesh"]
    if enable_density:
        args += [b"--density"]

    c_argv = <char**> malloc(sizeof(char*) * len(args))
    for idx, s in enumerate(args):
        c_argv[idx] = s

    try:
        PoissonReconLib(len(args), c_argv)
    finally:
        free(c_argv)


    face_cols, vertex_cols = 3, 3
    face_rows = int_data.size() // face_cols
    vertex_rows = double_data.size() // vertex_cols

    cdef int *ptr_faces = &int_data[0]
    cdef double *ptr_vertices = &double_data[0]

    faces = np.zeros((face_rows*face_cols,), dtype=np.int32 )
    vertices = np.zeros((vertex_rows*vertex_cols,), dtype=np.float64)

    for i in range(face_rows*face_cols):
        faces[i] = ptr_faces[i]

    for i in range(vertex_rows*vertex_cols):
        vertices[i] = ptr_vertices[i]

    int_data.clear()
    double_data.clear()

    return faces.reshape(face_rows,face_cols), vertices.reshape(vertex_rows,vertex_cols)
