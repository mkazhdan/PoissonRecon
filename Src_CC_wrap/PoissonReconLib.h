//##########################################################################
//#                                                                        #
//#               CLOUDCOMPARE WRAPPER: PoissonReconLib                    #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 or later of the License.      #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#               COPYRIGHT: Daniel Girardeau-Montaut                      #
//#                                                                        #
//##########################################################################

#ifndef CC_POISSON_RECON_LIB_WRAPPER
#define CC_POISSON_RECON_LIB_WRAPPER

//! Wrapper to use PoissonRecon (Kazhdan et. al) as a library
class PoissonReconLib
{
public:

	//! Algorithm parameters
	struct Parameters
	{
		//! Default initializer
		Parameters();

		//! The maximum depth of the tree that will be used for surface reconstruction
		/** Running at depth d corresponds to solving on a 2^d x 2^d x 2^d.
			Note that since the reconstructor adapts the octree to the sampling density,
			the specified reconstruction depth is only an upper bound.
		**/
		int depth;

		//! The depth beyond which the octree will be adapted.
		/** At coarser depths, the octree will be complete, containing all 2^d x 2^d x 2^d nodes.
		**/
		int fullDepth;
		
		//! The depth up to which a conjugate-gradients solver will be used to solve the linear system.
		/** Beyond this depth Gauss-Seidel relaxation will be used.
		**/
		int cgDepth;

		//! The ratio between the diameter of the cube used for reconstruction and the diameter of the samples' bounding cube.
		/** Specifies the factor of the bounding cube that the input samples should fit into.
		**/
		float scale;

		//! The minimum number of sample points that should fall within an octree node as the octree construction is adapted to sampling density.
		/** This parameter specifies the minimum number of points that should fall within an octree node.
			For noise-free samples, small values in the range [1.0 - 5.0] can be used. For more noisy samples, larger values
			in the range [15.0 - 20.0] may be needed to provide a smoother, noise-reduced, reconstruction.
		**/
		float samplesPerNode;

		//! The importance that interpolation of the point samples is given in the formulation of the screened Poisson equation.
		/** The results of the original (unscreened) Poisson Reconstruction can be obtained by setting this value to 0.
		**/
		float pointWeight;

		//! The (maximum if CG) number of solver iterations
		int iters;

		//! This parameter specifies the number of threads across which the solver should be parallelized
		int threads;

		//! If this flag is enabled, the size of a sample's normals is used as a confidence value, affecting the sample's constribution to the reconstruction process
		bool confidence;
		//! If this flag is enabled, the sampling density is written out with the vertices
		bool density;

		//! Pull factor for color interpolation
		float colorInterp;

		//DGM: the parameters below are not documented in PoissonRecon

		bool showResidual;
		int kernelDepth;
		int maxSolveDepth;
		
		enum BoundaryType { FREE, DIRICHLET, NEUMANN };
		BoundaryType boundary;

		//DGM: the below parameters are hidden in PoissonRecon

		//! If this flag is enabled, the isosurface extraction does not add a planar polygon's barycenter in order to ensure that the output mesh is manifold
		bool nonManifold;
		//! This flag specifies the accuracy cut-off to be used for CG
		float cgAccuracy;
		//! This flag specifies the exponent scale for the adaptive weighting
		int adaptiveExp;

		float width = 0;
		bool linear_fit = false;
	};

	template <typename Real>
	class ICloud
	{
	public:
		virtual size_t size() const = 0;
		virtual bool hasNormals() const = 0;
		virtual bool hasColors() const = 0;
		virtual void getPoint(size_t index, Real* coords) const = 0;
		virtual void getNormal(size_t index, Real* coords) const = 0;
		virtual void getColor(size_t index, Real* rgb) const = 0;
	};

	template <typename Real>
	class IMesh
	{
	public:
		virtual void addVertex(const Real* coords) = 0;
		virtual void addNormal(const Real* coords) = 0;
		virtual void addColor(const Real* rgb) = 0;
		virtual void addDensity(double d) = 0;
		virtual void addTriangle(size_t i1, size_t i2, size_t i3) = 0;
	};


	static bool Reconstruct(const Parameters& params,
							const PoissonReconLib::ICloud<float>& inCloud,
							PoissonReconLib::IMesh<float>& ouMesh);
};

#endif // CC_POISSON_RECON_LIB_6_11_WRAPPER
