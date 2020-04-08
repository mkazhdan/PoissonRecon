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

#include "PoissonReconLib.h"

//PoissonRecon
#include "../Src/FEMTree.h"

#include "PointData.h"

#include <cassert>

namespace {
	// The order of the B-Spline used to splat in data for color interpolation
	constexpr int DATA_DEGREE = 0;
	// The order of the B-Spline used to splat in the weights for density estimation
	constexpr int WEIGHT_DEGREE = 2;
	// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
	constexpr int NORMAL_DEGREE = 2;
	// The default finite-element degree
	constexpr int DEFAULT_FEM_DEGREE = 1;
	// The dimension of the system
	constexpr int DIMENSION = 3;

	inline float ComputeNorm(const float vec[3])
	{
		return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	}
	
	inline double ComputeNorm(const double vec[3])
	{
		return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	}
}

PoissonReconLib::Parameters::Parameters()
{
#ifdef WITH_OPENMP
	threads = omp_get_num_procs();
#endif
}

template <typename _Real>
class Vertex : public PointData<_Real>
{
public:

	typedef _Real Real;

	Vertex(const Point<Real, 3>& point)
		: PointData<Real>()
		, point(point)
		, w(0)
	{}

	Vertex(const Point<Real, 3>& point, const PointData<Real>& data, double _w = 0.0)
		: PointData<Real>(data.normal, data.color)
		, point(point)
		, w(_w)
	{}

	Vertex() : Vertex(Point<Real, 3>(0, 0, 0))
	{}

	Vertex& operator *= (Real s)
	{
		PointData<Real>::operator *= (s);
		point *= s;
		w *= s;

		return *this;
	}

	Vertex& operator /= (Real s)
	{
		PointData<Real>::operator *= (1 / s);
		point /= s;
		w /= s;
		return *this;
	}

	Vertex& operator+=(const Vertex& p)
	{
		PointData<Real>::operator += (p);
		point += p.point;
		w += p.w;

		return *this;
	}

public:
	Point<Real, 3> point;
	double w;
};

template <typename Real>
class PointStream : public InputPointStreamWithData<Real, DIMENSION, PointData<Real> >
{
public:
	PointStream(const PoissonReconLib::ICloud<Real>& _cloud)
		: cloud(_cloud), xform(nullptr), currentIndex(0)
	{}

	void reset(void) override
	{
		currentIndex = 0;
	}

	bool nextPoint(Point<Real, 3>& p, PointData<Real>& d) override
	{
		if (currentIndex >= cloud.size())
		{
			return false;
		}
		cloud.getPoint(currentIndex, p.coords);

		if (xform != nullptr)
		{
			p = (*xform) * p;
		}

		if (cloud.hasNormals())
		{
			cloud.getNormal(currentIndex, d.normal);
		}
		else
		{
			d.normal[0] = d.normal[1] = d.normal[2];
		}

		if (cloud.hasColors())
		{
			cloud.getColor(currentIndex, d.color);
		}
		else
		{
			d.color[0] = d.color[1] = d.color[2];
		}

		currentIndex++;
		return true;
	}

public:
	const PoissonReconLib::ICloud<Real>& cloud;
	XForm<Real, 4>* xform;
	size_t currentIndex;
};

template <unsigned int Dim, class Real>
struct FEMTreeProfiler {
	FEMTree<Dim, Real>& tree;
	double t;

	FEMTreeProfiler(FEMTree<Dim, Real>& t) : tree(t) {}
	void start(void) {
		t = Time(), FEMTree<Dim, Real>::ResetLocalMemoryUsage();
	}
	void dumpOutput(const char* header) const {
		FEMTree<Dim, Real>::MemoryUsage();
		//if (header) {
		//	utility::LogDebug("{} {} (s), {} (MB) / {} (MB) / {} (MB)", header,
		//		Time() - t,
		//		FEMTree<Dim, Real>::LocalMemoryUsage(),
		//		FEMTree<Dim, Real>::MaxMemoryUsage(),
		//		MemoryInfo::PeakMemoryUsageMB());
		//}
		//else {
		//	utility::LogDebug("{} (s), {} (MB) / {} (MB) / {} (MB)", Time() - t,
		//		FEMTree<Dim, Real>::LocalMemoryUsage(),
		//		FEMTree<Dim, Real>::MaxMemoryUsage(),
		//		MemoryInfo::PeakMemoryUsageMB());
		//}
	}
};

template <class Real, unsigned int Dim>
XForm<Real, Dim + 1> GetBoundingBoxXForm(Point<Real, Dim> min,
	Point<Real, Dim> max,
	Real scaleFactor) {
	Point<Real, Dim> center = (max + min) / 2;
	Real scale = max[0] - min[0];
	for (unsigned int d = 1; d < Dim; d++) {
		scale = std::max<Real>(scale, max[d] - min[d]);
	}
	scale *= scaleFactor;
	for (unsigned int i = 0; i < Dim; i++) {
		center[i] -= scale / 2;
	}
	XForm<Real, Dim + 1> tXForm = XForm<Real, Dim + 1>::Identity(),
		sXForm = XForm<Real, Dim + 1>::Identity();
	for (unsigned int i = 0; i < Dim; i++) {
		sXForm(i, i) = static_cast<Real>(1. / scale), tXForm(Dim, i) = -center[i];
	}
	return sXForm * tXForm;
}

template <class Real, unsigned int Dim>
XForm<Real, Dim + 1> GetBoundingBoxXForm(	Point<Real, Dim> min,
											Point<Real, Dim> max,
											Real width,
											Real scaleFactor,
											int& depth)
{
	// Get the target resolution (along the largest dimension)
	Real resolution = (max[0] - min[0]) / width;
	for (unsigned int d = 1; d < Dim; d++)
	{
		resolution = std::max<Real>(resolution, (max[d] - min[d]) / width);
	}
	resolution *= scaleFactor;
	depth = 0;
	while ((1 << depth) < resolution)
	{
		depth++;
	}

	Point<Real, Dim> center = (max + min) / 2;
	Real scale = (1 << depth) * width;

	for (unsigned int i = 0; i < Dim; i++)
	{
		center[i] -= scale / 2;
	}
	XForm<Real, Dim + 1> tXForm = XForm<Real, Dim + 1>::Identity();
	XForm<Real, Dim + 1> sXForm = XForm<Real, Dim + 1>::Identity();
	for (unsigned int i = 0; i < Dim; i++)
	{
		sXForm(i, i) = static_cast<Real>(1.0 / scale);
		tXForm(Dim, i) = -center[i];
	}
	return sXForm * tXForm;
}

template <class Real, unsigned int Dim>
XForm<Real, Dim + 1> GetPointXForm(	InputPointStream<Real, Dim>& stream,
									Real width,
									Real scaleFactor,
									int& depth)
{
	Point<Real, Dim> min, max;
	stream.boundingBox(min, max);
	return GetBoundingBoxXForm(min, max, width, scaleFactor, depth);
}

template <class Real, unsigned int Dim>
XForm<Real, Dim + 1> GetPointXForm(	InputPointStream<Real, Dim>& stream,
									Real scaleFactor)
{
	Point<Real, Dim> min, max;
	stream.boundingBox(min, max);
	return GetBoundingBoxXForm(min, max, scaleFactor);
}

template <unsigned int Dim, typename Real>
struct ConstraintDual
{
	Real target, weight;
	ConstraintDual(Real t, Real w) : target(t), weight(w) {}
	CumulativeDerivativeValues<Real, Dim, 0> operator()(const Point<Real, Dim>& p) const
	{
		return CumulativeDerivativeValues<Real, Dim, 0>(target * weight);
	};
};

template <unsigned int Dim, typename Real>
struct SystemDual
{
	SystemDual(Real w) : weight(w) {}
	CumulativeDerivativeValues<Real, Dim, 0> operator()(const Point<Real, Dim>& p,
														const CumulativeDerivativeValues<Real, Dim, 0>& dValues) const
	{
		return dValues * weight;
	};
	
	CumulativeDerivativeValues<double, Dim, 0> operator()(	const Point<Real, Dim>& p,
															const CumulativeDerivativeValues<double, Dim, 0>& dValues) const
	{
		return dValues * weight;
	};

	Real weight;
};

template <unsigned int Dim>
struct SystemDual<Dim, double>
{
	typedef double Real;

	SystemDual(Real w) : weight(w) {}
	CumulativeDerivativeValues<Real, Dim, 0> operator()(	const Point<Real, Dim>& p,
															const CumulativeDerivativeValues<Real, Dim, 0>& dValues) const
	{
		return dValues * weight;
	};

	Real weight;
};

template <typename Vertex, typename Real, typename SetVertexFunction, unsigned int... FEMSigs, typename... SampleData>
void ExtractMesh(	const PoissonReconLib::Parameters& params,
					UIntPack<FEMSigs...>,
					std::tuple<SampleData...>,
					FEMTree<sizeof...(FEMSigs), Real>& tree,
					const DenseNodeData<Real, UIntPack<FEMSigs...>>& solution,
					Real isoValue,
					const std::vector<typename FEMTree<sizeof...(FEMSigs), Real>::PointSample>* samples,
					std::vector< PointData<Real> >* sampleData,
					const typename FEMTree<sizeof...(FEMSigs), Real>::template DensityEstimator<WEIGHT_DEGREE>* density,
					const SetVertexFunction& SetVertex,
					XForm<Real, sizeof...(FEMSigs) + 1> iXForm,
					PoissonReconLib::IMesh<Real>& out_mesh)
{
	static const int Dim = sizeof...(FEMSigs);
	typedef UIntPack<FEMSigs...> Sigs;
	static const unsigned int DataSig = FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;

	const bool non_manifold = true;
	const bool polygon_mesh = false;

	CoredVectorMeshData<Vertex, node_index_type> mesh;
	
	if (samples && sampleData)
	{
		typedef typename FEMTree<Dim, Real>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;

		SparseNodeData< ProjectiveData<PointData<Real>, Real>, IsotropicUIntPack<Dim, DataSig>> _sampleData =
			tree.template setMultiDepthDataField<DataSig, false>(	*samples,
																	*sampleData,
																	(DensityEstimator*)nullptr);

		for (const RegularTreeNode<Dim, FEMTreeNodeData, depth_and_offset_type>* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
		{
			ProjectiveData<PointData<Real>, Real>* color = _sampleData(n);
			if (color)
				(*color) *= static_cast<Real>(pow(params.colorPullFactor, tree.depth(n)));
		}
		
		IsoSurfaceExtractor<Dim, Real, Vertex>::template Extract< PointData<Real> >(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, &_sampleData, solution, isoValue, mesh, SetVertex, !params.linearFit, !non_manifold, polygon_mesh, false);
	}
	else
	{
		IsoSurfaceExtractor<Dim, Real, Vertex>::template Extract< PointData<Real> >(Sigs(), UIntPack<WEIGHT_DEGREE>(), UIntPack<DataSig>(), tree, density, nullptr, solution, isoValue, mesh, SetVertex, !params.linearFit, !non_manifold, polygon_mesh, false);
	}

	mesh.resetIterator();

	for (size_t vidx = 0; vidx < mesh.outOfCorePointCount(); ++vidx)
	{
		Vertex v;
		mesh.nextOutOfCorePoint(v);
		v.point = iXForm * v.point;

		out_mesh.addVertex(v.point.coords);
		if (sampleData)
		{
			//out_mesh.addNormal(v.normal);
			out_mesh.addColor(v.color);
		}
		if (params.density)
		{
			out_mesh.addDensity(v.w);
		}
	}
	
	for (size_t tidx = 0; tidx < mesh.polygonCount(); ++tidx)
	{
		std::vector<CoredVertexIndex<node_index_type>> triangle;
		mesh.nextPolygon(triangle);
		if (triangle.size() == 3)
		{
			out_mesh.addTriangle(triangle[0].idx, triangle[1].idx, triangle[2].idx);
		}
		else
		{
			assert(false);
		}
	}
}

template <class Real, typename... SampleData, unsigned int... FEMSigs>
static bool Execute(PointStream<Real>& pointStream,
					PoissonReconLib::IMesh<Real>& out_mesh,
					const PoissonReconLib::Parameters& params,
					UIntPack<FEMSigs...> )
{
	static const int Dim = sizeof...(FEMSigs);
	typedef UIntPack<FEMSigs...> Sigs;
	typedef UIntPack<FEMSignature<FEMSigs>::Degree...> Degrees;
	typedef UIntPack<FEMDegreeAndBType<NORMAL_DEGREE, DerivativeBoundary<FEMSignature<FEMSigs>::BType, 1>::BType>::Signature...> NormalSigs;
	typedef typename FEMTree<Dim, Real>::template DensityEstimator<WEIGHT_DEGREE> DensityEstimator;
	typedef typename FEMTree<Dim, Real>::template InterpolationInfo<Real, 0> InterpolationInfo;

	// Compute scaling transformation (and optionally the depth)
	int depth = params.depth;
	XForm<Real, Dim + 1> xForm = XForm<Real, Dim + 1>::Identity();
	{
		if (params.finestCellWidth > 0)
		{
			Real scaleFactor = static_cast<Real>(params.scale > 0 ? params.scale : 1.0);
			xForm = GetPointXForm<Real, Dim>(pointStream, params.finestCellWidth, scaleFactor, depth) * xForm; //warning: depth may change!
		}
		else if (params.scale > 0)
		{
			xForm = GetPointXForm<Real, Dim>(pointStream, static_cast<Real>(params.scale)) * xForm;
		}
		pointStream.xform = &xForm;
	}

	if (depth < 2)
	{
		//depth should be greater than 2
		assert(false);
		return false;
	}

	//default parameters
	const int solve_depth = depth;
	const int full_depth = 5;
	const bool exact_interpolation = false;
	const Real target_value = static_cast<Real>(0.5);

	// Read in the samples (and color data)
	FEMTree<Dim, Real> tree(MEMORY_ALLOCATOR_BLOCK_SIZE);
	typedef std::vector<typename FEMTree<Dim, Real>::PointSample> SampleSet;
	typedef std::vector< PointData<Real> > SampleDataSet;
	std::unique_ptr<SampleSet> samples;
	std::unique_ptr<SampleDataSet> sampleData;
	try
	{
		samples.reset(new SampleSet);
		sampleData.reset(new SampleDataSet);

		if (params.normalConfidence > 0)
		{
			auto ProcessDataWithConfidence = [&](const Point<Real, Dim>& p, PointData<Real>& d)
			{
				Real l = ComputeNorm(d.normal);
				if (std::isnan(l) || l == 0)
					return static_cast<Real>(-1.0);

				return static_cast<Real>(pow(l, params.normalConfidence));
			};

			FEMTreeInitializer<Dim, Real>::template Initialize< PointData<Real> >(tree.spaceRoot(), pointStream, depth, *samples, *sampleData, true, tree.nodeAllocators[0], tree.initializer(), ProcessDataWithConfidence);
		}
		else
		{
			auto ProcessData = [](const Point<Real, Dim>& p, PointData<Real>& d)
			{
				Real l = ComputeNorm(d.normal);
				if (std::isnan(l) || l == 0)
					return static_cast<Real>(-1.0);

				d.normal[0] /= l;
				d.normal[1] /= l;
				d.normal[2] /= l;
				return static_cast<Real>(1.0);
			};

			FEMTreeInitializer<Dim, Real>::template Initialize< PointData<Real> >(tree.spaceRoot(), pointStream, solve_depth, *samples, *sampleData, true, tree.nodeAllocators[0], tree.initializer(), ProcessData);
		}
	}
	catch (std::exception e)
	{
		return false;
	}

	DenseNodeData<Real, Sigs> solution;
	std::unique_ptr<DensityEstimator> density;
	SparseNodeData<Point<Real, Dim>, NormalSigs>* normalInfo = nullptr;
	Real pointWeightSum = 0;
	{
		tree.resetNodeIndices();

		// Get the kernel density estimator
		{
			int kernelDepth = solve_depth - 2;
			assert(kernelDepth >= 0);

			density.reset(tree.template setDensityEstimator<WEIGHT_DEGREE>(*samples, kernelDepth, params.samplesPerNode, 1));
		}

		// Transform the Hermite samples into a vector field
		{
			normalInfo = new SparseNodeData<Point<Real, Dim>, NormalSigs>();
			
			if (params.normalConfidenceBias > 0)
			{
				std::function<bool(PointData<Real>, Point<Real, Dim>&, Real&)> ConversionAndBiasFunction = [&](PointData<Real> in, Point<Real, Dim>& out, Real& bias)
				{
					// Point<Real, Dim> n = in.template data<0>();
					Point<Real, Dim> n(in.normal[0], in.normal[1], in.normal[2]);
					Real l = static_cast<Real>(Length(n));
					// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
					if (l == 0)
						return false;

					out = n / l;
					bias = static_cast<Real>(log(l) * params.normalConfidenceBias / log(1 << (Dim - 1)));

					return true;
				};

				*normalInfo = tree.setDataField(NormalSigs(), *samples, *sampleData, density.get(), pointWeightSum, ConversionAndBiasFunction);
			}
			else
			{
				std::function<bool(PointData<Real>, Point<Real, Dim>&)> ConversionFunction = [](PointData<Real> in, Point<Real, Dim>& out)
				{
					Point<Real, Dim> n(in.normal[0], in.normal[1], in.normal[2]);
					Real l = static_cast<Real>(Length(n));
					// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
					if (l == 0)
						return false;

					out = n / l;
					return true;
				};

				*normalInfo = tree.setDataField(NormalSigs(), *samples, *sampleData, density.get(), pointWeightSum, ConversionFunction);
			}

			auto InvertNormal = [&](unsigned int, size_t i)
			{
				(*normalInfo)[i] *= static_cast<Real>(-1.0);
			};
			
			ThreadPool::Parallel_for(0, normalInfo->size(), InvertNormal);
		}

		if (!params.density)
		{
			density.reset();
		}
		if (!params.withColors || params.colorPullFactor == 0)
		{
			sampleData.reset();
		}

		// Trim the tree and prepare for multigrid
		{
			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
			
			tree.template finalizeForMultigrid<MAX_DEGREE>( full_depth,
															typename FEMTree<Dim, Real>::template HasNormalDataFunctor<NormalSigs>(*normalInfo),
															normalInfo,
															density.get() );
		}

		// Add the FEM constraints
		DenseNodeData<Real, Sigs> constraints;
		{
			constraints = tree.initDenseNodeData(Sigs());
			typename FEMIntegrator::template Constraint<Sigs, IsotropicUIntPack<Dim, 1>, NormalSigs, IsotropicUIntPack<Dim, 0>, Dim> F;
			unsigned int derivatives2[Dim];
			
			for (unsigned int d = 0; d < Dim; d++)
				derivatives2[d] = 0;
			
			typedef IsotropicUIntPack<Dim, 1> Derivatives1;
			typedef IsotropicUIntPack<Dim, 0> Derivatives2;
			for (unsigned int d = 0; d < Dim; d++)
			{
				unsigned int derivatives1[Dim];
				for (unsigned int dd = 0; dd < Dim; dd++)
					derivatives1[dd] = (dd == d ? 1 : 0);
				
				F.weights[d][TensorDerivatives<Derivatives1>::Index(derivatives1)][TensorDerivatives<Derivatives2>::Index(derivatives2)] = 1;
			}
			
			tree.addFEMConstraints(F, *normalInfo, constraints, solve_depth);
		}

		// Free up the normal info
		if (normalInfo)
		{
			delete normalInfo;
			normalInfo = nullptr;
		}

		// Add the interpolation constraints
		InterpolationInfo* iInfo = nullptr;
		if (params.pointWeight > 0)
		{
			if (exact_interpolation)
			{
				iInfo = FEMTree<Dim, Real>::template InitializeExactPointInterpolationInfo<Real, 0>(		tree,
																											*samples,
																											ConstraintDual<Dim, Real>(target_value, static_cast<Real>(params.pointWeight) * pointWeightSum),
																											SystemDual<Dim, Real>(static_cast<Real>(params.pointWeight) * pointWeightSum),
																											true,
																											0);
			}
			else
			{
				iInfo = FEMTree<Dim, Real>::template InitializeApproximatePointInterpolationInfo<Real, 0>(	tree,
																											*samples,
																											ConstraintDual<Dim, Real>(target_value, static_cast<Real>(params.pointWeight) * pointWeightSum),
																											SystemDual<Dim, Real>(static_cast<Real>(params.pointWeight) * pointWeightSum),
																											true,
																											1);
			}
			tree.addInterpolationConstraints(constraints, solve_depth, *iInfo);
		}

		// Solve the linear system
		{
			typename FEMTree<Dim, Real>::SolverInfo sInfo;
			{
				sInfo.cgDepth = 0;
				sInfo.cascadic = true;
				sInfo.vCycles = 1;
				sInfo.iters = params.iters;
				sInfo.cgAccuracy = params.cgAccuracy;
				sInfo.verbose = false;
				sInfo.showResidual = false;
				sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE;
				sInfo.sliceBlockSize = 1;
				sInfo.baseDepth = params.baseDepth;
				sInfo.baseVCycles = params.baseVCycles;
			}
			typename FEMIntegrator::template System<Sigs, IsotropicUIntPack<Dim, 1> > F({ 0.0, 1.0 });
			
			solution = tree.solveSystem(Sigs(), F, constraints, solve_depth,sInfo, iInfo);
		}

		// Free up the interpolation info
		if (iInfo)
		{
			delete iInfo;
			iInfo = nullptr;
		}
	}

	Real isoValue = 0;
	{
		double valueSum = 0, weightSum = 0;
		typename FEMTree<Dim, Real>::template MultiThreadedEvaluator<Sigs, 0> evaluator(&tree, solution);
		
		std::vector<double> valueSums(ThreadPool::NumThreads(), 0);
		std::vector<double> weightSums(ThreadPool::NumThreads(), 0);

		auto func = [&](unsigned int thread, size_t j)
		{
			const ProjectiveData<Point<Real, Dim>, Real>& sample = (*samples)[j].sample;
			if (sample.weight > 0)
			{
				weightSums[thread] += sample.weight;
				valueSums[thread] += evaluator.values(sample.data / sample.weight, thread, (*samples)[j].node)[0] * sample.weight;
			}
		};
		
		ThreadPool::Parallel_for( 0, samples->size(), func);
		
		for (size_t t = 0; t < valueSums.size(); t++)
		{
			valueSum += valueSums[t];
			weightSum += weightSums[t];
		}
		
		isoValue = static_cast<Real>(valueSum / weightSum);

		if (!params.withColors || params.colorPullFactor == 0)
		{
			samples.reset();
		}
	}

	auto SetVertex = [] (Vertex<Real>& v, Point<Real, Dim> p, double w, PointData<Real> d)
	{
		v = Vertex<Real>(p, d, w);
	};

	ExtractMesh<Vertex<Real>, Real>(params,
									UIntPack<FEMSigs...>(),
									std::tuple<SampleData...>(),
									tree,
									solution,
									isoValue,
									samples.get(),
									sampleData.get(),
									density.get(),
									SetVertex,
									xForm.inverse(),
									out_mesh);

	return true;
}


bool PoissonReconLib::Reconstruct(	const Parameters& params,
									const ICloud<float>& inCloud,
									IMesh<float>& outMesh )
{
	if (!inCloud.hasNormals())
	{
		//we need normals
		return false;
	}

#ifdef WITH_OPENMP
	ThreadPool::Init((ThreadPool::ParallelType)(int)ThreadPool::OPEN_MP, std::thread::hardware_concurrency());
#else
	ThreadPool::Init((ThreadPool::ParallelType)(int)ThreadPool::THREAD_POOL, std::thread::hardware_concurrency());
#endif

	PointStream<float> pointStream(inCloud);

	bool success = false;

	switch (params.boundary)
	{
	case Parameters::FREE:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_FREE>::Signature> FEMSigsFree;
		success = Execute<float>(pointStream, outMesh, params, FEMSigsFree());
		break;
	case Parameters::DIRICHLET:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_DIRICHLET>::Signature> FEMSigsDirichlet;
		success = Execute<float>(pointStream, outMesh, params, FEMSigsDirichlet());
		break;
	case Parameters::NEUMANN:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_NEUMANN>::Signature> FEMSigsNeumann;
		success = Execute<float>(pointStream, outMesh, params, FEMSigsNeumann());
		break;
	default:
		assert(false);
		break;
	}

	ThreadPool::Terminate();

	return success;
}

bool PoissonReconLib::Reconstruct(	const Parameters& params,
									const ICloud<double>& inCloud,
									IMesh<double>& outMesh )
{
	if (!inCloud.hasNormals())
	{
		//we need normals
		return false;
	}

#ifdef WITH_OPENMP
	ThreadPool::Init((ThreadPool::ParallelType)(int)ThreadPool::OPEN_MP, std::thread::hardware_concurrency());
#else
	ThreadPool::Init((ThreadPool::ParallelType)(int)ThreadPool::THREAD_POOL, std::thread::hardware_concurrency());
#endif

	PointStream<double> pointStream(inCloud);

	bool success = false;

	switch (params.boundary)
	{
	case Parameters::FREE:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_FREE>::Signature> FEMSigsFree;
		success = Execute<double>(pointStream, outMesh, params, FEMSigsFree());
		break;
	case Parameters::DIRICHLET:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_DIRICHLET>::Signature> FEMSigsDirichlet;
		success = Execute<double>(pointStream, outMesh, params, FEMSigsDirichlet());
		break;
	case Parameters::NEUMANN:
		typedef IsotropicUIntPack<DIMENSION, FEMDegreeAndBType<DEFAULT_FEM_DEGREE, BOUNDARY_NEUMANN>::Signature> FEMSigsNeumann;
		success = Execute<double>(pointStream, outMesh, params, FEMSigsNeumann());
		break;
	default:
		assert(false);
		break;
	}

	ThreadPool::Terminate();

	return success;
}
