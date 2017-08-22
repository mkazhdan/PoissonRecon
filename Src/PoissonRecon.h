/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#undef FAST_COMPILE
#undef ARRAY_DEBUG
#define BRUNO_LEVY_FIX
#define FOR_RELEASE

//ABELL - Necessary?
/**
#if defined( _WIN32 ) || defined( _WIN64 )
#include <Windows.h>
#include <Psapi.h>
#endif // _WIN32 || _WIN64
**/

#include <iostream>
#include <vector>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include <memory>

#include "MultiGridOctreeData.h"
#include "Mesh.h"
#include "point_source/TransformedPointSource.h"

constexpr BoundaryType BType = BOUNDARY_NEUMANN;
constexpr int Degree = 2;

void DumpOutput( const char* format , ... );
void DumpOutput2( std::vector< char* >& comments , const char* format , ... );
#include "MultiGridOctreeData.h"

#define DEFAULT_FULL_DEPTH 5

#define XSTR(x) STR(x)
#define STR(x) #x
/**
#if DEFAULT_FULL_DEPTH
#pragma message ( "[WARNING] Setting default full depth to " XSTR(DEFAULT_FULL_DEPTH) )
#endif // DEFAULT_FULL_DEPTH
**/

#include <stdarg.h>

template<typename Real>
XForm4x4<Real> GetPointXForm(PointSource& source, Real scaleFactor)
{
	Point3D<double> min , max;
	source.boundingBox(min, max);
	Point3D<double> center = ( max + min ) / 2;
	double scale = std::max( max[0]-min[0],
        std::max(max[1]-min[1], max[2]-min[2]));
	scale *= scaleFactor;
	for( int i=0 ; i<3 ; i++ )
        center[i] -= scale/2;

	XForm4x4<Real> tXForm = XForm4x4<Real>::Identity();
    XForm4x4<Real> sXForm = XForm4x4<Real>::Identity();
	for( int i=0 ; i<3 ; i++ )
    {
        sXForm(i,i) = (Real)(1./scale );
        tXForm(3,i) = -(Real)center[i];
    }
	return sXForm * tXForm;
}

template<typename Real>
using ProjData = ProjectiveData<Point3D<Real>, Real>;

template<typename Real>
using DataSample = ProjData<Real>;

template<typename Real>
using DataSampleVec = std::vector<DataSample<Real>>;

template<typename Real>
using OctreeSample = typename Octree<Real>::PointSample;

template <typename Real>
using OctreeSampleVec = std::vector<OctreeSample<Real>>;

template<typename Real>
using DensityEstimator = typename Octree<Real>::DensityEstimator;

template<typename Real>
using InterpolationInfo =
    typename Octree<Real>::template InterpolationInfo<false>;

template<typename Real>
struct PoissonOpts
{
    int m_threads;
    int m_voxelDepth;
    bool m_primalVoxel;
    bool m_verbose;
    bool m_confidence;
    bool m_hasColor;
    Real m_color;
    std::string m_voxelFilename;
    std::string m_xformFilename;
    Real m_samplesPerNode;
    int m_depth;
    int m_cgDepth;
    int m_iterations;
    bool m_showResidual;
    Real m_lowResIterMult;
    int m_kernelDepth;
    int m_solveDepth; // MaxSolveDepth );
    Real m_solverAccuracy;
    Real m_pointWeight;
    int m_adaptExponent; // 1
    bool m_density;
    bool m_linearFit;
    bool m_nonManifold;
    bool m_polygonMesh;
    int m_fullDepth;
    Real m_scale;

    PoissonOpts() : m_threads(omp_get_num_procs()), m_voxelDepth(-1),
        m_primalVoxel(false), m_verbose(false), m_confidence(false),
        m_color(16), m_samplesPerNode(1.5F), m_depth(8), m_cgDepth(0),
        m_iterations(8), m_showResidual(false), m_lowResIterMult(1),
        m_kernelDepth(0), m_solveDepth(0), m_solverAccuracy(1e-3F),
        m_pointWeight(4), m_adaptExponent(1), m_density(false),
        m_linearFit(false), m_nonManifold(false), m_polygonMesh(false),
        m_fullDepth(5), m_scale(1.1F)
    {}

    void dump()
    {
        std::cerr << "Threads = " << m_threads << "!\n";
        std::cerr << "Voxel depth = " << m_voxelDepth << "!\n";
        std::cerr << "Primal voxel = " << m_primalVoxel << "!\n";
        std::cerr << "Verbose = " << m_verbose << "!\n";
        std::cerr << "Confidence = " << m_confidence << "!\n";
        std::cerr << "Has color = " << m_hasColor << "!\n";
        std::cerr << "Color = " << m_color << "!\n";
        std::cerr << "Voxel filename = " << m_voxelFilename << "!\n";
        std::cerr << "Xform filename = " << m_xformFilename << "!\n";
        std::cerr << "Samples per node = " << m_samplesPerNode << "!\n";
        std::cerr << "Depth = " << m_depth << "!\n";
        std::cerr << "CG depth = " << m_cgDepth << "!\n";
        std::cerr << "Iterations = " << m_iterations << "!\n";
        std::cerr << "Show residual = " << m_showResidual << "!\n";
        std::cerr << "Low res multiplier = " << m_lowResIterMult << "!\n";
        std::cerr << "Kernel depth = " << m_kernelDepth << "!\n";
        std::cerr << "Solve depth = " << m_solveDepth << "!\n";
        std::cerr << "Solver accuracy = " << m_solverAccuracy << "!\n";
        std::cerr << "Point weight = " << m_pointWeight << "!\n";
        std::cerr << "Adapt exponent = " << m_adaptExponent << "!\n";
        std::cerr << "Density = " << m_density << "!\n";
        std::cerr << "Linear fit = " << m_linearFit << "!\n";
        std::cerr << "Non-manifold = " << m_nonManifold << "!\n";
        std::cerr << "Poly mesh = " << m_polygonMesh << "!\n";
        std::cerr << "Full depth = " << m_fullDepth << "!\n";
        std::cerr << "Scale = " << m_scale << "!\n";
    }
};

class Profiler
{
public:
    void start()
    {}
    void dumpOutput(const std::string& s)
        { std::cerr << s << std::endl; }
    void dumpOutput2(std::vector<std::string>& comments, const std::string& s)
    {
        comments.push_back(s);
        dumpOutput(s);
    }
};

template<class Real>
class PoissonRecon
{
public:
    PoissonRecon(PoissonOpts<Real> opts, PointSource& pointSource) :
        m_opts(opts), m_pointSource(pointSource),
        m_xForm(XForm4x4<Real>::Identity()),
        m_samples(new OctreeSampleVec<Real>), m_isoValue(0)
    {}
    void execute();
    void evaluate();
    void extractMesh(Kazhdan::Mesh& mesh);
    std::vector<std::string> comments() const
        { return m_comments; }
    XForm4x4<Real> inverseTransform() const
        { return m_iXForm; }

private:
    void readData();
    void calcDensity();
    void calcNormalData();
    void readXForm(const std::string& filename);
    void trim();
    void addFEMConstraints();
    void addInterpolationConstraints();
    void solve();
    void writeVoxels();

    template<typename Vertex>
    void writeSurface(Kazhdan::Mesh& mesh);

private:
    PoissonOpts<Real> m_opts;
    PointSource& m_pointSource;
    XForm4x4<Real> m_xForm;
    XForm4x4<Real> m_iXForm; // Inverse transform.
    Octree<Real> m_tree;
    Profiler m_profiler;  // OctreeProfiler<Real>(tree);
    DataSampleVec<Real> *m_sampleData;
    OctreeSampleVec<Real> *m_samples;
    Real m_isoValue;
    double m_startTime;
    DensityEstimator<Real> *m_density;
    SparseNodeData<Point3D<Real> > m_normalInfo;
    Real m_pointWeightSum;
    DenseNodeData<Real> m_constraints;
    DenseNodeData<Real> m_solution;
    InterpolationInfo<Real> *m_interp;
    std::vector<std::string> m_comments;
};


inline void writeVoxelValues(FILE *fp, size_t count, float *values)
{
    fwrite(values, sizeof(float), count * count * count, fp);
}

inline void writeVoxelValues(FILE *fp, size_t count, double *values)
{
    while (count--)
    {
        float f = (float) *values++;
        fwrite(&f, sizeof(float), 1, fp);
    }
}


template <typename Real>
void PoissonRecon<Real>::writeVoxels()
{
    if (m_opts.m_voxelFilename.empty())
        return;
    FILE* fp = fopen(m_opts.m_voxelFilename.data(), "wb" );
    if ( !fp )
        std::cerr << "Failed to open voxel file for writing: " <<
            m_opts.m_voxelFilename << std::endl;
    else
    {
        int res = 0;
        Pointer( Real ) values =
            m_tree.template voxelEvaluate<Real, Degree, BType>(m_solution, res,
            m_isoValue, m_opts.m_voxelDepth, m_opts.m_primalVoxel);
        fwrite( &res , sizeof(int) , 1 , fp );
        writeVoxelValues(fp, res * res * res, values);
        fclose( fp );
        DeletePointer( values );
    }
}

template <typename Real>
void PoissonRecon<Real>::extractMesh(Kazhdan::Mesh& mesh)
{
    if (m_opts.m_density && m_opts.m_hasColor)
    {
        writeSurface<PlyColorAndValueVertex<float>>(mesh);
    }
    else if (m_opts.m_density)
    {
        writeSurface<PlyValueVertex<float>>(mesh);
    }
    else if (m_opts.m_hasColor)
    {
        writeSurface<PlyColorVertex<float>>(mesh);
    }
    else
    {
        writeSurface<PlyVertex<float>>(mesh);
    }
}

template <typename Real>
template<typename Vertex>
void PoissonRecon<Real>::writeSurface(Kazhdan::Mesh& mesh)
{
    using ColorData = SparseNodeData<ProjectiveData<Point3D<Real>, Real> >;

    m_profiler.start();
    ColorData colorData = m_tree.template setDataField<DATA_DEGREE, false,
        WEIGHT_DEGREE>(
        *m_samples, *m_sampleData, (DensityEstimator<Real> *)nullptr);
    for (const OctNode<TreeNodeData>* n = m_tree.tree().nextNode(); n;
        n = m_tree.tree().nextNode(n))
    {
        ProjectiveData<Point3D<Real>, Real> *color = colorData(n);
        if (color)
            *color *= pow(m_opts.m_color, m_tree.depth(n));
    }

    m_tree.template getMCIsoSurface<Degree, BType, WEIGHT_DEGREE, DATA_DEGREE, Vertex>
        (m_density, &colorData, m_solution, m_isoValue, mesh,
         !m_opts.m_linearFit, !m_opts.m_nonManifold, m_opts.m_polygonMesh);
}

template<typename Real>
void PoissonRecon<Real>::readXForm(const std::string& filename)
{
    if (filename.empty())
        return;

    FILE* fp = fopen(filename.data(), "r" );
    if( !fp )
    {
        fprintf( stderr , "[WARNING] Could not read x-form from: %s\n" ,
            filename.data());
        return;
    }
    for( int i=0 ; i<4 ; i++ )
        for( int j=0 ; j<4 ; j++ )
        {
            float f;
            if ( fscanf( fp , " %f " , &f )!= 1 )
            {
                fprintf( stderr , "[ERROR] Execute: Failed to read xform\n" );
                exit( 0 );
            }
            m_xForm(i,j) = (Real)f;
        }
    fclose( fp );
}

template<typename Real>
int loadOctTree(Octree<Real>& tree, XForm4x4<Real>& xForm, PointSource& source,
    int depth, bool confidence, OctreeSampleVec<Real> *samples,
    DataSampleVec<Real> *sampleData)
{
    try
    {
        ColorPointSource& colorSource =
            dynamic_cast<ColorPointSource &>(source);

        ColorTransformedPointSource xsource(xForm, colorSource);
        return tree.template init< Point3D< Real >>(xsource, depth,
            confidence, *samples, sampleData);
    }
    catch (std::bad_cast)
    {}

    TransformedPointSource xsource(xForm, source);
    return tree.template init<Point3D<Real>>(xsource, depth,
        confidence, *samples, sampleData);
}

template<typename Real>
void PoissonRecon<Real>::readData()
{
    m_profiler.start();

    m_sampleData = new DataSampleVec<Real>();
    TransformedPointSource xPointSource(m_xForm, m_pointSource);
    m_xForm = GetPointXForm(xPointSource, m_opts.m_scale) * m_xForm;
    m_iXForm = m_xForm.inverse();
    int pointCount = loadOctTree(m_tree, m_xForm, m_pointSource, m_opts.m_depth,
        m_opts.m_confidence, m_samples, m_sampleData);
#pragma omp parallel for num_threads( m_opts.m_threads)
    for( int i=0 ; i<(int)m_samples->size() ; i++ )
        (*m_samples)[i].sample.data.n *= (Real)-1;

    m_profiler.dumpOutput2(m_comments , "# Read input into tree:");
    m_tree.resetNodeIndices();
}

template<typename Real>
void PoissonRecon<Real>::calcDensity()
{
    // Get the kernel density estimator [If discarding, compute anew.
    // Otherwise, compute once.]
    m_profiler.start();
    m_density = m_tree.template setDensityEstimator<WEIGHT_DEGREE>(*m_samples,
        m_opts.m_kernelDepth, m_opts.m_samplesPerNode);
	m_profiler.dumpOutput2(m_comments , "#   Got kernel density:");
}

template<typename Real>
void PoissonRecon<Real>::calcNormalData()
{
    // Transform the Hermite samples into a vector field.
    m_profiler.start();
    m_normalInfo = m_tree.template setNormalField< NORMAL_DEGREE,
        WEIGHT_DEGREE >(*m_samples, *m_density , m_pointWeightSum ,
        BType==BOUNDARY_NEUMANN);
    m_profiler.dumpOutput2(m_comments , "#     Got normal field:" );
}

template<typename Real>
void PoissonRecon<Real>::trim()
{
    // Trim the tree and prepare for multigrid
    m_profiler.start();
    std::vector<int> indexMap;

    constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ?
        NORMAL_DEGREE : Degree;
    m_tree.template inalizeForBroodedMultigrid<MAX_DEGREE, Degree, BType>
        (m_opts.m_fullDepth, typename Octree<Real>::template
        HasNormalDataFunctor<NORMAL_DEGREE>( m_normalInfo ) , &indexMap);

    m_normalInfo.remapIndices(indexMap);
    if (m_density)
        m_density->remapIndices( indexMap );
    m_profiler.dumpOutput2(m_comments , "#       Finalized tree:" );
}

template<typename Real>
void PoissonRecon<Real>::addFEMConstraints()
{
    m_profiler.start();

    // This just initializes a vector of data to 0.  Blarf.
    m_constraints = m_tree.initDenseNodeData();
//ABELL - This is really slow!
    m_tree.template addFEMConstraints<Degree, BType, NORMAL_DEGREE, BType>(
        FEMVFConstraintFunctor<NORMAL_DEGREE, BType, Degree, BType >( 1., 0.),
            m_normalInfo, m_constraints , m_opts.m_solveDepth);
    m_profiler.dumpOutput2(m_comments , "#  Set FEM constraints:");
}

template<typename Real>
void PoissonRecon<Real>::addInterpolationConstraints()
{
    m_profiler.start();
    if (m_opts.m_pointWeight > 0)
    {
        m_interp = new InterpolationInfo<Real>(m_tree, *m_samples, .5,
            m_opts.m_adaptExponent, m_pointWeightSum * m_opts.m_pointWeight, 0);
        m_tree.template addInterpolationConstraints<Degree, BType>(*m_interp,
            m_constraints, m_opts.m_solveDepth);
    }
    m_profiler.dumpOutput2(m_comments, "#Set point constraints:");
}

template<typename Real>
void PoissonRecon<Real>::solve()
{
    m_profiler.start();
    typename Octree<Real>::SolverInfo solverInfo;
    solverInfo.cgDepth = m_opts.m_cgDepth;
    solverInfo.iters = m_opts.m_iterations;
    solverInfo.cgAccuracy = m_opts.m_solverAccuracy;
    solverInfo.verbose = m_opts.m_verbose;
    solverInfo.showResidual = m_opts.m_showResidual;
    solverInfo.lowResIterMultiplier =
        std::max((Real)1.0, m_opts.m_lowResIterMult);

    m_solution = m_tree.template solveSystem<Degree, BType>(
        FEMSystemFunctor<Degree, BType>(0, 1, 0), m_interp, m_constraints,
        m_opts.m_solveDepth, solverInfo);
}

template<typename Real>
void PoissonRecon<Real>::execute()
{
    readXForm(m_opts.m_xformFilename);
    m_comments.push_back("Running Screened Poisson Reconstruction "
        "(Version 9.01)");
//    m_debug.dump(m_comments.back());
    m_tree.threads = m_opts.m_threads;
	OctNode< TreeNodeData >::SetAllocator( MEMORY_ALLOCATOR_BLOCK_SIZE );
    readData();

    calcDensity();
    calcNormalData();
    trim();
    addFEMConstraints();
    addInterpolationConstraints();
    solve();
}

template<typename Real>
void PoissonRecon<Real>::evaluate()
{
    m_profiler.start();
    Real valueSum = 0;
    Real weightSum = 0;
    typename Octree<Real>::template MultiThreadedEvaluator<Degree, BType>
        evaluator(&m_tree, m_solution, m_opts.m_threads);
    for (auto & s : *m_samples)
    {
        ProjectiveData<OrientedPoint3D<Real>, Real>& sample = s.sample;
        Real w = sample.weight;
        if (w > 0)
        {
            weightSum += w;
            valueSum += evaluator.value(sample.data.p / sample.weight,
                omp_get_thread_num(), s.node)  * w;
        }
    }
    m_isoValue = valueSum / weightSum;

    m_profiler.dumpOutput("Got average:");
}

