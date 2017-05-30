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

#include "CmdLineParser.h"
#include "MultiGridOctreeData.h"
#include "MyTime.h"
#include "PointStream.h"

constexpr BoundaryType BType = BOUNDARY_NEUMANN;
constexpr int Degree = 2;

void DumpOutput( const char* format , ... );
void DumpOutput2( std::vector< char* >& comments , const char* format , ... );
#include "MultiGridOctreeData.h"

#define DEFAULT_FULL_DEPTH 5

#define XSTR(x) STR(x)
#define STR(x) #x
#if DEFAULT_FULL_DEPTH
#pragma message ( "[WARNING] Setting default full depth to " XSTR(DEFAULT_FULL_DEPTH) )
#endif // DEFAULT_FULL_DEPTH

#include <stdarg.h>
char* outputFile=NULL;
int echoStdout=0;
void DumpOutput( const char* format , ... )
{
	if( outputFile )
	{
		FILE* fp = fopen( outputFile , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
}
void DumpOutput2(std::vector<std::string>& comments , const char* format , ... )
{
	if( outputFile )
	{
		FILE* fp = fopen( outputFile , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
    char str[1024];
	va_list args;
	va_start( args , format );
	vsprintf( str , format , args );
	va_end( args );
	if( str[strlen(str)-1]=='\n' )
        str[strlen(str)-1] = 0;
    comments.push_back(str);
}


cmdLineString
	In( "in" ) ,
	Out( "out" ) ,
	VoxelGrid( "voxel" ) ,
	XForm( "xForm" );

cmdLineReadable
#if defined( _WIN32 ) || defined( _WIN64 )
	Performance( "performance" ) ,
#endif // _WIN32 || _WIN64
	ShowResidual( "showResidual" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	Confidence( "confidence" ) ,
	NormalWeights( "nWeights" ) ,
	NonManifold( "nonManifold" ) ,
	ASCII( "ascii" ) ,
	Density( "density" ) ,
	LinearFit( "linearFit" ) ,
	PrimalVoxel( "primalVoxel" ) ,
#ifndef FAST_COMPILE
	Double( "double" ) ,
#endif // !FAST_COMPILE
	Verbose( "verbose" );

cmdLineInt
#ifndef FAST_COMPILE
	//Degree( "degree" , 2 ) ,
#endif // !FAST_COMPILE
	Depth( "depth" , 8 ) ,
	CGDepth( "cgDepth" , 0 ) ,
	KernelDepth( "kernelDepth" ) ,
	AdaptiveExponent( "adaptiveExp" , 1 ) ,
	Iters( "iters" , 8 ) ,
	VoxelDepth( "voxelDepth" , -1 ) ,
	FullDepth( "fullDepth" , DEFAULT_FULL_DEPTH ) ,
#ifndef FAST_COMPILE
//	BType( "bType" , BOUNDARY_NEUMANN+1 ) ,
#endif // !FAST_COMPILE
	MaxSolveDepth( "maxSolveDepth" ) ,
	Threads( "threads" , omp_get_num_procs() );

cmdLineFloat
	Color( "color" , 16.f ) ,
	SamplesPerNode( "samplesPerNode" , 1.5f ) ,
	Scale( "scale" , 1.1f ) ,
	CGSolverAccuracy( "cgAccuracy" , float(1e-3) ) ,
	LowResIterMultiplier( "iterMultiplier" , 1.f ) ,
	PointWeight( "pointWeight" , 4.f );


cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
//	&Degree ,
&Double ,
//     &BType ,
#endif // !FAST_COMPILE
	&In , &Depth , &Out , &XForm ,
	&Scale , &Verbose , &CGSolverAccuracy , &NoComments , &LowResIterMultiplier ,
	&KernelDepth , &SamplesPerNode , &Confidence , &NormalWeights , &NonManifold , &PolygonMesh , &ASCII , &ShowResidual , &VoxelDepth ,
	&PointWeight , &VoxelGrid , &Threads , &MaxSolveDepth ,
	&AdaptiveExponent ,
	&Density ,
	&FullDepth ,
	&CGDepth , &Iters ,
	&Color ,
	&LinearFit ,
	&PrimalVoxel ,
#if defined( _WIN32 ) || defined( _WIN64 )
	&Performance ,
#endif // _WIN32 || _WIN64
};


void ShowUsage(char* ex)
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>\n" , In.name );

	printf( "\t[--%s <ouput triangle mesh>]\n" , Out.name );

	printf( "\t[--%s <ouput voxel grid>]\n" , VoxelGrid.name );

#ifndef FAST_COMPILE
//	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );

//	printf( "\t[--%s <boundary type>=%d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , BoundaryNames[i] );
#endif // FAST_COMPILE

	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name , Depth.value );

	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );

	printf( "\t[--%s <minimum number of samples per node>=%f]\n" , SamplesPerNode.name, SamplesPerNode.value );

	printf( "\t[--%s <interpolation weight>=%.3e]\n" , PointWeight.name , PointWeight.value );

	printf( "\t[--%s]\n" , Confidence.name );

	printf( "\t[--%s]\n" , NormalWeights.name );

#ifndef FOR_RELEASE
	printf( "\t[--%s <adaptive weighting exponent>=%d]\n", AdaptiveExponent.name , AdaptiveExponent.value );
#endif // !FOR_RELEASE

	printf( "\t[--%s <iterations>=%d]\n" , Iters.name , Iters.value );

#ifndef FOR_RELEASE
	printf( "\t[--%s <low-resolution iteration multiplier>=%f]\n" , LowResIterMultiplier.name , LowResIterMultiplier.value );
#endif // FOR_RELEASE

	printf( "\t[--%s <conjugate-gradients depth>=%d]\n" , CGDepth.name , CGDepth.value );

#ifndef FOR_RELEASE
	printf( "\t[--%s <conjugate-gradients solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
#endif // !FOR_RELEASE

	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );

	printf( "\t[--%s <depth at which to extract the voxel grid>=<%s>]\n" , VoxelDepth.name , Depth.name );

	printf( "\t[--%s]\n" , PrimalVoxel.name );

	printf( "\t[--%s <pull factor>]\n" , Color.name );

	printf( "\t[--%s]\n" , Density.name );

	printf( "\t[--%s]\n" , LinearFit.name );

	printf( "\t[--%s]\n" , PolygonMesh.name);

#ifndef FOR_RELEASE
	printf( "\t[--%s]\n" , NonManifold.name );
#endif // !FOR_RELEASE

#ifdef _OPENMP
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
#endif // _OPENMP

	printf( "\t[--%s]\n" , Verbose.name );

#ifndef FOR_RELEASE
#if defined( _WIN32 ) || defined( _WIN64 )
	printf( "\t[--%s]\n" , Performance.name );
#endif // _WIN32 || _WIN64
#endif // !FOR_RELEASE

#ifndef FOR_RELEASE
	printf( "\t[--%s]\n" , ASCII.name );

	printf( "\t[--%s]\n" , NoComments.name );

#ifndef FAST_COMPILE
	printf( "\t[--%s]\n" , Double.name );
#endif // FAST_COMPILE
#endif // !FOR_RELEASE
}

template< class Real >
struct ColorInfo
{
	static Point3D< Real > ReadASCII( FILE* fp )
	{
		Point3D< unsigned char > c;
		if( fscanf( fp , " %c %c %c " , &c[0] , &c[1] , &c[2] )!=3 ) fprintf( stderr , "[ERROR] Failed to read color\n" ) , exit( 0 );
		return Point3D< Real >( (Real)c[0] , (Real)c[1] , (Real)c[2] );
	};
	static bool ValidPlyProperties( const bool* props ){ return ( props[0] || props[3] ) && ( props[1] || props[4] ) && ( props[2] || props[5] ); }
	const static PlyProperty PlyProperties[];
};


template<>
const PlyProperty ColorInfo< float >::PlyProperties[] =
{
	{ "r"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "g"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "b"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "red"   , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "green" , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "blue"  , PLY_UCHAR , PLY_FLOAT , int( offsetof( Point3D< float > , coords[2] ) ) , 0 , 0 , 0 , 0 }
};

template<>
const PlyProperty ColorInfo< double >::PlyProperties[] =
{
	{ "r"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "g"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "b"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "red"   , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "green" , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "blue"  , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Point3D< double > , coords[2] ) ) , 0 , 0 , 0 , 0 }
};

#if defined( _WIN32 ) || defined( _WIN64 )
double PeakMemoryUsageMB( void )
{
	HANDLE h = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	return GetProcessMemoryInfo( h , &pmc , sizeof(pmc) ) ? ( (double)pmc.PeakWorkingSetSize )/(1<<20) : 0;
}
#endif // _WIN32 || _WIN64


template< class Real >
struct OctreeProfiler
{
	Octree< Real >& tree;
	double t;

	OctreeProfiler( Octree< Real >& t ) : tree(t) { ; }
	void start( void ){ t = Time() , tree.resetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}
	void dumpOutput( const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) DumpOutput( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         DumpOutput(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) DumpOutput( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         DumpOutput(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}

	void dumpOutput2( std::vector<std::string >& comments , const char* header ) const
	{
		tree.memoryUsage();
#if defined( _WIN32 ) || defined( _WIN64 )
		if( header ) DumpOutput2( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
		else         DumpOutput2( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() , PeakMemoryUsageMB() );
#else // !_WIN32 && !_WIN64
		if( header ) DumpOutput2( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
		else         DumpOutput2( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , tree.localMemoryUsage() , tree.maxMemoryUsage() );
#endif // _WIN32 || _WIN64
	}
};

template< class Real >
XForm4x4< Real > GetPointXForm( OrientedPointStream< Real >& stream , Real scaleFactor )
{
	Point3D< Real > min , max;
	stream.boundingBox( min , max );
	Point3D< Real > center = ( max + min ) / 2;
	Real scale = std::max< Real >( max[0]-min[0] , std::max< Real >( max[1]-min[1] , max[2]-min[2] ) );
	scale *= scaleFactor;
	for( int i=0 ; i<3 ; i++ ) center[i] -= scale/2;
	XForm4x4< Real > tXForm = XForm4x4< Real >::Identity() , sXForm = XForm4x4< Real >::Identity();
	for( int i=0 ; i<3 ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(3,i) = -center[i];
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


template <typename Real>
using PointStream = OrientedPointStream<Real>;

template <typename Real>
using PointStreamWithData = OrientedPointStreamWithData<Real , Point3D<Real>>;

template <typename Real>
using XPointStream = TransformedOrientedPointStream<Real>;

template <typename Real>
using XPointStreamWithData = TransformedOrientedPointStreamWithData< Real,
    Point3D<Real>>;


template<typename Real>
struct PoissonOpts
{
    int m_threads;
    int m_voxelDepth;
    bool m_primalVoxel;
    bool m_verbose;
    bool m_confidence;
	//if( Verbose.set ) echoStdout=1;
    //bool hasColor(Color.set && Color.value > 0);
    bool m_hasColor;
    Real m_color;
    std::string m_voxelFilename;
    std::string m_xformFilename;
    std::string m_inputFilename; // Input.value
    std::string m_outputFilename; // Output.value
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
    bool m_ascii;

    PoissonOpts() : m_threads(1), m_voxelDepth(-1), m_primalVoxel(false),
        m_verbose(false), m_confidence(false), m_color(16),
        m_samplesPerNode(1.5), m_depth(8), m_cgDepth(0), m_iterations(8),
        m_showResidual(false), m_lowResIterMult(1), m_kernelDepth(0),
        m_solveDepth(0), m_solverAccuracy(1e-3), m_pointWeight(4),
        m_adaptExponent(1), m_density(false), m_linearFit(false),
        m_nonManifold(false), m_polygonMesh(false), m_ascii(false)
    {}
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
    PoissonRecon(PoissonOpts<Real> opts) : m_opts(opts),
        m_xForm(XForm4x4<Real>::Identity()),
        m_samples(new OctreeSampleVec<Real>), m_isoValue(0)
    {}
    void execute();
    void evaluate();
    void writeOutput();

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
    void writePly();
    PointStream<Real> *createPointStream();

    template<typename Vertex>
    void writeSurface(CoredFileMeshData<Vertex>& mesh);

private:
    PoissonOpts<Real> m_opts;
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


void writeVoxelValues(FILE *fp, size_t count, float *values)
{
    fwrite(values, sizeof(float), count * count * count, fp);
}

void writeVoxelValues(FILE *fp, size_t count, double *values)
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
std::cerr << "Writing voxels!\n";
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

template<typename Real>
void PoissonRecon<Real>::writePly()
{
    if (m_opts.m_density && m_opts.m_hasColor)
    {
std::cerr << "Color value vertex!\n";
        CoredFileMeshData<PlyColorAndValueVertex<float>> mesh;
        writeSurface(mesh);
    }
    else if (m_opts.m_density)
    {
std::cerr << "value vertex!\n";
        CoredFileMeshData<PlyValueVertex<float>> mesh;
        writeSurface(mesh);
    }
    else if (m_opts.m_hasColor)
    {
        CoredFileMeshData<PlyColorVertex<float>> mesh;
        writeSurface(mesh);
    }
    else
    {
        CoredFileMeshData<PlyVertex<float>> mesh;
        writeSurface(mesh);
    }
}

template<typename Real>
template<typename Vertex>
void PoissonRecon<Real>::writeSurface(CoredFileMeshData<Vertex>& mesh)
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

    std::cerr << "color data size = " << colorData.size() << "!\n";
    std::cerr << "Solution size = " << m_solution.size() << "!\n";
    std::cerr << "Linear fit/non-manifold/polygon mesh = " <<
        m_opts.m_linearFit << "/" << m_opts.m_nonManifold << "/" <<
        m_opts.m_polygonMesh << "!\n";
    std::cerr << "ISO value = " << m_isoValue << "!\n";
    m_tree.template getMCIsoSurface<Degree, BType, WEIGHT_DEGREE, DATA_DEGREE>
        (m_density, &colorData, m_solution, m_isoValue, mesh,
         !m_opts.m_linearFit, !m_opts.m_nonManifold, m_opts.m_polygonMesh);
    std::cerr << "Vertices / Polygons: " <<
        (mesh.outOfCorePointCount() + mesh.inCorePoints.size()) <<
        " / " << mesh.polygonCount();
    std::cerr << "Mesh OOC/IC/POLY = " << mesh.outOfCorePointCount() << "/" <<
        mesh.inCorePoints.size() << "/" << mesh.polygonCount() << "!\n";
    m_profiler.dumpOutput2(m_comments, std::string("#   Got ") +
        (m_opts.m_polygonMesh ? "polygons:" : "triangles:"));

    std::unique_ptr<char *> buf(new char *[m_comments.size()]);
    for (size_t i = 0; i < m_comments.size(); ++i)
        *(buf.get() + i) = (char *) m_comments[i].data();
    PlyWritePolygons((char *)m_opts.m_outputFilename.data(), &mesh,
        (m_opts.m_ascii ? PLY_ASCII : PLY_BINARY_NATIVE), buf.get(),
        m_comments.size(), m_iXForm);
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

void dumpParams(std::vector<std::string>& comments)
{
	int paramNum = sizeof(params)/sizeof(params[0]);
	char str[1024];

	for( int i=0 ; i<paramNum ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) )
                DumpOutput2( comments , "\t--%s %s\n" , params[i]->name , str );
			else
                DumpOutput2( comments , "\t--%s\n" , params[i]->name );
		}
}

template<typename Real>
PointStream<Real> *PoissonRecon<Real>::createPointStream()
{
    const std::string& f = m_opts.m_inputFilename;
    std::cerr << "Filename = " << f << "!\n";
    std::string ext = GetFileExtension(f);
    for (auto& c : ext)
        c = tolower(c);
    std::cerr << "Extension = " << ext << "!\n";

    PointStream<Real> *pointStream(nullptr);
    if (m_opts.m_hasColor)
    {
        m_sampleData = new DataSampleVec<Real>();
        if (ext == "bnpts")
            pointStream = new BinaryOrientedPointStreamWithData<Real,
                Point3D< Real > , float , Point3D< unsigned char > >(f.data());
        else if (ext == "ply")
            pointStream = new PLYOrientedPointStreamWithData<Real ,
                Point3D< Real > >(f.data(), ColorInfo< Real >::PlyProperties,
                6 , ColorInfo< Real >::ValidPlyProperties);
        else
            pointStream = new  ASCIIOrientedPointStreamWithData<Real ,
                Point3D< Real > >(f.data(), ColorInfo< Real >::ReadASCII );
    }
    else
    {
        if (ext == "bnpts")
            pointStream = new BinaryOrientedPointStream< Real, float>(f.data());
        else if (ext == "ply")
            pointStream = new PLYOrientedPointStream<Real>(f.data());
        else
            pointStream = new  ASCIIOrientedPointStream<Real>(f.data());
    }
    return pointStream;
}

template<typename Real>
int loadOctTree(Octree<Real>& tree, XForm4x4<Real>& xForm,
    PointStream<Real> *pointStream, int depth, bool confidence,
    OctreeSampleVec<Real> *samples, DataSampleVec<Real> *sampleData)
{
    int pointCount;

    if( sampleData )
    {
        XPointStreamWithData<Real> _pointStream(xForm,
            (PointStreamWithData<Real>&)*pointStream );
        pointCount = tree.template init< Point3D< Real >>(_pointStream, depth,
            confidence, *samples, sampleData);
    }
    else
    {
        XPointStream<Real> _pointStream(xForm , *pointStream);
        pointCount = tree.template init<Point3D<Real>>( _pointStream, depth,
            confidence, *samples, sampleData);
    }
    return pointCount;
}

template<typename Real>
void PoissonRecon<Real>::readData()
{
    m_profiler.start();
    PointStream<Real> *pointStream = createPointStream();
    XPointStream<Real> _pointStream(m_xForm , *pointStream);
    m_xForm = GetPointXForm( _pointStream , (Real)Scale.value ) * m_xForm;
    m_iXForm = m_xForm.inverse();
    int pointCount = loadOctTree(m_tree, m_xForm, pointStream, m_opts.m_depth,
        m_opts.m_confidence, m_samples, m_sampleData);
#pragma omp parallel for num_threads( Threads.value )
    for( int i=0 ; i<(int)m_samples->size() ; i++ )
        (*m_samples)[i].sample.data.n *= (Real)-1;
    //Ugh
    delete pointStream;

//ABELL
/**
    dump( "Input Points / Samples: " + std::to_string(pointCount) + " / " +
        std::to_string(m_samples->size()) + "\n");
**/
    m_profiler.dumpOutput2(m_comments , "# Read input into tree:");
    m_tree.resetNodeIndices();
}

template<typename Real>
void PoissonRecon<Real>::calcDensity()
{
    // Get the kernel density estimator [If discarding, compute anew.
    // Otherwise, compute once.]
    m_profiler.start();
    std::cerr << "Kernel depth/samples per node = " <<
        m_opts.m_kernelDepth << "/" << m_opts.m_samplesPerNode << "!\n";
    m_density = m_tree.template setDensityEstimator<WEIGHT_DEGREE>(*m_samples,
        m_opts.m_kernelDepth, m_opts.m_samplesPerNode);
    std::cerr << "Samples size = " << m_samples->size() << "!\n";
    std::cerr << "Density size = " << m_density->size() << "!\n";
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
        (FullDepth.value, typename Octree<Real>::template
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

    // Degree doesn't seem to be used here.
    int solveDepth = MaxSolveDepth.value;
    // This just initializes a vector of data to 0.  Blarf.
    m_constraints = m_tree.template initDenseNodeData();
std::cerr << "Constraints size = " << m_constraints.size() << "!\n";
std::cerr << "Normal info size = " << m_normalInfo.size() << "!\n";
std::cerr << "Solve depth = " << m_opts.m_solveDepth << "!\n";
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

std::cerr << "CG depth = " << m_opts.m_cgDepth << "!\n";
std::cerr << "iters = " << m_opts.m_iterations << "!\n";
std::cerr << "cgAccuracy = " << m_opts.m_solverAccuracy << "!\n";
std::cerr << "verbose = " << m_opts.m_verbose << "!\n";
std::cerr << "show residual = " << m_opts.m_showResidual << "!\n";
std::cerr << "LowResIter = " << solverInfo.lowResIterMultiplier << "!\n";

    m_solution = m_tree.template solveSystem<Degree, BType>(
        FEMSystemFunctor<Degree, BType>(0, 1, 0), m_interp, m_constraints,
        m_opts.m_solveDepth, solverInfo);
std::cerr << "Interp size = " << m_interp->iData.size() << "!\n";
std::cerr << "Interp value/grad = " << m_interp->valueWeight << "/" <<
    m_interp->gradientWeight << "!\n";
std::cerr << "Solution.size = " << m_solution.size() << "!\n";
std::cerr << "Contraints.size = " << m_constraints.size() << "!\n";
std::cerr << "Solve depth = " << m_opts.m_solveDepth << "!\n";
}

template<typename Real>
void PoissonRecon<Real>::execute()
{
//ABELL
//    Reset();
    readXForm(m_opts.m_xformFilename);
    m_comments.push_back("Running Screened Poisson Reconstruction "
        "(Version 9.01)");
//ABELL
//    m_debug.dump(m_comments.back());
    m_tree.threads = m_opts.m_threads;
	OctNode< TreeNodeData >::SetAllocator( MEMORY_ALLOCATOR_BLOCK_SIZE );
    readData();

	Real pointWeightSum;
    calcDensity();
    calcNormalData();
    trim();
    addFEMConstraints();
    addInterpolationConstraints();
//ABELL
    std::cerr << "Leaf Nodes / Active Nodes / GhostNodes: " +
        std::to_string(m_tree.leaves()) + " / " +
        std::to_string(m_tree.nodes()) + " / " +
        std::to_string(m_tree.ghostNodes()) + "\n";
/**
    DumpOutput("Leaf Nodes / Active Nodes / GhostNodes: " +
        std::to_string(m_tree.leaves()) + " / " +
        std::to_string(m_tree.nodes()) + " / " +
        std::to_string(m_tree.ghostNodes()) + "\n");
    DumpOutput("Memory usage:
**/
    solve();
}

template<typename Real>
void PoissonRecon<Real>::evaluate()
{
    m_profiler.start();
    double valueSum = 0;
    double weightSum = 0;
    typename Octree<Real>::template MultiThreadedEvaluator<Degree, BType>
        evaluator(&m_tree, m_solution, m_opts.m_threads);
    std::cerr << "Samples size = " << m_samples->size() << "!\n";
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
    std::cerr << "Iso-Value: " << m_isoValue << "!\n";
}

template<typename Real>
void PoissonRecon<Real>::writeOutput()
{
    writeVoxels();
    writePly();
}

int main(int argc, char *argv[])
{
    using Real = float;

	cmdLineParse(argc-1, &argv[1] ,sizeof(params)/sizeof(params[0]),
        params , 1);

    PoissonOpts<Real> opts;
//Right now we only care about some options.

    if (!In.set)
    {
        std::cerr << "Must supply input file!\n";
//        ShowUsage(argv[0]);
        return 0;
    }
    if (!Out.set)
    {
        std::cerr << "Must supply output file!\n";
        return 0;
    }
    opts.m_inputFilename = In.value;
    opts.m_outputFilename = Out.value;
    if (Depth.set)
        opts.m_depth = Depth.value;
    if (Color.set)
    {
        opts.m_color = Color.value;
        opts.m_hasColor = (Color.value > 0);
    }
    if (Density.set)
        opts.m_density = true;
    if (PointWeight.set)
        opts.m_pointWeight = PointWeight.value;

//ABELL
//    m_startTime = Time());
	if( !MaxSolveDepth.set )
        MaxSolveDepth.value = Depth.value;

    opts.m_solveDepth = MaxSolveDepth.value;
    opts.m_kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
	if (opts.m_kernelDepth > opts.m_depth)
	{
		std::cerr << "[WARNING] kernelDepth can't be greater than depth: " <<
            opts.m_kernelDepth << " <= " << opts.m_depth << "." << std::endl;
		opts.m_kernelDepth = opts.m_depth;
	}

    std::cerr << "Opts filename = " << opts.m_inputFilename << "!\n";
    PoissonRecon<Real> recon(opts);

    recon.execute();
    recon.evaluate();
    recon.writeOutput();

    return 0;
}

/**
template< class Real , int Degree , BoundaryType BType , class Vertex >
int _Execute( int argc , char* argv[] )
{
	typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;

	Reset< Real >();


		// Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
	    SparseNodeData< Point3D< Real > >* normalInfo = NULL;
		{
			profiler.start();
			normalInfo = new SparseNodeData< Point3D< Real > >();
			*normalInfo = tree.template setNormalField< NORMAL_DEGREE >( *samples , *density , pointWeightSum , BType==BOUNDARY_NEUMANN );
			profiler.dumpOutput2( comments , "#     Got normal field:" );
		}

		if( !Density.set ) delete density , density = NULL;

		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			std::vector< int > indexMap;

			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
			tree.template inalizeForBroodedMultigrid< MAX_DEGREE , Degree , BType >( FullDepth.value , typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >( *normalInfo ) , &indexMap );

			if( normalInfo ) normalInfo->remapIndices( indexMap );
			if( density ) density->remapIndices( indexMap );
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}

		DenseNodeData< Real > constraints;
		int solveDepth = MaxSolveDepth.value;
		// Add the FEM constraints
		{
			profiler.start();
			constraints = tree.template initDenseNodeData( );
			tree.template addFEMConstraints< Degree , BType , NORMAL_DEGREE , BType >( FEMVFConstraintFunctor< NORMAL_DEGREE , BType , Degree , BType >( 1. , 0. ) , *normalInfo , constraints , solveDepth );
			profiler.dumpOutput2( comments , "#  Set FEM constraints:" );
		}

		// Free up the normal info [If we don't need it for subseequent iterations.]
		delete normalInfo , normalInfo = NULL;

		InterpolationInfo* iInfo = NULL;
		// Add the interpolation constraints
		if( PointWeight.value>0 )
		{
			profiler.start();
	        Real targetValue = (Real)0.5;
			iInfo = new InterpolationInfo( tree , *samples , targetValue , AdaptiveExponent.value , (Real)PointWeight.value * pointWeightSum , (Real)0 );
			tree.template addInterpolationConstraints< Degree , BType >( *iInfo , constraints , solveDepth );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
		}

		DumpOutput( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );
		DumpOutput( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename Octree< Real >::SolverInfo solverInfo;
			solverInfo.cgDepth = CGDepth.value , solverInfo.iters = Iters.value , solverInfo.cgAccuracy = CGSolverAccuracy.value , solverInfo.verbose = Verbose.set , solverInfo.showResidual = ShowResidual.set , solverInfo.lowResIterMultiplier = std::max< double >( 1. , LowResIterMultiplier.value );
			solution = tree.template solveSystem< Degree , BType >( FEMSystemFunctor< Degree , BType >( 0 , 1. , 0 ) , iInfo , constraints , solveDepth , solverInfo );
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	CoredFileMeshData< Vertex > mesh;

	{
		profiler.start();
		double valueSum = 0 , weightSum = 0;
		typename Octree< Real >::template MultiThreadedEvaluator< Degree , BType > evaluator( &tree , solution , Threads.value );
#pragma omp parallel for num_threads( Threads.value ) reduction( + : valueSum , weightSum )
		for( int j=0 ; j<samples->size() ; j++ )
		{
			ProjectiveData< OrientedPoint3D< Real > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if( w>0 ) weightSum += w , valueSum += evaluator.value( sample.data.p / sample.weight , omp_get_thread_num() , (*samples)[j].node ) * w;
		}
		isoValue = (Real)( valueSum / weightSum );
		if( !( Color.set && Color.value>0 ) && samples ) delete samples , samples = NULL;
		profiler.dumpOutput( "Got average:" );
		DumpOutput( "Iso-Value: %e\n" , isoValue );
	}

	if( VoxelGrid.set )
	{
		profiler.start();
        writeVoxels<Real, Degree, BType>(tree, solution, VoxelGrid.value,
            isoValue, VoxelDepth.value, PrimalVoxel.set);
		profiler.dumpOutput( "Got voxel grid:" );
    }

	if( Out.set )
	{
		profiler.start();
		SparseNodeData< ProjectiveData< Point3D< Real > , Real > >* colorData = NULL;
		if( sampleData )
		{
			colorData = new SparseNodeData< ProjectiveData< Point3D< Real > , Real > >();
			*colorData = tree.template setDataField< DATA_DEGREE , false >( *samples , *sampleData , (DensityEstimator*)NULL );
			delete sampleData , sampleData = NULL;
			for( const OctNode< TreeNodeData >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
			{
				ProjectiveData< Point3D< Real > , Real >* clr = (*colorData)( n );
				if( clr ) (*clr) *= (Real)pow( Color.value , tree.depth( n ) );
			}
		}
		tree.template getMCIsoSurface< Degree , BType , WEIGHT_DEGREE , DATA_DEGREE >( density , colorData , solution , isoValue , mesh , !LinearFit.set , !NonManifold.set , PolygonMesh.set );
		DumpOutput( "Vertices / Polygons: %d / %d\n" , mesh.outOfCorePointCount()+mesh.inCorePoints.size() , mesh.polygonCount() );
		if( PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
		else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

		if( colorData ) delete colorData , colorData = NULL;

        int plyType(ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE);
        char **commentData(NoComments.set ? NULL : &comments[0]);
        int numComments(NoComments.set ? 0 : (int)comments.size());

        PlyWritePolygons( Out.value , &mesh , plyType, commentData, numComments, iXForm );
	}
	if( density ) delete density , density = NULL;
	DumpOutput2( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , tree.maxMemoryUsage() );

	return 1;
}

#if defined( _WIN32 ) || defined( _WIN64 )
inline double to_seconds( const FILETIME& ft )
{
	const double low_to_sec=100e-9; // 100 nanoseconds
	const double high_to_sec=low_to_sec*4294967296.0;
	return ft.dwLowDateTime*low_to_sec+ft.dwHighDateTime*high_to_sec;
}
#endif // _WIN32 || _WIN64
**/

/**
#ifndef FAST_COMPILE
template< class Real , class Vertex >
int Execute( int argc , char* argv[] )
{
	switch( BType.value )
	{
		case BOUNDARY_FREE+1:
		{
			switch( Degree.value )
			{
			case 1: return _Execute< Real , 1 , BOUNDARY_FREE , Vertex >( argc , argv );
			case 2: return _Execute< Real , 2 , BOUNDARY_FREE , Vertex >( argc , argv );
			case 3: return _Execute< Real , 3 , BOUNDARY_FREE , Vertex >( argc , argv );
			case 4: return _Execute< Real , 4 , BOUNDARY_FREE , Vertex >( argc , argv );
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) ; return EXIT_FAILURE;
			}
		}
		case BOUNDARY_NEUMANN+1:
		{
			switch( Degree.value )
			{
			case 1: return _Execute< Real , 1 , BOUNDARY_NEUMANN , Vertex >( argc , argv );
			case 2: return _Execute< Real , 2 , BOUNDARY_NEUMANN , Vertex >( argc , argv );
			case 3: return _Execute< Real , 3 , BOUNDARY_NEUMANN , Vertex >( argc , argv );
			case 4: return _Execute< Real , 4 , BOUNDARY_NEUMANN , Vertex >( argc , argv );
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) ; return EXIT_FAILURE;
			}
		}
		case BOUNDARY_DIRICHLET+1:
		{
			switch( Degree.value )
			{
			case 1: return _Execute< Real , 1 , BOUNDARY_DIRICHLET , Vertex >( argc , argv );
			case 2: return _Execute< Real , 2 , BOUNDARY_DIRICHLET , Vertex >( argc , argv );
			case 3: return _Execute< Real , 3 , BOUNDARY_DIRICHLET , Vertex >( argc , argv );
			case 4: return _Execute< Real , 4 , BOUNDARY_DIRICHLET , Vertex >( argc , argv );
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) ; return EXIT_FAILURE;
			}
		}
		default: fprintf( stderr , "[ERROR] Not a valid boundary type: %d\n" , BType.value ) ; return EXIT_FAILURE;
	}
}
#endif // !FAST_COMPILE

int main( int argc , char* argv[] )
{
#ifdef ARRAY_DEBUG
	fprintf( stderr , "[WARNING] Running in array debugging mode\n" );
#endif // ARRAY_DEBUG
#if defined( WIN32 ) && defined( MAX_MEMORY_GB )
	if( MAX_MEMORY_GB>0 )
	{
		SIZE_T peakMemory = 1;
		peakMemory <<= 30;
		peakMemory *= MAX_MEMORY_GB;
		printf( "Limiting memory usage to %.2f GB\n" , float( peakMemory>>30 ) );
		HANDLE h = CreateJobObject( NULL , NULL );
		AssignProcessToJobObject( h , GetCurrentProcess() );

		JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
		jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_MEMORY;
		jeli.JobMemoryLimit = peakMemory;
		if( !SetInformationJobObject( h , JobObjectExtendedLimitInformation , &jeli , sizeof( jeli ) ) )
			fprintf( stderr , "Failed to set memory limit\n" );
	}
#endif // defined( WIN32 ) && defined( MAX_MEMORY_GB )
	double t = Time();

	cmdLineParse( argc-1 , &argv[1] , sizeof(params)/sizeof(cmdLineReadable*) , params , 1 );
#ifdef FAST_COMPILE
	static const int Degree = 2;
	static const BoundaryType BType = BOUNDARY_NEUMANN;
	fprintf( stderr , "[WARNING] Compiling for degree-%d, boundary-%s, single-precision _only_\n" , Degree , BoundaryNames[ BType ] );
	if( Density.set )
		if( Color.set && Color.value>0 ) return _Execute< float , Degree , BType , PlyColorAndValueVertex< float > >( argc , argv );
		else                             return _Execute< float , Degree , BType , PlyValueVertex< float > >( argc , argv );
	else
		if( Color.set && Color.value>0 ) return _Execute< float , Degree , BType , PlyColorVertex< float > >( argc , argv );
		else                             return _Execute< float , Degree , BType , PlyVertex< float > >( argc , argv );
#else // !FAST_COMPILE
	{
		if( Density.set )
			if( Color.set && Color.value>0 )
				if( Double.set ) Execute< double , PlyColorAndValueVertex< float > >( argc , argv );
				else             Execute< float  , PlyColorAndValueVertex< float > >( argc , argv );
			else
				if( Double.set ) Execute< double , PlyValueVertex< float > >( argc , argv );
				else             Execute< float  , PlyValueVertex< float > >( argc , argv );
		else
			if( Color.set && Color.value>0 )
				if( Double.set ) Execute< double , PlyColorVertex< float > >( argc , argv );
				else             Execute< float  , PlyColorVertex< float > >( argc , argv );
			else
				if( Double.set ) Execute< double , PlyVertex< float > >( argc , argv );
				else             Execute< float  , PlyVertex< float > >( argc , argv );
	}
#endif // FAST_COMPILE
#if defined( _WIN32 ) || defined( _WIN64 )
	if( Performance.set )
	{
		HANDLE cur_thread=GetCurrentThread();
		FILETIME tcreat, texit, tkernel, tuser;
		if( GetThreadTimes( cur_thread , &tcreat , &texit , &tkernel , &tuser ) )
			printf( "Time (Wall/User/Kernel): %.2f / %.2f / %.2f\n" , Time()-t , to_seconds( tuser ) , to_seconds( tkernel ) );
		else printf( "Time: %.2f\n" , Time()-t );
		HANDLE h = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		if( GetProcessMemoryInfo( h , &pmc , sizeof(pmc) ) ) printf( "Peak Memory (MB): %d\n" , (int)( pmc.PeakWorkingSetSize>>20 ) );
	}
#endif // _WIN32 || _WIN64
	return EXIT_SUCCESS;
}
**/
