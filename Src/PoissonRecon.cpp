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
/**/
#if defined( _WIN32 ) || defined( _WIN64 )
#include <Windows.h>
#include <Psapi.h>
#endif // _WIN32 || _WIN64
/**/

#include "CmdLineParser.h"
#include "PoissonRecon.h"
#include "point_source/AsciiPointSource.h"
#include "point_source/BinaryPointSource.h"
#include "point_source/PlyPointSource.h"
#include "MemMesh.h"

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

using MeshPtr = std::unique_ptr<Kazhdan::Mesh>;

template <typename Vertex, typename Real>
void writeOutput(PoissonRecon<Real>& recon, Kazhdan::Mesh& mesh)
{
    std::cerr << "Vertices / Polygons: " << mesh.pointCount() << " / " <<
        mesh.polygonCount() << std::endl;

    std::vector<std::string> comments(recon.comments());
    std::unique_ptr<char *> buf(new char *[comments.size()]);
    for (size_t i = 0; i < comments.size(); ++i)
        *(buf.get() + i) = (char *)comments[i].data();
    PlyWritePolygons<Vertex, Real>(Out.value, mesh,
        (ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE), buf.get(), comments.size(),
        recon.inverseTransform());
}

template <typename Real>
void extractMesh(PoissonRecon<Real>& recon, bool color, bool density)
{
    MeshPtr mesh;
    if (color && density)
    {
        mesh.reset(new Kazhdan::CompleteMemMesh);
        recon.extractMesh(*mesh);
        writeOutput<PlyColorAndValueVertex<float>, Real>(recon, *mesh);
    }
    else if (color)
    {
        mesh.reset(new Kazhdan::ColorMemMesh);
        recon.extractMesh(*mesh);
        writeOutput<PlyColorVertex<float>, Real>(recon, *mesh);
    }
    else if (density)
    {
        mesh.reset(new Kazhdan::DensityMemMesh);
        recon.extractMesh(*mesh);
        writeOutput<PlyValueVertex<float>, Real>(recon, *mesh);
    }
    else
    {
        mesh.reset(new Kazhdan::MemMesh);
        recon.extractMesh(*mesh);
        writeOutput<PlyVertex<float>, Real>(recon, *mesh);
    }
}

PointSourcePtr createPointSource(const std::string& filename, bool color)
{
    std::string ext = GetFileExtension(filename);
    for (auto& c : ext)
        c = tolower(c);

    PointSource *pointSource;
    if (color)
    {
        if (ext == "bnpts")
            pointSource = new ColorBinaryPointSource(filename);
        else if (ext == "ply")
            pointSource = new ColorPLYPointSource(filename);
        else
            pointSource = new ColorASCIIPointSource(filename);
    }
    else
    {
        if (ext == "bnpts")
            pointSource = new BinaryPointSource(filename);
        else if (ext == "ply")
            pointSource = new PLYPointSource(filename);
        else
            pointSource = new ASCIIPointSource(filename);
    }
    return PointSourcePtr(pointSource);
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
        ShowUsage(argv[0]);
        return 0;
    }
    if (!Out.set)
    {
        std::cerr << "Must supply output file!\n";
        return 0;
    }
    if (Depth.set)
        opts.m_depth = Depth.value;
    if (Color.set)
    {
        opts.m_color = Color.value;
        opts.m_hasColor = (Color.value > 0);
    }
    opts.m_density = Density.set;
    if (PointWeight.set)
        opts.m_pointWeight = PointWeight.value;
    opts.m_primalVoxel = PrimalVoxel.set;
    opts.m_verbose = Verbose.set;
    opts.m_confidence = Confidence.set;
    opts.m_showResidual = ShowResidual.set;
    opts.m_linearFit = LinearFit.set;
    opts.m_nonManifold = NonManifold.set;
    opts.m_polygonMesh = PolygonMesh.set;

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

    PointSourcePtr pointSource(createPointSource(In.value, opts.m_hasColor));

    PoissonRecon<Real> recon(opts, *pointSource);
    recon.execute();
    recon.evaluate();

    extractMesh(recon, opts.m_hasColor, opts.m_density);

    return 0;
}

