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

#include "PreProcessor.h"

#undef USE_DOUBLE								// If enabled, double-precesion is used
#define DATA_DEGREE 0							// The order of the B-Spline used to splat in data for color interpolation
												// This can be changed to zero if more interpolatory performance is desired.
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2							// The order of the B-Spline used to splat int the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 2					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN	// The default finite-element boundary type
#define DEFAULT_DIMENSION 3						// The dimension of the system

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
#include "Ply.h"
#include "VertexFactory.h"
#include "Image.h"
#include "RegularGrid.h"

MessageWriter messageWriter;

double BaseSSDWeights[] = { 5e+1f , 5e-4f , 1e-5f }; 

enum NormalType
{
	NORMALS_NONE ,
	NORMALS_SAMPLES ,
	NORMALS_GRADIENTS ,
	NORMALS_COUNT
};
const char* NormalsNames[] = { "none" , "samples" , "gradients" };

cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" ) ,
	TempDir( "tempDir" ) ,
	Grid( "grid" ) ,	
	Tree( "tree" ) ,
	Transform( "xForm" );

cmdLineReadable
	Performance( "performance" ) ,
	ShowResidual( "showResidual" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	ASCII( "ascii" ) ,
	Density( "density" ) ,
	NonLinearFit( "nonLinearFit" ) ,
	PrimalGrid( "primalGrid" ) ,
	ExactInterpolation( "exact" ) ,
	Colors( "colors" ) ,
	InCore( "inCore" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
#ifndef FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // !FAST_COMPILE
	Depth( "depth" , 8 ) ,
	KernelDepth( "kernelDepth" ) ,
	SolveDepth( "solveDepth" ) ,
	Iters( "iters" , 8 ) ,
	FullDepth( "fullDepth" , 5 ) ,
	BaseDepth( "baseDepth" ) ,
	BaseVCycles( "baseVCycles" , 4 ) ,
#ifndef FAST_COMPILE
	BType( "bType" , DEFAULT_FEM_BOUNDARY+1 ) ,
#endif // !FAST_COMPILE
	Normals( "normals" , NORMALS_NONE ) ,
	MaxMemoryGB( "maxMemory" , 0 ) ,
#ifdef _OPENMP
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
#else // !_OPENMP
	ParallelType( "parallel" , (int)ThreadPool::THREAD_POOL ) ,
#endif // _OPENMP
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "chunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)numMaxThreads);

cmdLineParameter< float >
	DataX( "data" , 32.f ) ,
	SamplesPerNode( "samplesPerNode" , 1.5f ) ,
	Scale( "scale" , 1.1f ) ,
	Width( "width" , 0.f ) ,
	Confidence( "confidence" , 0.f ) ,
	ConfidenceBias( "confidenceBias" , 0.f ) ,
	LowDepthCutOff( "lowDepthCutOff" , 0.f ) ,
	CGSolverAccuracy( "cgAccuracy" , 1e-3f ) ,
	ValueWeight   (    "valueWeight" , 1.f ) ,
	GradientWeight( "gradientWeight" , 1.f ) ,
	BiLapWeight   (    "biLapWeight" , 1.f );


cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
	&Degree , &BType ,
#endif // !FAST_COMPILE
	&SolveDepth ,
	&In , &Depth , &Out , &Transform ,
	&Width ,
	&Scale , &Verbose , &CGSolverAccuracy , &NoComments ,
	&KernelDepth , &SamplesPerNode , &Confidence , &NonManifold , &PolygonMesh , &ASCII , &ShowResidual ,
	&ConfidenceBias ,
	&ValueWeight , &GradientWeight , &BiLapWeight ,
	&Grid , &Threads ,
	&Tree ,
	&Density ,
	&FullDepth ,
	&BaseDepth , &BaseVCycles ,
	&Iters ,
	&DataX ,
	&Colors ,
	&Normals ,
	&NonLinearFit ,
	&PrimalGrid ,
	&TempDir ,
	&ExactInterpolation ,
	&Performance ,
	&MaxMemoryGB ,
	&InCore ,
	&ParallelType ,
	&ScheduleType ,
	&ThreadChunkSize ,
	&LowDepthCutOff ,
	NULL
};

void ShowUsage(char* ex)
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>\n" , In.name );
	printf( "\t[--%s <ouput triangle mesh>]\n" , Out.name );
	printf( "\t[--%s <ouput grid>]\n" , Grid.name );
	printf( "\t[--%s <ouput fem tree>]\n" , Tree.name );
#ifndef FAST_COMPILE
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
	printf( "\t[--%s <boundary type>=%d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , BoundaryNames[i] );
#endif // !FAST_COMPILE
	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s <maximum solution depth>=%d]\n" , SolveDepth.name , SolveDepth.value );
	printf( "\t[--%s <grid width>]\n" , Width.name );
	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <coarse MG solver depth>]\n" , BaseDepth.name );
	printf( "\t[--%s <coarse MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );
	printf( "\t[--%s <minimum number of samples per node>=%f]\n" , SamplesPerNode.name, SamplesPerNode.value );
	printf( "\t[--%s <zero-crossing weight>=%.3e]\n" , ValueWeight.name , ValueWeight.value );
	printf( "\t[--%s <gradient weight>=%.3e]\n" , GradientWeight.name , GradientWeight.value );
	printf( "\t[--%s <bi-laplacian weight>=%.3e]\n" , BiLapWeight.name , BiLapWeight.value );
	printf( "\t[--%s <iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s]\n" , ExactInterpolation.name );
	printf( "\t[--%s <pull factor>=%f]\n" , DataX.name , DataX.value );
	printf( "\t[--%s]\n" , Colors.name );
	printf( "\t[--%s <normal type>=%d]\n" , Normals.name , Normals.value );
	for( int i=0 ; i<NORMALS_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , NormalsNames[i] );
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <parallel type>=%d]\n" , ParallelType.name , ParallelType.value );
	for( size_t i=0 ; i<ThreadPool::ParallelNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ParallelNames[i].c_str() );
	printf( "\t[--%s <schedue type>=%d]\n" , ScheduleType.name , ScheduleType.value );
	for( size_t i=0 ; i<ThreadPool::ScheduleNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ScheduleNames[i].c_str() );
	printf( "\t[--%s <thread chunk size>=%d]\n" , ThreadChunkSize.name , ThreadChunkSize.value );
	printf( "\t[--%s <normal confidence exponent>=%f]\n" , Confidence.name , Confidence.value );
	printf( "\t[--%s <normal confidence bias exponent>=%f]\n" , ConfidenceBias.name , ConfidenceBias.value );
	printf( "\t[--%s <low depth cut-off>=%f]\n" , LowDepthCutOff.name , LowDepthCutOff.value );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s <cg solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s]\n" , Performance.name );
	printf( "\t[--%s]\n" , Density.name );
	printf( "\t[--%s]\n" , NonLinearFit.name );
	printf( "\t[--%s]\n" , PrimalGrid.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , NoComments.name );
	printf( "\t[--%s]\n" , TempDir.name );
	printf( "\t[--%s]\n" , InCore.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

double Weight( double v , double start , double end )
{
	v = ( v - start ) / ( end - start );
	if     ( v<0 ) return 1.;
	else if( v>1 ) return 0.;
	else
	{
		// P(x) = a x^3 + b x^2 + c x + d
		//		P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
		// =>	d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
		// =>	c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
		// =>	a = 2 , b = -3 , c = 0 , d = 1
		// =>	P(x) = 2 x^3 - 3 x^2 + 1
		return 2. * v * v * v - 3. * v * v + 1.;
	}
}

template< unsigned int Dim , class Real >
struct FEMTreeProfiler
{
	double t;

	void start( void ){ t = Time() , FEMTree< Dim , Real >::ResetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput2( std::vector< std::string >& comments , const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %d (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
};

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real scaleFactor )
{
	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = max[0] - min[0];
	for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
	scale *= scaleFactor;
	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real width , Real scaleFactor , int& depth )
{
	// Get the target resolution (along the largest dimension)
	Real resolution = ( max[0]-min[0] ) / width;
	for( int d=1 ; d<Dim ; d++ ) resolution = std::max< Real >( resolution , ( max[d]-min[d] ) / width );
	resolution *= scaleFactor;
	depth = 0;
	while( (1<<depth)<resolution ) depth++;

	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = (1<<depth) * width;

	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}

template< typename Real , unsigned int Dim , typename AuxData >
using InputOrientedPointStreamInfo = typename FEMTreeInitializer< Dim , Real >::template InputPointStream< VectorTypeUnion< Real , typename VertexFactory::NormalFactory< Real , Dim >::VertexType , AuxData > >;

template< typename Real , unsigned int Dim , typename AuxData >
using InputOrientedPointStream = typename InputOrientedPointStreamInfo< Real , Dim , AuxData >::StreamType;

template< class Real , unsigned int Dim , typename AuxData >
XForm< Real , Dim+1 > GetPointXForm( InputOrientedPointStream< Real , Dim , AuxData > &stream , Real width , Real scaleFactor , int& depth )
{
	Point< Real , Dim > min , max;
	InputOrientedPointStreamInfo< Real , Dim , AuxData >::BoundingBox( stream , min , max );
	return GetBoundingBoxXForm( min , max , width , scaleFactor , depth );
}
template< class Real , unsigned int Dim , typename AuxData >
XForm< Real , Dim+1 > GetPointXForm( InputOrientedPointStream< Real , Dim , AuxData > &stream , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	InputOrientedPointStreamInfo< Real , Dim , AuxData >::BoundingBox( stream , min , max );
	return GetBoundingBoxXForm( min , max , scaleFactor );
}

template< unsigned int Dim , typename Real , typename TotalPointSampleData >
struct ConstraintDual
{
	Real target , vWeight , gWeight;
	ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const TotalPointSampleData& data ) const 
	{
		Point< Real , Dim > n = data.template get<0>();
		CumulativeDerivativeValues< Real , Dim , 1 > cdv;
		cdv[0] = target*vWeight;
		for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
		return cdv;
	}
};
template< unsigned int Dim , typename Real , typename TotalPointSampleData >
struct SystemDual
{
	CumulativeDerivativeValues< Real , Dim , 1 > weight;
	SystemDual( Real v , Real g )
	{
		weight[0] = v;
		for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
	}
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	}
	CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	};
};
template< unsigned int Dim , class TotalPointSampleData >
struct SystemDual< Dim , double , TotalPointSampleData >
{
	typedef double Real;
	CumulativeDerivativeValues< Real , Dim , 1 > weight;
	SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	}
};

template< typename Real , typename SetVertexFunction , typename InputSampleDataType , typename VertexFactory , unsigned int ... FEMSigs >
void ExtractMesh
(
	UIntPack< FEMSigs ... > ,
	FEMTree< sizeof ... ( FEMSigs ) , Real >& tree ,
	const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution ,
	Real isoValue ,
	const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample > *samples ,
	std::vector< InputSampleDataType > *sampleData ,
	const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE > *density ,
	const VertexFactory &vertexFactory ,
	const InputSampleDataType &zeroInputSampleDataType ,
	SetVertexFunction SetVertex ,
	std::vector< std::string > &comments ,
	XForm< Real , sizeof...(FEMSigs)+1 > unitCubeToModel
)
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef typename VertexFactory::VertexType Vertex;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;

	FEMTreeProfiler< Dim , Real > profiler;

	char tempHeader[1024];
	{
		char tempPath[1024];
		tempPath[0] = 0;
		if( TempDir.set ) strcpy( tempPath , TempDir.value );
		else SetTempDirectory( tempPath , sizeof(tempPath) );
		if( strlen(tempPath)==0 ) sprintf( tempPath , ".%c" , FileSeparator );
		if( tempPath[ strlen( tempPath )-1 ]==FileSeparator ) sprintf( tempHeader , "%sPR_" , tempPath );
		else                                                  sprintf( tempHeader , "%s%cPR_" , tempPath , FileSeparator );
	}

	CoredMeshData< Vertex , node_index_type > *mesh;
	if( InCore.set ) mesh = new CoredVectorMeshData< Vertex , node_index_type >();
	else             mesh = new CoredFileMeshData< node_index_type , VertexFactory >( vertexFactory , tempHeader );
	profiler.start();
	typename IsoSurfaceExtractor< Dim , Real , Vertex >::IsoStats isoStats;
	if( sampleData )
	{
		SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setExtrapolatedDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
		for( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
		{
			ProjectiveData< InputSampleDataType , Real >* clr = _sampleData( n );
			if( clr ) (*clr) *= (Real)pow( DataX.value , tree.depth( n ) );
		}
		isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , &_sampleData , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , NonLinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
	}
#if defined( __GNUC__ ) && __GNUC__ < 5
#ifdef SHOW_WARNINGS
#warning "you've got me gcc version<5"
#endif // SHOW_WARNINGS
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , (SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , NonLinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , NonLinearFit.set , Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
	messageWriter( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)( mesh->outOfCoreVertexNum()+mesh->inCoreVertices.size() ) , (unsigned long long)mesh->polygonNum() );

	std::string isoStatsString = isoStats.toString() + std::string( "\n" );
	messageWriter( isoStatsString.c_str() );
	if( PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
	else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

	std::vector< std::string > noComments;
	typename VertexFactory::Transform unitCubeToModelTransform( unitCubeToModel );
	PLY::WritePolygons< VertexFactory , node_index_type , Real , Dim >( Out.value , vertexFactory , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , unitCubeToModelTransform );

	delete mesh;
}

template< typename Real , unsigned int Dim >
void WriteGrid( const char *fileName , ConstPointer( Real ) values , unsigned int res , XForm< Real , Dim+1 > voxelToModel , bool verbose )
{
	char *ext = GetFileExtension( fileName );

	if( Dim==2 && ImageWriter::ValidExtension( ext ) )
	{
		unsigned int totalResolution = 1;
		for( int d=0 ; d<Dim ; d++ ) totalResolution *= res;

		// Compute average
		Real avg = 0;
		std::vector< Real > avgs( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , totalResolution , [&]( unsigned int thread , size_t i ){ avgs[thread] += values[i]; } );
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) avg += avgs[t];
		avg /= (Real)totalResolution;

		// Compute standard deviation
		Real std = 0;
		std::vector< Real > stds( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , totalResolution , [&]( unsigned int thread , size_t i ){ stds[thread] += ( values[i] - avg ) * ( values[i] - avg ); } );
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) std += stds[t];
		std = (Real)sqrt( std / totalResolution );

		if( verbose )
		{
			printf( "Grid to image: [%.2f,%.2f] -> [0,255]\n" , avg - 2*std , avg + 2*std );
			printf( "Transform:\n" );
			for( int i=0 ; i<Dim+1 ; i++ )
			{
				printf( "\t" );
				for( int j=0 ; j<Dim+1 ; j++ ) printf( " %f" , voxelToModel(j,i) );
				printf( "\n" );
			}
		}

		unsigned char *pixels = new unsigned char[ totalResolution*3 ];
		ThreadPool::Parallel_for( 0 , totalResolution , [&]( unsigned int , size_t i )
		{
			Real v = (Real)std::min< Real >( (Real)1. , std::max< Real >( (Real)-1. , ( values[i] - avg ) / (2*std ) ) );
			v = (Real)( ( v + 1. ) / 2. * 256. );
			unsigned char color = (unsigned char )std::min< Real >( (Real)255. , std::max< Real >( (Real)0. , v ) );
			for( int c=0 ; c<3 ; c++ ) pixels[i*3+c ] = color;
		}
		);
		ImageWriter::Write( fileName , pixels , res , res , 3 );
		delete[] pixels;
	}
	else if( !strcasecmp( ext , "iso" ) )
	{
		FILE *fp = fopen( fileName , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: " , fileName );
		int r = (int)res;
		fwrite( &r , sizeof(int) , 1 , fp );
		size_t count = 1;
		for( unsigned int d=0 ; d<Dim ; d++ ) count *= res;
		fwrite( values , sizeof(Real) , count , fp );
		fclose( fp );
	}
	else
	{
		unsigned int _res[Dim];
		for( int d=0 ; d<Dim ; d++ ) _res[d] = res;
		RegularGrid< Real , Dim >::Write( fileName , _res , values , voxelToModel );
	}
	delete[] ext;
}

template< class Real , typename AuxDataFactory , unsigned int ... FEMSigs >
void Execute( UIntPack< FEMSigs ... > , const AuxDataFactory &auxDataFactory )
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
	typedef UIntPack< FEMDegreeAndBType< NORMAL_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;
	static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 1 > InterpolationInfo;
	using namespace VertexFactory;

	// The factory for constructing an input sample
	typedef Factory< Real , PositionFactory< Real , Dim > , Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > > InputSampleFactory;

	// The factory for constructing an input sample's data
	typedef Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > InputSampleDataFactory;

	// The input point stream information: First piece of data is the normal; the remainder is the auxiliary data
	typedef InputOrientedPointStreamInfo< Real , Dim , typename AuxDataFactory::VertexType > InputPointStreamInfo;

	// The type of the input sample
	typedef typename InputPointStreamInfo::PointAndDataType InputSampleType;

	// The type of the input sample's data
	typedef typename InputPointStreamInfo::DataType InputSampleDataType;

	typedef            InputDataStream< InputSampleType >  InputPointStream;
	typedef TransformedInputDataStream< InputSampleType > XInputPointStream;

	InputSampleFactory inputSampleFactory( PositionFactory< Real , Dim >() , InputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory ) );
	InputSampleDataFactory inputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory );

	typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;
	std::vector< std::string > comments;
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "** Running SSD Reconstruction (Version %s) **\n" , VERSION );
	messageWriter( comments , "************************************************\n" );
	messageWriter( comments , "************************************************\n" );
	if( !Threads.set ) messageWriter( comments , "Running with %d threads\n" , Threads.value );

	bool needNormalData = DataX.value>0 && Normals.value;
	bool needAuxData = DataX.value>0 && auxDataFactory.bufferSize();

	ThreadPool::Init( (ThreadPool::ParallelType)ParallelType.value , Threads.value );

	XForm< Real , Dim+1 > modelToUnitCube , unitCubeToModel;
	if( Transform.set )
	{
		FILE* fp = fopen( Transform.value , "r" );
		if( !fp )
		{
			WARN( "Could not read x-form from: " , Transform.value );
			modelToUnitCube = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<Dim+1 ; i++ ) for( int j=0 ; j<Dim+1 ; j++ )
			{
				float f;
				if( fscanf( fp , " %f " , &f )!=1 ) ERROR_OUT( "Failed to read xform" );
				modelToUnitCube(i,j) = (Real)f;
			}
			fclose( fp );
		}
	}
	else modelToUnitCube = XForm< Real , Dim+1 >::Identity();

	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}

	double startTime = Time();
	Real isoValue = 0;

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	FEMTreeProfiler< Dim , Real > profiler;

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --" , Depth.name , " and --" , Width.name , " set, ignoring --" , Width.name );
		Width.value = 0;
	}

	size_t pointCount;

	Real pointWeightSum;
	std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
	std::vector< InputSampleDataType >* sampleData = NULL;
	DensityEstimator* density = NULL;
	SparseNodeData< Point< Real , Dim > , NormalSigs >* normalInfo = NULL;
	Real targetValue = (Real)0.;

	// Read in the samples (and color data)
	{
		profiler.start();
		InputPointStream* pointStream;
		char* ext = GetFileExtension( In.value );
		sampleData = new std::vector< InputSampleDataType >();
		std::vector< InputSampleType > inCorePoints;
		if( InCore.set )
		{
			InputPointStream *_pointStream;
			if     ( !strcasecmp( ext , "bnpts" ) ) _pointStream = new BinaryInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else if( !strcasecmp( ext , "ply"   ) ) _pointStream = new    PLYInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else                                    _pointStream = new  ASCIIInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			InputSampleType p;
			while( _pointStream->next( p ) ) inCorePoints.push_back( p );
			delete _pointStream;

			pointStream = new MemoryInputDataStream< InputSampleType >( inCorePoints.size() , &inCorePoints[0] );		}
		else
		{
			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
			else                                    pointStream = new  ASCIIInputDataStream< InputSampleFactory >( In.value , inputSampleFactory );
		}
		delete[] ext;

		typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
		auto XFormFunctor = [&]( InputSampleType &p ){ p = _modelToUnitCube( p ); };
		XInputPointStream _pointStream( XFormFunctor , *pointStream );
		if( Width.value>0 ) modelToUnitCube = GetPointXForm< Real , Dim , typename AuxDataFactory::VertexType >( _pointStream , Width.value , (Real)( Scale.value>0 ? Scale.value : 1. ) , Depth.value ) * modelToUnitCube;
		else                modelToUnitCube = Scale.value>0 ? GetPointXForm< Real , Dim , typename AuxDataFactory::VertexType >( _pointStream , (Real)Scale.value ) * modelToUnitCube : modelToUnitCube;

		if( !SolveDepth.set ) SolveDepth.value = Depth.value;
		if( SolveDepth.value>Depth.value )
		{
			WARN( "Solution depth cannot exceed system depth: " , SolveDepth.value , " <= " , Depth.value );
			SolveDepth.value = Depth.value;
		}

		{
			typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
			auto XFormFunctor = [&]( InputSampleType &p ){ p = _modelToUnitCube( p ); };
			XInputPointStream _pointStream( XFormFunctor , *pointStream );
			auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
			{
				Real l = (Real)Length( d.template get<0>() );
				if( !l || !std::isfinite( l ) ) return (Real)-1.;
				return (Real)pow( l , Confidence.value );
			};
			auto ProcessData = []( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
			{
				Real l = (Real)Length( d.template get<0>() );
				if( !l || !std::isfinite( l ) ) return (Real)-1.;
				d.template get<0>() /= l;
				return (Real)1.;
			};
			if( Confidence.value>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( tree.spaceRoot() , _pointStream , Depth.value , *samples , *sampleData , true , tree.nodeAllocators.size() ? tree.nodeAllocators[0] : NULL , tree.initializer() , ProcessDataWithConfidence );
			else                     pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( tree.spaceRoot() , _pointStream , Depth.value , *samples , *sampleData , true , tree.nodeAllocators.size() ? tree.nodeAllocators[0] : NULL , tree.initializer() , ProcessData );
		}
		unitCubeToModel = modelToUnitCube.inverse();
		delete pointStream;

		messageWriter( "Input Points / Samples: %llu / %llu\n" , pointCount , (unsigned long long)samples->size() );
		profiler.dumpOutput2( comments , "# Read input into tree:" );
	}
	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
	if( kernelDepth>Depth.value )
	{
		WARN( KernelDepth.name , " can't be greater than " , Depth.name , ": " , KernelDepth.value , " <= " , Depth.value );
		kernelDepth = Depth.value;
	}

	DenseNodeData< Real , Sigs > solution;
	{
		DenseNodeData< Real , Sigs > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = Depth.value;

		tree.resetNodeIndices( 0 );

		// Get the kernel density estimator
		{
			profiler.start();
			density = tree.template setDensityEstimator< 1 , WEIGHT_DEGREE >( *samples , kernelDepth , SamplesPerNode.value );
			profiler.dumpOutput2( comments , "#   Got kernel density:" );
		}

		// Transform the Hermite samples into a vector field
		{
			profiler.start();
			normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
			std::function< bool ( InputSampleDataType , Point< Real , Dim >& ) > ConversionFunction = []( InputSampleDataType in , Point< Real , Dim > &out )
			{
				Point< Real , Dim > n = in.template get<0>();
				Real l = (Real)Length( n );
				// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
				if( !l ) return false;
				out = n / l;
				return true;
			};
			std::function< bool ( InputSampleDataType , Point< Real , Dim >& , Real & ) > ConversionAndBiasFunction = []( InputSampleDataType in , Point< Real , Dim > &out , Real &bias )
			{
				Point< Real , Dim > n = in.template get<0>();
				Real l = (Real)Length( n );
				// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
				if( !l ) return false;
				out = n / l;
				bias = (Real)( log( l ) * ConfidenceBias.value / log( 1<<(Dim-1) ) );
				return true;
			};

			if( ConfidenceBias.value>0 ) *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , BaseDepth.value , Depth.value , (Real)LowDepthCutOff.value , pointWeightSum , ConversionAndBiasFunction );
			else                         *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , BaseDepth.value , Depth.value , (Real)LowDepthCutOff.value , pointWeightSum , ConversionFunction );
			profiler.dumpOutput2( comments , "#     Got normal field:" );
			messageWriter( "Point weight / Estimated Measure: %g / %g\n" , pointWeightSum , pointCount*pointWeightSum );
		}

		if( !Density.set ) delete density , density = NULL;

		// Add the interpolation constraints
		if( ValueWeight.value>0 || GradientWeight.value>0 )
		{
			profiler.start();
			if( ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , InputSampleDataType , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , InputSampleDataType >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , InputSampleDataType >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , false );
			else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , InputSampleDataType , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , InputSampleDataType >( targetValue , (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum  ) , SystemDual< Dim , Real , InputSampleDataType >( (Real)ValueWeight.value * pointWeightSum , (Real)GradientWeight.value * pointWeightSum ) , true , 1 );
			profiler.dumpOutput2( comments , "#Initialized point interpolation constraints:" );
		}

		// Trim the tree and prepare for multigrid
		{
			profiler.start();
			constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
			tree.template finalizeForMultigrid< MAX_DEGREE , Degrees::Max() >( BaseDepth.value , FullDepth.value , typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs >( *normalInfo ) , []( const FEMTreeNode * ){ return false; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , density ) );
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}

		// Free up the normal info [If we don't need it for subsequent iterations.]
		if( normalInfo ) delete normalInfo , normalInfo = NULL;

		// Add the interpolation constraints
		if( ValueWeight.value>0 || GradientWeight.value>0 )
		{
			profiler.start();
			constraints = tree.initDenseNodeData( Sigs() );
			tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
			if( !needNormalData && !needAuxData ) delete sampleData , sampleData = NULL;
		}

		messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() );
		messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = 0 , sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.iters = Iters.value , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.baseVCycles = BaseVCycles.value;
			typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 2 > > F( { 0. , 0. , (double)BiLapWeight.value } );
			solution = tree.solveSystem( Sigs() , F , constraints , SolveDepth.value , sInfo , std::make_tuple( iInfo ) );
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	{
		profiler.start();
		double valueSum = 0 , weightSum = 0;
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &tree , solution );
		std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
		{
			ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
		}
		);
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
		isoValue = (Real)( valueSum / weightSum );
		if( !needNormalData && !needAuxData ) delete samples , samples = NULL;
		profiler.dumpOutput( "Got average:" );
		messageWriter( "Iso-Value: %e = %g / %g\n" , isoValue , valueSum , weightSum );
	}
	if( Tree.set )
	{
		FILE* fp = fopen( Tree.value , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: " , Tree.value );
		FEMTree< Dim , Real >::WriteParameter( fp );
		DenseNodeData< Real , Sigs >::WriteSignatures( fp );
		tree.write( fp , modelToUnitCube );
		solution.write( fp );
		fclose( fp );
	}

	if( Grid.set )
	{
		int res = 0;
		profiler.start();
		Pointer( Real ) values = tree.template regularGridEvaluate< true >( solution , res , -1 , PrimalGrid.set );
		size_t resolution = 1;
		profiler.dumpOutput( "Got grid:" );
		XForm< Real , Dim+1 > voxelToUnitCube = XForm< Real , Dim+1 >::Identity();
		if( PrimalGrid.set ) for( int d=0 ; d<Dim ; d++ ) voxelToUnitCube( d , d ) = (Real)( 1. / (res-1) );
		else                 for( int d=0 ; d<Dim ; d++ ) voxelToUnitCube( d , d ) = (Real)( 1. / res ) , voxelToUnitCube( Dim , d ) = (Real)( 0.5 / res );
		WriteGrid< Real , DEFAULT_DIMENSION >( Grid.value , values , res , unitCubeToModel * voxelToUnitCube , Verbose.set );
		DeletePointer( values );
	}

	if( Out.set )
	{
		if( Normals.value )
		{
			if( Density.set )
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
				if( Normals.value==NORMALS_SAMPLES )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<0>() , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
				else if( Normals.value==NORMALS_GRADIENTS )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = -g/(1<<Depth.value) , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
			}
			else
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , auxDataFactory );
				if( Normals.value==NORMALS_SAMPLES )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = d.template get<0>() , v.template get<2>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
				else if( Normals.value==NORMALS_GRADIENTS )
				{
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = -g/(1<<Depth.value) , v.template get<2>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
			}
		}
		else
		{
			if( Density.set )
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
				auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = w , v.template get<2>() = d.template get<1>(); };
				ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
			}
			else
			{
				typedef Factory< Real , PositionFactory< Real , Dim > , AuxDataFactory > VertexFactory;
				VertexFactory vertexFactory( PositionFactory< Real , Dim >() , auxDataFactory );
				auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<1>(); };
				ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
			}
		}
		if( sampleData ){ delete sampleData ; sampleData = NULL; }
	}
	if( density ) delete density , density = NULL;
	messageWriter( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , FEMTree< Dim , Real >::MaxMemoryUsage() );
}

#ifndef FAST_COMPILE
template< unsigned int Dim , class Real , BoundaryType BType , typename AuxDataFactory >
void Execute( const AuxDataFactory &auxDataFactory )
{
	switch( Degree.value )
	{
//		case 1: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BType >::Signature >() , auxDataFactory );
		case 2: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BType >::Signature >() , auxDataFactory );
		case 3: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BType >::Signature >() , auxDataFactory );
//		case 4: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BType >::Signature >() , auxDataFactory );
		default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
	}
}

template< unsigned int Dim , class Real , typename AuxDataFactory >
void Execute( const AuxDataFactory &auxDataFactory )
{
	switch( BType.value )
	{
		case BOUNDARY_FREE+1:      return Execute< Dim , Real , BOUNDARY_FREE      >( auxDataFactory );
		case BOUNDARY_NEUMANN+1:   return Execute< Dim , Real , BOUNDARY_NEUMANN   >( auxDataFactory );
		case BOUNDARY_DIRICHLET+1: return Execute< Dim , Real , BOUNDARY_DIRICHLET >( auxDataFactory );
		default: ERROR_OUT( "Not a valid boundary type: " , BType.value );
	}
}
#endif // !FAST_COMPILE

int main( int argc , char* argv[] )
{
	Timer timer;
#ifdef USE_SEG_FAULT_HANDLER
	WARN( "using seg-fault handler" );
	StackTracer::exec = argv[0];
	signal( SIGSEGV , SignalHandler );
#endif // USE_SEG_FAULT_HANDLER
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG

	cmdLineParse( argc-1 , &argv[1] , params );
	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}
	if( GradientWeight.value<=0 ) ERROR_OUT( "Gradient weight must be positive: " , GradientWeight.value , "> 0" );
	if( BiLapWeight.value<=0 ) ERROR_OUT( "Bi-Laplacian weight must be positive: " , BiLapWeight.value , " > 0" );
	if( !BaseDepth.set ) BaseDepth.value = FullDepth.value;
	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}

	ValueWeight.value    *= (float)BaseSSDWeights[0];
	GradientWeight.value *= (float)BaseSSDWeights[1];
	BiLapWeight.value    *= (float)BaseSSDWeights[2];

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

#ifdef FAST_COMPILE
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
	typedef IsotropicUIntPack< DEFAULT_DIMENSION , FEMDegreeAndBType< Degree , BType >::Signature > FEMSigs;
	WARN( "Compiled for degree-" , Degree , ", boundary-" , BoundaryNames[ BType ] , ", " , sizeof(Real)==4 ? "single" : "double" , "-precision _only_" );

	char *ext = GetFileExtension( In.value );
	if( !strcasecmp( ext , "ply" ) )
	{
		typedef VertexFactory::Factory< Real , typename VertexFactory::PositionFactory< Real , DEFAULT_DIMENSION > , typename VertexFactory::NormalFactory< Real , DEFAULT_DIMENSION > > Factory;
		Factory factory;
		bool *readFlags = new bool[ factory.plyReadNum() ];
		std::vector< PlyProperty > unprocessedProperties;
		PLY::ReadVertexHeader( In.value , factory , readFlags , unprocessedProperties );
		if( !factory.plyValidReadProperties<0>( readFlags ) ) ERROR_OUT( "Ply file does not contain positions" );
		if( !factory.plyValidReadProperties<1>( readFlags ) ) ERROR_OUT( "Ply file does not contain normals" );
		delete[] readFlags;

		if( unprocessedProperties.size() ) Execute< Real >( FEMSigs() , VertexFactory::DynamicFactory< Real >( unprocessedProperties ) );
		else                               Execute< Real >( FEMSigs() , VertexFactory::EmptyFactory< Real >() );
	}
	else
	{
		if( Colors.set ) Execute< Real >( FEMSigs() , VertexFactory::RGBColorFactory< Real >() );
		else             Execute< Real >( FEMSigs() , VertexFactory::EmptyFactory< Real >() );
	}
	delete[] ext;
#else // !FAST_COMPILE
	char *ext = GetFileExtension( In.value );
	if( !strcasecmp( ext , "ply" ) )
	{
		typedef VertexFactory::Factory< Real , typename VertexFactory::PositionFactory< Real , DEFAULT_DIMENSION > , typename VertexFactory::NormalFactory< Real , DEFAULT_DIMENSION > > Factory;
		Factory factory;
		bool *readFlags = new bool[ factory.plyReadNum() ];
		std::vector< PlyProperty > unprocessedProperties;
		PLY::ReadVertexHeader( In.value , factory , readFlags , unprocessedProperties );
		if( !factory.plyValidReadProperties<0>( readFlags ) ) ERROR_OUT( "Ply file does not contain positions" );
		if( !factory.plyValidReadProperties<1>( readFlags ) ) ERROR_OUT( "Ply file does not contain normals" );
		delete[] readFlags;

		if( unprocessedProperties.size() ) Execute< DEFAULT_DIMENSION , Real >( VertexFactory::DynamicFactory< Real >( unprocessedProperties ) );
		else                               Execute< DEFAULT_DIMENSION , Real >( VertexFactory::EmptyFactory< Real >() );
	}
	else
	{
		if( Colors.set ) Execute< DEFAULT_DIMENSION , Real >( VertexFactory::RGBColorFactory< Real >() );
		else             Execute< DEFAULT_DIMENSION , Real >( VertexFactory::EmptyFactory< Real >() );
	}
	delete[] ext;
#endif // FAST_COMPILE
	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}
	return EXIT_SUCCESS;
}
