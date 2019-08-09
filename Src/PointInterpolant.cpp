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

//#undef USE_DOUBLE								// If enabled, double-precesion is used
#define USE_DOUBLE								// If enabled, double-precesion is used
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define DEFAULT_FEM_DEGREE 2					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_FREE		// The default finite-element boundary type
#define DEFAULT_DIMENSION 2						// The dimension of the system

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
#include "Ply.h"
#include "PointStreamData.h"
#include "Image.h"

MessageWriter messageWriter;

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
	PrimalGrid( "primalGrid" ) ,
	ExactInterpolation( "exact" ) ,
	InCore( "inCore" ) ,
	NoValueConstraints( "noValues" ) ,
	UseGradientConstraints( "useGradients" ) ,
	NoComments( "noComments" ) ,
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	NonLinearFit( "nonLinearFit" ) ,
	ASCII( "ascii" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
#ifndef FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // !FAST_COMPILE
	Depth( "depth" , 8 ) ,
	Iters( "iters" , 8 ) ,
	FullDepth( "fullDepth" , 5 ) ,
	BaseDepth( "baseDepth" , 5 ) ,
	BaseVCycles( "baseVCycles" , 4 ) ,
#ifndef FAST_COMPILE
	BType( "bType" , DEFAULT_FEM_BOUNDARY+1 ) ,
	Dimension( "dim" , DEFAULT_DIMENSION ) ,
#endif // !FAST_COMPILE
	MaxMemoryGB( "maxMemory" , 0 ) ,
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "chunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() );

cmdLineParameter< float >
	Scale( "scale" , 1.1f ) ,
	Width( "width" , 0.f ) ,
	CGSolverAccuracy( "cgAccuracy" , 1e-3f ) ,
	IsoValue( "iso" , 0.f ) ,
	ValueWeight   (    "valueWeight" , 1000.f ) ,
	GradientWeight( "gradientWeight" , 1.f ) ,
	LapWeight     (      "lapWeight" , 0.f ) ,
	BiLapWeight   (    "biLapWeight" , 1.f );

cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
	&Degree , &BType , &Dimension ,
#endif // !FAST_COMPILE
	&In , &Out , &Depth , &Transform ,
	&Width ,
	&Scale , &Verbose , &CGSolverAccuracy , &NoComments ,
	&NonManifold , &PolygonMesh , &ASCII , &ShowResidual ,
	&ValueWeight , &GradientWeight ,
	&LapWeight , &BiLapWeight ,
	&Grid , &Threads ,
	&Tree ,
	&FullDepth ,
	&BaseDepth , &BaseVCycles ,
	&Iters ,
	&IsoValue ,
	&PrimalGrid ,
	&ExactInterpolation ,
	&Performance ,
	&MaxMemoryGB ,
	&InCore ,
	&ParallelType ,
	&ScheduleType ,
	&ThreadChunkSize ,
	&NoValueConstraints ,
	&UseGradientConstraints ,
	&NonLinearFit ,
	NULL
};

void ShowUsage(char* ex)
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>\n" , In.name );
	printf( "\t[--%s <ouput mesh>]\n" , Out.name );
	printf( "\t[--%s <ouput grid>]\n" , Grid.name );
	printf( "\t[--%s <ouput fem tree>]\n" , Tree.name );
#ifndef FAST_COMPILE
	printf( "\t[--%s <dimension>=%d]\n" , Dimension.name , Dimension.value );
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
	printf( "\t[--%s <boundary type>=%d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , BoundaryNames[i] );
#endif // !FAST_COMPILE
	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s <grid width>]\n" , Width.name );
	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <coarse MG solver depth>=%d]\n" , BaseDepth.name , BaseDepth.value );
	printf( "\t[--%s <coarse MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );
	printf( "\t[--%s <zero-crossing weight>=%.3e]\n" , ValueWeight.name , ValueWeight.value );
	printf( "\t[--%s <gradient weight>=%.3e]\n" , GradientWeight.name , GradientWeight.value );
	printf( "\t[--%s <laplacian weight>=%.3e]\n" , LapWeight.name , LapWeight.value );
	printf( "\t[--%s <bi-laplacian weight>=%.3e]\n" , BiLapWeight.name , BiLapWeight.value );
	printf( "\t[--%s <iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s]\n" , ExactInterpolation.name );
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <parallel type>=%d]\n" , ParallelType.name , ParallelType.value );
	for( size_t i=0 ; i<ThreadPool::ParallelNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ParallelNames[i].c_str() );
	printf( "\t[--%s <schedue type>=%d]\n" , ScheduleType.name , ScheduleType.value );
	for( size_t i=0 ; i<ThreadPool::ScheduleNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ScheduleNames[i].c_str() );
	printf( "\t[--%s <thread chunk size>=%d]\n" , ThreadChunkSize.name , ThreadChunkSize.value );
	printf( "\t[--%s <cg solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s <iso-value>=%f]\n" , IsoValue.name , IsoValue.value );
	printf( "\t[--%s]\n" , NoValueConstraints.name );
	printf( "\t[--%s]\n" , UseGradientConstraints.name );
	printf( "\t[--%s]\n" , Performance.name );
	printf( "\t[--%s]\n" , PrimalGrid.name );
	printf( "\t[--%s]\n" , NoComments.name );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , NonLinearFit.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , InCore.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< unsigned int Dim , class Real >
struct FEMTreeProfiler
{
	FEMTree< Dim , Real >& tree;
	double t;

	FEMTreeProfiler( FEMTree< Dim , Real >& t ) : tree(t) { ; }
	void start( void ){ t = Time() , FEMTree< Dim , Real >::ResetLocalMemoryUsage(); }
	void print( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) printf( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         printf(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput( const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter(    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
	}
	void dumpOutput2( std::vector< std::string >& comments , const char* header ) const
	{
		FEMTree< Dim , Real >::MemoryUsage();
		if( header ) messageWriter( comments , "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" , header , Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
		else         messageWriter( comments ,    "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n" ,          Time()-t , FEMTree< Dim , Real >::LocalMemoryUsage() , FEMTree< Dim , Real >::MaxMemoryUsage() , MemoryInfo::PeakMemoryUsageMB() );
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

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real width , Real scaleFactor , int& depth )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , width , scaleFactor , depth );
}
template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetPointXForm( InputPointStream< Real , Dim >& stream , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	stream.boundingBox( min , max );
	return GetBoundingBoxXForm( min , max , scaleFactor );
}

template< unsigned int Dim , typename Real , typename TotalPointSampleData > struct ValueAndGradientFromSample;

template< unsigned int Dim , typename Real >
struct ValueAndGradientFromSample< Dim , Real , MultiPointStreamData< Real , PointStreamValue< Real > , PointStreamNormal< Real , Dim > > >
{
	typedef MultiPointStreamData< Real , PointStreamValue< Real > , PointStreamNormal< Real , Dim > > TotalPointSampleData;
	std::pair< Real , Point< Real , Dim > > operator()( TotalPointSampleData d ) const { return std::pair< Real , Point< Real , Dim > >( d.template data<0>() , d.template data<1>() ); }
};

template< unsigned int Dim , typename Real >
struct ValueAndGradientFromSample< Dim , Real , MultiPointStreamData< Real , PointStreamValue< Real > > >
{
	typedef MultiPointStreamData< Real , PointStreamValue< Real > > TotalPointSampleData;
	std::pair< Real , Point< Real , Dim > > operator()( TotalPointSampleData d ) const { return std::pair< Real , Point< Real , Dim > >( d.template data<0>() , Point< Real , Dim >() ); }
};

template< unsigned int Dim , typename Real >
struct ValueAndGradientFromSample< Dim , Real , MultiPointStreamData< Real , PointStreamNormal< Real , Dim > > >
{
	typedef MultiPointStreamData< Real , PointStreamNormal< Real , Dim > > TotalPointSampleData;
	std::pair< Real , Point< Real , Dim > > operator()( TotalPointSampleData d ) const { return std::pair< Real , Point< Real , Dim > >( (Real)0 , d.template data<0>() ); }
};


template< unsigned int Dim , typename Real , typename TotalPointSampleData > struct ConstraintDual;

template< unsigned int Dim , typename Real >
struct ConstraintDual< Dim , Real , MultiPointStreamData< Real , PointStreamValue< Real > , PointStreamNormal< Real , Dim > > >
{
	typedef MultiPointStreamData< Real , PointStreamValue< Real > , PointStreamNormal< Real , Dim > > TotalPointSampleData;
	Real vWeight , gWeight;
	ConstraintDual( Real v , Real g ) : vWeight(v) , gWeight(g) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const TotalPointSampleData& data ) const 
	{
		Real value = data.template data<0>();
		Point< Real , Dim > gradient = data.template data<1>();
		CumulativeDerivativeValues< Real , Dim , 1 > cdv;
		cdv[0] = value*vWeight;
		for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = gradient[d]*gWeight;
		return cdv;
	}
};
template< unsigned int Dim , typename Real >
struct ConstraintDual< Dim , Real , MultiPointStreamData< Real , PointStreamValue< Real > > >
{
	typedef MultiPointStreamData< Real , PointStreamValue< Real > > TotalPointSampleData;
	Real vWeight , gWeight;
	ConstraintDual( Real v , Real g ) : vWeight(v) , gWeight(g) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const TotalPointSampleData& data ) const 
	{
		Real value = data.template data<0>();
		CumulativeDerivativeValues< Real , Dim , 1 > cdv;
		cdv[0] = value*vWeight;
		return cdv;
	}
};
template< unsigned int Dim , typename Real >
struct ConstraintDual< Dim , Real , MultiPointStreamData< Real , PointStreamNormal< Real , Dim > > >
{
	typedef MultiPointStreamData< Real , PointStreamNormal< Real , Dim > > TotalPointSampleData;
	Real vWeight , gWeight;
	ConstraintDual( Real v , Real g ) : vWeight(v) , gWeight(g) { }
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const TotalPointSampleData& data ) const 
	{
		Point< Real , Dim > gradient = data.template data<0>();
		CumulativeDerivativeValues< Real , Dim , 1 > cdv;
		for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = gradient[d]*gWeight;
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
	SystemDual( Real v , Real g )
	{
		weight[0] = v;
		for( unsigned int d=0 ; d<Dim ; d++ ) weight[1+d] = g;
	}
	CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const TotalPointSampleData& data , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
	{
		return dValues * weight;
	}
};

template< typename Vertex , typename Real , unsigned int ... FEMSigs , typename TotalPointSampleData >
void ExtractMesh( UIntPack< FEMSigs ... > , FEMTree< sizeof ... ( FEMSigs ) , Real >& tree , const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution , Real isoValue , const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample >* samples , std::function< void ( Vertex& , Point< Real , sizeof ... ( FEMSigs ) > , Real , TotalPointSampleData ) > SetVertex , std::vector< std::string > &comments , XForm< Real , sizeof...(FEMSigs)+1 > iXForm )
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	static const unsigned int DataSig = FEMDegreeAndBType< WEIGHT_DEGREE , BOUNDARY_FREE >::Signature;

	FEMTreeProfiler< Dim , Real > profiler( tree );

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
	else             mesh = new CoredFileMeshData< Vertex , node_index_type >( tempHeader );
	profiler.start();
	typename IsoSurfaceExtractor< Dim , Real , Vertex >::IsoStats isoStats;
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
	isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , (typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE >*)NULL , (SparseNodeData< ProjectiveData< TotalPointSampleData , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#else // !__GNUC__ || __GNUC__ >=5
	isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< TotalPointSampleData >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , (typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE >*)NULL , NULL , solution , isoValue , *mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , false );
#endif // __GNUC__ || __GNUC__ < 4
	messageWriter( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)( mesh->outOfCorePointCount()+mesh->inCorePoints.size() ) , (unsigned long long)mesh->polygonCount() );
	std::string isoStatsString = isoStats.toString() + std::string( "\n" );
	messageWriter( isoStatsString.c_str() );
	if( PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
	else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

	std::vector< std::string > noComments;
	if( !PlyWritePolygons< Vertex , node_index_type , Real , Dim >( Out.value , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , iXForm ) )
		ERROR_OUT( "Could not write mesh to: " , Out.value );
	delete mesh;
}

template< typename Real , unsigned int Dim >
void WriteGrid( ConstPointer( Real ) values , int res , const char *fileName )
{
	int resolution = 1;
	for( int d=0 ; d<Dim ; d++ ) resolution *= res;

	char *ext = GetFileExtension( fileName );

	if( Dim==2 && ImageWriter::ValidExtension( ext ) )
	{
		Real avg = 0;
		std::vector< Real > avgs( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , resolution , [&]( unsigned int thread , size_t i ){ avgs[thread] += values[i]; } );
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) avg += avgs[t];
		avg /= (Real)resolution;

		Real std = 0;
		std::vector< Real > stds( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , resolution , [&]( unsigned int thread , size_t i ){ stds[thread] += ( values[i] - avg ) * ( values[i] - avg ); } );
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) std += stds[t];
		std = (Real)sqrt( std / resolution );

		if( Verbose.set ) printf( "Grid to image: [%.2f,%.2f] -> [0,255]\n" , avg - 2*std , avg + 2*std );

		unsigned char *pixels = new unsigned char[ resolution*3 ];
		ThreadPool::Parallel_for( 0 , resolution , [&]( unsigned int , size_t i )
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
	else
	{

		FILE *fp = fopen( fileName , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open grid file for writing: " , fileName );
		else
		{
			fwrite( &res , sizeof(int) , 1 , fp );
			if( typeid(Real)==typeid(float) ) fwrite( values , sizeof(float) , resolution , fp );
			else
			{
				float *fValues = new float[resolution];
				for( int i=0 ; i<resolution ; i++ ) fValues[i] = float( values[i] );
				fwrite( fValues , sizeof(float) , resolution , fp );
				delete[] fValues;
			}
			fclose( fp );
		}
	}
	delete[] ext;
}


template< class Real , typename TotalPointSampleData , unsigned int ... FEMSigs >
void Execute( int argc , char* argv[] , UIntPack< FEMSigs ... > )
{
	static const int Dim = sizeof ... ( FEMSigs );
	typedef UIntPack< FEMSigs ... > Sigs;
	typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
	typedef UIntPack< FEMDegreeAndBType< WEIGHT_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > DataSigs;
	typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 1 > InterpolationInfo;
	typedef InputPointStreamWithData< Real , Dim , TotalPointSampleData > InputPointStream;
	typedef TransformedInputPointStreamWithData< Real , Dim , TotalPointSampleData > XInputPointStream;
	std::vector< std::string > comments;
	messageWriter( comments , "***********************************************\n" );
	messageWriter( comments , "***********************************************\n" );
	messageWriter( comments , "** Running Point Interpolant (Version %s) **\n" , VERSION );
	messageWriter( comments , "***********************************************\n" );
	messageWriter( comments , "***********************************************\n" );
	if( !Threads.set ) messageWriter( comments , "Running with %d threads\n" , Threads.value );

	ThreadPool::Init( (ThreadPool::ParallelType)ParallelType.value , Threads.value );

	XForm< Real , Dim+1 > xForm , iXForm;
	if( Transform.set )
	{
		FILE* fp = fopen( Transform.value , "r" );
		if( !fp )
		{
			WARN( "Could not read x-form from: " , Transform.value );
			xForm = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<Dim+1 ; i++ ) for( int j=0 ; j<Dim+1 ; j++ )
			{
				float f;
				if( fscanf( fp , " %f " , &f )!=1 ) ERROR_OUT( "Failed to read xform" );
				xForm(i,j) = (Real)f;
			}
			fclose( fp );
		}
	}
	else xForm = XForm< Real , Dim+1 >::Identity();

	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}

	double startTime = Time();

	FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
	FEMTreeProfiler< Dim , Real > profiler( tree );

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --" , Depth.name , " and --" , Width.name , " set, ignoring --" , Width.name );
		Width.value = 0;
	}

	size_t pointCount;

	std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
	std::vector< TotalPointSampleData >* sampleData = NULL;

	// Read in the samples
	{
		profiler.start();
		InputPointStream* pointStream;
		char* ext = GetFileExtension( In.value );
		sampleData = new std::vector< TotalPointSampleData >();
		std::vector< std::pair< Point< Real , Dim > , TotalPointSampleData > > inCorePoints;
		if( InCore.set )
		{
			InputPointStream *_pointStream;
			if     ( !strcasecmp( ext , "bnpts" ) ) _pointStream = new BinaryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadBinary );
			else if( !strcasecmp( ext , "ply"   ) ) _pointStream = new    PLYInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::PlyReadProperties() , TotalPointSampleData::PlyReadNum , TotalPointSampleData::ValidPlyReadProperties );
			else                                    _pointStream = new  ASCIIInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadASCII );
			Point< Real , Dim > p;
			TotalPointSampleData d;
			while( _pointStream->nextPoint( p , d ) ) inCorePoints.push_back( std::pair< Point< Real , Dim > , TotalPointSampleData >( p , d ) );
			delete _pointStream;

			pointStream = new MemoryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( inCorePoints.size() , &inCorePoints[0] );
		}
		else
		{
			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadBinary );
			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::PlyReadProperties() , TotalPointSampleData::PlyReadNum , TotalPointSampleData::ValidPlyReadProperties );
			else                                    pointStream = new  ASCIIInputPointStreamWithData< Real , Dim , TotalPointSampleData >( In.value , TotalPointSampleData::ReadASCII );
		}
		delete[] ext;
		typename TotalPointSampleData::Transform _xForm( xForm );
		XInputPointStream _pointStream( [&]( Point< Real , Dim >& p , TotalPointSampleData& d ){ p = xForm*p , d = _xForm(d); } , *pointStream );
		if( Width.value>0 ) xForm = GetPointXForm< Real , Dim >( _pointStream , Width.value , (Real)( Scale.value>0 ? Scale.value : 1. ) , Depth.value ) * xForm;
		else                xForm = Scale.value>0 ? GetPointXForm< Real , Dim >( _pointStream , (Real)Scale.value ) * xForm : xForm;
		{
			typename TotalPointSampleData::Transform _xForm( xForm );
			XInputPointStream _pointStream( [&]( Point< Real , Dim >& p , TotalPointSampleData& d ){ p = xForm*p , d = _xForm(d); } , *pointStream );
			auto ProcessData = []( const Point< Real , Dim >& p , TotalPointSampleData& d ){ return (Real)1.; };
			pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< TotalPointSampleData >( tree.spaceRoot() , _pointStream , Depth.value , *samples , *sampleData , true , tree.nodeAllocators.size() ? tree.nodeAllocators[0] : NULL , tree.initializer() , ProcessData );
		}
		iXForm = xForm.inverse();
		delete pointStream;

		messageWriter( "Input Points / Samples: %llu / %llu\n" , pointCount , (unsigned long long)samples->size() );
		profiler.dumpOutput2( comments , "# Read input into tree:" );
	}

	DenseNodeData< Real , Sigs > solution;
	{
		DenseNodeData< Real , Sigs > constraints;
		InterpolationInfo* iInfo = NULL;
		int solveDepth = Depth.value;

		tree.resetNodeIndices();

		// Prepare for multigrid
		{
			profiler.start();
			tree.template finalizeForMultigrid< Degrees::Max() >( FullDepth.value , []( const typename FEMTree< Dim , Real >::FEMTreeNode * ){ return true; } );
			profiler.dumpOutput2( comments , "#       Finalized tree:" );
		}

		// Add the interpolation constraints
		{
			profiler.start();
			if( ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value , (Real)GradientWeight.value ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value , (Real)GradientWeight.value ) , true , false );
			else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , TotalPointSampleData , 1 >( tree , *samples , GetPointer( *sampleData ) , ConstraintDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value , (Real)GradientWeight.value ) , SystemDual< Dim , Real , TotalPointSampleData >( (Real)ValueWeight.value , (Real)GradientWeight.value ) , true , 1 );
			constraints = tree.initDenseNodeData( Sigs() );
			tree.addInterpolationConstraints( constraints , solveDepth , *iInfo );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
		}

		messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() );
		messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

		// Solve the linear system
		{
			profiler.start();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = 0 , sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.iters = Iters.value , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
			typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 2 > > F( { 0. , (double)LapWeight.value , (double)BiLapWeight.value } );
			solution = tree.solveSystem( Sigs() , F , constraints , solveDepth , sInfo , iInfo );
			profiler.dumpOutput2( comments , "# Linear system solved:" );
			if( iInfo ) delete iInfo , iInfo = NULL;
		}
	}

	if( Verbose.set )
	{
		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 1 > evaluator( &tree , solution );
		std::pair< double , double > valueStat(0,0) , gradientStat(0,0);
		std::vector< std::pair< double , double > > valueStats( ThreadPool::NumThreads() , std::pair< double , double >(0,0) ) , gradientStats( ThreadPool::NumThreads() , std::pair< double , double >(0,0) );
		ValueAndGradientFromSample< Dim , Real , TotalPointSampleData > valueAndGradientFromSample;
		ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
		{
			ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
			Real w = sample.weight;
			if( w>0 )
			{
				CumulativeDerivativeValues< Real , Dim , 1 > values = evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node );
				Real value = values[0];
				Point< Real , Dim > gradient;
				for( int d=0 ; d<Dim ; d++ ) gradient[d] = values[d+1];
				std::pair< Real , Point< Real , Dim > > valueAndGradient = valueAndGradientFromSample( (*sampleData)[j] / w );
				valueStats[ thread ].first += ( value - valueAndGradient.first ) * ( value - valueAndGradient.first ) * w;
				valueStats[ thread ].second += ( value * value + valueAndGradient.first * valueAndGradient.first ) * w;
				gradientStats[ thread ].first += Point< Real , Dim >::SquareNorm( gradient - valueAndGradient.second ) * w;
				gradientStats[ thread ].second += ( Point< Real , Dim >::SquareNorm( gradient ) + Point< Real , Dim >::SquareNorm( valueAndGradient.second ) ) * w;
			}
		}
		);
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) valueStat.first += valueStats[t].first , valueStat.second += valueStats[t].second , gradientStat.first += gradientStats[t].first , gradientStat.second += gradientStats[t].second;
		if( ValueWeight.value>0 && GradientWeight.value>0 ) messageWriter( "Value / Gradient Error: %g / %g\n" , (Real)sqrt( valueStat.first / valueStat.second ) , (Real)sqrt( gradientStat.first / gradientStat.second ) );
		else if( ValueWeight.value>0 ) messageWriter( "Value Error: %g\n" , (Real)sqrt( valueStat.first / valueStat.second ) );
		else if( GradientWeight.value>0 ) messageWriter( "Gradient Error: %g\n" , (Real)sqrt( gradientStat.first / gradientStat.second ) );
	}

	delete samples , samples = NULL;
	delete sampleData , sampleData = NULL;


	if( Tree.set )
	{
		FILE* fp = fopen( Tree.value , "wb" );
		if( !fp ) ERROR_OUT( "Failed to open file for writing: " , Tree.value );
		FEMTree< Dim , Real >::WriteParameter( fp );
		DenseNodeData< Real , Sigs >::WriteSignatures( fp );
		tree.write( fp , xForm );
		solution.write( fp );
		fclose( fp );
	}

	if( Grid.set )
	{
		int res = 0;
		profiler.start();
		Pointer( Real ) values = tree.template regularGridEvaluate< true >( solution , res , -1 , PrimalGrid.set );
		size_t resolution = 1;
		for( int d=0 ; d<Dim ; d++ ) resolution *= res;
		profiler.dumpOutput( "Got grid:" );
		WriteGrid< Real , Dim >( values , res , Grid.value );
		DeletePointer( values );
		if( Verbose.set )
		{
			printf( "Transform:\n" );
			for( int i=0 ; i<Dim+1 ; i++ )
			{
				printf( "\t" );
				for( int j=0 ; j<Dim+1 ; j++ ) printf( " %f" , iXForm(j,i) );
				printf( "\n" );
			}
		}
	}

	if( Out.set )
	{
		typedef PlyVertex< Real , Dim > Vertex;
		std::function< void ( Vertex& , Point< Real , Dim > , Real , TotalPointSampleData ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real , TotalPointSampleData ){ v.point = p; };
		ExtractMesh< Vertex >( UIntPack< FEMSigs ... >() , tree , solution , (Real)IsoValue.value , samples , SetVertex , comments , iXForm );
	}

	messageWriter( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , FEMTree< Dim , Real >::MaxMemoryUsage() );
}

template< class Real , unsigned int ... FEMSigs >
void Execute( int argc , char* argv[] , UIntPack< FEMSigs ... > )
{
	static const int Dim = sizeof ... ( FEMSigs );
	if     ( !UseGradientConstraints.set ) Execute< Real , MultiPointStreamData< Real , PointStreamValue< Real >                                   > >( argc , argv , UIntPack< FEMSigs ... >() );
	else if( NoValueConstraints.set )      Execute< Real , MultiPointStreamData< Real ,                            PointStreamNormal< Real , Dim > > >( argc , argv , UIntPack< FEMSigs ... >() );
	else                                   Execute< Real , MultiPointStreamData< Real , PointStreamValue< Real > , PointStreamNormal< Real , Dim > > >( argc , argv , UIntPack< FEMSigs ... >() );
}

#ifndef FAST_COMPILE
template< unsigned int Dim , class Real >
void Execute( int argc , char* argv[] )
{
	switch( BType.value )
	{
	case BOUNDARY_FREE+1:
	{
		switch( Degree.value )
		{
//			case 1: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BOUNDARY_FREE >::Signature >() );
			case 2: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >() );
			case 3: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >() );
//			case 4: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 1 - 3 are supported" );
		}
	}
	case BOUNDARY_NEUMANN+1:
	{
		switch( Degree.value )
		{
//			case 1: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BOUNDARY_NEUMANN >::Signature >() );
			case 2: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_NEUMANN >::Signature >() );
			case 3: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_NEUMANN >::Signature >() );
//			case 4: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 1 - 3 are supported" );
		}
	}
	case BOUNDARY_DIRICHLET+1:
	{
		switch( Degree.value )
		{
//			case 1: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BOUNDARY_DIRICHLET >::Signature >() );
			case 2: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BOUNDARY_DIRICHLET >::Signature >() );
			case 3: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >() );
//			case 4: return Execute< Real >( argc , argv , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BOUNDARY_DIRICHLET >::Signature >() );
			default: ERROR_OUT( "Only B-Splines of degree 1 - 3 are supported" );
		}
	}
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
#if 0
	if( !In.set || !Out.set ) ERROR_OUT( "Need input and output" );
	unsigned int width , height;
	unsigned char *pixels = ImageReader::ReadColor( In.value , width , height );
	FILE *fp = fopen( Out.value , "wb" );
	if( !fp ) ERROR_OUT( "Failed to open file for reading: %s" , Out.value );
	for( int i=0 ; i<10000 ; i++ )
	{
		int x = rand() % width , y = rand() % height;
		double gray = (double)( pixels[ 3*(y*width+x) + 0 ] + pixels[ 3*(y*width+x) + 1 ] + pixels[ 3*(y*width+x) + 2 ] ) / ( 255. * 3 );
		fprintf( fp , "%d %d  %f\n" , x , y , gray );
	}
	fclose( fp );
#else

	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}
	if( NoValueConstraints.set ) ValueWeight.value = 0;
	if( !UseGradientConstraints.set ) GradientWeight.value = 0;

	if( ValueWeight.value<0 ) ERROR_OUT( "Value weight must be non-negative: " , ValueWeight.value , "> 0" );
	if( GradientWeight.value<0 ) ERROR_OUT( "Gradient weight must be non-negative: " , GradientWeight.value , "> 0" );
	if( !ValueWeight.value && !GradientWeight.value ) ERROR_OUT( "Either value or gradient weight must be positive" );

	if( LapWeight.value<0 ) ERROR_OUT( "Laplacian weight must be non-negative: " , LapWeight.value , " > 0" );
	if( BiLapWeight.value<0 ) ERROR_OUT( "Bi-Laplacian weight must be non-negative: " , BiLapWeight.value , " > 0" );
	if( !LapWeight.value && !BiLapWeight.value ) ERROR_OUT( "Eiter Laplacian or bi-Laplacian weight must be positive" );

	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

#ifdef FAST_COMPILE
	static const int Dimension = DIMENSION;
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
	typedef IsotropicUIntPack< Dimension , FEMDegreeAndBType< Degree , BType >::Signature > FEMSigs;
	WARN( "Compiled for degree-" , Degree , ", boundary-" , BoundaryNames[ BType ] , ", " , sizeof(Real)==4 ? "single" : "double" , "-precision _only_" );
	Execute< Real >( argc , argv , FEMSigs() );
#else // !FAST_COMPILE
	if( Dimension.value==2 ) Execute< 2 , Real >( argc , argv );
	else if( Dimension.value==3 ) Execute< 3 , Real >( argc , argv );
	else ERROR_OUT( "Only Degrees 2 and 3 are supported" );
#endif // FAST_COMPILE
	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}
#endif
	return EXIT_SUCCESS;
}
