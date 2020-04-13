/*
Copyright (c) 2016, Michael Kazhdan
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

#include "Mesh/PoissonRecon/PreProcessor.h"

#undef USE_DOUBLE				// If enabled, double-precesion is used
#define DIMENSION 3				// The dimension of the system
#define DEFAULT_FEM_DEGREE 1	// The default finite-element degree

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <functional>
#include "Mesh/PoissonRecon/MyMiscellany.h"
#include "Mesh/PoissonRecon/CmdLineParser.h"
#include "Mesh/PoissonRecon/PPolynomial.h"
#include "Mesh/PoissonRecon/FEMTree.h"

MessageWriter messageWriter;

cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" ) ,
	InXForm( "inXForm" ) ,
	OutXForm( "outXForm" );

cmdLineReadable
	Performance( "performance" ) ,
	ShowResidual( "showResidual" ) ,
	ExactInterpolation( "exact" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
#ifndef FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
#endif // !FAST_COMPILE
	GSIterations( "iters" , 8 ) ,
	Depth( "depth" , 8 ) ,
	FullDepth( "fullDepth" , 5 ) ,
	BaseDepth( "baseDepth" , 0 ) ,
	BaseVCycles( "baseVCycles" , 1 ) ,
	MaxMemoryGB( "maxMemory" , 0 ) ,
#ifdef _OPENMP
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
#else // !_OPENMP
	ParallelType( "parallel" , (int)ThreadPool::THREAD_POOL ) ,
#endif // _OPENMP
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "chunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() );

cmdLineParameter< float >
	Scale( "scale" , 2.f ) ,
	CGSolverAccuracy( "cgAccuracy" , float(1e-3) ) ,
	DiffusionTime( "diffusion" , 0.0005f ) ,
	WeightScale( "wScl" , 0.125f ) ,
	WeightExponent( "wExp" , 6.f ) ,
	ValueWeight( "valueWeight" , 1e-2f );

cmdLineReadable* params[] =
{
#ifndef FAST_COMPILE
	&Degree ,
#endif // !FAST_COMPILE
	&In , &Out , &Depth , &InXForm , &OutXForm ,
	&Scale , &Verbose , &CGSolverAccuracy ,
	&ShowResidual ,
	&ValueWeight , &DiffusionTime ,
	&Threads ,
	&FullDepth ,
	&GSIterations ,
	&WeightScale , &WeightExponent ,
	&BaseDepth , &BaseVCycles ,
	&Performance ,
	&ExactInterpolation ,
	&MaxMemoryGB ,
	&ParallelType , 
	&ScheduleType , 
	&ThreadChunkSize ,
	NULL
};


void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output EDT solution>]\n" , Out.name );
#ifndef FAST_COMPILE
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
#endif // !FAST_COMPILE
	printf( "\t[--%s <maximum reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s <full depth>=%d]\n" , FullDepth.name , FullDepth.value );
	printf( "\t[--%s <coarse MG solver depth>=%d]\n" , BaseDepth.name , BaseDepth.value );
	printf( "\t[--%s <coarse MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );
	printf( "\t[--%s <diffusion time>=%.3e]\n" , DiffusionTime.name , DiffusionTime.value );
	printf( "\t[--%s <value interpolation weight>=%.3e]\n" , ValueWeight.name , ValueWeight.value );
	printf( "\t[--%s <iterations>=%d]\n" , GSIterations.name , GSIterations.value );
	printf( "\t[--%s]\n" , ExactInterpolation.name );
#ifdef _OPENMP
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
#endif // _OPENMP
	printf( "\t[--%s <cg solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <successive under-relaxation weight>=%f]\n" , WeightScale.name , WeightScale.value );
	printf( "\t[--%s <successive under-relaxation exponent>=%f]\n" , WeightExponent.name , WeightExponent.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s]\n" , Performance.name );
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
XForm< Real , Dim+1 > GetPointXForm( const std::vector< Point< Real , Dim > >& vertices , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	min = max = vertices[0];
	for( int i=0 ; i<vertices.size() ; i++ ) for( int j=0 ; j<Dim ; j++ ) min[j] = std::min< Real >( min[j] , vertices[i][j] ) , max[j] = std::max< Real >( max[j] , vertices[i][j] );
	Point< Real , Dim > center = ( max + min ) / 2;

	Real scale = max[0]-min[0];
	for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
	scale *= scaleFactor;
	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	return sXForm * tXForm;
}

template< class Real , unsigned int Dim >
void Print( const XForm< Real , Dim >& xForm )
{
	for( int j=0 ; j<Dim ; j++ )
	{
		for( int i=0 ; i<Dim ; i++ ) printf( " %f" , xForm(i,j) );
		printf( "\n" );
	}
}

template< unsigned int Dim , class Real >
struct ConstraintDual
{
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p ) const { return CumulativeDerivativeValues< Real , Dim , 0 >( ); }
};
template< unsigned int Dim , class Real >
struct SystemDual
{
	Real weight;
	SystemDual( Real w ) : weight(w){ }
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( Point< Real , Dim > p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues*weight; }
	CumulativeDerivativeValues< double , Dim , 0 > operator()( Point< Real , Dim > p , const CumulativeDerivativeValues< double , Dim , 0 >& dValues ) const { return dValues * weight; };
};
template< unsigned int Dim >
struct SystemDual< Dim , double >
{
	typedef double Real;
	Real weight;
	SystemDual( Real w ) : weight(w){ }
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( Point< Real , Dim > p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues*weight; }
};

template< unsigned int Dim , class Real , unsigned int FEMSig >
void _Execute( int argc , char* argv[] )
{
	ThreadPool::Init( (ThreadPool::ParallelType)ParallelType.value , Threads.value );
	static const unsigned int Degree = FEMSignature< FEMSig >::Degree;
	typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;
	std::vector< std::string > comments;
	messageWriter( comments , "*****************************************\n" );
	messageWriter( comments , "*****************************************\n" );
	messageWriter( comments , "** Running EDT in Heat (Version %s) **\n" , VERSION );
	messageWriter( comments , "*****************************************\n" );
	messageWriter( comments , "*****************************************\n" );
	if( !Threads.set ) messageWriter( comments , "Running with %d threads\n" , Threads.value );


	XForm< Real , Dim+1 > xForm , iXForm;
	if( InXForm.set )
	{
		FILE* fp = fopen( InXForm.value , "r" );
		if( !fp )
		{
			WARN( "Could not open file for reading x-form: " , InXForm.value );
			xForm = XForm< Real , Dim+1 >::Identity();
		}
		else
		{
			for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ )
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
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return;
	}
	
	std::vector< NodeAndPointSample< Dim , Real > > geometrySamples;
	std::vector< NodeAndPointSample< Dim , Real > > heatPositions;
	std::vector< Point< Real , Dim > > heatGradients;

	// Read the mesh into the tree
	{
		profiler.start();
		// Read the mesh
		std::vector< Point< Real , Dim > > vertices;
		std::vector< TriangleIndex< node_index_type > > triangles;
		{
			int file_type;
			std::vector< PlyVertex< float , Dim > > _vertices;
			std::vector< std::vector< int > > _polygons;
			std::vector< std::string > comments;
			if( !PlyReadPolygons( In.value , _vertices , _polygons , PlyVertex< float , Dim >::PlyReadProperties() , PlyVertex< float , Dim >::PlyReadNum , file_type , comments ) )
				ERROR_OUT( "Failed to read ply file: " , In.value );
			vertices.resize( _vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ ) for( int d=0 ; d<Dim ; d++ ) vertices[i][d] = _vertices[i].point[d];
			triangles.resize( _polygons.size() );
			for( int i=0 ; i<triangles.size() ; i++ ) for( int j=0 ; j<Dim ; j++ ) triangles[i][j] = _polygons[i][j];
		}
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = xForm * vertices[i];
		XForm< Real , Dim+1 > _xForm = GetPointXForm< Real , Dim >( vertices , (Real)Scale.value );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = _xForm * vertices[i];
		xForm = _xForm * xForm;
		FEMTreeInitializer< Dim , Real >::Initialize( tree.spaceRoot() , vertices , triangles , Depth.value , geometrySamples , true , tree.nodeAllocators , tree.initializer() );
		iXForm = xForm.inverse();
		if( OutXForm.set )
		{
			FILE* fp = fopen( OutXForm.value , "w" );
			if( !fp ) WARN( "Could not open file for writing x-form: %s" );
			else
			{
				for( int i=0 ; i<Dim+1 ; i++ )
				{
					for( int j=0 ; j<Dim+1 ; j++ ) fprintf( fp , " %f" , (float)iXForm(i,j) );
					fprintf( fp , "\n" );
				}
				fclose( fp );
			}
		}

		double area = 0;
		std::vector< double > areas( ThreadPool::NumThreads() , 0 );
		ThreadPool::Parallel_for( 0 , triangles.size() , [&]( unsigned int thread , size_t i )
		{
			Simplex< Real , Dim , Dim-1 > s;
			for( int k=0 ; k<Dim ; k++ ) for( int j=0 ; j<Dim ; j++ ) s[k][j] = vertices[ triangles[i][k] ][j];
			Real a2 = s.squareMeasure();
			if( a2>0 ) areas[thread] += sqrt(a2) / 2;
		}
		);
		for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) area += areas[t];
		messageWriter( "Input Vertices / Triangle / Samples / Area: %llu / %llu / %llu / %g\n" , (unsigned long long)vertices.size() , (unsigned long long)triangles.size() , (unsigned long long)geometrySamples.size() , area );
		profiler.dumpOutput2( comments , "# Read input into tree:" );
	}

	// Thicken the tree around the mesh
	{
		profiler.start();
		typename FEMTree< Dim , Real >::FEMTreeNode** nodes = new typename FEMTree< Dim , Real >::FEMTreeNode*[ geometrySamples.size() ];
		for( int i=0 ; i<geometrySamples.size() ; i++ ) nodes[i] = geometrySamples[i].node;
		tree.template thicken< Degree >( nodes , (int)geometrySamples.size() );
		profiler.dumpOutput2( comments , "#       Thickened tree:" );
		delete[] nodes;
	}

	// Finalize the topology of the tree
	{
		profiler.start();
		tree.template finalizeForMultigrid< Degree >( FullDepth.value , typename FEMTree< Dim , Real >::TrivialHasDataFunctor() );
		profiler.dumpOutput2( comments , "#       Finalized tree:" );
	}

	messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes: %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() );
	messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );

	SparseNodeData< Point< Real , Dim+1 > , IsotropicUIntPack< Dim , FEMTrivialSignature > > leafValues;
	const double GradientCutOff = 0;

	// Compute the heat solution
	DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > > heatSolution;
	DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > > constraints;

	// Add the FEM constraints
	{
		profiler.start();
		constraints = tree.initDenseNodeData( IsotropicUIntPack< Dim , FEMSig >() );
		DenseNodeData< Point< Real , 1 > , IsotropicUIntPack< Dim , FEMTrivialSignature > > _constraints( tree.nodesSize() );
		for( int i=0 ; i<geometrySamples.size() ; i++ ) _constraints[ geometrySamples[i].node ][0] = geometrySamples[i].sample.weight * ( 1<<(Depth.value*Dim) );
		typename FEMIntegrator::template ScalarConstraint< IsotropicUIntPack< Dim , FEMSig > , IsotropicUIntPack< Dim , 0 > , IsotropicUIntPack< Dim , FEMTrivialSignature > , IsotropicUIntPack< Dim , 0 > > F( {1.} );		tree.addFEMConstraints( F , _constraints , constraints , Depth.value );
		profiler.dumpOutput2( comments , "# Set heat constraints:" );
	}

	// Solve the linear system
	{
		profiler.start();
		typename FEMTree< Dim , Real >::SolverInfo sInfo;
		sInfo.cgDepth = 0 , sInfo.cascadic = false , sInfo.iters = GSIterations.value , sInfo.vCycles = 1 , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
		sInfo.useSupportWeights = true;
		sInfo.sorRestrictionFunction  = [&]( Real w , Real ){ return ( Real )( WeightScale.value * pow( w , WeightExponent.value ) ); };
		{
			typename FEMIntegrator::template System< IsotropicUIntPack< Dim , FEMSig > , IsotropicUIntPack< Dim , 1 > > F( { 1. , (double)DiffusionTime.value } );
			heatSolution = tree.solveSystem( IsotropicUIntPack< Dim , FEMSig >() , F , constraints , Depth.value , sInfo );
		}
		sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
		profiler.dumpOutput2( comments , "#   Heat system solved:" );
	}

	// Evaluate the gradients at the leaves
	{
		profiler.start();

		typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< IsotropicUIntPack< Dim , FEMSig > , 0 > evaluator( &tree , heatSolution );
		typedef typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template ConstNeighbors< IsotropicUIntPack< Dim , 3 > > OneRingNeighbors;
		typedef typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template ConstNeighborKey< IsotropicUIntPack< Dim , 1 > , IsotropicUIntPack< Dim , 1 > > OneRingNeighborKey;
		std::vector< OneRingNeighborKey > oneRingNeighborKeys( ThreadPool::NumThreads() );
		int treeDepth = tree.tree().maxDepth();
		for( int i=0 ; i<oneRingNeighborKeys.size() ; i++ ) oneRingNeighborKeys[i].set( treeDepth );
		DenseNodeData< Real , IsotropicUIntPack< Dim , FEMTrivialSignature > > leafCenterValues = tree.initDenseNodeData( IsotropicUIntPack< Dim , FEMTrivialSignature >() );

		ThreadPool::Parallel_for( tree.nodesBegin(0) , tree.nodesEnd(Depth.value) , [&]( unsigned int thread , size_t i )
		{
			if( tree.isValidSpaceNode( tree.node((node_index_type)i) ) )
			{
				Point< Real , Dim > center ; Real width;
				tree.centerAndWidth( (node_index_type)i , center , width );
				leafCenterValues[i] = evaluator.values( center , thread )[0];
			}
		}
		);

		auto CenterGradient = [&] ( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* leaf , int thread )
		{
			int d , off[Dim] ; Point< Real , Dim > p ; Real width , _width = (Real)1./(1<<Depth.value);
			tree.depthAndOffset( leaf , d , off ) , tree.centerAndWidth( leaf->nodeData.nodeIndex , p , width );
			int res = 1<<d , _res = 1<<Depth.value;
			Point< Real , Dim > g;
			unsigned int index1[Dim] , index2[Dim];
			for( int dd=0 ; dd<Dim ; dd++ ) index1[dd] = index2[dd] = 1;
			const OneRingNeighbors& neighbors = oneRingNeighborKeys[thread].getNeighbors( leaf );
			for( int c=0 ; c<Dim ; c++ )
			{
				Real value1 , value2;
				if( off[c]-1>=0  ) index1[c] = 0;
				if( off[c]+1<res ) index2[c] = 2;
				const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* node1 = neighbors.neighbors().data[ GetWindowIndex( IsotropicUIntPack< Dim , 3 >() , index1 ) ];
				const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* node2 = neighbors.neighbors().data[ GetWindowIndex( IsotropicUIntPack< Dim , 3 >() , index2 ) ];
				if( d==Depth.value && tree.isValidSpaceNode( node2 ) ) value2 = leafCenterValues[ node2->nodeData.nodeIndex ];
				else
				{
					Point< Real , Dim > delta;
					delta[c] = ( (int)index2[c]-1 ) * _width;
					value2 = evaluator.values( p+delta , thread )[0];
				}
				if( d==Depth.value && tree.isValidSpaceNode( node1 ) ) value1 = leafCenterValues[ node1->nodeData.nodeIndex ];
				else
				{
					Point< Real , Dim > delta;
					delta[c] = ( (int)index1[c]-1 ) * _width;
					value1 = evaluator.values( p+delta , thread )[0];
				}
				
				g[c] = ( value2 - value1 ) / ( (Real)( index2[c] - index1[c] ) );

				index1[c] = index2[c] = 1;
			}

			return g * _res;
		};

		for( node_index_type i=tree.nodesBegin(0) ; i<tree.nodesEnd(Depth.value) ; i++ ) if( tree.isValidSpaceNode( tree.node(i) ) && !tree.isValidSpaceNode( tree.node(i)->children ) )
		{
			RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* leaf = ( RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* )tree.node(i);
			leafValues[leaf] *= 0;
		}

		ThreadPool::Parallel_for( tree.nodesBegin(0) , tree.nodesEnd(Depth.value) , [&]( unsigned int thread , size_t i  )
		{
			if( tree.isValidSpaceNode( tree.node((node_index_type)i) ) && !tree.isValidSpaceNode( tree.node((node_index_type)i)->children ) )
			{
				RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* leaf = ( RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* )tree.node((node_index_type)i);
				Point< Real , Dim > g = CenterGradient( leaf , thread );
				Real len = (Real)Length( g );
				if( len>GradientCutOff ) g /= len;
				Point< Real , Dim+1 >* leafValue = leafValues(leaf);
				if( leafValue ) for( int d=0 ; d<Dim ; d++ ) (*leafValue)[d+1] = -g[d];
				else ERROR_OUT( "Leaf value doesn't exist" );
			}
		}
		);
		profiler.dumpOutput2( comments , "#  Evaluated gradients:" );
	}


	// Compute the EDT
	{
		// Evaluate the gradients at the center of the leaf nodes
		DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > > edtSolution , constraints;
		InterpolationInfo *valueInfo = NULL;

		// Add the FEM constraints
		{
			profiler.start();
			constraints = tree.initDenseNodeData( IsotropicUIntPack< Dim , FEMSig >() );
			typename FEMIntegrator::template Constraint< IsotropicUIntPack< Dim , FEMSig > , IsotropicUIntPack< Dim , 1 > , IsotropicUIntPack< Dim , FEMTrivialSignature > , IsotropicUIntPack< Dim , 0 > , Dim+1 > F;
			typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
			typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
			unsigned int derivatives2[Dim];
			for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
			for( int d=0 ; d<Dim ; d++ )
			{
				unsigned int derivatives1[Dim];
				for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ? 1 : 0;
				F.weights[d+1][TensorDerivatives< Derivatives1 >::Index( derivatives1 )][ TensorDerivatives< Derivatives2 >::Index( derivatives2 )] = 1.;
			}
			tree.addFEMConstraints( F , leafValues , constraints , Depth.value );
			profiler.dumpOutput2( comments , "#  Set EDT constraints:" );
		}

		// Add the interpolation constraints
		if( ValueWeight.value>0 )
		{
			profiler.start();
			if( ExactInterpolation.set ) valueInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 >( tree , geometrySamples , ConstraintDual< Dim , Real >() , SystemDual< Dim , Real >( std::max< Real >( 0 , (Real)ValueWeight.value ) ) , true , false );
			else                         valueInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 >( tree , geometrySamples , ConstraintDual< Dim , Real >() , SystemDual< Dim , Real >( std::max< Real >( 0 , (Real)ValueWeight.value ) ) , true , 0 );
			tree.addInterpolationConstraints( constraints , Depth.value , *valueInfo );
			profiler.dumpOutput2( comments , "#Set point constraints:" );
		}

		// Solve the linear system
		{
			profiler.start();
			typename FEMTree< Dim , Real >::SolverInfo sInfo;
			sInfo.cgDepth = 0 , sInfo.cascadic = true , sInfo.vCycles = 1 , sInfo.cgAccuracy = CGSolverAccuracy.value , sInfo.verbose = Verbose.set , sInfo.showResidual = ShowResidual.set , sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , sInfo.sliceBlockSize = 1;
			sInfo.iters = GSIterations.value;
			sInfo.baseDepth = BaseDepth.value , sInfo.baseVCycles = BaseVCycles.value;
			sInfo.useSupportWeights = true;
			sInfo.sorRestrictionFunction  = [&]( Real w , Real ){ return (Real)( WeightScale.value * pow( w , WeightExponent.value ) ); }; 
			typename FEMIntegrator::template System< IsotropicUIntPack< Dim , FEMSig > , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
			edtSolution = tree.solveSystem( IsotropicUIntPack< Dim , FEMSig >() , F , constraints , Depth.value , sInfo , valueInfo );
			profiler.dumpOutput2( comments , "#    EDT system solved:" );
		}
		if( valueInfo ) delete valueInfo , valueInfo = NULL;

		{
			auto GetAverageValueAndError = [&]( const FEMTree< Dim , Real >* tree , const DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > >& coefficients , double& average , double& error )
			{
				double errorSum = 0 , valueSum = 0 , weightSum = 0;
				typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< IsotropicUIntPack< Dim , FEMSig > , 0 > evaluator( tree , coefficients );
				std::vector< double > errorSums( ThreadPool::NumThreads() , 0 ) , valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
				ThreadPool::Parallel_for( 0 , geometrySamples.size() , [&]( unsigned int thread , size_t j )
				{
					ProjectiveData< Point< Real , Dim > , Real >& sample = geometrySamples[j].sample;
					Real w = sample.weight;
					Real value = evaluator.values( sample.data / sample.weight , thread , geometrySamples[j].node )[0];
					errorSums[thread] += value * value * w;
					valueSums[thread] += value * w;
					weightSums[thread] += w;
				}
				);
				for( unsigned int t=0 ; t<ThreadPool::NumThreads() ; t++ ) errorSum += errorSums[t] , valueSum += valueSums[t] , weightSum += weightSums[t];
				average = valueSum / weightSum , error = sqrt( errorSum / weightSum );
			};
			double average , error;
			GetAverageValueAndError( &tree , edtSolution , average , error );
			if( Verbose.set ) printf( "Interpolation average / error: %g / %g\n" , average , error );
			ThreadPool::Parallel_for( tree.nodesBegin(0) , tree.nodesEnd(0) , [&]( unsigned int , size_t i ){ edtSolution[i] -= (Real)average; } );
		}

		if( Out.set )
		{
			FILE* fp = fopen( Out.value , "wb" );
			if( !fp ) ERROR_OUT( "Failed to open file for writing: " , Out.value );
			FEMTree< Dim , Real >::WriteParameter( fp );
			DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > >::WriteSignatures( fp );
			tree.write( fp , xForm );
			edtSolution.write( fp );
			fclose( fp );
		}
	}
}

#ifndef FAST_COMPILE
template< unsigned int Dim , class Real >
void Execute( int argc , char* argv[] )
{
	switch( Degree.value )
	{
		case 1: return _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_FREE >::Signature >( argc , argv );
		case 2: return _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >( argc , argv );
		case 3: return _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >( argc , argv );
		case 4: return _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >( argc , argv );
		default: ERROR_OUT( "Only B-Splines of degree 1 - 4 are supported" );
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
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
	messageWriter.echoSTDOUT = Verbose.set;

#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE

#ifdef FAST_COMPILE
	static const int Degree = DEFAULT_FEM_DEGREE;
	static const BoundaryType BType = BOUNDARY_FREE;

	WARN( "Compiled for degree-" , Degree , ", boundary-" , BoundaryNames[ BType ] , ", " , sizeof(Real)==4 ? "single" : "double" , "-precision _only_" );
	if( BaseDepth.value>FullDepth.value )
	{
		if( BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , BaseDepth.value , " <= " , FullDepth.value );
		BaseDepth.value = FullDepth.value;
	}
	_Execute< DIMENSION , Real , FEMDegreeAndBType< Degree , BType >::Signature >( argc , argv );
#else // !FAST_COMPILE
	Execute< DIMENSION , Real >( argc , argv );
#endif // FAST_COMPILE
	if( Performance.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}
	return EXIT_SUCCESS;
}
