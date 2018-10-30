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

#undef ARRAY_DEBUG

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

cmdLineParameter< char* >
	In( "in" ) ,
	OutMesh( "mesh" ) ,
	OutGrid( "grid" );

cmdLineReadable
	PolygonMesh( "polygonMesh" ) ,
	NonManifold( "nonManifold" ) ,
	FlipOrientation( "flip" ) ,
	ASCII( "ascii" ) ,
	NonLinearFit( "nonLinearFit" ) ,
	PrimalGrid( "primalGrid" ) ,
	Verbose( "verbose" );

cmdLineParameter< int >
	Threads( "threads" , omp_get_num_procs() );

cmdLineParameter< float >
	IsoValue( "iso" , 0.f );

cmdLineReadable* params[] =
{
	&In , 
	&OutMesh , &NonManifold , &PolygonMesh , &FlipOrientation , &ASCII , &NonLinearFit , &IsoValue ,
	&OutGrid , &PrimalGrid ,
	&Threads ,
	&Verbose , 
	NULL
};


void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input tree>\n" , In.name );
	printf( "\t[--%s <ouput triangle mesh>]\n" , OutMesh.name );
	printf( "\t[--%s <ouput grid>]\n" , OutGrid.name );
#ifdef _OPENMP
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
#endif // _OPENMP
	printf( "\t[--%s <iso-value for extraction>=%f]\n" , IsoValue.name , IsoValue.value );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s]\n" , NonLinearFit.name );
	printf( "\t[--%s]\n" , FlipOrientation.name );
	printf( "\t[--%s]\n" , PrimalGrid.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< unsigned int Dim , class Real , unsigned int FEMSig >
void _Execute( const FEMTree< Dim , Real >* tree , FILE* fp )
{
	static const unsigned int Degree = FEMSignature< FEMSig >::Degree;
	DenseNodeData< Real , IsotropicUIntPack< Dim , FEMSig > > coefficients;

	coefficients.read( fp );

	// Output the grid
	if( OutGrid.set )
	{
		FILE* _fp = fopen( OutGrid.value , "wb" );
		if( !_fp ) fprintf( stderr , "Failed to open grid file for writing: %s\n" , OutGrid.value );
		else
		{
			int res = 0;
			double t = Time();
			Pointer( Real ) values = tree->template regularGridEvaluate< true >( coefficients , res , -1 , PrimalGrid.set );
			if( Verbose.set ) printf( "Got grid: %.2f(s)\n" , Time()-t );

			int cells = 1;
			for( int d=0 ; d<Dim ; d++ ) cells *= res;
			Pointer( float ) fValues = AllocPointer< float >( cells );

			Real min , max;
			min = max = values[0];
			for( int i=0 ; i<cells ; i++ ) min = std::min< Real >( min , values[i] ) , max = std::max< Real >( max , values[i] );

			int strides[Dim] , _strides[Dim] ; strides[0] = _strides[0] = 1;
			int idx[Dim+1] , _idx[Dim+1] ; idx[0] = _idx[0] = 0;
			for( int d=1 ; d<Dim ; d++ ) strides[d] = strides[d-1] * res , _strides[d] = _strides[d-1] * res;
			WindowLoop< Dim >::Run
			(
				0 , res ,
				[&]( int d , int i ){ _idx[d+1] = _idx[d] + _strides[d] * i , idx[d+1] = idx[d] + strides[d] * i; } ,
				[&]( void ){ fValues[ _idx[Dim] ] = (float)values[ idx[Dim] ]; }
			);
			DeletePointer( values );

			fwrite( &res , sizeof(int) , 1 , _fp );
			fwrite( fValues , sizeof(float) , cells , _fp );
			fclose( _fp );
			FreePointer( fValues );
		}
	}

	// Output the mesh
	if( OutMesh.set )
	{
		double t = Time();
		typedef PlyVertex< Real , Dim > Vertex;
		CoredFileMeshData< Vertex > mesh;
		std::function< void ( Vertex& , Point< Real , Dim > , Real , Real ) > SetVertex = []( Vertex& v , Point< Real , Dim > p , Real , Real ){ v.point = p; };
#if defined( __GNUC__ ) && __GNUC__ < 5
		#warning "you've got me gcc version<5"
			static const unsigned int DataSig = FEMDegreeAndBType< 0 , BOUNDARY_FREE >::Signature;
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , ( SparseNodeData< ProjectiveData< Real , Real > , IsotropicUIntPack< Dim , DataSig > > * )NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#else // !__GNUC__ || __GNUC__ >=5
		IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< Real >( IsotropicUIntPack< Dim , FEMSig >() , UIntPack< 0 >() , UIntPack< FEMTrivialSignature >() , *tree , ( typename FEMTree< Dim , Real >::template DensityEstimator< 0 >* )NULL , NULL , coefficients , IsoValue.value , mesh , SetVertex , NonLinearFit.set , !NonManifold.set , PolygonMesh.set , FlipOrientation.set );
#endif // __GNUC__ || __GNUC__ < 4

		if( Verbose.set ) printf( "Got iso-surface: %.2f(s)\n" , Time()-t );
		if( Verbose.set ) printf( "Vertices / Polygons: %d / %d\n" , (int)( mesh.outOfCorePointCount()+mesh.inCorePoints.size() ) , (int)mesh.polygonCount() );

		std::vector< std::string > comments;
		if( !PlyWritePolygons< Vertex , Real , Dim >( OutMesh.value , &mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , comments , XForm< Real , Dim+1 >::Identity() ) )
			fprintf( stderr , "[ERROR] Could not write mesh to: %s\n" , OutMesh.value ) , exit( 0 );
	}
}


template< unsigned int Dim , class Real >
int Execute( FILE* fp , int degree , BoundaryType bType )
{
	FEMTree< Dim , Real > tree( fp , MEMORY_ALLOCATOR_BLOCK_SIZE );

	if( Verbose.set ) printf( "Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n" , (int)tree.leaves() , (int)tree.nodes() , (int)tree.ghostNodes() );

	switch( bType )
	{
	case BOUNDARY_FREE:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_FREE >::Signature >( &tree , fp ) ; break;
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) , exit( 0 );
		}
	}
	break;
	case BOUNDARY_NEUMANN:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN >::Signature >( &tree , fp ) ; break;
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) , exit( 0 );
		}
	}
	break;
	case BOUNDARY_DIRICHLET:
	{
		switch( degree )
		{
			case 1: _Execute< Dim , Real , FEMDegreeAndBType< 1 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 2: _Execute< Dim , Real , FEMDegreeAndBType< 2 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 3: _Execute< Dim , Real , FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			case 4: _Execute< Dim , Real , FEMDegreeAndBType< 4 , BOUNDARY_DIRICHLET >::Signature >( &tree , fp ) ; break;
			default: fprintf( stderr , "[ERROR] Only B-Splines of degree 1 - 4 are supported" ) , exit( 0 );
		}
	}
	break;
	default: fprintf( stderr , "[ERROR] Not a valid boundary type: %d\n" , bType ) , exit( 0 );
	}
	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
#ifdef ARRAY_DEBUG
	fprintf( stderr , "[WARNING] Array debugging enabled\n" );
#endif // ARRAY_DEBUG
	cmdLineParse( argc-1 , &argv[1] , params );
	omp_set_num_threads( Threads.value > 1 ? Threads.value : 1 );
	if( Verbose.set )
	{
		printf( "**************************************************\n" );
		printf( "**************************************************\n" );
		printf( "** Running Octree Visualization (Version %s) **\n" , VERSION );
		printf( "**************************************************\n" );
		printf( "**************************************************\n" );
	}

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	FILE* fp = fopen( In.value , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , In.value ) , exit( 0 );
	FEMTreeRealType realType ; int degree ; BoundaryType bType;
	int dimension;
	ReadFEMTreeParameter( fp , realType , dimension );
	{
		unsigned int dim = dimension;
		unsigned int* sigs = ReadDenseNodeDataSignatures( fp , dim );
		if( dimension!=dim ) fprintf( stderr , "[ERROR] Octree and node data dimensions don't math: %d != %d\n" , dimension , dim ) , exit( 0 );
		for( unsigned int d=1 ; d<dim ; d++ ) if( sigs[0]!=sigs[d] )
		{
			fprintf( stderr , "[ERROR] Anisotropic signatures:\n" );
			for( unsigned int dd=0 ; dd<dim ; dd++ ) printf( "\t%d] %d %s\n" , dd , FEMSignatureDegree( sigs[dd] ) , BoundaryNames[ FEMSignatureBType( sigs[dd] ) ] );
			exit( 0 );
		}
		degree = FEMSignatureDegree( sigs[0] );
		bType = FEMSignatureBType( sigs[0] );
		delete[] sigs;
	}
	if( Verbose.set ) printf( "%d-dimension , %s-precision , degree-%d , %s-boundary\n" , dimension , FEMTreeRealNames[ realType ] , degree , BoundaryNames[ bType ] );

	switch( dimension )
	{
	case 2:
		switch( realType )
		{
			case FEM_TREE_REAL_FLOAT:  Execute< 2 , float  >( fp , degree , bType ) ; break;
			case FEM_TREE_REAL_DOUBLE: Execute< 2 , double >( fp , degree , bType ) ; break;
			default: fprintf( stderr , "[ERROR] Unrecognized real type: %d\n" , realType ) , exit( 0 );
		}
		break;
	case 3:
		switch( realType )
		{
			case FEM_TREE_REAL_FLOAT:  Execute< 3 , float  >( fp , degree , bType ) ; break;
			case FEM_TREE_REAL_DOUBLE: Execute< 3 , double >( fp , degree , bType ) ; break;
			default: fprintf( stderr , "[ERROR] Unrecognized real type: %d\n" , realType ) , exit( 0 );
		}
		break;
	default: fprintf( stderr , "[ERROR] Only dimensions 1-4 supported\n" ) , exit( 0 );
	}

	fclose( fp );
	return EXIT_SUCCESS;
}
