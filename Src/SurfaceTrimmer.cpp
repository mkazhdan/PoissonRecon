/*
Copyright (c) 2013, Michael Kazhdan
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

#define DEFAULT_DIMENSION 3

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include "FEMTree.h"
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "MAT.h"
#include "Geometry.h"
#include "Ply.h"
#include "VertexFactory.h"

MessageWriter messageWriter;


cmdLineParameter< char* >
	In( "in" ) ,
	Out( "out" );
cmdLineParameter< int >
	Smooth( "smooth" , 5 );
cmdLineParameter< float >
	Trim( "trim" ) ,
	IslandAreaRatio( "aRatio" , 0.001f );
cmdLineReadable
	PolygonMesh( "polygonMesh" ) ,
	Long( "long" ) ,
	ASCII( "ascii" ) ,
	Verbose( "verbose" );


cmdLineReadable* params[] =
{
	&In , &Out , &Trim , &PolygonMesh , &Smooth , &IslandAreaRatio , &Verbose , &Long , &ASCII ,
	NULL
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t --%s <trimming value>\n" , Trim.name );
	printf( "\t[--%s <ouput polygon mesh>]\n" , Out.name );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Smooth.name , Smooth.value );
	printf( "\t[--%s <relative area of islands>=%f]\n" , IslandAreaRatio.name , IslandAreaRatio.value );
	printf( "\t[--%s]\n" , PolygonMesh.name );
	printf( "\t[--%s]\n" , Long.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

template< typename Real , unsigned int Dim , typename ... AuxData >
using ValuedPointData = VectorTypeUnion< Real , Point< Real , Dim > , Real , AuxData ... >;

template< typename Index >
struct EdgeKey
{
	Index key1 , key2;
	EdgeKey( Index k1=0 , Index k2=0 ) : key1(k1) , key2(k2) {}
	bool operator == ( const EdgeKey &key ) const  { return key1==key.key1 && key2==key.key2; }
#if 1
	struct Hasher{ size_t operator()( const EdgeKey &key ) const { return (size_t)( key.key1 * key.key2 ); } };
#else
	struct Hasher{ size_t operator()( const EdgeKey &key ) const { return key.key1 ^ key.key2; } };
#endif
};

template< typename Real , unsigned int Dim ,  typename ... AuxData >
ValuedPointData< Real , Dim , AuxData ... > InterpolateVertices( const ValuedPointData< Real , Dim , AuxData ... >& v1 , const ValuedPointData< Real , Dim , AuxData ... >& v2 , Real value )
{
	if( v1.template get<1>()==v2.template get<1>() ) return (v1+v2)/Real(2.);
	Real dx = ( v1.template get<1>()-value ) / ( v1.template get<1>()-v2.template get<1>() );
	return v1 * (Real)(1.-dx) + v2*dx;
}

template< typename Real , unsigned int Dim , typename Index , typename ... AuxData >
void SmoothValues( std::vector< ValuedPointData< Real , Dim , AuxData ... > >& vertices , const std::vector< std::vector< Index > >& polygons )
{
	std::vector< int > count( vertices.size() );
	std::vector< Real > sums( vertices.size() , 0 );
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		int sz = int(polygons[i].size());
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
			Index v1 = polygons[i][j1] , v2 = polygons[i][j2];
			count[v1]++ , count[v2]++;
			sums[v1] += vertices[v2].template get<1>() , sums[v2] += vertices[v1].template get<1>();
		}
	}
	for( size_t i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() = ( sums[i] + vertices[i].template get<1>() ) / ( count[i] + 1 );
}

template< typename Real , unsigned int Dim , typename Index , typename ... AuxData >
void SplitPolygon
(
	const std::vector< Index >& polygon ,
	std::vector< ValuedPointData< Real , Dim , AuxData ... > >& vertices ,
	std::vector< std::vector< Index > >* ltPolygons , std::vector< std::vector< Index > >* gtPolygons ,
	std::vector< bool >* ltFlags , std::vector< bool >* gtFlags ,
	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >& vertexTable,
	Real trimValue
)
{
	int sz = int( polygon.size() );
	std::vector< bool > gt( sz );
	int gtCount = 0;
	for( int j=0 ; j<sz ; j++ )
	{
		gt[j] = ( vertices[ polygon[j] ].template get<1>()>trimValue );
		if( gt[j] ) gtCount++;
	}
	if     ( gtCount==sz ){ if( gtPolygons ) gtPolygons->push_back( polygon ) ; if( gtFlags ) gtFlags->push_back( false ); }
	else if( gtCount==0  ){ if( ltPolygons ) ltPolygons->push_back( polygon ) ; if( ltFlags ) ltFlags->push_back( false ); }
	else
	{
		int start;
		for( start=0 ; start<sz ; start++ ) if( gt[start] && !gt[(start+sz-1)%sz] ) break;

		bool gtFlag = true;
		std::vector< Index > poly;

		// Add the initial vertex
		{
			int j1 = (start+int(sz)-1)%sz , j2 = start;
			Index v1 = polygon[j1] , v2 = polygon[j2] , vIdx;
			typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = vertexTable.find( EdgeKey< Index >(v1,v2) );
			if( iter==vertexTable.end() )
			{
				vertexTable[ EdgeKey< Index >(v1,v2) ] = vIdx = (Index)vertices.size();
				vertices.push_back( InterpolateVertices( vertices[v1] , vertices[v2] , trimValue ) );
			}
			else vIdx = iter->second;
			poly.push_back( vIdx );
		}

		for( int _j=0  ; _j<=sz ; _j++ )
		{
			int j1 = (_j+start+sz-1)%sz , j2 = (_j+start)%sz;
			Index v1 = polygon[j1] , v2 = polygon[j2];
			if( gt[j2]==gtFlag ) poly.push_back( v2 );
			else
			{
				Index vIdx;
				typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = vertexTable.find( EdgeKey< Index >(v1,v2) );
				if( iter==vertexTable.end() )
				{
					vertexTable[ EdgeKey< Index >(v1,v2) ] = vIdx = (Index)vertices.size();
					vertices.push_back( InterpolateVertices( vertices[v1] , vertices[v2] , trimValue ) );
				}
				else vIdx = iter->second;
				poly.push_back( vIdx );
				if( gtFlag ){ if( gtPolygons ) gtPolygons->push_back( poly ) ; if( ltFlags ) ltFlags->push_back( true ); }
				else        { if( ltPolygons ) ltPolygons->push_back( poly ) ; if( gtFlags ) gtFlags->push_back( true ); }
				poly.clear() , poly.push_back( vIdx ) , poly.push_back( v2 );
				gtFlag = !gtFlag;
			}
		}
	}
}

template< class Real , unsigned int Dim , typename Index , class Vertex >
void Triangulate( const std::vector< Vertex >& vertices , const std::vector< std::vector< Index > >& polygons , std::vector< std::vector< Index > >& triangles )
{
	triangles.clear();
	for( size_t i=0 ; i<polygons.size() ; i++ )
		if( polygons[i].size()>3 )
		{
			std::vector< Point< Real , Dim > > _vertices( polygons[i].size() );
			for( int j=0 ; j<int( polygons[i].size() ) ; j++ ) _vertices[j] = vertices[ polygons[i][j] ].template get<0>();
			std::vector< TriangleIndex< Index > > _triangles = MinimalAreaTriangulation< Index , Real , Dim >( ( ConstPointer( Point< Real , Dim > ) )GetPointer( _vertices ) , _vertices.size() );

			// Add the triangles to the mesh
			size_t idx = triangles.size();
			triangles.resize( idx+_triangles.size() );
			for( int j=0 ; j<int(_triangles.size()) ; j++ )
			{
				triangles[idx+j].resize(3);
				for( int k=0 ; k<3 ; k++ ) triangles[idx+j][k] = polygons[i][ _triangles[j].idx[k] ];
			}
		}
		else if( polygons[i].size()==3 ) triangles.push_back( polygons[i] );
}

template< class Real , unsigned int Dim , typename Index , class Vertex >
double PolygonArea( const std::vector< Vertex >& vertices , const std::vector< Index >& polygon )
{
	auto Area =[]( Point< Real , Dim > v1 , Point< Real , Dim > v2 , Point< Real , Dim > v3 )
	{
		Point< Real , Dim > v[] = { v2-v1 , v3-v1 };
		XForm< Real , 2 > Mass;
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) Mass(i,j) = Point< Real , Dim >::Dot( v[i] , v[j] );
		return (Real)( sqrt( Mass.determinant() ) / 2. );
	};

	if( polygon.size()<3 ) return 0.;
	else if( polygon.size()==3 ) return Area( vertices[polygon[0]].template get<0>() , vertices[polygon[1]].template get<0>() , vertices[polygon[2]].template get<0>() );
	else
	{
		Point< Real , DEFAULT_DIMENSION > center;
		for( size_t i=0 ; i<polygon.size() ; i++ ) center += vertices[ polygon[i] ].template get<0>();
		center /= Real( polygon.size() );
		double area = 0;
		for( size_t i=0 ; i<polygon.size() ; i++ ) area += Area( center , vertices[ polygon[i] ].template get<0>() , vertices[ polygon[ (i+1)%polygon.size() ] ].template get<0>() );
		return area;
	}
}

template< typename Index , class Vertex >
void RemoveHangingVertices( std::vector< Vertex >& vertices , std::vector< std::vector< Index > >& polygons )
{
	std::unordered_map< Index, Index > vMap;
	std::vector< bool > vertexFlags( vertices.size() , false );
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) vertexFlags[ polygons[i][j] ] = true;
	Index vCount = 0;
	for( Index i=0 ; i<(Index)vertices.size() ; i++ ) if( vertexFlags[i] ) vMap[i] = vCount++;
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) polygons[i][j] = vMap[ polygons[i][j] ];

	std::vector< Vertex > _vertices( vCount );
	for( Index i=0 ; i<(Index)vertices.size() ; i++ ) if( vertexFlags[i] ) _vertices[ vMap[i] ] = vertices[i];
	vertices = _vertices;
}

template< typename Index >
void SetConnectedComponents( const std::vector< std::vector< Index > >& polygons , std::vector< std::vector< Index > >& components )
{
	std::vector< Index > polygonRoots( polygons.size() );
	for( size_t i=0 ; i<polygons.size() ; i++ ) polygonRoots[i] = (Index)i;
	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher > edgeTable;
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		int sz = int( polygons[i].size() );
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
			Index v1 = polygons[i][j1] , v2 = polygons[i][j2];
			EdgeKey< Index > eKey = EdgeKey< Index >(v1,v2);
			typename std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher >::iterator iter = edgeTable.find(eKey);
			if( iter==edgeTable.end() ) edgeTable[ eKey ] = (Index)i;
			else
			{
				Index p = iter->second;
				while( polygonRoots[p]!=p )
				{
					Index temp = polygonRoots[p];
					polygonRoots[p] = (Index)i;
					p = temp;
				}
				polygonRoots[p] = (Index)i;
			}
		}
	}
	for( size_t i=0 ; i<polygonRoots.size() ; i++ )
	{
		Index p = (Index)i;
		while( polygonRoots[p]!=p ) p = polygonRoots[p];
		Index root = p;
		p = (Index)i;
		while( polygonRoots[p]!=p )
		{
			Index temp = polygonRoots[p];
			polygonRoots[p] = root;
			p = temp;
		}
	}
	int cCount = 0;
	std::unordered_map< Index , Index > vMap;
	for( Index i=0 ; i<(Index)polygonRoots.size() ; i++ ) if( polygonRoots[i]==i ) vMap[i] = cCount++;
	components.resize( cCount );
	for( Index i=0 ; i<(Index)polygonRoots.size() ; i++ ) components[ vMap[ polygonRoots[i] ] ].push_back(i);
}

template< typename Real , unsigned int Dim , typename Index , typename ... AuxDataFactories >
int Execute( AuxDataFactories ... auxDataFactories )
{
	typedef VertexFactory::Factory< Real , typename VertexFactory::PositionFactory< Real , Dim > , typename VertexFactory::ValueFactory< Real > , AuxDataFactories ... > Factory;
	typedef typename Factory::VertexType Vertex;
	typename VertexFactory::PositionFactory< Real , Dim > pFactory;
	typename VertexFactory::ValueFactory< Real > vFactory;
	Factory factory( pFactory , vFactory , auxDataFactories ... );
	Real min , max;

	std::vector< Vertex > vertices;
	std::vector< std::vector< Index > > polygons;

	int ft;
	std::vector< std::string > comments;
	PLY::ReadPolygons< Factory , Index >( In.value , factory , vertices , polygons , ft , comments );

	for( int i=0 ; i<Smooth.value ; i++ ) SmoothValues< Real , Dim , Index >( vertices , polygons );
	min = max = vertices[0].template get<1>();
	for( size_t i=0 ; i<vertices.size() ; i++ ) min = std::min< Real >( min , vertices[i].template get<1>() ) , max = std::max< Real >( max , vertices[i].template get<1>() );

	std::unordered_map< EdgeKey< Index > , Index , typename EdgeKey< Index >::Hasher > vertexTable;
	std::vector< std::vector< Index > > ltPolygons , gtPolygons;
	std::vector< bool > ltFlags , gtFlags;

	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "** Running Surface Trimmer (Version %s) **\n" , VERSION );
	messageWriter( comments , "*********************************************\n" );
	messageWriter( comments , "*********************************************\n" );
	char str[1024];
	for( int i=0 ; params[i] ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( str );
			if( strlen( str ) ) messageWriter( comments , "\t--%s %s\n" , params[i]->name , str );
			else                messageWriter( comments , "\t--%s\n" , params[i]->name );
		}
	if( Verbose.set ) printf( "Value Range: [%f,%f]\n" , min , max );

	double t=Time();
	for( size_t i=0 ; i<polygons.size() ; i++ ) SplitPolygon( polygons[i] , vertices , &ltPolygons , &gtPolygons , &ltFlags , &gtFlags , vertexTable , Trim.value );
	if( IslandAreaRatio.value>0 )
	{
		std::vector< std::vector< Index > > _ltPolygons , _gtPolygons;
		std::vector< std::vector< Index > > ltComponents , gtComponents;
		SetConnectedComponents( ltPolygons , ltComponents );
		SetConnectedComponents( gtPolygons , gtComponents );
		std::vector< double > ltAreas( ltComponents.size() , 0. ) , gtAreas( gtComponents.size() , 0. );
		std::vector< bool > ltComponentFlags( ltComponents.size() , false ) , gtComponentFlags( gtComponents.size() , false );
		double area = 0.;
		for( size_t i=0 ; i<ltComponents.size() ; i++ )
		{
			for( size_t j=0 ; j<ltComponents[i].size() ; j++ )
			{
				ltAreas[i] += PolygonArea< Real , Dim , Index , Vertex >( vertices , ltPolygons[ ltComponents[i][j] ] );
				ltComponentFlags[i] = ( ltComponentFlags[i] || ltFlags[ ltComponents[i][j] ] );
			}
			area += ltAreas[i];
		}
		for( size_t i=0 ; i<gtComponents.size() ; i++ )
		{
			for( size_t j=0 ; j<gtComponents[i].size() ; j++ )
			{
				gtAreas[i] += PolygonArea< Real , Dim , Index , Vertex >( vertices , gtPolygons[ gtComponents[i][j] ] );
				gtComponentFlags[i] = ( gtComponentFlags[i] || gtFlags[ gtComponents[i][j] ] );
			}
			area += gtAreas[i];
		}
		for( size_t i=0 ; i<ltComponents.size() ; i++ )
		{
			if( ltAreas[i]<area*IslandAreaRatio.value && ltComponentFlags[i] ) for( size_t j=0 ; j<ltComponents[i].size() ; j++ ) _gtPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
			else                                                               for( size_t j=0 ; j<ltComponents[i].size() ; j++ ) _ltPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
		}
		for( size_t i=0 ; i<gtComponents.size() ; i++ )
		{
			if( gtAreas[i]<area*IslandAreaRatio.value && gtComponentFlags[i] ) for( size_t j=0 ; j<gtComponents[i].size() ; j++ ) _ltPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
			else                                                               for( size_t j=0 ; j<gtComponents[i].size() ; j++ ) _gtPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
		}
		ltPolygons = _ltPolygons , gtPolygons = _gtPolygons;
	}
	if( !PolygonMesh.set )
	{
		{
			std::vector< std::vector< Index > > polys = ltPolygons;
			Triangulate< Real , Dim , Index , Vertex >( vertices , ltPolygons , polys ) , ltPolygons = polys;
		}
		{
			std::vector< std::vector< Index > > polys = gtPolygons;
			Triangulate< Real , Dim , Index , Vertex >( vertices , gtPolygons , polys ) , gtPolygons = polys;
		}
	}

	RemoveHangingVertices( vertices , gtPolygons );
	char comment[1024];
	sprintf( comment , "#Trimmed In: %9.1f (s)" , Time()-t );
	comments.push_back( comment );
	if( Out.set ) PLY::WritePolygons( Out.value , factory , vertices , gtPolygons , ASCII.set ? PLY_ASCII : ft , comments );

	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , &argv[1] , params );
	messageWriter.echoSTDOUT = Verbose.set;

	if( !In.set || !Trim.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	typedef float Real;
	static constexpr unsigned int Dim = DEFAULT_DIMENSION;
	typedef VertexFactory::Factory< Real , typename VertexFactory::PositionFactory< Real , Dim > , typename VertexFactory::ValueFactory< Real > > Factory;
	Factory factory;
	bool *readFlags = new bool[ factory.plyReadNum() ];
	std::vector< PlyProperty > unprocessedProperties;
	PLY::ReadVertexHeader( In.value , factory , readFlags , unprocessedProperties );
	if( !factory.plyValidReadProperties<0>( readFlags ) ) ERROR_OUT( "Ply file does not contain positions" );
	if( !factory.plyValidReadProperties<1>( readFlags ) ) ERROR_OUT( "Ply file does not contain values" );
	delete[] readFlags;

	if( Long.set ) return Execute< Real , Dim , long long >( VertexFactory::DynamicFactory< Real >( unprocessedProperties ) );
	else           return Execute< Real , Dim , int       >( VertexFactory::DynamicFactory< Real >( unprocessedProperties ) );

}
