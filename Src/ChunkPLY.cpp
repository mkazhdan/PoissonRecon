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

#define NEW_CHUNKS
#define DISABLE_PARALLELIZATION

// -80,-130,40
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <sstream>
#include <unordered_map>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "Geometry.h"
#include "Ply.h"
#include "PointStream.h"
#include "PointStreamData.h"

cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameter< float > Width( "width" , -1.f ) , PadRadius( "radius" , 0.f );
cmdLineParameterArray< float , 6 > BoundingBox( "bBox" );
cmdLineReadable ASCII( "ascii" ) , Verbose( "verbose" );

cmdLineReadable* params[] = { &In , &Out , &Width , &PadRadius , &ASCII , &Verbose , &BoundingBox , NULL };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input ply file (vertices or polygons)>\n" , In.name );
	printf( "\t[--%s <output ply file name/header>]\n" , Out.name );
	printf( "\t[--%s <chunk width>=%f]\n" , Width.name , Width.value );
	printf( "\t[--%s <padding radius (as a fraction of the width)>=%f]\n" , PadRadius.name , PadRadius.value );
	printf( "\t[--%s <minx miny minz maxx maxy maxz>]\n" , BoundingBox.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

void PrintBoundingBox( Point< float , 3 > min , Point< float , 3 > max )
{
	printf( "[" );
	for( unsigned int d=0 ; d<3 ; d++ ) printf( " %f" , min[d] );
	printf( " ] [" );
	for( unsigned int d=0 ; d<3 ; d++ ) printf( " %f" , max[d] );
	printf( " ]" );
}

template< typename ... VertexData >
void WritePoints( const char *fileName , int ft , const std::vector< PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > > &vertices , const std::vector< std::string > &comments )
{
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;
	typedef MultiPointStreamData< float , VertexData ... > StreamData;

	PLYOutputPointStreamWithData< float , 3 , StreamData > pointStream( fileName , vertices.size() , ft , StreamData::PlyWriteProperties() , StreamData::PlyWriteNum );
	for( size_t i=0 ; i<vertices.size() ; i++ ) pointStream.nextPoint( vertices[i].point , vertices[i].data );
}


template< typename ... VertexData >
void WriteMesh( const char *fileName , int ft , const std::vector< PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > > &vertices , const std::vector< std::vector< long long > > &polygons , const std::vector< std::string > &comments )
{
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;

	if( vertices.size()>std::numeric_limits< int >::max() )
	{
		if( vertices.size()>std::numeric_limits< unsigned int >::max() ) ERROR_OUT( "more vertices than can be indexed by an unsigned int: %llu" , (unsigned long long)vertices.size() );
		WARN( "more vertices than can be indexed by an int, using unsigned int instead: %llu" , (unsigned long long)vertices.size() );
		std::vector< std::vector< unsigned int > > outPolygons;
		outPolygons.resize( polygons.size() );
		for( size_t i=0 ; i<polygons.size() ; i++ )
		{
			outPolygons[i].resize( polygons[i].size() );
			for( int j=0 ; j<polygons[i].size() ; j++ ) outPolygons[i][j] = (unsigned int)polygons[i][j];
		}
		if( !PlyWritePolygons< Vertex >( fileName , vertices , outPolygons , Vertex::PlyWriteProperties() , Vertex::PlyWriteNum , ft , comments ) ) ERROR_OUT( "Could not write mesh to: " , fileName );
	}
	else
	{
		std::vector< std::vector< int > > outPolygons;
		outPolygons.resize( polygons.size() );
		for( size_t i=0 ; i<polygons.size() ; i++ )
		{
			outPolygons[i].resize( polygons[i].size() );
			for( int j=0 ; j<polygons[i].size() ; j++ ) outPolygons[i][j] = (int)polygons[i][j];
		}
		if( !PlyWritePolygons< Vertex >( fileName , vertices , outPolygons , Vertex::PlyWriteProperties() , Vertex::PlyWriteNum , ft , comments ) ) ERROR_OUT( "Could not write mesh to: " , fileName );
	}
}

template< typename Vertex >
void GetBoundingBox( const std::vector< Vertex > &vertices , Point< float , 3 > &min , Point< float , 3 > &max )
{
	min = max = vertices[0].point;
	for( size_t i=0 ; i<vertices.size() ; i++ ) for( unsigned int d=0 ; d<3 ; d++ )
	{
		min[d] = std::min< float >( min[d] , vertices[i].point[d] );
		max[d] = std::max< float >( max[d] , vertices[i].point[d] );
	}
}

template< typename Vertex >
void GetSubPoints( const std::vector< Vertex > &vertices , Point< float , 3 > min , Point< float , 3 > max , std::vector< Vertex > &subVertices )
{
	subVertices.resize( 0 );

	for( size_t i=0 ; i<vertices.size() ; i++ )
	{
		bool inside = true;
		for( unsigned int d=0 ; d<3 ; d++ ) if( vertices[i].point[d]<min[d] || vertices[i].point[d]>=max[d] ) inside = false;
		if( inside ) subVertices.push_back( vertices[i] );
	}
}

template< typename Vertex >
void GetSubVertices( const std::vector< Vertex > &vertices , std::vector< std::vector< long long > > &polygons , std::vector< Vertex > &subVertices )
{
	subVertices.resize( 0 );
	long long count = 0;
	std::unordered_map< long long , long long > vMap;
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ )
	{
		auto iter = vMap.find( polygons[i][j] );
		if( iter==vMap.end() ) vMap[ polygons[i][j] ] = count++;
	}

	subVertices.resize( vMap.size() );

	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ )
	{
		long long oldIdx = polygons[i][j];
		long long newIdx = vMap[ oldIdx ];
		subVertices[ newIdx ] = vertices[ oldIdx ];
		polygons[i][j] = newIdx;
	}
}

template< typename Vertex >
void GetSubMesh( const std::vector< Vertex > &vertices , const std::vector< std::vector< long long > > &polygons , Point< float , 3 > min , Point< float , 3 > max , std::vector< Vertex > &subVertices , std::vector< std::vector< long long > > &subPolygons )
{
	subPolygons.resize( 0 );

	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		Point< float , 3 > center;
		for( size_t j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
		center /= (float)polygons[i].size();
		bool inside = true;
		for( unsigned int d=0 ; d<3 ; d++ ) if( center[d]<min[d] || center[d]>=max[d] ) inside = false;
		if( inside ) subPolygons.push_back( polygons[i] );
	}

	GetSubVertices( vertices , subPolygons , subVertices );
}

template< typename ... VertexData >
void Execute( void )
{
	Timer timer;
	double t = Time();
	typedef PlyVertexWithData< float , 3 , MultiPointStreamData< float , VertexData ... > > Vertex;
	std::vector< Vertex > vertices;
	std::vector< std::vector< long long > > polygons;

	int ft;
	std::vector< std::string > comments;
	PlyReadPolygons< Vertex >( In.value , vertices , polygons , Vertex::PlyReadProperties() , Vertex::PlyReadNum , ft , comments );
	printf( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)vertices.size() , (unsigned long long)polygons.size() );
	printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
	printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

	Point< float , 3 > min , max;
	GetBoundingBox( vertices , min , max );
	PrintBoundingBox( min , max ) ; printf( "\n" );

	float width = Width.value;


	if( BoundingBox.set )
	{
		Point< float , 3 > min( BoundingBox.values[0] , BoundingBox.values[1] , BoundingBox.values[2] );
		Point< float , 3 > max( BoundingBox.values[4] , BoundingBox.values[5] , BoundingBox.values[6] );
		auto InBoundingBox = [&]( Point< float , 3 > p )
		{
			return
				p[0]>=min[0] && p[0]<max[0] &&
				p[1]>=min[1] && p[1]<max[1] &&
				p[2]>=min[2] && p[2]<max[2];
		};
		if( polygons.size() )
		{
			std::vector< std::vector< long long > > _polygons;

			Timer timer;
#ifdef NEW_CHUNKS
			{
				size_t polygonCount = 0;
				for( size_t i=0 ; i<polygons.size() ; i++ )
				{
					Point< float , 3 > center;
					for( int j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
					center /= polygons[i].size();
					if( InBoundingBox( center ) ) polygonCount++;
				}
				_polygons.reserve( polygonCount );
			}
#endif // NEW_CHUNKS
			for( size_t i=0 ; i<polygons.size() ; i++ )
			{
				Point< float , 3 > center;
				for( int j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
				center /= polygons[i].size();
				if( InBoundingBox( center ) ) _polygons.push_back( polygons[i] );
			}
			printf( "\tChunked polygons:\n" );
			printf( "\t\tTime (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
			printf( "\t\tPeak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

			if( Out.set )
				if( _polygons.size() )
				{
					std::vector< Vertex > _vertices;
					GetSubVertices( vertices , _polygons , _vertices );

					if( Verbose.set )
					{
						printf( "\t\t%s\n" , Out.value );
						printf( "\t\t\tVertices / Polygons: %llu / %llu\n" , (unsigned long long)_vertices.size() , (unsigned long long)_polygons.size() );
					}

					WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , _vertices , _polygons , comments );
				}
				else WARN( "no polygons in bounding box" );
		}
		else
		{
			std::vector< Vertex > _vertices;

			Timer timer;
#ifdef NEW_CHUNKS
			{
				size_t vertexCount = 0;
				for( size_t i=0 ; i<vertices.size() ; i++ ) if( InBoundingBox( vertices[i].point ) ) vertexCount++;
				_vertices.reserve( vertexCount );
			}
#endif // NEW_CHUNKS

			for( size_t i=0 ; i<vertices.size() ; i++ ) if( InBoundingBox( vertices[i].point ) ) _vertices.push_back( vertices[i] );
			printf( "\tChunked vertices:\n" );
			printf( "\t\tTime (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
			printf( "\t\tPeak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

			if( Out.set )
				if( _vertices.size() )
				{
					if( Verbose.set )
					{
						printf( "\t\t%s\n" , Out.value );
						printf( "\t\t\tPoints: %llu\n" , (unsigned long long)_vertices.size() );
					}

					WritePoints( Out.value , ASCII.set ? PLY_ASCII : ft , _vertices , comments );
				}
				else WARN( "no vertices in bounding box" );
		}
	}
	else if( width>0 )
	{
		float radius = PadRadius.value * width;
		size_t vCount=0 , pCount=0;
		for( unsigned int d=0 ; d<3 ; d++ ) min[d] -= width/10000.f , max[d] += width/10000.f;
		int begin[] = { (int)floor( min[0]/width ) , (int)floor( min[1]/width ) , (int)floor( min[2]/width ) };
		int end  [] = { (int)ceil ( max[0]/width ) , (int)ceil ( max[1]/width ) , (int)ceil ( max[2]/width ) };
		int size [] = { end[0]-begin[0] , end[1]-begin[1] , end[2]-begin[2] };
		struct Range{ int begin[3] , end[3]; };
		auto SetRange = [&]( Point< float , 3 > p , Range &range )
		{
			for( int d=0 ; d<3 ; d++ )
			{
				range.begin[d] = (int)floor( (p[d]-radius)/width ) , range.end[d] = (int)ceil( (p[d]+radius)/width );
				if( range.begin[d]==range.end[d] ) range.end[d]++;
			}
		};
		auto Index1D = [&]( int x , int y , int z )
		{
			x -= begin[0] , y -= begin[1] , z -= begin[2];
			return x + y*size[0] + z*size[0]*size[1];
		};
		auto Index3D = [&]( int idx , int &x , int &y , int &z )
		{
			x = idx % size[0];
			idx /= size[0];
			y = idx % size[1];
			idx /= size[1];
			z = idx % size[2];
			x += begin[0] , y += begin[1] , z += begin[2];
		};

		if( polygons.size() )
		{
			std::vector< std::vector< std::vector< long long > > > _polygons( size[0]*size[1]*size[2] );
			Range range;

			Timer timer;
#ifdef NEW_CHUNKS
			{
				std::vector< size_t > polygonCounts( size[0]*size[1]*size[2] , 0 );
				for( size_t i=0 ; i<polygons.size() ; i++ )
				{
					Point< float , 3 > center;
					for( int j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
					center /= polygons[i].size();
					SetRange( center , range );
					for( int x=range.begin[0] ; x<range.end[0] ; x++ ) for( int y=range.begin[1] ; y<range.end[1] ; y++ ) for( int z=range.begin[2] ; z<range.end[2] ; z++ )
						polygonCounts[ Index1D(x,y,z) ]++;
				}
				for( size_t i=0 ; i<polygonCounts.size() ; i++ ) _polygons[i].reserve( polygonCounts[i] );
			}
#endif // NEW_CHUNKS
			for( size_t i=0 ; i<polygons.size() ; i++ )
			{
				Point< float , 3 > center;
				for( int j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
				center /= polygons[i].size();
				SetRange( center , range );
				for( int x=range.begin[0] ; x<range.end[0] ; x++ ) for( int y=range.begin[1] ; y<range.end[1] ; y++ ) for( int z=range.begin[2] ; z<range.end[2] ; z++ )
					_polygons[ Index1D(x,y,z) ].push_back( polygons[i] );
			}
			printf( "\tChunked polygons:\n" );
			printf( "\t\tTime (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
			printf( "\t\tPeak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

			if( Out.set )
			{
#ifdef DISABLE_PARALLELIZATION
#else // !DISABLE_PARALLELIZATION
#pragma omp parallel for
#endif // DISABLE_PARALLELIZATION
				for( int i=0 ; i<_polygons.size() ; i++ ) if( _polygons[i].size() )
				{
					std::vector< Vertex > _vertices;
					GetSubVertices( vertices , _polygons[i] , _vertices );
					std::stringstream stream;
					int x , y , z;
					Index3D( i , x , y , z );

					stream << Out.value << "." << x << "." << y << "." << z << ".ply";


					Point< float , 3 > min , max;
					min = Point< float , 3 >( x+0 , y+0 , z+0 ) * width;
					max = Point< float , 3 >( x+1 , y+1 , z+1 ) * width;
					if( Verbose.set )
					{
						static std::mutex mutex;
						std::lock_guard< std::mutex > lock( mutex );
						printf( "\t\t%s\n" , stream.str().c_str() );
						printf( "\t\t\tVertices / Polygons: %llu / %llu\n" , (unsigned long long)_vertices.size() , (unsigned long long)_polygons[i].size() );
						printf( "\t\t\t" ) ; PrintBoundingBox( min , max ) ; printf( "\n" );
					}

					WriteMesh( stream.str().c_str() , ASCII.set ? PLY_ASCII : ft , _vertices , _polygons[i] , comments );
					vCount += _vertices.size() , pCount += _polygons[i].size();
				}
			}
		}
		else
		{
			std::vector< std::vector< Vertex > > _vertices( size[0]*size[1]*size[2] );
			Range range;

			Timer timer;
#ifdef NEW_CHUNKS
			{
				std::vector< size_t > vertexCounts( size[0]*size[1]*size[2] , 0 );
				for( size_t i=0 ; i<vertices.size() ; i++ )
				{
					SetRange( vertices[i].point , range );
					for( int x=range.begin[0] ; x<range.end[0] ; x++ ) for( int y=range.begin[1] ; y<range.end[1] ; y++ ) for( int z=range.begin[2] ; z<range.end[2] ; z++ )
						vertexCounts[ Index1D(x,y,z) ]++;
				}
				for( size_t i=0 ; i<vertexCounts.size() ; i++ ) _vertices[i].reserve( vertexCounts[i] );
			}
#endif // NEW_CHUNKS

			for( size_t i=0 ; i<vertices.size() ; i++ )
			{
				SetRange( vertices[i].point , range );
				for( int x=range.begin[0] ; x<range.end[0] ; x++ ) for( int y=range.begin[1] ; y<range.end[1] ; y++ ) for( int z=range.begin[2] ; z<range.end[2] ; z++ )
					_vertices[ Index1D(x,y,z) ].push_back( vertices[i] );
			}
			printf( "\tChunked vertices:\n" );
			printf( "\t\tTime (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
			printf( "\t\tPeak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

			if( Out.set )
			{
#ifdef DISABLE_PARALLELIZATION
#else // !DISABLE_PARALLELIZATION
#pragma omp parallel for
#endif  // DISABLE_PARALLELIZATION
				for( int i=0 ; i<_vertices.size() ; i++ ) if( _vertices[i].size() )
				{
					std::stringstream stream;
					int x , y , z;
					Index3D( i , x , y , z );
					stream << Out.value << "." << x << "." << y << "." << z << ".ply";
					Point< float , 3 > min , max;
					min = Point< float , 3 >( x+0 , y+0 , z+0 ) * width;
					max = Point< float , 3 >( x+1 , y+1 , z+1 ) * width;

					if( Verbose.set )
					{
						static std::mutex mutex;
						std::lock_guard< std::mutex > lock( mutex );
						printf( "\t\t%s\n" , stream.str().c_str() );
						printf( "\t\t\tPoints: %llu\n" , (unsigned long long)_vertices[i].size() );
						printf( "\t\t\t" ) ; PrintBoundingBox( min , max ) ; printf( "\n" );
					}

					WritePoints( stream.str().c_str() , ASCII.set ? PLY_ASCII : ft , _vertices[i] , comments );
					vCount += _vertices[i].size();
				}
			}
		}
		if( !radius )
		{
			if( polygons.size() )
			{
				if( pCount!=polygons.size() ) WARN( "polygon counts don't match: " , polygons.size() , " != " , pCount );
			}
			else
			{
				if( vCount!=vertices.size() ) WARN( "vertex counts don't match:" , vertices.size() , " != " , vCount );
			}
		}
	}
	else
	{
		if( Out.set )
		{
			if( polygons.size() ) WriteMesh( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , polygons , comments );
			else WritePoints( Out.value , ASCII.set ? PLY_ASCII : ft , vertices , comments );
		}
	}
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , &argv[1] , params );
#ifdef USE_SEG_FAULT_HANDLER
	WARN( "using seg-fault handler" );
	StackTracer::exec = argv[0];
	signal( SIGSEGV , SignalHandler );
#endif // USE_SEG_FAULT_HANDLER
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	Timer timer;

	typedef MultiPointStreamData< float , PointStreamValue< float > , PointStreamNormal< float , 3 > , PointStreamColor< float > > VertexData;
	typedef PlyVertexWithData< float , 3 , VertexData > Vertex;
	bool readFlags[ Vertex::PlyReadNum ];
	if( !PlyReadHeader( In.value , Vertex::PlyReadProperties() , Vertex::PlyReadNum , readFlags ) ) ERROR_OUT( "Failed to read ply header: " , In.value );

	bool hasValue  = VertexData::ValidPlyReadProperties< 0 >( readFlags + 3 );
	bool hasNormal = VertexData::ValidPlyReadProperties< 1 >( readFlags + 3 );
	bool hasColor  = VertexData::ValidPlyReadProperties< 2 >( readFlags + 3 );

	if( hasValue )
		if( hasColor )
			if( hasNormal ) Execute< PointStreamValue< float > , PointStreamNormal< float , 3 > , PointStreamColor< float > >();
			else            Execute< PointStreamValue< float > ,                                  PointStreamColor< float > >();
		else
			if( hasNormal ) Execute< PointStreamValue< float > , PointStreamNormal< float , 3 >                             >();
			else            Execute< PointStreamValue< float >                                                              >();
	else
		if( hasColor )
			if( hasNormal ) Execute< PointStreamNormal< float , 3 > , PointStreamColor< float > >();
			else            Execute<                                  PointStreamColor< float > >();
		else
			if( hasNormal ) Execute< PointStreamNormal< float , 3 >                             >();
			else            Execute<                                                            >();
	printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
	printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );

	return EXIT_SUCCESS;
}
