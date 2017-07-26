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

#undef ARRAY_DEBUG

/**
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include <algorithm>
#include "Geometry.h"
#include "Ply.h"
#include "MyTime.h"
**/

#include "MAT.h"

long long EdgeKey( int key1 , int key2 )
{
	if( key1<key2 ) return ( ( (long long)key1 )<<32 ) | ( (long long)key2 );
	else            return ( ( (long long)key2 )<<32 ) | ( (long long)key1 );
}

template< class Real >
inline Point3D< Real > CrossProduct( Point3D< Real > p1 , Point3D< Real > p2 )
{
    return Point3D< Real >( p1[1]*p2[2]-p1[2]*p2[1] ,
        p1[2]*p2[0]-p1[0]*p2[2] , p1[0]*p1[1]-p1[1]*p2[0] );
}

template< class Real >
double TriangleArea( Point3D< Real > v1 , Point3D< Real > v2 ,
    Point3D< Real > v3 )
{
	Point3D< Real > n = CrossProduct( v2-v1 , v3-v1 );
	return sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] ) / 2.;
}

template <typename Vertex, typename Real>
class SurfaceTrimmer
{
public:
    typedef std::vector<int> Polygon;

private:
    float m_trim;
    float m_islandRatio;
    bool m_polygonMesh;
    std::vector<Vertex>& m_vertices;
    std::vector<Polygon>& m_polygons;

public:
    SurfaceTrimmer(std::vector<Vertex>& vertices,
        std::vector<Polygon>& polygons, Real trim, Real islandRatio = .001,
        bool polygonMesh = false) :
        m_trim(trim), m_islandRatio(islandRatio), m_polygonMesh(polygonMesh),
        m_vertices(vertices), m_polygons(polygons)
    {}
    int Execute();
    void SmoothValues();

private:
    Vertex InterpolateVertices( const Vertex& v1 , const Vertex& v2 ,
        Real value );
    void SplitPolygon(const Polygon& polygon, std::vector<Polygon>* ltPolygons,
        std::vector<Polygon>* gtPolygons , std::vector< bool >* ltFlags,
        std::vector< bool >* gtFlags ,
        std::unordered_map< long long, int >& vertexTable, Real trimValue);
    void Triangulate( const std::vector<Polygon>& polygons ,
        std::vector< std::vector< int > >& triangles );
    double PolygonArea( const Polygon& polygon );
    void RemoveHangingVertices( std::vector<Polygon>& polygons );
    void SetConnectedComponents( const std::vector<Polygon>& polygons ,
        std::vector< std::vector< int > >& components );
    double calcAreas(std::vector<Polygon>& polygons,
        std::vector<std::vector<int>>& components, std::vector<bool>& flags,
        std::vector<double>& areas);
    void reSortPolygons(std::vector<Polygon>& ltPolygons,
        std::vector<std::vector<int>>& ltComponents, std::vector<bool>& ltFlags,
        std::vector<double>& ltAreas,
        std::vector<Polygon>& gtPolygons,
        std::vector<std::vector<int>>& gtComponents, std::vector<bool>& gtFlags,
        std::vector<double> gtAreas, double islandArea);
};

template <typename Vertex, typename Real>
Vertex SurfaceTrimmer<Vertex, Real>::InterpolateVertices( const Vertex& v1 ,
    const Vertex& v2 , Real value )
{
	typename Vertex::Wrapper _v1(v1) , _v2(v2);
	if( _v1.value==_v2.value ) return Vertex( (_v1+_v2)/Real(2.) );

	Real dx = ( _v1.value-value ) / ( _v1.value-_v2.value );
	return Vertex( _v1*(1.f-dx) + _v2*dx );
}

template <typename Vertex, typename Real>
void SurfaceTrimmer<Vertex, Real>::SmoothValues()
{
	std::vector< int > count( m_vertices.size() );
	std::vector< Real > sums( m_vertices.size() , 0 );
	for( size_t i=0 ; i<m_polygons.size() ; i++ )
	{
		int sz = int(m_polygons[i].size());
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
			int v1 = m_polygons[i][j1] , v2 = m_polygons[i][j2];
			count[v1]++ , count[v2]++;
			sums[v1] += m_vertices[v2].value , sums[v2] += m_vertices[v1].value;
		}
	}
	for( size_t i=0 ; i<m_vertices.size() ; i++ )
        m_vertices[i].value = ( sums[i] + m_vertices[i].value ) / ( count[i] + 1 );
}

template<typename Vertex , typename Real>
void SurfaceTrimmer<Vertex, Real>::SplitPolygon
	(
	const Polygon& polygon ,
	std::vector<Polygon>* ltPolygons , std::vector<Polygon>* gtPolygons ,
	std::vector< bool >* ltFlags , std::vector< bool >* gtFlags ,
	std::unordered_map< long long, int >& vertexTable,
	Real trimValue
	)
{
	int sz = int( polygon.size() );
	std::vector< bool > gt( sz );
	int gtCount = 0;
	for( int j=0 ; j<sz ; j++ )
	{
		gt[j] = ( m_vertices[ polygon[j] ].value>trimValue );
		if( gt[j] ) gtCount++;
	}
	if     ( gtCount==sz ){ if( gtPolygons ) gtPolygons->push_back( polygon ) ; if( gtFlags ) gtFlags->push_back( false ); }
	else if( gtCount==0  ){ if( ltPolygons ) ltPolygons->push_back( polygon ) ; if( ltFlags ) ltFlags->push_back( false ); }
	else
	{
		int start;
		for( start=0 ; start<sz ; start++ ) if( gt[start] && !gt[(start+sz-1)%sz] ) break;

		bool gtFlag = true;
		std::vector< int > poly;

		// Add the initial vertex
		{
			int j1 = (start+int(sz)-1)%sz , j2 = start;
			int v1 = polygon[j1] , v2 = polygon[j2];
			int vIdx;
			std::unordered_map< long long, int >::iterator iter = vertexTable.find(EdgeKey(v1, v2));
			if( iter==vertexTable.end() )
			{
				vertexTable[ EdgeKey( v1 , v2 ) ] = vIdx = int( m_vertices.size() );
				m_vertices.push_back( InterpolateVertices( m_vertices[v1] , m_vertices[v2] , trimValue ) );
			}
			else vIdx = iter->second;
			poly.push_back( vIdx );
		}

		for( int _j=0  ; _j<=sz ; _j++ )
		{
			int j1 = (_j+start+sz-1)%sz , j2 = (_j+start)%sz;
			int v1 = polygon[j1] , v2 = polygon[j2];
			if( gt[j2]==gtFlag ) poly.push_back( v2 );
			else
			{
				int vIdx;
				std::unordered_map< long long, int >::iterator iter = vertexTable.find(EdgeKey(v1, v2));
				if( iter==vertexTable.end() )
				{
					vertexTable[ EdgeKey( v1 , v2 ) ] = vIdx = int( m_vertices.size() );
					m_vertices.push_back( InterpolateVertices( m_vertices[v1] , m_vertices[v2] , trimValue ) );
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

template<typename Vertex, typename Real>
void SurfaceTrimmer<Vertex, Real>::Triangulate(
    const std::vector<Polygon>& polygons,
    std::vector< std::vector< int > >& triangles )
{
	triangles.clear();
	for( size_t i=0 ; i<polygons.size() ; i++ )
		if( polygons.size()>3 )
		{
			MinimalAreaTriangulation< Real > mat;
			std::vector< Point3D< Real > > _vertices( polygons[i].size() );
			std::vector< TriangleIndex > _triangles;
			for( int j=0 ; j<int( polygons[i].size() ) ; j++ )
                _vertices[j] = m_vertices[ polygons[i][j] ].point;
			mat.GetTriangulation( _vertices , _triangles );

			// Add the triangles to the mesh
			size_t idx = triangles.size();
			triangles.resize( idx+_triangles.size() );
			for( int j=0 ; j<int(_triangles.size()) ; j++ )
			{
				triangles[idx+j].resize(3);
				for( int k=0 ; k<3 ; k++ )
                    triangles[idx+j][k] = polygons[i][ _triangles[j].idx[k] ];
			}
		}
		else if( polygons[i].size()==3 ) triangles.push_back( polygons[i] );
}

template<typename Vertex, typename Real>
double SurfaceTrimmer<Vertex, Real>::PolygonArea( const Polygon& polygon )
{
	if( polygon.size()<3 ) return 0.;
	else if( polygon.size()==3 ) return TriangleArea( m_vertices[polygon[0]].point , m_vertices[polygon[1]].point , m_vertices[polygon[2]].point );
	else
	{
		Point3D< Real > center;
		for( size_t i=0 ; i<polygon.size() ; i++ ) center += m_vertices[ polygon[i] ].point;
		center /= Real( polygon.size() );
		double area = 0;
		for( size_t i=0 ; i<polygon.size() ; i++ ) area += TriangleArea( center , m_vertices[ polygon[i] ].point , m_vertices[ polygon[ (i+1)%polygon.size() ] ].point );
		return area;
	}
}

template<typename Vertex, typename Real>
void SurfaceTrimmer<Vertex, Real>::RemoveHangingVertices(
    std::vector<Polygon>& polygons )
{
	std::unordered_map< int, int > vMap;
	std::vector< bool > vertexFlags( m_vertices.size() , false );
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) vertexFlags[ polygons[i][j] ] = true;
	int vCount = 0;
	for( int i=0 ; i<int(m_vertices.size()) ; i++ ) if( vertexFlags[i] ) vMap[i] = vCount++;
	for( size_t i=0 ; i<polygons.size() ; i++ ) for( size_t j=0 ; j<polygons[i].size() ; j++ ) polygons[i][j] = vMap[ polygons[i][j] ];

	std::vector< Vertex > _vertices( vCount );
	for( int i=0 ; i<int(m_vertices.size()) ; i++ )
        if( vertexFlags[i] ) _vertices[ vMap[i] ] = m_vertices[i];
	m_vertices = _vertices;
}

template<typename Vertex, typename Real>
void SurfaceTrimmer<Vertex, Real>::SetConnectedComponents(
    const std::vector<Polygon>& polygons ,
    std::vector< std::vector< int > >& components )
{
	std::vector< int > polygonRoots( polygons.size() );
	for( size_t i=0 ; i<polygons.size() ; i++ ) polygonRoots[i] = int(i);
	std::unordered_map< long long, int > edgeTable;
	for( size_t i=0 ; i<polygons.size() ; i++ )
	{
		int sz = int( polygons[i].size() );
		for( int j=0 ; j<sz ; j++ )
		{
			int j1 = j , j2 = (j+1)%sz;
			int v1 = polygons[i][j1] , v2 = polygons[i][j2];
			long long eKey = EdgeKey( v1 , v2 );
			std::unordered_map< long long, int >::iterator iter = edgeTable.find(eKey);
			if( iter==edgeTable.end() ) edgeTable[ eKey ] = int(i);
			else
			{
				int p = iter->second;
				while( polygonRoots[p]!=p )
				{
					int temp = polygonRoots[p];
					polygonRoots[p] = int(i);
					p = temp;
				}
				polygonRoots[p] = int(i);
			}
		}
	}
	for( size_t i=0 ; i<polygonRoots.size() ; i++ )
	{
		int p = int(i);
		while( polygonRoots[p]!=p ) p = polygonRoots[p];
		int root = p;
		p = int(i);
		while( polygonRoots[p]!=p )
		{
			int temp = polygonRoots[p];
			polygonRoots[p] = root;
			p = temp;
		}
	}
	int cCount = 0;
	std::unordered_map< int , int > vMap;
	for( int i= 0 ; i<int(polygonRoots.size()) ; i++ ) if( polygonRoots[i]==i ) vMap[i] = cCount++;
	components.resize( cCount );
	for( int i=0 ; i<int(polygonRoots.size()) ; i++ ) components[ vMap[ polygonRoots[i] ] ].push_back(i);
}

template <typename Vertex, typename Real>
double SurfaceTrimmer<Vertex, Real>::calcAreas(std::vector<Polygon>& polygons,
    std::vector<std::vector<int>>& components, std::vector<bool>& flags,
    std::vector<double>& areas)
{
    double area = 0;
    for (size_t i = 0; i < components.size() ; i++)
    {
        for (size_t j=0 ; j<components[i].size() ; j++)
        {
            areas[i] += PolygonArea(polygons[ components[i][j] ] );
            flags[i] = ( flags[i] || flags[ components[i][j] ] );
        }
        area += areas[i];
    }
    return area;
}

template<typename Vertex, typename Real>
void SurfaceTrimmer<Vertex, Real>::reSortPolygons(
    std::vector<Polygon>& ltPolygons,
    std::vector<std::vector<int>>& ltComponents, std::vector<bool>& ltFlags,
    std::vector<double>& ltAreas,
    std::vector<Polygon>& gtPolygons,
    std::vector<std::vector<int>>& gtComponents, std::vector<bool>& gtFlags,
    std::vector<double> gtAreas, double islandArea)
{
    std::vector<Polygon> _ltPolygons , _gtPolygons;

    for( size_t i=0 ; i<ltComponents.size() ; i++ )
    {
        if( ltAreas[i] < islandArea && ltFlags[i] )
            for( size_t j=0 ; j<ltComponents[i].size() ; j++ )
                _gtPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
        else
            for( size_t j=0 ; j<ltComponents[i].size() ; j++ )
                _ltPolygons.push_back( ltPolygons[ ltComponents[i][j] ] );
    }

    for (size_t i=0 ; i<gtComponents.size() ; i++ )
    {
        if( gtAreas[i] < islandArea && gtFlags[i] )
            for( size_t j=0 ; j<gtComponents[i].size() ; j++ )
                _ltPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
        else
            for( size_t j=0 ; j<gtComponents[i].size() ; j++ )
                _gtPolygons.push_back( gtPolygons[ gtComponents[i][j] ] );
    }
    ltPolygons = _ltPolygons;
    gtPolygons = _gtPolygons;
}

template<typename Vertex, typename Real>
int SurfaceTrimmer<Vertex, Real>::Execute()
{
	std::unordered_map<long long, int> vertexTable;
	std::vector<Polygon> ltPolygons, gtPolygons;
	std::vector<bool> ltFlags , gtFlags;

    for (Polygon& poly : m_polygons)
	    SplitPolygon( poly, &ltPolygons , &gtPolygons , &ltFlags , &gtFlags ,
            vertexTable , m_trim);

	if(m_islandRatio >0)
	{
		std::vector<std::vector<int>> ltComponents , gtComponents;
		SetConnectedComponents( ltPolygons , ltComponents );
		SetConnectedComponents( gtPolygons , gtComponents );

		std::vector< double > ltAreas( ltComponents.size() , 0. );
        std::vector<double> gtAreas( gtComponents.size() , 0. );
		std::vector<bool> ltComponentFlags( ltComponents.size() , false );
        std::vector<bool> gtComponentFlags( gtComponents.size() , false );

        double area = 0;
        area += calcAreas(ltPolygons, ltComponents, ltFlags, ltAreas);
        area += calcAreas(gtPolygons, gtComponents, gtFlags, gtAreas);

        reSortPolygons(ltPolygons, ltComponents, ltFlags, ltAreas,
            gtPolygons, gtComponents, gtFlags, gtAreas,
            area * m_islandRatio);
    }

	if( !m_polygonMesh)
	{
        std::vector<Polygon> polys = ltPolygons;
        Triangulate(ltPolygons , polys );
        ltPolygons = polys;

        polys = gtPolygons;
        Triangulate(gtPolygons , polys );
        gtPolygons = polys;
    }

	RemoveHangingVertices(gtPolygons );
    m_polygons = gtPolygons;

	return EXIT_SUCCESS;
}
