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
// [COMMENTS]
// -- Throughout the code, should make a distinction between indices and offsets
// -- Make an instance of _evaluate that samples the finite-elements correctly (specifically, to handle the boundaries)
// -- Make functions like depthAndOffset parity dependent (ideally all "depth"s should be relative to the B-Slpline resolution
// -- Make all points relative to the unit-cube, regardless of degree parity
// -- It's possible that for odd degrees, the iso-surfacing will fail because the leaves in the SortedTreeNodes do not form a partition of space
// -- [MAYBE] Treat normal field as a sum of delta functions, rather than a smoothed signal (again, so that high degrees aren't forced to generate smooth reconstructions)
// -- [MAYBE] Make the degree of the B-Spline with which the normals are splatted independent of the degree of the FEM system. (This way, higher degree systems aren't forced to generate smoother normal fields.)
// -- [MAYBE] Remove the isValidFEM/isValidSpace functions since the octree supports all degrees/boundary types (up to the max degree for which finalizedBrooded... was called)

// [TODO]
// -- Currently, the implementation assumes that the boundary constraints are the same for vector fields and scalar fields
// -- Modify the setting of the flags so that only the subset of the broods that are needed 

#ifndef MULTI_GRID_OCTREE_DATA_INCLUDED
#define MULTI_GRID_OCTREE_DATA_INCLUDED

#define NEW_CODE
#define FAST_SET_UP				// If enabled, kernel density estimation is done aglomeratively

#define POINT_DATA_RES 0		// Specifies the resolution of the subgrid storing points with each voxel (0==1 but is faster)

#define DATA_DEGREE 0			// The order of the B-Spline used to splat in data for color interpolation
#define WEIGHT_DEGREE 2			// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2			// The order of the B-Spline used to splat int the normals for constructing the Laplacian constraints
//#define MAX_MEMORY_GB 15		// The maximum memory the application is allowed to use
#define MAX_MEMORY_GB 0

#include <unordered_map>
#include <omp.h>
#include "BSplineData.h"
#include "PointStream.h"
#include "Geometry.h"
#include "Octree.h"
#include "SparseMatrix.h"

#ifndef _OPENMP
int omp_get_num_procs( void ){ return 1; }
int omp_get_thread_num( void ){ return 0; }
#endif // _OPENMP

#define DERIVATIVES( Degree ) ( ( Degree>1 ) ? 2 : ( Degree==1 ? 1 : 0 ) )

class TreeNodeData
{
public:
	enum
	{
		SPACE_FLAG = 1 ,
		FEM_FLAG = 2 ,
		GHOST_FLAG = 1<<7
	};
	int nodeIndex;
	char flags;

	void setGhostFlag( bool f ){ if( f ) flags |= GHOST_FLAG ; else flags &= ~GHOST_FLAG; }
	bool getGhostFlag( void ) const { return ( flags & GHOST_FLAG )!=0; }
	TreeNodeData( void );
	~TreeNodeData( void );
};

class VertexData
{
	typedef OctNode< TreeNodeData > TreeOctNode;
public:
	static const int VERTEX_COORDINATE_SHIFT = ( sizeof( long long ) * 8 ) / 3;
	static long long   EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth , int index[DIMENSION] );
	static long long   EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth );
	static long long   FaceIndex( const TreeOctNode* node , int fIndex , int maxDepth,int index[DIMENSION] );
	static long long   FaceIndex( const TreeOctNode* node , int fIndex , int maxDepth );
	static long long CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int index[DIMENSION] );
	static long long CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth );
	static long long CenterIndex( const TreeOctNode* node , int maxDepth , int index[DIMENSION] );
	static long long CenterIndex( const TreeOctNode* node , int maxDepth );
	static long long CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int index[DIMENSION] );
	static long long CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int index[DIMENSION] );
	static long long CornerIndexKey( const int index[DIMENSION] );
};

// This class stores the octree nodes, sorted by depth and then by z-slice.
// To support primal representations, the initializer takes a function that
// determines if a node should be included/indexed in the sorted list.
// [NOTE] Indexing of nodes is _GLOBAL_
class SortedTreeNodes
{
	typedef OctNode< TreeNodeData > TreeOctNode;
protected:
	Pointer( Pointer( int ) ) _sliceStart;
	int _levels;
public:
	Pointer( TreeOctNode* ) treeNodes;
	int begin( int depth ) const{ return _sliceStart[depth][0]; }
	int   end( int depth ) const{ return _sliceStart[depth][(size_t)1<<depth]; }
	int begin( int depth , int slice ) const{ return _sliceStart[depth][slice  ]  ; }
	int   end( int depth , int slice ) const{ if(depth<0||depth>=_levels||slice<0||slice>=(1<<depth)) printf( "uh oh\n" ) ; return _sliceStart[depth][slice+1]; }
	int size( void ) const { return _sliceStart[_levels-1][(size_t)1<<(_levels-1)]; }
	int size( int depth ) const { if(depth<0||depth>=_levels) printf( "uhoh\n" ); return _sliceStart[depth][(size_t)1<<depth] - _sliceStart[depth][0]; }
	int size( int depth , int slice ) const { return _sliceStart[depth][slice+1] - _sliceStart[depth][slice]; }
	int levels( void ) const { return _levels; }

	SortedTreeNodes( void );
	~SortedTreeNodes( void );
	void set( TreeOctNode& root , std::vector< int >* map );
	void set( TreeOctNode& root );

	template< int Indices >
	struct  _Indices
	{
		int idx[Indices];
		_Indices( void ){ memset( idx , -1 , sizeof( int ) * Indices ); }
		int& operator[] ( int i ) { return idx[i]; }
		const int& operator[] ( int i ) const { return idx[i]; }
	};
	typedef _Indices< Square::CORNERS > SquareCornerIndices;
	typedef _Indices< Square::EDGES > SquareEdgeIndices;
	typedef _Indices< Square::FACES > SquareFaceIndices;

	struct SliceTableData
	{
		Pointer( SquareCornerIndices ) cTable;
		Pointer( SquareEdgeIndices   ) eTable;
		Pointer( SquareFaceIndices   ) fTable;
		int cCount , eCount , fCount , nodeOffset , nodeCount;
		SliceTableData( void ){ fCount = eCount = cCount = 0 , cTable = NullPointer( SquareCornerIndices ) , eTable = NullPointer( SquareEdgeIndices ) , fTable = NullPointer( SquareFaceIndices ) , _cMap = _eMap = _fMap = NullPointer( int ); }
		~SliceTableData( void ){ clear(); }
#ifdef BRUNO_LEVY_FIX
		void clear( void ){ DeletePointer( cTable ) ; DeletePointer( eTable ) ; DeletePointer( fTable ) ; DeletePointer( _cMap ) ; DeletePointer( _eMap ) ; DeletePointer( _fMap ) ; fCount = eCount = cCount = 0; }
#else // !BRUNO_LEVY_FIX
		void clear( void ){ DeletePointer( cTable ) ; DeletePointer( eTable ) ; DeletePointer( fTable ) ; fCount = eCount = cCount = 0; }
#endif // BRUNO_LEVY_FIX
		SquareCornerIndices& cornerIndices( const TreeOctNode* node );
		SquareCornerIndices& cornerIndices( int idx );
		const SquareCornerIndices& cornerIndices( const TreeOctNode* node ) const;
		const SquareCornerIndices& cornerIndices( int idx ) const;
		SquareEdgeIndices& edgeIndices( const TreeOctNode* node );
		SquareEdgeIndices& edgeIndices( int idx );
		const SquareEdgeIndices& edgeIndices( const TreeOctNode* node ) const;
		const SquareEdgeIndices& edgeIndices( int idx ) const;
		SquareFaceIndices& faceIndices( const TreeOctNode* node );
		SquareFaceIndices& faceIndices( int idx );
		const SquareFaceIndices& faceIndices( const TreeOctNode* node ) const;
		const SquareFaceIndices& faceIndices( int idx ) const;
	protected:
		Pointer( int ) _cMap;
		Pointer( int ) _eMap;
		Pointer( int ) _fMap;
		friend class SortedTreeNodes;
	};
	struct XSliceTableData
	{
		Pointer( SquareCornerIndices ) eTable;
		Pointer( SquareEdgeIndices ) fTable;
		int fCount , eCount , nodeOffset , nodeCount;
		XSliceTableData( void ){ fCount = eCount = 0 , eTable = NullPointer( SquareCornerIndices ) , fTable = NullPointer( SquareEdgeIndices ) , _eMap = _fMap = NullPointer( int ); }
		~XSliceTableData( void ){ clear(); }
#ifdef BRUNO_LEVY_FIX
		void clear( void ) { DeletePointer( fTable ) ; DeletePointer( eTable ) ; DeletePointer( _eMap ) ; DeletePointer( _fMap ) ; fCount = eCount = 0; }
#else // !BRUNO_LEVY_FIX
		void clear( void ) { DeletePointer( fTable ) ; DeletePointer( eTable ) ; fCount = eCount = 0; }
#endif // BRUNO_LEVY_FIX
		SquareCornerIndices& edgeIndices( const TreeOctNode* node );
		SquareCornerIndices& edgeIndices( int idx );
		const SquareCornerIndices& edgeIndices( const TreeOctNode* node ) const;
		const SquareCornerIndices& edgeIndices( int idx ) const;
		SquareEdgeIndices& faceIndices( const TreeOctNode* node );
		SquareEdgeIndices& faceIndices( int idx );
		const SquareEdgeIndices& faceIndices( const TreeOctNode* node ) const;
		const SquareEdgeIndices& faceIndices( int idx ) const;
	protected:
		Pointer( int ) _eMap;
		Pointer( int ) _fMap;
		friend class SortedTreeNodes;
	};
	void setSliceTableData (  SliceTableData& sData , int depth , int offset , int threads ) const;
	void setXSliceTableData( XSliceTableData& sData , int depth , int offset , int threads ) const;
};

template< int Degree >
struct PointSupportKey : public OctNode< TreeNodeData >::NeighborKey< BSplineSupportSizes< Degree >::SupportEnd , -BSplineSupportSizes< Degree >::SupportStart >
{
	static const int LeftRadius  =  BSplineSupportSizes< Degree >::SupportEnd;
	static const int RightRadius = -BSplineSupportSizes< Degree >::SupportStart;
	static const int Size = LeftRadius + RightRadius + 1;
};
template< int Degree >
struct ConstPointSupportKey : public OctNode< TreeNodeData >::ConstNeighborKey< BSplineSupportSizes< Degree >::SupportEnd , -BSplineSupportSizes< Degree >::SupportStart >
{
	static const int LeftRadius  =  BSplineSupportSizes< Degree >::SupportEnd;
	static const int RightRadius = -BSplineSupportSizes< Degree >::SupportStart;
	static const int Size = LeftRadius + RightRadius + 1;
};

template< class Real , bool HasGradients >
struct SinglePointData
{
	Point3D< Real > position;
	Real weight;
	Real value , _value;
	SinglePointData  operator +  ( const SinglePointData& p ) const { return SinglePointData( position + p.position , value + p.value , weight + p.weight ); }
	SinglePointData& operator += ( const SinglePointData& p ){ position += p.position ; weight += p.weight , value += p.value ; return *this; }
	SinglePointData  operator *  ( Real s ) const { return SinglePointData( position*s , weight*s , value*s ); }
	SinglePointData& operator *= ( Real s ){ position *= s , weight *= s , value *= s ; return *this; }
	SinglePointData  operator /  ( Real s ) const { return SinglePointData( position/s , weight/s , value/s ); }
	SinglePointData& operator /= ( Real s ){ position /= s , weight /= s , value /= s ; return *this; }
	SinglePointData( void ) : position( Point3D< Real >() ) , weight(0) , value(0) , _value(0) { ; }
	SinglePointData( Point3D< Real > p , Real v , Real w ) { position = p , value = v , weight = w , _value = (Real)0; }
};
template< class Real >
struct SinglePointData< Real , true > : public SinglePointData< Real , false >
{
	using SinglePointData< Real , false >::position;
	using SinglePointData< Real , false >::weight;
	using SinglePointData< Real , false >::value;
	using SinglePointData< Real , false >::_value;
	Point3D< Real > gradient , _gradient;
	SinglePointData  operator +  ( const SinglePointData& p ) const { return SinglePointData( position + p.position , weight + p.weight , value + p.value , gradient + p.gradient ); }
	SinglePointData& operator += ( const SinglePointData& p ){ position += p.position , weight += p.weight , value += p.value , gradient += p.gradient ; return *this; }
	SinglePointData  operator *  ( Real s ) const { return SinglePointData( position*s , weight*s , value*s , gradient*s ); }
	SinglePointData& operator *= ( Real s ){ position *= s , weight *= s , value *= s , gradient *= s ; return *this; }
	SinglePointData  operator /  ( Real s ) const { return SinglePointData( position/s , weight/s , value/s , gradient/s ); }
	SinglePointData& operator /= ( Real s ){ position /= s , weight /= s , value /= s , gradient /= s ; return *this; }
	SinglePointData( void ) : SinglePointData< Real , false >() , gradient( Point3D< Real >() ) , _gradient( Point3D< Real >() ) { ; }
	SinglePointData( Point3D< Real > p , Real v , Point3D< Real > g , Real w ) : SinglePointData< Real , false >( p , v , w ) { gradient = g , _gradient = Point3D< Real >(); }
};

#if POINT_DATA_RES
template< class Real , bool HasGradients >
struct PointData
{
	static const int RES = POINT_DATA_RES;
	static const int SAMPLES = RES * RES * RES;

	SinglePointData< Real , HasGradients > points[SAMPLES];
	SinglePointData< Real , HasGradients >& operator[] ( int idx ) { return points[idx]; }
	const SinglePointData< Real , HasGradients >& operator[] ( int idx ) const { return points[idx]; }

	static void SetIndices( Point3D< Real > p , Point3D< Real > c , Real w , int x[3] )
	{
		for( int d=0 ; d<3 ; d++ ) x[d] = std::max< int >( 0 , std::min< int >( RES-1 , int( floor( ( p[d]-( c[d]-w/2 ) ) / w * RES ) ) ) );
	}

	void addPoint( SinglePointData< Real , HasGradients > p , Point3D< Real > center , Real width  )
	{
		int x[3];
		SetIndices( p.position , center , width , x );
		points[ x[0]+x[1]*RES+x[2]*RES*RES ] += p;
	}

	PointData  operator +  ( const PointData& p ) const { PointData _p ; for( int c=0 ; c<SAMPLES ;  c++ ) _p.points[c] = points[c] + _p.points[c] ; return _p; }
	PointData& operator += ( const PointData& p ){ for( int c=0 ; c<SAMPLES ; c++ ) points[c] += p.points[c] ; return *this; }
	PointData  operator *  ( Real s ) const { PointData _p ; for( int c=0 ; c<SAMPLES ;  c++ ) _p.points[c] = points[c] * s ; return _p; }
	PointData& operator *= ( Real s ){ for( int c=0 ; c<SAMPLES ; c++ ) points[c] *= s ; return *this; }
	PointData  operator /  ( Real s ) const { PointData _p ; for( int c=0 ; c<SAMPLES ;  c++ ) _p.points[c] = points[c] / s ; return _p; }
	PointData& operator /= ( Real s ){ for( int c=0 ; c<SAMPLES ; c++ ) points[c] /= s ; return *this; }
};
#else // !POINT_DATA_RES
template< class Real , bool HasGradients > using PointData = SinglePointData< Real , HasGradients >;
#endif // POINT_DATA_RES

template< class Data , int Degree >
struct SparseNodeData
{
	size_t size( void ) const { return _data.size(); }
	const Data& operator[] ( int idx ) const { return _data[idx]; }
	Data& operator[] ( int idx ) { return _data[idx]; }
	void reserve( size_t sz ){ if( sz>_indices.size() ) _indices.resize( sz , -1 ); }
	Data* operator()( const OctNode< TreeNodeData >* node ){ return ( node->nodeData.nodeIndex<0 || node->nodeData.nodeIndex>=(int)_indices.size() || _indices[ node->nodeData.nodeIndex ]<0 ) ? NULL : &_data[ _indices[ node->nodeData.nodeIndex ] ]; }
	const Data* operator()( const OctNode< TreeNodeData >* node ) const { return ( node->nodeData.nodeIndex<0 || node->nodeData.nodeIndex>=(int)_indices.size() || _indices[ node->nodeData.nodeIndex ]<0 ) ? NULL : &_data[ _indices[ node->nodeData.nodeIndex ] ]; }
	Data& operator[]( const OctNode< TreeNodeData >* node )
	{
		if( node->nodeData.nodeIndex>=(int)_indices.size() ) _indices.resize( node->nodeData.nodeIndex+1 , -1 );
		if( _indices[ node->nodeData.nodeIndex ]==-1 )
		{
			_indices[ node->nodeData.nodeIndex ] = (int)_data.size();
			_data.push_back( Data() );
		}
		return _data[ _indices[ node->nodeData.nodeIndex ] ];
	}
	void remapIndices( const std::vector< int >& map )
	{
		std::vector< int > temp = _indices;
		_indices.resize( map.size() );
		for( size_t i=0 ; i<map.size() ; i++ )
			if( map[i]<(int)temp.size() ) _indices[i] = temp[ map[i] ];
			else                          _indices[i] = -1;
	}
	template< class _Data , int _Degree > friend struct SparseNodeData;
	template< class _Data , int _Degree >
	void init( const SparseNodeData< _Data , _Degree >& snd ){ _indices = snd._indices , _data.resize( snd._data.size() ); }
	void remove( const OctNode< TreeNodeData >* node ){ if( node->nodeData.nodeIndex<(int)_indices.size() && node->nodeData.nodeIndex>=0 ) _indices[ node->nodeData.nodeIndex ] = -1; }
protected:
	std::vector< int > _indices;
	std::vector< Data > _data;
};
template< class Data , int Degree >
struct DenseNodeData
{
	DenseNodeData( void ){ _data = NullPointer( Data ) ; _sz = 0; }
	DenseNodeData( size_t sz ){ _sz = sz ; if( sz ) _data = NewPointer< Data >( sz ) ; else _data = NullPointer( Data ); }
	DenseNodeData( const DenseNodeData&  d ) : DenseNodeData() { _resize( d._sz ) ; if( _sz ) memcpy( _data , d._data , sizeof(Data) * _sz ); }
	DenseNodeData(       DenseNodeData&& d ){ _data = d._data , _sz = d._sz ; d._data = NullPointer( Data ) , d._sz = 0; }
	DenseNodeData& operator = ( const DenseNodeData&  d ){ _resize( d._sz ) ; if( _sz ) memcpy( _data , d._data , sizeof(Data) * _sz ) ; return *this; }
	DenseNodeData& operator = (       DenseNodeData&& d ){ size_t __sz = _sz ; Pointer( Data ) __data = _data ; _data = d._data , _sz = d._sz ; d._data = __data , d._sz = __sz ; return *this; }
	~DenseNodeData( void ){ DeletePointer( _data ) ; _sz = 0; }

	Data& operator[] ( int idx ) { return _data[idx]; }
	const Data& operator[] ( int idx ) const { return _data[idx]; }
	size_t size( void ) const { return _sz; }
	Data& operator[]( const OctNode< TreeNodeData >* node ) { return _data[ node->nodeData.nodeIndex ]; }
	Data* operator()( const OctNode< TreeNodeData >* node ) { return ( node==NULL || node->nodeData.nodeIndex>=(int)_sz ) ? NULL : &_data[ node->nodeData.nodeIndex ]; }
	const Data* operator()( const OctNode< TreeNodeData >* node ) const { return ( node==NULL || node->nodeData.nodeIndex>=(int)_sz ) ? NULL : &_data[ node->nodeData.nodeIndex ]; }
	int index( const OctNode< TreeNodeData >* node ) const { return ( !node || node->nodeData.nodeIndex<0 || node->nodeData.nodeIndex>=(int)_data.size() ) ? -1 : node->nodeData.nodeIndex; }
protected:
	size_t _sz;
	void _resize( size_t sz ){ DeletePointer( _data ) ; if( sz ) _data = NewPointer< Data >( sz ) ; else _data = NullPointer( Data ) ; _sz = sz; }
	Pointer( Data ) _data;
};

// This is may be necessary in case the memory usage is larger than what fits on the stack
template< class C , int N > struct Stencil
{
	Stencil( void ){ _values = NewPointer< C >( N * N * N ); }
	~Stencil( void ){ DeletePointer( _values ); }
	C& operator()( int i , int j , int k ){ return _values[ i*N*N + j*N + k ]; }
	const C& operator()( int i , int j , int k ) const { return _values[ i*N*N + j*N + k ]; }
protected:
	Pointer( C ) _values;
};

template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 >
class SystemCoefficients
{
	typedef typename BSplineIntegrationData< Degree1 , BType1 , Degree2 , BType2 >::FunctionIntegrator FunctionIntegrator;
	static const int OverlapSize  = BSplineOverlapSizes< Degree1 , Degree2 >::OverlapSize;
	static const int OverlapStart = BSplineOverlapSizes< Degree1 , Degree2 >::OverlapStart;
	static const int OverlapEnd   = BSplineOverlapSizes< Degree1 , Degree2 >::OverlapEnd;
public:
	typedef typename BSplineIntegrationData< Degree1 , BType1 , Degree2 , BType2 >::FunctionIntegrator::template      Integrator< DERIVATIVES( Degree1 ) , DERIVATIVES( Degree2 ) >      Integrator;
	typedef typename BSplineIntegrationData< Degree1 , BType1 , Degree2 , BType2 >::FunctionIntegrator::template ChildIntegrator< DERIVATIVES( Degree1 ) , DERIVATIVES( Degree2 ) > ChildIntegrator;

	// The FEMSystemFunctor is a class that takes an object of type Integrator/ChildIntegrator, as well as a pair of indices of octree nodes
	// and returns the corresponding system coefficient.
	template< class _FEMSystemFunctor > static void SetCentralSystemStencil ( const _FEMSystemFunctor& F , const      Integrator& integrator , Stencil< double , OverlapSize >& stencil           );
	template< class _FEMSystemFunctor > static void SetCentralSystemStencils( const _FEMSystemFunctor& F , const ChildIntegrator& integrator , Stencil< double , OverlapSize >  stencils[2][2][2] );
	template< bool Reverse , class _FEMSystemFunctor > static void SetCentralConstraintStencil ( const _FEMSystemFunctor& F , const      Integrator& integrator , Stencil<          double   , OverlapSize >& stencil           );
	template< bool Reverse , class _FEMSystemFunctor > static void SetCentralConstraintStencils( const _FEMSystemFunctor& F , const ChildIntegrator& integrator , Stencil<          double   , OverlapSize >  stencils[2][2][2] );
	template< bool Reverse , class _FEMSystemFunctor > static void SetCentralConstraintStencil ( const _FEMSystemFunctor& F , const      Integrator& integrator , Stencil< Point3D< double > , OverlapSize >& stencil           );
	template< bool Reverse , class _FEMSystemFunctor > static void SetCentralConstraintStencils( const _FEMSystemFunctor& F , const ChildIntegrator& integrator , Stencil< Point3D< double > , OverlapSize >  stencils[2][2][2] );
};

template< int FEMDegree , BoundaryType BType >
struct FEMSystemFunctor
{
	double massWeight , lapWeight , biLapWeight;
	FEMSystemFunctor( double mWeight=0 , double lWeight=0 , double bWeight=0 ) : massWeight( mWeight ) , lapWeight( lWeight ) , biLapWeight( bWeight ) { ; }
	double integrate( const typename SystemCoefficients< FEMDegree , BType , FEMDegree , BType >::     Integrator& integrator , const int off1[] , const int off2[] ) const { return _integrate( integrator , off1 , off2 ); }
	double integrate( const typename SystemCoefficients< FEMDegree , BType , FEMDegree , BType >::ChildIntegrator& integrator , const int off1[] , const int off2[] ) const { return _integrate( integrator , off1 , off2 ); }
	bool vanishesOnConstants( void ) const { return massWeight==0; }
protected:
	template< class I > double _integrate( const I& integrator , const int off1[] , const int off2[] ) const;
};
template< int SFDegree , BoundaryType SFBType , int FEMDegree , BoundaryType FEMBType >
struct FEMSFConstraintFunctor
{
	double massWeight , lapWeight , biLapWeight;
	FEMSFConstraintFunctor( double mWeight=0 , double lWeight=0 , double bWeight=0 ) : massWeight( mWeight ) , lapWeight( lWeight ) , biLapWeight( bWeight ) { ; }
	template< bool Reverse >
	double integrate( const typename SystemCoefficients< Reverse ? FEMDegree : SFDegree , Reverse ? FEMBType : SFBType , Reverse ? SFDegree : FEMDegree , Reverse ? SFBType : FEMBType >::     Integrator& integrator , const int off1[] , const int off2[] ) const { return _integrate< Reverse >( integrator , off1 , off2 ); }
	template< bool Reverse >
	double integrate( const typename SystemCoefficients< Reverse ? FEMDegree : SFDegree , Reverse ? FEMBType : SFBType , Reverse ? SFDegree : FEMDegree , Reverse ? SFBType : FEMBType >::ChildIntegrator& integrator , const int off1[] , const int off2[] ) const { return _integrate< Reverse >( integrator , off1 , off2 ); }
protected:
	template< bool Reverse , class I > double _integrate( const I& integrator , const int off1[] , const int off[2] ) const;
};
template< int VFDegree , BoundaryType VFBType , int FEMDegree , BoundaryType FEMBType >
struct FEMVFConstraintFunctor
{
	double lapWeight , biLapWeight;
	FEMVFConstraintFunctor( double lWeight=0 , double bWeight=0 ) : lapWeight( lWeight ) , biLapWeight( bWeight ) { ; }
	template< bool Reverse >
	Point3D< double > integrate( const typename SystemCoefficients< Reverse ? FEMDegree : VFDegree , Reverse ? FEMBType : VFBType , Reverse ? VFDegree : FEMDegree , Reverse ? VFBType : FEMBType >::     Integrator& integrator , const int off1[] , const int off2[] ) const { return _integrate< Reverse >( integrator , off1 , off2 ); }
	template< bool Reverse >
	Point3D< double > integrate( const typename SystemCoefficients< Reverse ? FEMDegree : VFDegree , Reverse ? FEMBType : VFBType , Reverse ? VFDegree : FEMDegree , Reverse ? VFBType : FEMBType >::ChildIntegrator& integrator , const int off1[] , const int off2[] ) const { return _integrate< Reverse >( integrator , off1 , off2 ); }
protected:
	template< bool Reverse , class I > Point3D< double > _integrate( const I& integrator , const int off1[] , const int off[2] ) const;
};

inline void SetGhostFlag( OctNode< TreeNodeData >* node , bool flag ){ if( node && node->parent ) node->parent->nodeData.setGhostFlag( flag ); }
inline bool GetGhostFlag( const OctNode< TreeNodeData >* node ){ return node==NULL || node->parent==NULL || node->parent->nodeData.getGhostFlag( ); }
inline bool IsActiveNode( const OctNode< TreeNodeData >* node ){ return !GetGhostFlag( node ); }

template< class Real >
class Octree
{
	typedef OctNode< TreeNodeData > TreeOctNode;
	static int _NodeCount;
	static void _NodeInitializer( TreeOctNode& node ){ node.nodeData.nodeIndex = _NodeCount++; }
public:
#if 0
	struct LocalDepth
	{
		LocalDepth( int d=0 ) : _d(d) { ; }
		operator int&()       { return _d; }
		operator int () const { return _d; }
	protected:
		int _d;
	};
	struct LocalOffset
	{
		LocalOffset( const int* off=NULL ){ if( off ) memcpy( _off , off , sizeof(_off) ) ; else memset( _off , 0 , sizeof( _off ) ); }
		operator        int*()       { return _off; }
		operator const  int*() const { return _off; }
	protected:
		int _off[3];
	};
#else
	typedef int LocalDepth;
	typedef int LocalOffset[3];
#endif

	static void ResetNodeCount( void ){ _NodeCount = 0 ; }
	static int NodeCount( void ){ return _NodeCount; }
	template< int FEMDegree , BoundaryType BType > void functionIndex( const TreeOctNode* node , int idx[3] ) const;

	struct PointSample{ const TreeOctNode* node ; ProjectiveData< OrientedPoint3D< Real > , Real > sample; };

	typedef typename TreeOctNode::     NeighborKey< 1 , 1 >      AdjacenctNodeKey;
	typedef typename TreeOctNode::ConstNeighborKey< 1 , 1 > ConstAdjacenctNodeKey;

	template< int FEMDegree , BoundaryType BType > bool isValidFEMNode( const TreeOctNode* node ) const;
	bool isValidSpaceNode( const TreeOctNode* node ) const;
	TreeOctNode* leaf( Point3D< Real > p );
	const TreeOctNode* leaf( Point3D< Real > p ) const;

	template< bool HasGradients >
	struct InterpolationInfo
	{
		SparseNodeData< PointData< Real , HasGradients > , 0 > iData;
		Real valueWeight , gradientWeight;
		InterpolationInfo( const class Octree< Real >& tree , const std::vector< PointSample >& samples , Real pointValue , int adaptiveExponent , Real v , Real g ) : valueWeight(v) , gradientWeight(g)
		{ iData = tree._densifyInterpolationInfo< HasGradients >( samples , pointValue , adaptiveExponent ); }
		PointData< Real , HasGradients >* operator()( const OctNode< TreeNodeData >* node ){ return iData(node); }
		const PointData< Real , HasGradients >* operator()( const OctNode< TreeNodeData >* node ) const { return iData(node); }
	};

	template< int DensityDegree > struct DensityEstimator : public SparseNodeData< Real , DensityDegree >
	{
		DensityEstimator( int kernelDepth ) : _kernelDepth( kernelDepth ){ ; }
		int kernelDepth( void ) const { return _kernelDepth; }
	protected:
		int _kernelDepth;
	};
protected:
	bool _isValidSpaceNode( const TreeOctNode* node ) const { return !GetGhostFlag( node ) && ( node->nodeData.flags & TreeNodeData::SPACE_FLAG ); }
	bool _isValidFEMNode( const TreeOctNode* node ) const { return !GetGhostFlag( node ) && ( node->nodeData.flags & TreeNodeData::FEM_FLAG ); }

	TreeOctNode* _tree;
	TreeOctNode* _spaceRoot;
	SortedTreeNodes _sNodes;
	LocalDepth _fullDepth , _maxDepth;

	static bool _InBounds( Point3D< Real > p );

	int _depthOffset;
	int _localToGlobal( LocalDepth d ) const { return d + _depthOffset; }
	LocalDepth _localDepth( const TreeOctNode* node ) const { return node->depth() - _depthOffset; }
	LocalDepth _localMaxDepth( const TreeOctNode* tree ) const { return tree->maxDepth() - _depthOffset; }
	int _localInset( LocalDepth d ) const { return _depthOffset<=1 ? 0 : 1<<( d + _depthOffset - 1 ); }
	void _localDepthAndOffset( const TreeOctNode* node , LocalDepth& d , LocalOffset& off ) const
	{
		node->depthAndOffset( d , off ) ; d -= _depthOffset;
		int inset = _localInset( d );
		off[0] -= inset , off[1] -= inset , off[2] -= inset;
	}
	template< int FEMDegree , BoundaryType BType > static int _BSplineBegin( LocalDepth depth ){ return BSplineEvaluationData< FEMDegree , BType >::Begin( depth ); }
	template< int FEMDegree , BoundaryType BType > static int _BSplineEnd  ( LocalDepth depth ){ return BSplineEvaluationData< FEMDegree , BType >::End  ( depth ); }
	template< int FEMDegree , BoundaryType BType >
	bool _outOfBounds( const TreeOctNode* node ) const
	{
		if( !node ) return true;
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		return d<0 || BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[0] ) || BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[1] ) || BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[2] );
	}
	int _sNodesBegin( LocalDepth d ) const { return _sNodes.begin( _localToGlobal( d ) ); }
	int _sNodesEnd  ( LocalDepth d ) const { return _sNodes.end  ( _localToGlobal( d ) ); }
	int _sNodesSize ( LocalDepth d ) const { return _sNodes.size ( _localToGlobal( d ) ); }
	int _sNodesBegin( LocalDepth d , int slice ) const { return _sNodes.begin( _localToGlobal( d ) , slice + _localInset( d ) ); }
	int _sNodesEnd  ( LocalDepth d , int slice ) const { return _sNodes.end  ( _localToGlobal( d ) , slice + _localInset( d ) ); }
	int _sNodesSize ( LocalDepth d , int slice ) const { return _sNodes.size ( _localToGlobal( d ) , slice + _localInset( d ) ); }

	template< int FEMDegree > static bool _IsInteriorlySupported( LocalDepth depth , const LocalOffset off )
	{
		if( depth>=0 )
		{
			int begin , end;
			BSplineSupportSizes< FEMDegree >::InteriorSupportedSpan( depth , begin , end );
			return ( off[0]>=begin && off[0]<end && off[1]>=begin && off[1]<end && off[2]>=begin && off[2]<end );
		}
		else return false;
	}
	template< int FEMDegree > bool _isInteriorlySupported( const TreeOctNode* node ) const
	{
		if( !node ) return false;
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		return _IsInteriorlySupported< FEMDegree >( d , off );
	}
	template< int FEMDegree1 , int FEMDegree2 > static bool _IsInteriorlyOverlapped( LocalDepth depth , const LocalOffset off )
	{
		if( depth>=0 )
		{
			int begin , end;
			BSplineIntegrationData< FEMDegree1 , BOUNDARY_NEUMANN , FEMDegree2 , BOUNDARY_NEUMANN >::InteriorOverlappedSpan( depth , begin , end );
			return ( off[0]>=begin && off[0]<end && off[1]>=begin && off[1]<end && off[2]>=begin && off[2]<end );
		}
		else return false;
	}
	template< int FEMDegree1 , int FEMDegree2 > bool _isInteriorlyOverlapped( const TreeOctNode* node ) const
	{
		if( !node ) return false;
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		return _IsInteriorlyOverlapped< FEMDegree1 , FEMDegree2 >( d , off );
	}
	void _startAndWidth( const TreeOctNode* node , Point3D< Real >& start , Real& width ) const
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		if( d>=0 ) width = Real( 1.0 / (1<<  d ) );
		else       width = Real( 1.0 * (1<<(-d)) );
		for( int dd=0 ; dd<DIMENSION ; dd++ ) start[dd] = Real( off[dd] ) * width;
	}
	void _centerAndWidth( const TreeOctNode* node , Point3D< Real >& center , Real& width ) const
	{
		int d , off[3];
		_localDepthAndOffset( node , d , off );
		width = Real( 1.0 / (1<<d) );
		for( int dd=0 ; dd<DIMENSION ; dd++ ) center[dd] = Real( off[dd] + 0.5 ) * width;
	}
	int _childIndex( const TreeOctNode* node , Point3D< Real > p ) const
	{
		Point3D< Real > c ; Real w;
		_centerAndWidth( node , c , w );
		return ( p[0]<c[0] ? 0 : 1 ) | ( p[1]<c[1] ? 0 : 2 ) | ( p[2]<c[2] ? 0 : 4 );
	}

	template< int Degree , BoundaryType BType > void _setFullDepth( TreeOctNode* node , LocalDepth depth ) const;
	template< int Degree , BoundaryType BType > void _setFullDepth( LocalDepth depth );

	template< int LeftRadius , int RightRadius >
	static typename TreeOctNode::ConstNeighbors< LeftRadius + RightRadius + 1 >& _neighbors( TreeOctNode::ConstNeighborKey< LeftRadius , RightRadius >& key , const TreeOctNode* node ){ return key.neighbors[ node->depth() ]; }
	template< int LeftRadius , int RightRadius >
	static typename TreeOctNode::Neighbors< LeftRadius + RightRadius + 1 >& _neighbors( TreeOctNode::NeighborKey< LeftRadius , RightRadius >& key , const TreeOctNode* node ){ return key.neighbors[ node->depth() ]; }
	template< int LeftRadius , int RightRadius >
	static const typename TreeOctNode::template Neighbors< LeftRadius + RightRadius + 1 >& _neighbors( const typename TreeOctNode::template NeighborKey< LeftRadius , RightRadius >& key , const TreeOctNode* node ){ return key.neighbors[ node->depth() ]; }
	template< int LeftRadius , int RightRadius >
	static const typename TreeOctNode::template ConstNeighbors< LeftRadius + RightRadius + 1 >& _neighbors( const typename TreeOctNode::template ConstNeighborKey< LeftRadius , RightRadius >& key , const TreeOctNode* node ){ return key.neighbors[ node->depth() ]; }

public:
	LocalDepth depth( const TreeOctNode* node ) const { return _localDepth( node ); }
	void depthAndOffset( const TreeOctNode* node , LocalDepth& depth , LocalOffset& offset ) const { _localDepthAndOffset( node , depth , offset ); }

	int nodesBegin( LocalDepth d ) const { return _sNodes.begin( _localToGlobal( d ) ); }
	int nodesEnd  ( LocalDepth d ) const { return _sNodes.end  ( _localToGlobal( d ) ); }
	int nodesSize ( LocalDepth d ) const { return _sNodes.size ( _localToGlobal( d ) ); }
	int nodesBegin( LocalDepth d , int slice ) const { return _sNodes.begin( _localToGlobal( d ) , slice + _localInset( d ) ); }
	int nodesEnd  ( LocalDepth d , int slice ) const { return _sNodes.end  ( _localToGlobal( d ) , slice + _localInset( d ) ); }
	int nodesSize ( LocalDepth d , int slice ) const { return _sNodes.size ( _localToGlobal( d ) , slice + _localInset( d ) ); }
	const TreeOctNode* node( int idx ) const { return _sNodes.treeNodes[idx]; }
protected:

	////////////////////////////////////
	// System construction code       //
	// MultiGridOctreeData.System.inl //
	////////////////////////////////////
	template< int FEMDegree >
	void _setMultiColorIndices( int start , int end , std::vector< std::vector< int > >& indices ) const;
	struct _SolverStats
	{
		double evaluateTime , systemTime , solveTime;
		double bNorm2 , inRNorm2 , outRNorm2;
	};
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	int _solveSystemGS( const FEMSystemFunctor& F , const BSplineData< FEMDegree , BType >& bsData , InterpolationInfo< HasGradients >* interpolationInfo , LocalDepth depth , DenseNodeData< Real , FEMDegree >& solution , DenseNodeData< Real , FEMDegree >& constraints , DenseNodeData< Real , FEMDegree >& metSolutionConstraints , int iters , bool coarseToFine , _SolverStats& stats , bool computeNorms );
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	int _solveSystemCG( const FEMSystemFunctor& F , const BSplineData< FEMDegree , BType >& bsData , InterpolationInfo< HasGradients >* interpolationInfo , LocalDepth depth , DenseNodeData< Real , FEMDegree >& solution , DenseNodeData< Real , FEMDegree >& constraints , DenseNodeData< Real , FEMDegree >& metSolutionConstraints , int iters , bool coarseToFine , _SolverStats& stats , bool computeNorms , double accuracy );
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	int _setMatrixRow( const FEMSystemFunctor& F , const InterpolationInfo< HasGradients >* interpolationInfo , const typename TreeOctNode::Neighbors< BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& neighbors , Pointer( MatrixEntry< Real > ) row , int offset , const typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template Integrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& integrator , const Stencil< double , BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& stencil , const BSplineData< FEMDegree , BType >& bsData ) const;
	template< int FEMDegree , BoundaryType BType >
	int _getMatrixRowSize( const typename TreeOctNode::Neighbors< BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& neighbors ) const;

	template< int FEMDegree1 , int FEMDegree2 > static void _SetParentOverlapBounds( const TreeOctNode* node , int& startX , int& endX , int& startY , int& endY , int& startZ , int& endZ );
	// Updates the constraints @(depth) based on the solution coefficients @(depth-1)

	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	void _updateConstraintsFromCoarser( const FEMSystemFunctor& F , const InterpolationInfo< HasGradients >* interpolationInfo , const typename TreeOctNode::Neighbors< BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& neighbors , const typename TreeOctNode::Neighbors< BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& pNeighbors , TreeOctNode* node , DenseNodeData< Real , FEMDegree >& constraints , const DenseNodeData< Real , FEMDegree >& metSolution , const typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template ChildIntegrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& childIntegrator , const Stencil< double , BSplineOverlapSizes< FEMDegree , FEMDegree >::OverlapSize >& stencil , const BSplineData< FEMDegree , BType >& bsData ) const;

	// evaluate the points @(depth) using coefficients @(depth-1)
	template< int FEMDegree , BoundaryType BType , bool HasGradients >
	void _setPointValuesFromCoarser( InterpolationInfo< HasGradients >& interpolationInfo , LocalDepth highDepth , const BSplineData< FEMDegree , BType >& bsData , const DenseNodeData< Real , FEMDegree >& upSampledCoefficients );

	// Updates the cumulative integral constraints @(depth-1) based on the change in solution coefficients @(depth)
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor >
	void _updateCumulativeIntegralConstraintsFromFiner( const FEMSystemFunctor& F , 
		const BSplineData< FEMDegree , BType >& bsData , LocalDepth highDepth , const DenseNodeData< Real , FEMDegree >& fineSolution , DenseNodeData< Real , FEMDegree >& cumulativeConstraints ) const;
	// Updates the cumulative interpolation constraints @(depth-1) based on the change in solution coefficient @(depth)
	template< int FEMDegree , BoundaryType BType , bool HasGradients >
	void _updateCumulativeInterpolationConstraintsFromFiner( const InterpolationInfo< HasGradients >& interpolationInfo ,
		const BSplineData< FEMDegree , BType >& bsData , LocalDepth highDepth , const DenseNodeData< Real , FEMDegree >& fineSolution , DenseNodeData< Real , FEMDegree >& cumulativeConstraints ) const;

	template< int FEMDegree , BoundaryType BType >
	Real _coarserFunctionValue( Point3D< Real > p , const PointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const BSplineData< FEMDegree , BType >& bsData , const DenseNodeData< Real , FEMDegree >& upSampledCoefficients ) const;
	template< int FEMDegree , BoundaryType BType >
	Point3D< Real > _coarserFunctionGradient( Point3D< Real > p , const PointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const BSplineData< FEMDegree , BType >& bsData , const DenseNodeData< Real , FEMDegree >& upSampledCoefficients ) const;
	template< int FEMDegree , BoundaryType BType >
	Real   _finerFunctionValue( Point3D< Real > p , const PointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const BSplineData< FEMDegree , BType >& bsData , const DenseNodeData< Real , FEMDegree >& coefficients ) const;
	template< int FEMDegree , BoundaryType BType >
	Point3D< Real >   _finerFunctionGradient( Point3D< Real > p , const PointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const BSplineData< FEMDegree , BType >& bsData , const DenseNodeData< Real , FEMDegree >& coefficients ) const;
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	int _getSliceMatrixAndUpdateConstraints( const FEMSystemFunctor& F , const InterpolationInfo< HasGradients >* interpolationInfo , SparseMatrix< Real >& matrix , DenseNodeData< Real , FEMDegree >& constraints , typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template Integrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& integrator , typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template ChildIntegrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& childIntegrator , const BSplineData< FEMDegree , BType >& bsData , LocalDepth depth , int slice , const DenseNodeData< Real , FEMDegree >& metSolution , bool coarseToFine );
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	int _getMatrixAndUpdateConstraints( const FEMSystemFunctor& F , const InterpolationInfo< HasGradients >* interpolationInfo , SparseMatrix< Real >& matrix , DenseNodeData< Real , FEMDegree >& constraints , typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template Integrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& integrator , typename BSplineIntegrationData< FEMDegree , BType , FEMDegree , BType >::FunctionIntegrator::template ChildIntegrator< DERIVATIVES( FEMDegree ) , DERIVATIVES( FEMDegree ) >& childIntegrator , const BSplineData< FEMDegree , BType >& bsData , LocalDepth depth , const DenseNodeData< Real , FEMDegree >& metSolution , bool coarseToFine );

	// Down samples constraints @(depth) to constraints @(depth-1)
	template< class C , int FEMDegree , BoundaryType BType > void _downSample( LocalDepth highDepth , DenseNodeData< C , FEMDegree >& constraints ) const;
	// Up samples coefficients @(depth-1) to coefficients @(depth)
	template< class C , int FEMDegree , BoundaryType BType > void _upSample( LocalDepth highDepth , DenseNodeData< C , FEMDegree >& coefficients ) const;
	template< class C , int FEMDegree , BoundaryType BType > static void _UpSample( LocalDepth highDepth , ConstPointer( C ) lowCoefficients , Pointer( C ) highCoefficients , int threads );
public:
	template< class C , int FEMDegree , BoundaryType BType > DenseNodeData< C , FEMDegree > coarseCoefficients( const  DenseNodeData< C , FEMDegree >& coefficients ) const;
	template< class C , int FEMDegree , BoundaryType BType > DenseNodeData< C , FEMDegree > coarseCoefficients( const SparseNodeData< C , FEMDegree >& coefficients ) const;
protected:

	/////////////////////////////////////////////
	// Code for splatting point-sample data    //
	// MultiGridOctreeData.WeightedSamples.inl //
	/////////////////////////////////////////////
	template< int WeightDegree >
	void _addWeightContribution( DensityEstimator< WeightDegree >& densityWeights , TreeOctNode* node , Point3D< Real > position , PointSupportKey< WeightDegree >& weightKey , Real weight=Real(1.0) );
	template< int WeightDegree , class PointSupportKey >
	Real _getSamplesPerNode( const DensityEstimator< WeightDegree >& densityWeights , const TreeOctNode* node , Point3D< Real > position , PointSupportKey& weightKey ) const;
	template< int WeightDegree , class PointSupportKey >
	void _getSampleDepthAndWeight( const DensityEstimator< WeightDegree >& densityWeights , const TreeOctNode* node , Point3D< Real > position , PointSupportKey& weightKey , Real& depth , Real& weight ) const;
	template< int WeightDegree , class PointSupportKey >
	void _getSampleDepthAndWeight( const DensityEstimator< WeightDegree >& densityWeights , Point3D< Real > position , PointSupportKey& weightKey , Real& depth , Real& weight ) const;
	template< bool CreateNodes ,                    int DataDegree , class V > void      _splatPointData( TreeOctNode* node ,                                           Point3D< Real > point , V v , SparseNodeData< V , DataDegree >& data ,                                              PointSupportKey< DataDegree >& dataKey                                                   );
	template< bool CreateNodes , int WeightDegree , int DataDegree , class V > Real      _splatPointData( const DensityEstimator< WeightDegree >& densityWeights , Point3D< Real > point , V v , SparseNodeData< V , DataDegree >& data , PointSupportKey< WeightDegree >& weightKey , PointSupportKey< DataDegree >& dataKey , LocalDepth minDepth , LocalDepth maxDepth , int dim=DIMENSION );
	template< bool CreateNodes , int WeightDegree , int DataDegree , class V > Real _multiSplatPointData( const DensityEstimator< WeightDegree >* densityWeights , TreeOctNode* node , Point3D< Real > point , V v , SparseNodeData< V , DataDegree >& data , PointSupportKey< WeightDegree >& weightKey , PointSupportKey< DataDegree >& dataKey , int dim=DIMENSION );
	template< class V , int DataDegree , BoundaryType BType , class Coefficients > V _evaluate( const Coefficients& coefficients , Point3D< Real > p , const BSplineData< DataDegree , BType >& bsData , const ConstPointSupportKey< DataDegree >& dataKey ) const;
public:
	template< class V , int DataDegree , BoundaryType BType > Pointer( V ) voxelEvaluate( const DenseNodeData< V , DataDegree >& coefficients , int& res , Real isoValue=0.f , LocalDepth depth=-1 , bool primal=false );

	template< int NormalDegree >
	struct HasNormalDataFunctor
	{
		const SparseNodeData< Point3D< Real > , NormalDegree >& normalInfo;
		HasNormalDataFunctor( const SparseNodeData< Point3D< Real > , NormalDegree >& ni ) : normalInfo( ni ){ ; }
		bool operator() ( const TreeOctNode* node ) const
		{
			const Point3D< Real >* n = normalInfo( node );
			if( n )
			{
				const Point3D< Real >& normal = *n;
				if( normal[0]!=0 || normal[1]!=0 || normal[2]!=0 ) return true;
			}
			if( node->children ) for( int c=0 ; c<Cube::CORNERS ; c++ ) if( (*this)( node->children + c ) ) return true;
			return false;
		}
	};
	struct TrivialHasDataFunctor{ bool operator() ( const TreeOctNode* node ) const{ return true; } };

	// [NOTE] The input/output for this method is pre-scaled by weight
	template< bool HasGradients > bool _setInterpolationInfoFromChildren( TreeOctNode* node , SparseNodeData< PointData< Real , HasGradients > , 0 >& iInfo ) const;
	template< bool HasGradients > SparseNodeData< PointData< Real , HasGradients > , 0 > _densifyInterpolationInfo( const std::vector< PointSample >& samples , Real pointValue , int adaptiveExponent ) const;

	template< int FEMDegree , BoundaryType BType > void _setValidityFlags( void );
	template< class HasDataFunctor > void _clipTree( const HasDataFunctor& f );

	template< int FEMDegree , BoundaryType BType > SparseNodeData<          Real   , 0 > leafValues   ( const DenseNodeData< Real , FEMDegree >& coefficients ) const;
	template< int FEMDegree , BoundaryType BType > SparseNodeData< Point3D< Real > , 0 > leafGradients( const DenseNodeData< Real , FEMDegree >& coefficients ) const;

	////////////////////////////////////
	// Evaluation Methods             //
	// MultiGridOctreeData.Evaluation //
	////////////////////////////////////
	static const int CHILDREN = Cube::CORNERS;
	template< int FEMDegree , BoundaryType BType >
	struct _Evaluator
	{
		typename BSplineEvaluationData< FEMDegree , BType >::Evaluator evaluator;
		typename BSplineEvaluationData< FEMDegree , BType >::ChildEvaluator childEvaluator;
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > cellStencil;
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > cellStencils  [CHILDREN];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > edgeStencil             [Cube::EDGES  ];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > edgeStencils  [CHILDREN][Cube::EDGES  ];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > faceStencil             [Cube::FACES  ];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > faceStencils  [CHILDREN][Cube::FACES  ];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > cornerStencil           [Cube::CORNERS];
		Stencil< double , BSplineSupportSizes< FEMDegree >::SupportSize > cornerStencils[CHILDREN][Cube::CORNERS];

		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dCellStencil;
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dCellStencils  [CHILDREN];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dEdgeStencil             [Cube::EDGES  ];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dEdgeStencils  [CHILDREN][Cube::EDGES  ];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dFaceStencil             [Cube::FACES  ];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dFaceStencils  [CHILDREN][Cube::FACES  ];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dCornerStencil           [Cube::CORNERS];
		Stencil< Point3D< double > , BSplineSupportSizes< FEMDegree >::SupportSize > dCornerStencils[CHILDREN][Cube::CORNERS];

		void set( LocalDepth depth );
		_Evaluator( void ){ _bsData = NULL; }
		~_Evaluator( void ){ if( _bsData ) delete _bsData , _bsData = NULL; }
	protected:
		BSplineData< FEMDegree , BType >* _bsData;
		friend Octree;
	};
	template< class V , int FEMDegree , BoundaryType BType >
	V _getCenterValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node ,                     const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< class V , int FEMDegree , BoundaryType BType >
	V _getCornerValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int corner        , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< class V , int FEMDegree , BoundaryType BType >
	V _getEdgeValue  ( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int edge          , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< class V , int FEMDegree , BoundaryType BType >
	V _getValue      ( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , Point3D< Real > p , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator ) const;

	template< int FEMDegree , BoundaryType BType >
	std::pair< Real , Point3D< Real > > _getCenterValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node ,                     const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< int FEMDegree , BoundaryType BType >
	std::pair< Real , Point3D< Real > > _getCornerValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int corner        , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< int FEMDegree , BoundaryType BType >
	std::pair< Real , Point3D< Real > > _getEdgeValueAndGradient  ( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int edge          , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const;
	template< int FEMDegree , BoundaryType BType >
	std::pair< Real , Point3D< Real > > _getValueAndGradient      ( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , Point3D< Real > p , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator ) const;

public:
	template< int Degree , BoundaryType BType >
	class MultiThreadedEvaluator
	{
		const Octree* _tree;
		int _threads;
		std::vector< ConstPointSupportKey< Degree > > _neighborKeys;
		_Evaluator< Degree , BType > _evaluator;
		const DenseNodeData< Real , Degree >& _coefficients;
		DenseNodeData< Real , Degree > _coarseCoefficients;
	public:
		MultiThreadedEvaluator( const Octree* tree , const DenseNodeData< Real , Degree >& coefficients , int threads=1 );
		Real value( Point3D< Real > p , int thread=0 , const TreeOctNode* node=NULL );
		std::pair< Real , Point3D< Real > > valueAndGradient( Point3D< Real > , int thread=0 , const TreeOctNode* node=NULL );
	};

	////////////////////////////////////////
	// Iso-Surfacing Methods              //
	// MultiGridOctreeData.IsoSurface.inl //
	////////////////////////////////////////
protected:
	struct _IsoEdge
	{
		long long edges[2];
		_IsoEdge( void ){ edges[0] = edges[1] = 0; }
		_IsoEdge( long long v1 , long long v2 ){ edges[0] = v1 , edges[1] = v2; }
		long long& operator[]( int idx ){ return edges[idx]; }
		const long long& operator[]( int idx ) const { return edges[idx]; }
	};
	struct _FaceEdges
	{
		_IsoEdge edges[2];
		int count;
	};
	template< class Vertex >
	struct _SliceValues
	{
		typename SortedTreeNodes::SliceTableData sliceData;
		Pointer( Real ) cornerValues ; Pointer( Point3D< Real > ) cornerGradients ; Pointer( char ) cornerSet;
		Pointer( long long ) edgeKeys ; Pointer( char ) edgeSet;
		Pointer( _FaceEdges ) faceEdges ; Pointer( char ) faceSet;
		Pointer( char ) mcIndices;
		std::unordered_map< long long, std::vector< _IsoEdge > > faceEdgeMap;
		std::unordered_map< long long, std::pair< int, Vertex > > edgeVertexMap;
		std::unordered_map< long long, long long > vertexPairMap;

		_SliceValues( void );
		~_SliceValues( void );
		void reset( bool nonLinearFit );
	protected:
		int _oldCCount , _oldECount , _oldFCount , _oldNCount;
	};
	template< class Vertex >
	struct _XSliceValues
	{
		typename SortedTreeNodes::XSliceTableData xSliceData;
		Pointer( long long ) edgeKeys ; Pointer( char ) edgeSet;
		Pointer( _FaceEdges ) faceEdges ; Pointer( char ) faceSet;
		std::unordered_map< long long, std::vector< _IsoEdge > > faceEdgeMap;
		std::unordered_map< long long, std::pair< int, Vertex > > edgeVertexMap;
		std::unordered_map< long long, long long > vertexPairMap;

		_XSliceValues( void );
		~_XSliceValues( void );
		void reset( void );
	protected:
		int _oldECount , _oldFCount;
	};
	template< class Vertex >
	struct _SlabValues
	{
	protected:
		_XSliceValues< Vertex > _xSliceValues[2];
		_SliceValues< Vertex > _sliceValues[2];
	public:
		_SliceValues< Vertex >& sliceValues( int idx ){ return _sliceValues[idx&1]; }
		const _SliceValues< Vertex >& sliceValues( int idx ) const { return _sliceValues[idx&1]; }
		_XSliceValues< Vertex >& xSliceValues( int idx ){ return _xSliceValues[idx&1]; }
		const _XSliceValues< Vertex >& xSliceValues( int idx ) const { return _xSliceValues[idx&1]; }
	};
	template< class Vertex , int FEMDegree , BoundaryType BType >
	void _setSliceIsoCorners( const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , Real isoValue , LocalDepth depth , int slice ,         std::vector< _SlabValues< Vertex > >& sValues , const _Evaluator< FEMDegree , BType >& evaluator , int threads );
	template< class Vertex , int FEMDegree , BoundaryType BType >
	void _setSliceIsoCorners( const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , Real isoValue , LocalDepth depth , int slice , int z , std::vector< _SlabValues< Vertex > >& sValues , const _Evaluator< FEMDegree , BType >& evaluator , int threads );
	template< int WeightDegree , int ColorDegree , BoundaryType BType , class Vertex >
	void _setSliceIsoVertices( const BSplineData< ColorDegree , BType >* colorBSData , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , Real isoValue , LocalDepth depth , int slice ,         int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues< Vertex > >& sValues , int threads );
	template< int WeightDegree , int ColorDegree , BoundaryType BType , class Vertex >
	void _setSliceIsoVertices( const BSplineData< ColorDegree , BType >* colorBSData , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , Real isoValue , LocalDepth depth , int slice , int z , int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues< Vertex > >& sValues , int threads );
	template< int WeightDegree , int ColorDegree , BoundaryType BType , class Vertex >
	void _setXSliceIsoVertices( const BSplineData< ColorDegree , BType >* colorBSData , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , Real isoValue , LocalDepth depth , int slab , int& vOffset , CoredMeshData< Vertex >& mesh , std::vector< _SlabValues< Vertex > >& sValues , int threads );
	template< class Vertex >
	void _setSliceIsoEdges( LocalDepth depth , int slice ,         std::vector< _SlabValues< Vertex > >& slabValues , int threads );
	template< class Vertex >
	void _setSliceIsoEdges( LocalDepth depth , int slice , int z , std::vector< _SlabValues< Vertex > >& slabValues , int threads );
	template< class Vertex >
	void _setXSliceIsoEdges( LocalDepth depth , int slice , std::vector< _SlabValues< Vertex > >& slabValues , int threads );
	template< class Vertex >
	void _copyFinerSliceIsoEdgeKeys( LocalDepth depth , int slice ,         std::vector< _SlabValues< Vertex > >& sValues , int threads );
	template< class Vertex >
	void _copyFinerSliceIsoEdgeKeys( LocalDepth depth , int slice , int z , std::vector< _SlabValues< Vertex > >& sValues , int threads );
	template< class Vertex >
	void _copyFinerXSliceIsoEdgeKeys( LocalDepth depth , int slab , std::vector< _SlabValues< Vertex > >& sValues , int threads );

	template< class Vertex >
	void _setIsoSurface( LocalDepth depth , int offset , const _SliceValues< Vertex >& bValues , const _SliceValues< Vertex >& fValues , const _XSliceValues< Vertex >& xValues , CoredMeshData< Vertex >& mesh , bool polygonMesh , bool addBarycenter , int& vOffset , int threads );

	template< class Vertex >
	static int _addIsoPolygons( CoredMeshData< Vertex >& mesh , std::vector< std::pair< int , Vertex > >& polygon , bool polygonMesh , bool addBarycenter , int& vOffset );

	template< int WeightDegree , int ColorDegree , BoundaryType BType , class Vertex >
	bool _getIsoVertex( const BSplineData< ColorDegree , BType >* colorBSData , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , Real isoValue , ConstPointSupportKey< WeightDegree >& weightKey , ConstPointSupportKey< ColorDegree >& colorKey , const TreeOctNode* node , int edgeIndex , int z , const _SliceValues< Vertex >& sValues , Vertex& vertex );
	template< int WeightDegree , int ColorDegree , BoundaryType BType , class Vertex >
	bool _getIsoVertex( const BSplineData< ColorDegree , BType >* colorBSData , const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , Real isoValue , ConstPointSupportKey< WeightDegree >& weightKey , ConstPointSupportKey< ColorDegree >& colorKey , const TreeOctNode* node , int cornerIndex , const _SliceValues< Vertex >& bValues , const _SliceValues< Vertex >& fValues , Vertex& vertex );

	void _init( TreeOctNode* node , LocalDepth maxDepth , bool (*Refine)( LocalDepth d , LocalOffset off ) );

	double _maxMemoryUsage , _localMemoryUsage;
public:
	int threads;
	double maxMemoryUsage( void ) const { return _maxMemoryUsage; }
	double localMemoryUsage( void ) const { return _localMemoryUsage; }
	void resetLocalMemoryUsage( void ){ _localMemoryUsage = 0; }
	double memoryUsage( void );

	Octree( void );

	void init( LocalDepth maxDepth , bool (*Refine)( LocalDepth d , LocalOffset off ) );
	template< class Data >
	int init( OrientedPointStream< Real >& pointStream , LocalDepth maxDepth , bool useConfidence , std::vector< PointSample >& samples , std::vector< ProjectiveData< Data , Real > >* sampleData );
	template< int DensityDegree >
	typename Octree::template DensityEstimator< DensityDegree >* setDensityEstimator( const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode );
	template< int NormalDegree , int DensityDegree >
	SparseNodeData< Point3D< Real > , NormalDegree > setNormalField( const std::vector< PointSample >& samples , const DensityEstimator< DensityDegree >& density , Real& pointWeightSum , bool forceNeumann );
	template< int DataDegree , bool CreateNodes , int DensityDegree , class Data >
	SparseNodeData< ProjectiveData< Data , Real > , DataDegree > setDataField( const std::vector< PointSample >& samples , std::vector< ProjectiveData< Data , Real > >& sampleData , const DensityEstimator< DensityDegree >* density );
	template< int MaxDegree , int FEMDegree , BoundaryType FEMBType , class HasDataFunctor > void inalizeForBroodedMultigrid( LocalDepth fullDepth , const HasDataFunctor& F , std::vector< int >* map=NULL );

	// Generate an empty set of constraints
	template< int FEMDegree > DenseNodeData< Real , FEMDegree > initDenseNodeData( void );

	// Add finite-elements constraints (derived from a sparse scalar field)
	template< int FEMDegree , BoundaryType FEMBType , int SFDegree , BoundaryType SFBType , class FEMSFConstraintFunctor > void addFEMConstraints( const FEMSFConstraintFunctor& F , const SparseNodeData< Real , SFDegree >& sfCoefficients , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth )
	{ return _addFEMConstraints< FEMDegree , FEMBType , SFDegree , SFBType , FEMSFConstraintFunctor , const SparseNodeData< Real   , SFDegree > , Real , double >( F , sfCoefficients , constraints , maxDepth ); }
	// Add finite-elements constraints (derived from a dense scalar field)
	template< int FEMDegree , BoundaryType FEMBType , int SFDegree , BoundaryType SFBType , class FEMSFConstraintFunctor > void addFEMConstraints( const FEMSFConstraintFunctor& F , const  DenseNodeData< Real , SFDegree >& sfCoefficients , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth )
	{ return _addFEMConstraints< FEMDegree , FEMBType , SFDegree , SFBType , FEMSFConstraintFunctor , const  DenseNodeData< Real   , SFDegree > , Real , double >( F , sfCoefficients , constraints , maxDepth ); }
	// Add finite-elements constraints (derived from a sparse vector field)
	template< int FEMDegree , BoundaryType FEMBType , int VFDegree , BoundaryType VFBType , class FEMVFConstraintFunctor > void addFEMConstraints( const FEMVFConstraintFunctor& F , const SparseNodeData< Point3D< Real > , VFDegree >& vfCoefficients , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth )
	{ return _addFEMConstraints< FEMDegree , FEMBType , VFDegree , VFBType , FEMVFConstraintFunctor , const SparseNodeData< Point3D< Real > , VFDegree > , Point3D< Real > , Point3D< double > >( F , vfCoefficients , constraints , maxDepth ); }
	// Add finite-elements constraints (derived from a dense vector field)
	template< int FEMDegree , BoundaryType FEMBType , int VFDegree , BoundaryType VFBType , class FEMVFConstraintFunctor > void addFEMConstraints( const FEMVFConstraintFunctor& F , const  DenseNodeData< Point3D< Real > , VFDegree >& vfCoefficients , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth )
	{ return _addFEMConstraints< FEMDegree , FEMBType , VFDegree , VFBType , FEMVFConstraintFunctor , const  DenseNodeData< Point3D< Real > , VFDegree > , Point3D< Real > , Point3D< double > >( F , vfCoefficients , constraints , maxDepth ); }
	// Add interpolation constraints
	template< int FEMDegree , BoundaryType FEMBType , bool HasGradients > void addInterpolationConstraints( const InterpolationInfo< HasGradients >& interpolationInfo , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth );

	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor > double dot( const DotFunctor& F , const SparseNodeData< Real , Degree1 >& coefficients1 , const SparseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , false >( F , (const InterpolationInfo< false >*)NULL , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor > double dot( const DotFunctor& F , const SparseNodeData< Real , Degree1 >& coefficients1 , const DenseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , false >( F , (const InterpolationInfo< false >*)NULL , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor > double dot( const DotFunctor& F , const DenseNodeData< Real , Degree1 >& coefficients1 , const SparseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , false >( F , (const InterpolationInfo< false >*)NULL , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor > double dot( const DotFunctor& F , const DenseNodeData< Real , Degree1 >& coefficients1 , const DenseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , false >( F , (const InterpolationInfo< false >*)NULL , coefficients1 , coefficients2 ); }

	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor , bool HasGradients > double dot( const DotFunctor& F , const InterpolationInfo< HasGradients >* iInfo , const SparseNodeData< Real , Degree1 >& coefficients1 , const SparseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , HasGradients >( F , iInfo , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor , bool HasGradients > double dot( const DotFunctor& F , const InterpolationInfo< HasGradients >* iInfo , const SparseNodeData< Real , Degree1 >& coefficients1 , const DenseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , HasGradients >( F , iInfo , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor , bool HasGradients > double dot( const DotFunctor& F , const InterpolationInfo< HasGradients >* iInfo , const DenseNodeData< Real , Degree1 >& coefficients1 , const SparseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , HasGradients >( F , iInfo , coefficients1 , coefficients2 ); }
	template< int Degree1 , BoundaryType BType1 , int Degree2 , BoundaryType BType2 , class DotFunctor , bool HasGradients > double dot( const DotFunctor& F , const InterpolationInfo< HasGradients >* iInfo , const DenseNodeData< Real , Degree1 >& coefficients1 , const DenseNodeData< Real , Degree2 >& coefficients2 ) const
	{ return _dot< Degree1 , BType1 , Degree2 , BType2 , DotFunctor , HasGradients >( F , iInfo , coefficients1 , coefficients2 ); }

	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	void setSystemMatrix( const FEMSystemFunctor& F , const InterpolationInfo< HasGradients >* interpolationInfo , LocalDepth depth , SparseMatrix< Real >& matrix ) const;

	// Solve the linear system
	struct SolverInfo
	{
		// How to solve
		LocalDepth cgDepth;
		int iters;
		double cgAccuracy , lowResIterMultiplier;
		// What to output
		bool verbose , showResidual;

		SolverInfo( void ) : cgDepth(0) , iters(1), cgAccuracy(0) , lowResIterMultiplier(0) , verbose(false) , showResidual(false) { ; }
	};
	template< int FEMDegree , BoundaryType BType , class FEMSystemFunctor , bool HasGradients >
	DenseNodeData< Real , FEMDegree > solveSystem( const FEMSystemFunctor& F , InterpolationInfo< HasGradients >* iData , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxSolveDepth , const SolverInfo& solverInfo );

	template< int FEMDegree , BoundaryType BType , int WeightDegree , int ColorDegree , class Vertex >
	void getMCIsoSurface( const DensityEstimator< WeightDegree >* densityWeights , const SparseNodeData< ProjectiveData< Point3D< Real > , Real > , ColorDegree >* colorData , const DenseNodeData< Real , FEMDegree >& solution , Real isoValue , CoredMeshData< Vertex >& mesh , bool nonLinearFit=true , bool addBarycenter=false , bool polygonMesh=false );


	const TreeOctNode& tree( void ) const{ return *_tree; }
	size_t leaves( void ) const { return _tree->leaves(); }
	size_t nodes( void ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( IsActiveNode( n ) ) count++ ; return count; }
	size_t ghostNodes( void ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( !IsActiveNode( n ) ) count++ ; return count; }
	inline size_t validSpaceNodes( void ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( isValidSpaceNode( n ) ) count++ ;  return count; }
	inline size_t validSpaceNodes( LocalDepth d ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( _localDepth(n)==d && isValidSpaceNode( n ) ) count++ ; return count; }
	template< int Degree , BoundaryType BType > size_t validFEMNodes( void ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( isValidFEMNode< Degree , BType >( n ) ) count++ ;  return count; }
	template< int Degree , BoundaryType BType > size_t validFEMNodes( LocalDepth d ) const { int count = 0 ; for( const TreeOctNode* n=_tree->nextNode() ; n ; n=_tree->nextNode( n ) ) if( _localDepth(n)==d && isValidFEMNode< Degree , BType >( n ) ) count++ ; return count; }
	LocalDepth depth( void ) const { return _localMaxDepth( _tree ); }
	void resetNodeIndices( void ){ _NodeCount = 0 ; for( TreeOctNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) _NodeInitializer( *node ) , node->nodeData.flags=0; }

protected:
	template< class D > static bool _IsZero( const D& d );
	template< class D > static Real _Dot( const D& d1 , const D& d2 );
	template< int FEMDegree , BoundaryType FEMBType , int CDegree , BoundaryType CBType , class FEMConstraintFunctor , class Coefficients , class D , class _D >
	void _addFEMConstraints( const FEMConstraintFunctor& F , const Coefficients& coefficients , DenseNodeData< Real , FEMDegree >& constraints , LocalDepth maxDepth );
	template< int FEMDegree1 , BoundaryType FEMBType1 , int FEMDegree2 , BoundaryType FEMBType2 , class DotFunctor , bool HasGradients , class Coefficients1 , class Coefficients2 >
	double _dot( const DotFunctor& F , const InterpolationInfo< HasGradients >* iInfo , const Coefficients1& coefficients1 , const Coefficients2& coefficients2 ) const;
};
template< class Real > int Octree< Real >::_NodeCount = 0;


template< class Real > void Reset( void ){ Octree< Real >::ResetNodeCount(); }



template< class Real >
template< int FEMDegree , BoundaryType BType>
void Octree< Real >::_Evaluator< FEMDegree , BType >::set( LocalDepth depth )
{
	static const int  LeftPointSupportRadius =  BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int RightPointSupportRadius = -BSplineSupportSizes< FEMDegree >::SupportStart;

	BSplineEvaluationData< FEMDegree , BType >::SetEvaluator( evaluator , depth );
	if( depth>0 ) BSplineEvaluationData< FEMDegree , BType >::SetChildEvaluator( childEvaluator , depth-1 );
	int center = ( 1<<depth )>>1;

	// First set the stencils for the current depth
	for( int x=-LeftPointSupportRadius ; x<=RightPointSupportRadius ; x++ ) for( int y=-LeftPointSupportRadius ; y<=RightPointSupportRadius ; y++ ) for( int z=-LeftPointSupportRadius ; z<=RightPointSupportRadius ; z++ )
	{
		int fIdx[] = { center+x , center+y , center+z };

		// The cell stencil
		{
			double vv[3] , dv[3];
			for( int dd=0 ; dd<DIMENSION ; dd++ )
			{
				vv[dd] = evaluator.centerValue( fIdx[dd] , center , false );
				dv[dd] = evaluator.centerValue( fIdx[dd] , center , true  );
			}
			cellStencil( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
			dCellStencil( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
		}

		//// The face stencil
		for( int f=0 ; f<Cube::FACES ; f++ )
		{
			int dir , off;
			Cube::FactorFaceIndex( f , dir , off );
			double vv[3] , dv[3];
			switch( dir )
			{
			case 0:
				vv[0] = evaluator.cornerValue( fIdx[0] , center+off , false );
				vv[1] = evaluator.centerValue( fIdx[1] , center     , false );
				vv[2] = evaluator.centerValue( fIdx[2] , center     , false );
				dv[0] = evaluator.cornerValue( fIdx[0] , center+off , true  );
				dv[1] = evaluator.centerValue( fIdx[1] , center     , true  );
				dv[2] = evaluator.centerValue( fIdx[2] , center     , true  );
				break;
			case 1:
				vv[0] = evaluator.centerValue( fIdx[0] , center     , false );
				vv[1] = evaluator.cornerValue( fIdx[1] , center+off , false );
				vv[2] = evaluator.centerValue( fIdx[2] , center     , false );
				dv[0] = evaluator.centerValue( fIdx[0] , center     , true  );
				dv[1] = evaluator.cornerValue( fIdx[1] , center+off , true  );
				dv[2] = evaluator.centerValue( fIdx[2] , center     , true  );
				break;
			case 2:
				vv[0] = evaluator.centerValue( fIdx[0] , center     , false );
				vv[1] = evaluator.centerValue( fIdx[1] , center     , false );
				vv[2] = evaluator.cornerValue( fIdx[2] , center+off , false );
				dv[0] = evaluator.centerValue( fIdx[0] , center     , true  );
				dv[1] = evaluator.centerValue( fIdx[1] , center     , true  );
				dv[2] = evaluator.cornerValue( fIdx[2] , center+off , true  );
				break;
			}
			faceStencil[f]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
			dFaceStencil[f]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
		}

		//// The edge stencil
		for( int e=0 ; e<Cube::EDGES ; e++ )
		{
			int orientation , i1 , i2;
			Cube::FactorEdgeIndex( e , orientation , i1 , i2 );
			double vv[3] , dv[3];
			switch( orientation )
			{
			case 0:
				vv[0] = evaluator.centerValue( fIdx[0] , center    , false );
				vv[1] = evaluator.cornerValue( fIdx[1] , center+i1 , false );
				vv[2] = evaluator.cornerValue( fIdx[2] , center+i2 , false );
				dv[0] = evaluator.centerValue( fIdx[0] , center    , true  );
				dv[1] = evaluator.cornerValue( fIdx[1] , center+i1 , true  );
				dv[2] = evaluator.cornerValue( fIdx[2] , center+i2 , true  );
				break;
			case 1:
				vv[0] = evaluator.cornerValue( fIdx[0] , center+i1 , false );
				vv[1] = evaluator.centerValue( fIdx[1] , center    , false );
				vv[2] = evaluator.cornerValue( fIdx[2] , center+i2 , false );
				dv[0] = evaluator.cornerValue( fIdx[0] , center+i1 , true  );
				dv[1] = evaluator.centerValue( fIdx[1] , center    , true  );
				dv[2] = evaluator.cornerValue( fIdx[2] , center+i2 , true  );
				break;
			case 2:
				vv[0] = evaluator.cornerValue( fIdx[0] , center+i1 , false );
				vv[1] = evaluator.cornerValue( fIdx[1] , center+i2 , false );
				vv[2] = evaluator.centerValue( fIdx[2] , center    , false );
				dv[0] = evaluator.cornerValue( fIdx[0] , center+i1 , true  );
				dv[1] = evaluator.cornerValue( fIdx[1] , center+i2 , true  );
				dv[2] = evaluator.centerValue( fIdx[2] , center    , true  );
				break;
			}
			edgeStencil[e]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
			dEdgeStencil[e]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
		}

		//// The corner stencil
		for( int c=0 ; c<Cube::CORNERS ; c++ )
		{
			int cx , cy  ,cz;
			Cube::FactorCornerIndex( c , cx , cy , cz );
			double vv[3] , dv[3];
			vv[0] = evaluator.cornerValue( fIdx[0] , center+cx , false );
			vv[1] = evaluator.cornerValue( fIdx[1] , center+cy , false );
			vv[2] = evaluator.cornerValue( fIdx[2] , center+cz , false );
			dv[0] = evaluator.cornerValue( fIdx[0] , center+cx , true  );
			dv[1] = evaluator.cornerValue( fIdx[1] , center+cy , true  );
			dv[2] = evaluator.cornerValue( fIdx[2] , center+cz , true  );
			cornerStencil[c]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
			dCornerStencil[c]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
		}
	}

	// Now set the stencils for the parents
	for( int child=0 ; child<CHILDREN ; child++ )
	{
		int childX , childY , childZ;
		Cube::FactorCornerIndex( child , childX , childY , childZ );
		for( int x=-LeftPointSupportRadius ; x<=RightPointSupportRadius ; x++ ) for( int y=-LeftPointSupportRadius ; y<=RightPointSupportRadius ; y++ ) for( int z=-LeftPointSupportRadius ; z<=RightPointSupportRadius ; z++ )
		{
			int fIdx[] = { center/2+x , center/2+y , center/2+z };

			//// The cell stencil
			{
				double vv[3] , dv[3];
				vv[0] = childEvaluator.centerValue( fIdx[0] , center+childX , false );
				vv[1] = childEvaluator.centerValue( fIdx[1] , center+childY , false );
				vv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ , false );
				dv[0] = childEvaluator.centerValue( fIdx[0] , center+childX , true  );
				dv[1] = childEvaluator.centerValue( fIdx[1] , center+childY , true  );
				dv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ , true  );
				cellStencils[child]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
				dCellStencils[child]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
			}

			//// The face stencil
			for( int f=0 ; f<Cube::FACES ; f++ )
			{
				int dir , off;
				Cube::FactorFaceIndex( f , dir , off );
				double vv[3] , dv[3];
				switch( dir )
				{
				case 0:
					vv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+off , false );
					vv[1] = childEvaluator.centerValue( fIdx[1] , center+childY     , false );
					vv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ     , false );
					dv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+off , true  );
					dv[1] = childEvaluator.centerValue( fIdx[1] , center+childY     , true  );
					dv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ     , true  );
					break;
				case 1:
					vv[0] = childEvaluator.centerValue( fIdx[0] , center+childX     , false );
					vv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+off , false );
					vv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ     , false );
					dv[0] = childEvaluator.centerValue( fIdx[0] , center+childX     , true  );
					dv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+off , true  );
					dv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ     , true  );
					break;
				case 2:
					vv[0] = childEvaluator.centerValue( fIdx[0] , center+childX     , false );
					vv[1] = childEvaluator.centerValue( fIdx[1] , center+childY     , false );
					vv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+off , false );
					dv[0] = childEvaluator.centerValue( fIdx[0] , center+childX     , true  );
					dv[1] = childEvaluator.centerValue( fIdx[1] , center+childY     , true  );
					dv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+off , true  );
					break;
				}
				faceStencils[child][f]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
				dFaceStencils[child][f]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
			}

			//// The edge stencil
			for( int e=0 ; e<Cube::EDGES ; e++ )
			{
				int orientation , i1 , i2;
				Cube::FactorEdgeIndex( e , orientation , i1 , i2 );
				double vv[3] , dv[3];
				switch( orientation )
				{
				case 0:
					vv[0] = childEvaluator.centerValue( fIdx[0] , center+childX    , false );
					vv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+i1 , false );
					vv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+i2 , false );
					dv[0] = childEvaluator.centerValue( fIdx[0] , center+childX    , true  );
					dv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+i1 , true  );
					dv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+i2 , true  );
					break;
				case 1:
					vv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+i1 , false );
					vv[1] = childEvaluator.centerValue( fIdx[1] , center+childY    , false );
					vv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+i2 , false );
					dv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+i1 , true  );
					dv[1] = childEvaluator.centerValue( fIdx[1] , center+childY    , true  );
					dv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+i2 , true  );
					break;
				case 2:
					vv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+i1 , false );
					vv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+i2 , false );
					vv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ    , false );
					dv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+i1 , true  );
					dv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+i2 , true  );
					dv[2] = childEvaluator.centerValue( fIdx[2] , center+childZ    , true  );
					break;
				}
				edgeStencils[child][e]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
				dEdgeStencils[child][e]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
			}

			//// The corner stencil
			for( int c=0 ; c<Cube::CORNERS ; c++ )
			{
				int cx , cy  ,cz;
				Cube::FactorCornerIndex( c , cx , cy , cz );
				double vv[3] , dv[3];
				vv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+cx , false );
				vv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+cy , false );
				vv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+cz , false );
				dv[0] = childEvaluator.cornerValue( fIdx[0] , center+childX+cx , true  );
				dv[1] = childEvaluator.cornerValue( fIdx[1] , center+childY+cy , true  );
				dv[2] = childEvaluator.cornerValue( fIdx[2] , center+childZ+cz , true  );
				cornerStencils[child][c]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = vv[0] * vv[1] * vv[2];
				dCornerStencils[child][c]( x+LeftPointSupportRadius , y+LeftPointSupportRadius , z+LeftPointSupportRadius ) = Point3D< double >( dv[0] * vv[1] * vv[2] , vv[0] * dv[1] * vv[2] , vv[0] * vv[1] * dv[2] );
			}
		}
	}
	if( _bsData ) delete _bsData;
	_bsData = new BSplineData< FEMDegree , BType >( depth );
}
template< class Real >
template< class V , int FEMDegree , BoundaryType BType >
V Octree< Real >::_getValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , Point3D< Real > p , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator ) const
{
	static const int SupportSize = BSplineSupportSizes< FEMDegree >::SupportSize;
	static const int  LeftSupportRadius = -BSplineSupportSizes< FEMDegree >::SupportStart;
	static const int RightSupportRadius =  BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int  LeftPointSupportRadius =   BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineSupportSizes< FEMDegree >::SupportStart;

	if( IsActiveNode( node->children ) ) fprintf( stderr , "[WARNING] getValue assumes leaf node\n" );
	V value(0);

	while( GetGhostFlag( node ) )
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* _n = neighbors.neighbors[i][j][k];

			if( _isValidFEMNode( _n ) )
			{
				int _pIdx[3];
				Point3D< Real > _s ; Real _w;
				_startAndWidth( _n , _s , _w );
				int _fIdx[3];
				functionIndex< FEMDegree , BType >( _n , _fIdx );
				for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
				value += 
					solution[ _n->nodeData.nodeIndex ] *
					(Real)
					(
						evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) *
						evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) *
						evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
						);
			}
		}
		node = node->parent;
	}

	LocalDepth d = _localDepth( node );

	for( int dd=0 ; dd<3 ; dd++ )
		if     ( p[dd]==0 ) p[dd] = (Real)(0.+1e-6);
		else if( p[dd]==1 ) p[dd] = (Real)(1.-1e-6);

		{
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* _n = neighbors.neighbors[i][j][k];
				if( _isValidFEMNode( _n ) )
				{
					int _pIdx[3];
					Point3D< Real > _s ; Real _w;
					_startAndWidth( _n , _s , _w );
					int _fIdx[3];
					functionIndex< FEMDegree , BType >( _n , _fIdx );
					for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
					value += 
						solution[ _n->nodeData.nodeIndex ] *
						(Real)
						(
							evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) *
							evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) *
							evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
							);
				}
			}
			if( d>0 )
			{
				const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
				for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
				{
					const TreeOctNode* _n = neighbors.neighbors[i][j][k];
					if( _isValidFEMNode( _n ) )
					{
						int _pIdx[3];
						Point3D< Real > _s ; Real _w;
						_startAndWidth( _n , _s , _w );
						int _fIdx[3];
						functionIndex< FEMDegree , BType >( _n , _fIdx );
						for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
						value += 
							coarseSolution[ _n->nodeData.nodeIndex ] *
							(Real)
							(
								evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) *
								evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) *
								evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
								);
					}
				}
			}
		}
		return value;
}
template< class Real >
template< int FEMDegree , BoundaryType BType >
std::pair< Real , Point3D< Real > > Octree< Real >::_getValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , Point3D< Real > p , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator ) const
{
	static const int SupportSize = BSplineSupportSizes< FEMDegree >::SupportSize;
	static const int  LeftSupportRadius = -BSplineSupportSizes< FEMDegree >::SupportStart;
	static const int RightSupportRadius =  BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int  LeftPointSupportRadius =   BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineSupportSizes< FEMDegree >::SupportStart;

	if( IsActiveNode( node->children ) ) fprintf( stderr , "[WARNING] _getValueAndGradient assumes leaf node\n" );
	Real value(0);
	Point3D< Real > gradient;

	while( GetGhostFlag( node ) )
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* _n = neighbors.neighbors[i][j][k];

			if( _isValidFEMNode( _n ) )
			{
				int _pIdx[3];
				Point3D< Real > _s; Real _w;
				_startAndWidth( _n , _s , _w );
				int _fIdx[3];
				functionIndex< FEMDegree , BType >( _n , _fIdx );
				for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
				value += 
					solution[ _n->nodeData.nodeIndex ] *
					(Real)
					(
						evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
					);
				gradient += 
					Point3D< Real >
					(
						evaluator._bsData->dBaseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
						evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->dBaseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
						evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->dBaseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
					) * solution[ _n->nodeData.nodeIndex ];
			}
		}
		node = node->parent;
	}


	LocalDepth d = _localDepth( node );

	for( int dd=0 ; dd<3 ; dd++ )
		if     ( p[dd]==0 ) p[dd] = (Real)(0.+1e-6);
		else if( p[dd]==1 ) p[dd] = (Real)(1.-1e-6);

		{
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* _n = neighbors.neighbors[i][j][k];

				if( _isValidFEMNode( _n ) )
				{
					int _pIdx[3];
					Point3D< Real > _s ; Real _w;
					_startAndWidth( _n , _s , _w );
					int _fIdx[3];
					functionIndex< FEMDegree , BType >( _n , _fIdx );
					for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
					value += 
						solution[ _n->nodeData.nodeIndex ] *
						(Real)
						(
							evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
						);
					gradient += 
						Point3D< Real >
						(
							evaluator._bsData->dBaseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
							evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->dBaseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
							evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->dBaseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
						) * solution[ _n->nodeData.nodeIndex ];
				}
			}
			if( d>0 )
			{
				const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
				for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
				{
					const TreeOctNode* _n = neighbors.neighbors[i][j][k];

					if( _isValidFEMNode( _n ) )
					{
						int _pIdx[3];
						Point3D< Real > _s ; Real _w;
						_startAndWidth( _n , _s , _w );
						int _fIdx[3];
						functionIndex< FEMDegree , BType >( _n , _fIdx );
						for( int dd=0 ; dd<3 ; dd++ ) _pIdx[dd] = std::max< int >( 0 , std::min< int >( SupportSize-1 , LeftSupportRadius + (int)floor( ( p[dd]-_s[dd] ) / _w ) ) );
						value += 
							coarseSolution[ _n->nodeData.nodeIndex ] *
							(Real)
							(
								evaluator._bsData->baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
							);
						gradient += 
							Point3D< Real >
							(
								evaluator._bsData->dBaseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
								evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData->dBaseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData-> baseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] ) ,
								evaluator._bsData-> baseBSplines[ _fIdx[0] ][ _pIdx[0] ]( p[0] ) * evaluator._bsData-> baseBSplines[ _fIdx[1] ][ _pIdx[1] ]( p[1] ) * evaluator._bsData->dBaseBSplines[ _fIdx[2] ][ _pIdx[2] ]( p[2] )
							) * coarseSolution[ _n->nodeData.nodeIndex ];
					}
				}
			}
		}
		return std::pair< Real , Point3D< Real > >( value , gradient );
}
template< class Real >
template< class V , int FEMDegree , BoundaryType BType >
V Octree< Real >::_getCenterValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineEvaluationData< FEMDegree , BType >::SupportSize;
	static const int  LeftPointSupportRadius =   BSplineEvaluationData< FEMDegree , BType >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineEvaluationData< FEMDegree , BType >::SupportStart;

	if( IsActiveNode( node->children ) ) fprintf( stderr , "[WARNING] getCenterValue assumes leaf node\n" );
	V value(0);
	LocalDepth d = _localDepth( node );

	if( isInterior )
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );
		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* n = neighbors.neighbors[i][j][k];
			if( IsActiveNode( n ) ) value += solution[ n->nodeData.nodeIndex ] * Real( evaluator.cellStencil( i , j , k ) );
		}
		if( d>0 )
		{
			int _corner = int( node - node->parent->children );
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* n = neighbors.neighbors[i][j][k];
				if( IsActiveNode( n ) ) value += coarseSolution[n->nodeData.nodeIndex] * Real( evaluator.cellStencils[_corner]( i , j , k ) );
			}
		}
	}
	else
	{
		LocalOffset cIdx;
		_localDepthAndOffset( node , d , cIdx );
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* n = neighbors.neighbors[i][j][k];

			if( _isValidFEMNode( n ) )
			{
				LocalDepth _d ; LocalOffset fIdx;
				_localDepthAndOffset( n , _d , fIdx );
				value +=
					solution[ n->nodeData.nodeIndex ] *
					Real(
						evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false ) *
						evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false ) *
						evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false )
					);
			}
		}
		if( d>0 )
		{
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* n = neighbors.neighbors[i][j][k];
				if( _isValidFEMNode( n ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( n , _d , fIdx );
					value +=
						coarseSolution[ n->nodeData.nodeIndex ] *
						Real(
							evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false ) *
							evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false ) *
							evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false )
						);
				}
			}
		}
	}
	return value;
}
template< class Real >
template< int FEMDegree , BoundaryType BType >
std::pair< Real , Point3D< Real > > Octree< Real >::_getCenterValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineEvaluationData< FEMDegree , BType >::SupportSize;
	static const int  LeftPointSupportRadius =   BSplineEvaluationData< FEMDegree , BType >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineEvaluationData< FEMDegree , BType >::SupportStart;

	if( IsActiveNode( node->children ) ) fprintf( stderr , "[WARNING] getCenterValueAndGradient assumes leaf node\n" );
	Real value(0);
	Point3D< Real > gradient;
	LocalDepth d = _localDepth( node );

	if( isInterior )
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );
		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* n = neighbors.neighbors[i][j][k];
			if( IsActiveNode( n ) )
			{
				value    +=          Real  ( evaluator. cellStencil( i , j , k ) ) * solution[ n->nodeData.nodeIndex ];
				gradient += Point3D< Real >( evaluator.dCellStencil( i , j , k ) ) * solution[ n->nodeData.nodeIndex ];
			}
		}
		if( d>0 )
		{
			int _corner = int( node - node->parent->children );
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* n = neighbors.neighbors[i][j][k];
				if( IsActiveNode( n ) )
				{
					value    +=          Real  ( evaluator. cellStencils[_corner]( i , j , k ) ) * coarseSolution[n->nodeData.nodeIndex];
					gradient += Point3D< Real >( evaluator.dCellStencils[_corner]( i , j , k ) ) * coarseSolution[n->nodeData.nodeIndex];
				}
			}
		}
	}
	else
	{
		LocalOffset cIdx;
		_localDepthAndOffset( node , d , cIdx );
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );

		for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
		{
			const TreeOctNode* n = neighbors.neighbors[i][j][k];

			if( _isValidFEMNode( n ) )
			{
				LocalDepth _d ; LocalOffset fIdx;
				_localDepthAndOffset( n , _d , fIdx );
				value +=
					Real
					(
						evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false )
					) * solution[ n->nodeData.nodeIndex ];
				gradient += 
					Point3D< Real >
					(
						evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , true  ) * evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false ) ,
						evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , true  ) * evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false ) ,
						evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , true  )
						) * solution[ n->nodeData.nodeIndex ];
			}
		}
		if( d>0 )
		{
			const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
			for( int i=0 ; i<SupportSize ; i++ ) for( int j=0 ; j<SupportSize ; j++ ) for( int k=0 ; k<SupportSize ; k++ )
			{
				const TreeOctNode* n = neighbors.neighbors[i][j][k];
				if( _isValidFEMNode( n ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( n , _d , fIdx );
					value +=
						Real
						(
							evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false )
						) * coarseSolution[ n->nodeData.nodeIndex ];
					gradient +=
						Point3D< Real >
						(
							evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , true  ) * evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false ) ,
							evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , true  ) * evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false ) ,
							evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false ) * evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false ) * evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , true  )
						) * coarseSolution[ n->nodeData.nodeIndex ];
				}
			}
		}
	}
	return std::pair< Real , Point3D< Real > >( value , gradient );
}
template< class Real >
template< class V , int FEMDegree , BoundaryType BType >
V Octree< Real >::_getEdgeValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int edge , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineEvaluationData< FEMDegree , BType >::SupportSize;
	static const int  LeftPointSupportRadius =  BSplineEvaluationData< FEMDegree , BType >::SupportEnd;
	static const int RightPointSupportRadius = -BSplineEvaluationData< FEMDegree , BType >::SupportStart;
	V value(0);
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );
	int startX = 0 , endX = SupportSize , startY = 0 , endY = SupportSize , startZ = 0 , endZ = SupportSize;
	int orientation , i1 , i2;
	Cube::FactorEdgeIndex( edge , orientation , i1 , i2 );
	switch( orientation )
	{
	case 0:
		cIdx[1] += i1 , cIdx[2] += i2;
		if( i1 ) startY++ ; else endY--;
		if( i2 ) startZ++ ; else endZ--;
		break;
	case 1:
		cIdx[0] += i1 , cIdx[2] += i2;
		if( i1 ) startX++ ; else endX--;
		if( i2 ) startZ++ ; else endZ--;
		break;
	case 2:
		cIdx[0] += i1 , cIdx[1] += i2;
		if( i1 ) startX++ ; else endX--;
		if( i2 ) startY++ ; else endY--;
		break;
	}

	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , d );
		for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		{
			const TreeOctNode* _node = neighbors.neighbors[x][y][z];
			if( _isValidFEMNode( _node ) )
			{
				if( isInterior ) value += solution[ _node->nodeData.nodeIndex ] * evaluator.edgeStencil[edge]( x , y , z );
				else
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					switch( orientation )
					{
					case 0:
						value +=
							solution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					case 1:
						value +=
							solution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					case 2:
						value +=
							solution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					}
				}
			}
		}
	}
	if( d>0 )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		// If the corner/child indices don't match, then the sample position is in the interior of the
		// coarser cell and so the full support resolution should be used.
		switch( orientation )
		{
		case 0:
			if( _cy!=i1 ) startY = 0 , endY = SupportSize;
			if( _cz!=i2 ) startZ = 0 , endZ = SupportSize;
			break;
		case 1:
			if( _cx!=i1 ) startX = 0 , endX = SupportSize;
			if( _cz!=i2 ) startZ = 0 , endZ = SupportSize;
			break;
		case 2:
			if( _cx!=i1 ) startX = 0 , endX = SupportSize;
			if( _cy!=i2 ) startY = 0 , endY = SupportSize;
			break;
		}
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
		for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		{
			const TreeOctNode* _node = neighbors.neighbors[x][y][z];
			if( _isValidFEMNode( _node ) )
			{
				if( isInterior ) value += coarseSolution[ _node->nodeData.nodeIndex ] * evaluator.edgeStencils[_corner][edge]( x , y , z );
				else
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					switch( orientation )
					{
					case 0:
						value +=
							coarseSolution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					case 1:
						value +=
							coarseSolution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					case 2:
						value +=
							coarseSolution[ _node->nodeData.nodeIndex ] *
							Real(
								evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
								evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
								evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false )
							);
						break;
					}
				}
			}
		}
	}
	return Real( value );
}
template< class Real >
template< int FEMDegree , BoundaryType BType >
std::pair< Real , Point3D< Real > > Octree< Real >::_getEdgeValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int edge , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineEvaluationData< FEMDegree , BType >::SupportSize;
	static const int  LeftPointSupportRadius =  BSplineEvaluationData< FEMDegree , BType >::SupportEnd;
	static const int RightPointSupportRadius = -BSplineEvaluationData< FEMDegree , BType >::SupportStart;
	double value = 0;
	Point3D< double > gradient;
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	int startX = 0 , endX = SupportSize , startY = 0 , endY = SupportSize , startZ = 0 , endZ = SupportSize;
	int orientation , i1 , i2;
	Cube::FactorEdgeIndex( edge , orientation , i1 , i2 );
	switch( orientation )
	{
	case 0:
		cIdx[1] += i1 , cIdx[2] += i2;
		if( i1 ) startY++ ; else endY--;
		if( i2 ) startZ++ ; else endZ--;
		break;
	case 1:
		cIdx[0] += i1 , cIdx[2] += i2;
		if( i1 ) startX++ ; else endX--;
		if( i2 ) startZ++ ; else endZ--;
		break;
	case 2:
		cIdx[0] += i1 , cIdx[1] += i2;
		if( i1 ) startX++ ; else endX--;
		if( i2 ) startY++ ; else endY--;
		break;
	}
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );
		for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		{
			const TreeOctNode* _node = neighbors.neighbors[x][y][z];
			if( _isValidFEMNode( _node ) )
			{
				if( isInterior )
				{
					value    += evaluator. edgeStencil[edge]( x , y , z ) * solution[ _node->nodeData.nodeIndex ];
					gradient += evaluator.dEdgeStencil[edge]( x , y , z ) * solution[ _node->nodeData.nodeIndex ];
				}
				else
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );

					double vv[3] , dv[3];
					switch( orientation )
					{
					case 0:
						vv[0] = evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.evaluator.centerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , true  );
						break;
					case 1:
						vv[0] = evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.evaluator.centerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , true  );
						break;
					case 2:
						vv[0] = evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.evaluator.centerValue( fIdx[2] , cIdx[2] , true  );
						break;
					}
					value += solution[ _node->nodeData.nodeIndex ] * vv[0] * vv[1] * vv[2];
					gradient += Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
				}
			}
		}
	}
	if( d>0 )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		// If the corner/child indices don't match, then the sample position is in the interior of the
		// coarser cell and so the full support resolution should be used.
		switch( orientation )
		{
		case 0:
			if( _cy!=i1 ) startY = 0 , endY = SupportSize;
			if( _cz!=i2 ) startZ = 0 , endZ = SupportSize;
			break;
		case 1:
			if( _cx!=i1 ) startX = 0 , endX = SupportSize;
			if( _cz!=i2 ) startZ = 0 , endZ = SupportSize;
			break;
		case 2:
			if( _cx!=i1 ) startX = 0 , endX = SupportSize;
			if( _cy!=i2 ) startY = 0 , endY = SupportSize;
			break;
		}
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
		for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
		{
			const TreeOctNode* _node = neighbors.neighbors[x][y][z];
			if( _isValidFEMNode( _node ) )
			{
				if( isInterior )
				{
					value    += evaluator. edgeStencils[_corner][edge]( x , y , z ) * coarseSolution[ _node->nodeData.nodeIndex ];
					gradient += evaluator.dEdgeStencils[_corner][edge]( x , y , z ) * coarseSolution[ _node->nodeData.nodeIndex ];
				}
				else
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					double vv[3] , dv[3];
					switch( orientation )
					{
					case 0:
						vv[0] = evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.childEvaluator.centerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , true  );
						break;
					case 1:
						vv[0] = evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.childEvaluator.centerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , true  );
						break;
					case 2:
						vv[0] = evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false );
						vv[1] = evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false );
						vv[2] = evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , false );
						dv[0] = evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , true  );
						dv[1] = evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , true  );
						dv[2] = evaluator.childEvaluator.centerValue( fIdx[2] , cIdx[2] , true  );
						break;
					}
					value += coarseSolution[ _node->nodeData.nodeIndex ] * vv[0] * vv[1] * vv[2];
					gradient += Point3D< double >( dv[0]*vv[1]*vv[2] , vv[0]*dv[1]*vv[2] , vv[0]*vv[1]*dv[2] ) * coarseSolution[ _node->nodeData.nodeIndex ];
				}
			}
		}
	}
	return std::pair< Real , Point3D< Real > >( Real( value ) , Point3D< Real >( gradient ) );
}

template< class Real >
template< class V , int FEMDegree , BoundaryType BType >
V Octree< Real >::_getCornerValue( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int corner , const DenseNodeData< V , FEMDegree >& solution , const DenseNodeData< V , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineSupportSizes< FEMDegree >::SupportSize;
	static const int  LeftPointSupportRadius =   BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineSupportSizes< FEMDegree >::SupportStart;

	V value(0);
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	int cx , cy , cz;
	int startX = 0 , endX = SupportSize , startY = 0 , endY = SupportSize , startZ = 0 , endZ = SupportSize;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	cIdx[0] += cx , cIdx[1] += cy , cIdx[2] += cz;
	{
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );
		if( cx==0 ) endX--;
		else      startX++;
		if( cy==0 ) endY--;
		else      startY++;
		if( cz==0 ) endZ--;
		else      startZ++;
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( IsActiveNode( _node ) ) value += solution[ _node->nodeData.nodeIndex ] * Real( evaluator.cornerStencil[corner]( x , y , z ) );
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _isValidFEMNode( _node ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					value +=
						solution[ _node->nodeData.nodeIndex ] *
						Real(
							evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
							evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
							evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false )
						);
				}
			}
	}
	if( d>0 )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		// If the corner/child indices don't match, then the sample position is in the interior of the
		// coarser cell and so the full support resolution should be used.
		if( cx!=_cx ) startX = 0 , endX = SupportSize;
		if( cy!=_cy ) startY = 0 , endY = SupportSize;
		if( cz!=_cz ) startZ = 0 , endZ = SupportSize;
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( IsActiveNode( _node ) ) value += coarseSolution[ _node->nodeData.nodeIndex ] * Real( evaluator.cornerStencils[_corner][corner]( x , y , z ) );
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _isValidFEMNode( _node ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					value +=
						coarseSolution[ _node->nodeData.nodeIndex ] *
						Real(
							evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false ) *
							evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false ) *
							evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false )
						);
				}
			}
	}
	return Real( value );
}
template< class Real >
template< int FEMDegree , BoundaryType BType >
std::pair< Real , Point3D< Real > > Octree< Real >::_getCornerValueAndGradient( const ConstPointSupportKey< FEMDegree >& neighborKey , const TreeOctNode* node , int corner , const DenseNodeData< Real , FEMDegree >& solution , const DenseNodeData< Real , FEMDegree >& coarseSolution , const _Evaluator< FEMDegree , BType >& evaluator , bool isInterior ) const
{
	static const int SupportSize = BSplineSupportSizes< FEMDegree >::SupportSize;
	static const int  LeftPointSupportRadius =   BSplineSupportSizes< FEMDegree >::SupportEnd;
	static const int RightPointSupportRadius = - BSplineSupportSizes< FEMDegree >::SupportStart;

	double value = 0;
	Point3D< double > gradient;
	LocalDepth d ; LocalOffset cIdx;
	_localDepthAndOffset( node , d , cIdx );

	int cx , cy , cz;
	int startX = 0 , endX = SupportSize , startY = 0 , endY = SupportSize , startZ = 0 , endZ = SupportSize;
	Cube::FactorCornerIndex( corner , cx , cy , cz );
	cIdx[0] += cx , cIdx[1] += cy , cIdx[2] += cz;
	{
		if( cx==0 ) endX--;
		else      startX++;
		if( cy==0 ) endY--;
		else      startY++;
		if( cz==0 ) endZ--;
		else      startZ++;
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node );
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( IsActiveNode( _node ) ) value += solution[ _node->nodeData.nodeIndex ] * evaluator.cornerStencil[corner]( x , y , z ) , gradient += evaluator.dCornerStencil[corner]( x , y , z ) * solution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _isValidFEMNode( _node ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					double v [] = { evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , false ) , evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , false ) , evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , false ) };
					double dv[] = { evaluator.evaluator.cornerValue( fIdx[0] , cIdx[0] , true  ) , evaluator.evaluator.cornerValue( fIdx[1] , cIdx[1] , true  ) , evaluator.evaluator.cornerValue( fIdx[2] , cIdx[2] , true  ) };
					value += solution[ _node->nodeData.nodeIndex ] * v[0] * v[1] * v[2];
					gradient += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * solution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	if( d>0 )
	{
		int _corner = int( node - node->parent->children );
		int _cx , _cy , _cz;
		Cube::FactorCornerIndex( _corner , _cx , _cy , _cz );
		if( cx!=_cx ) startX = 0 , endX = SupportSize;
		if( cy!=_cy ) startY = 0 , endY = SupportSize;
		if( cz!=_cz ) startZ = 0 , endZ = SupportSize;
		const typename TreeOctNode::ConstNeighbors< SupportSize >& neighbors = _neighbors< LeftPointSupportRadius , RightPointSupportRadius >( neighborKey , node->parent );
		if( isInterior )
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node=neighbors.neighbors[x][y][z];
				if( IsActiveNode( _node ) ) value += coarseSolution[ _node->nodeData.nodeIndex ] * evaluator.cornerStencils[_corner][corner]( x , y , z ) , gradient += evaluator.dCornerStencils[_corner][corner]( x , y , z ) * coarseSolution[ _node->nodeData.nodeIndex ];
			}
		else
			for( int x=startX ; x<endX ; x++ ) for( int y=startY ; y<endY ; y++ ) for( int z=startZ ; z<endZ ; z++ )
			{
				const TreeOctNode* _node = neighbors.neighbors[x][y][z];
				if( _isValidFEMNode( _node ) )
				{
					LocalDepth _d ; LocalOffset fIdx;
					_localDepthAndOffset( _node , _d , fIdx );
					double v [] = { evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , false ) , evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , false ) , evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , false ) };
					double dv[] = { evaluator.childEvaluator.cornerValue( fIdx[0] , cIdx[0] , true  ) , evaluator.childEvaluator.cornerValue( fIdx[1] , cIdx[1] , true  ) , evaluator.childEvaluator.cornerValue( fIdx[2] , cIdx[2] , true  ) };
					value += coarseSolution[ _node->nodeData.nodeIndex ] * v[0] * v[1] * v[2];
					gradient += Point3D< double >( dv[0]*v[1]*v[2] , v[0]*dv[1]*v[2] , v[0]*v[1]*dv[2] ) * coarseSolution[ _node->nodeData.nodeIndex ];
				}
			}
	}
	return std::pair< Real , Point3D< Real > >( Real( value ) , Point3D< Real >( gradient ) );
}
template< class Real >
template< int Degree , BoundaryType BType >
Octree< Real >::MultiThreadedEvaluator< Degree , BType >::MultiThreadedEvaluator( const Octree< Real >* tree , const DenseNodeData< Real , Degree >& coefficients , int threads ) : _coefficients( coefficients ) , _tree( tree )
{
	_threads = std::max< int >( 1 , threads );
	_neighborKeys.resize( _threads );
	_coarseCoefficients = _tree->template coarseCoefficients< Real , Degree , BType >( _coefficients );
	_evaluator.set( _tree->_maxDepth );
	for( int t=0 ; t<_threads ; t++ ) _neighborKeys[t].set( tree->_localToGlobal( _tree->_maxDepth ) );
}
template< class Real >
template< int Degree , BoundaryType BType >
Real Octree< Real >::MultiThreadedEvaluator< Degree , BType >::value( Point3D< Real > p , int thread , const TreeOctNode* node )
{
	if( !node ) node = _tree->leaf( p );
	ConstPointSupportKey< Degree >& nKey = _neighborKeys[thread];
	nKey.getNeighbors( node );
	return _tree->template _getValue< Real , Degree >( nKey , node , p , _coefficients , _coarseCoefficients , _evaluator );
}
template< class Real >
template< int Degree , BoundaryType BType >
std::pair< Real , Point3D< Real > > Octree< Real >::MultiThreadedEvaluator< Degree , BType >::valueAndGradient( Point3D< Real > p , int thread , const TreeOctNode* node )
{
	if( !node ) node = _tree->leaf( p );
	ConstPointSupportKey< Degree >& nKey = _neighborKeys[thread];
	nKey.getNeighbors( node );
	return _tree->template _getValueAndGradient< Degree >( nKey , node , p , _coefficients , _coarseCoefficients , _evaluator );
}

#ifdef FAST_SET_UP
#include <functional>
#endif // FAST_SET_UP
#include <cmath>
#include "PointStream.h"

#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
//#define MEMORY_ALLOCATOR_BLOCK_SIZE 0

const double MATRIX_ENTRY_EPSILON = 0;
const double EPSILON              = 1e-6;
const double ROUND_EPS            = 1e-5;

//////////////////
// TreeNodeData //
//////////////////
TreeNodeData::TreeNodeData( void ){ flags = 0; }
TreeNodeData::~TreeNodeData( void ) { }


////////////
// Octree //
////////////
template< class Real >
double Octree< Real >::memoryUsage( void )
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	_maxMemoryUsage = std::max< double >( mem , _maxMemoryUsage );
	_localMemoryUsage = std::max< double >( mem , _localMemoryUsage );
	return mem;
}

template< class Real > Octree< Real >::Octree( void ) : threads(1) , _maxMemoryUsage(0) , _localMemoryUsage(0)
{
	_tree = TreeOctNode::NewBrood( _NodeInitializer );
	_tree->initChildren( _NodeInitializer ) , _spaceRoot = _tree->children;
	_depthOffset = 1;
}

template< class Real >
template< int FEMDegree , BoundaryType BType >
void Octree< Real >::functionIndex( const TreeOctNode* node , int idx[3] ) const
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	for( int dd=0 ; dd<DIMENSION ; dd++ ) idx[dd] = BSplineData< FEMDegree , BType >::FunctionIndex( d , off[dd] );
}

template< class Real >
OctNode< TreeNodeData >* Octree< Real >::leaf( Point3D< Real > p )
{
	if( !_InBounds( p ) ) return NULL;
	Point3D< Real > center = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
	Real width = Real(1.0);
	TreeOctNode* node = _spaceRoot;
	while( node->children )
	{
		int cIndex = TreeOctNode::CornerIndex( center , p );
		node = node->children + cIndex;
		width /= 2;
		if( cIndex&1 ) center[0] += width/2;
		else           center[0] -= width/2;
		if( cIndex&2 ) center[1] += width/2;
		else           center[1] -= width/2;
		if( cIndex&4 ) center[2] += width/2;
		else           center[2] -= width/2;
	}
	return node;
}
template< class Real >
const OctNode< TreeNodeData >* Octree< Real >::leaf( Point3D< Real > p ) const
{
	if( !_InBounds( p ) ) return NULL;
	Point3D< Real > center = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
	Real width = Real(1.0);
	TreeOctNode* node = _spaceRoot;
	while( node->children )
	{
		int cIndex = TreeOctNode::CornerIndex( center , p );
		node = node->children + cIndex;
		width /= 2;
		if( cIndex&1 ) center[0] += width/2;
		else           center[0] -= width/2;
		if( cIndex&2 ) center[1] += width/2;
		else           center[1] -= width/2;
		if( cIndex&4 ) center[2] += width/2;
		else           center[2] -= width/2;
	}
	return node;
}
template< class Real > bool Octree< Real >::_InBounds( Point3D< Real > p ){ return p[0]>=Real(0.) && p[0]<=Real(1.0) && p[1]>=Real(0.) && p[1]<=Real(1.0) && p[2]>=Real(0.) && p[2]<=Real(1.0); }
template< class Real >
template< int FEMDegree , BoundaryType BType >
bool Octree< Real >::isValidFEMNode( const TreeOctNode* node ) const
{
	if( GetGhostFlag( node ) ) return false;
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	return !BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[0] ) && !BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[1] ) && !BSplineEvaluationData< FEMDegree , BType >::OutOfBounds( d , off[2] );
}
template< class Real >
bool Octree< Real >::isValidSpaceNode( const TreeOctNode* node ) const
{
	if( !node ) return false;
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	int res = 1<<d;
	return off[0]>=0 && off[0]<res && off[1]>=0 && off[1]<res && off[2]>=0 && off[2]<res;
}
template< class Real >
template< int Degree , BoundaryType BType >
void Octree< Real >::_setFullDepth( TreeOctNode* node , LocalDepth depth ) const
{
	bool refine = false;
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	if( d<depth )
		if( d<0 ) refine = true;
		else if( BType==BOUNDARY_FREE && !_outOfBounds< Degree , BType >( node ) ) refine = true;
		else if( !BSplineSupportSizes< Degree >::OutOfBounds( d , off[0] ) && !BSplineSupportSizes< Degree >::OutOfBounds( d , off[1] ) && !BSplineSupportSizes< Degree >::OutOfBounds( d , off[2] ) ) refine = true;
	if( refine )
	{
		if( !node->children ) node->initChildren( _NodeInitializer );
		for( int c=0 ; c<Cube::CORNERS ; c++ ) _setFullDepth< Degree , BType >( node->children+c , depth );
	}
}
template< class Real >
template< int Degree , BoundaryType BType >
void Octree< Real >::_setFullDepth( LocalDepth depth )
{
	if( !_tree->children ) _tree->initChildren( _NodeInitializer );
	for( int c=0 ; c<Cube::CORNERS ; c++ ) _setFullDepth< Degree , BType >( _tree->children+c , depth );
}

template< class Real , bool HasGradients >
struct _PointDataAccumulator_
{
#if POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , HasGradients >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Point3D< Real > center , Real width , Real weight );
#else // !POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , HasGradients >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Real weight );
#endif // POINT_DATA_RES
};
template< class Real >
struct _PointDataAccumulator_< Real , false >
{
#if POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , false >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Point3D< Real > center , Real width , Real weight ){ pData.addPoint( SinglePointData< Real , false >( position , value , weight ) , center , width ); }
#else // !POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , false >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Real weight ){ pData.position += position , pData.value += value , pData.weight += weight; }
#endif // POINT_DATA_RES
};
template< class Real >
struct _PointDataAccumulator_< Real , true >
{
#if POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , true >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Point3D< Real > center , Real width , Real weight ){ pData.addPoint( SinglePointData< Real , true >( position , value , gradient , weight ) , center , width ); }
#else // !POINT_DATA_RES
	static inline void _AddToPointData_( PointData< Real , true >& pData , Point3D< Real > position , Real value , Point3D< Real > gradient , Real weight ){ pData.position += position , pData.value += value , pData.gradient += gradient , pData.weight += weight; }
#endif // POINT_DATA_RES
};

template< class Real >
void Octree< Real >::_init( TreeOctNode* node , LocalDepth maxDepth , bool (*Refine)( LocalDepth , LocalOffset ) )
{
	if( _localDepth( node )<maxDepth )
	{
		LocalDepth d ; LocalOffset off;
		_localDepthAndOffset( node , d , off );
		if( Refine( d , off ) )
		{
			node->initChildren( _NodeInitializer );
			for( int c=0 ; c<Cube::CORNERS ; c++ ) _init( node->children + c , maxDepth , Refine );
		}
	}
}
template< class Real > void Octree< Real >::init( LocalDepth maxDepth , bool (*Refine)( LocalDepth , LocalOffset ) ){ _init( _spaceRoot , maxDepth , Refine ); }
template< class Real >
template< class Data >
int Octree< Real >::init( OrientedPointStream< Real >& pointStream , LocalDepth maxDepth , bool useConfidence , std::vector< PointSample >& samples , std::vector< ProjectiveData< Data , Real > >* sampleData )
{
	OrientedPointStreamWithData< Real , Data >& pointStreamWithData = ( OrientedPointStreamWithData< Real , Data >& )pointStream;

	// Add the point data
	int outOfBoundPoints = 0 , zeroLengthNormals = 0 , undefinedNormals = 0 , pointCount = 0;
	{
		std::vector< int > nodeToIndexMap;
		Point3D< Real > p , n;
		OrientedPoint3D< Real > _p;
		Data _d;
		while( ( sampleData ? pointStreamWithData.nextPoint( _p , _d ) : pointStream.nextPoint( _p ) ) )
		{
			p = Point3D< Real >(_p.p) , n = Point3D< Real >(_p.n);
			Real len = (Real)Length( n );
			if( !_InBounds(p) ){ outOfBoundPoints++ ; continue; }
			if( !len ){ zeroLengthNormals++ ; continue; }
			if( len!=len ){ undefinedNormals++ ; continue; }
			n /= len;
			Point3D< Real > center = Point3D< Real >( Real(0.5) , Real(0.5) , Real(0.5) );
			Real width = Real(1.0);
			TreeOctNode* temp = _spaceRoot;
			LocalDepth depth = _localDepth( temp );
			while( depth<maxDepth )
			{
				if( !temp->children ) temp->initChildren( _NodeInitializer );
				int cIndex = TreeOctNode::CornerIndex( center , p );
				temp = temp->children + cIndex;
				width /= 2;
				if( cIndex&1 ) center[0] += width/2;
				else           center[0] -= width/2;
				if( cIndex&2 ) center[1] += width/2;
				else           center[1] -= width/2;
				if( cIndex&4 ) center[2] += width/2;
				else           center[2] -= width/2;
				depth++;
			}
			Real weight = (Real)( useConfidence ? len : 1. );
			int nodeIndex = temp->nodeData.nodeIndex;
			if( nodeIndex>=nodeToIndexMap.size() ) nodeToIndexMap.resize( nodeIndex+1 , -1 );
			int idx = nodeToIndexMap[ nodeIndex ];
			if( idx==-1 )
			{
				idx = (int)samples.size();
				nodeToIndexMap[ nodeIndex ] = idx;
				samples.resize( idx+1 ) , samples[idx].node = temp;
				if( sampleData ) sampleData->resize( idx+1 );
			}
			samples[idx].sample += ProjectiveData< OrientedPoint3D< Real > , Real >( OrientedPoint3D< Real >( p * weight , n * weight ) , weight );
			if( sampleData ) (*sampleData)[ idx ] += ProjectiveData< Data , Real >( _d * weight , weight );
			pointCount++;
		}
		pointStream.reset();
	}
	if( outOfBoundPoints  ) fprintf( stderr , "[WARNING] Found out-of-bound points: %d\n" , outOfBoundPoints );
	if( zeroLengthNormals ) fprintf( stderr , "[WARNING] Found zero-length normals: %d\n" , zeroLengthNormals );
	if( undefinedNormals  ) fprintf( stderr , "[WARNING] Found undefined normals: %d\n" , undefinedNormals );

	memoryUsage();
	return pointCount;
}
template< class Real >
template< int DensityDegree >
typename Octree< Real >::template DensityEstimator< DensityDegree >* Octree< Real >::setDensityEstimator( const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode )
{
	LocalDepth maxDepth = _localMaxDepth( _tree );
	splatDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( splatDepth , maxDepth ) );
	DensityEstimator< DensityDegree >* _density = new DensityEstimator< DensityDegree >( splatDepth );
	DensityEstimator< DensityDegree >& density = *_density;
	PointSupportKey< DensityDegree > densityKey;
	densityKey.set( _localToGlobal( splatDepth ) );

#ifdef FAST_SET_UP
	std::vector< int > sampleMap( NodeCount() , -1 );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<samples.size() ; i++ ) if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = i;
	std::function< ProjectiveData< OrientedPoint3D< Real > , Real > ( TreeOctNode* ) > SetDensity = [&] ( TreeOctNode* node )
	{
		ProjectiveData< OrientedPoint3D< Real > , Real > sample;
		LocalDepth d = _localDepth( node );
		int idx = node->nodeData.nodeIndex;
		if( node->children )
			for( int c=0 ; c<Cube::CORNERS ; c++ )
			{
				ProjectiveData< OrientedPoint3D< Real > , Real > s = SetDensity( node->children + c );
				if( d<=splatDepth && s.weight>0 )
				{
					Point3D< Real > p = s.data.p / s.weight;
					Real w = s.weight / samplesPerNode;
					_addWeightContribution( density , node , p , densityKey , w );
				}
				sample += s;
			}
		else if( idx<sampleMap.size() && sampleMap[idx]!=-1 )
		{
			sample = samples[ sampleMap[ idx ] ].sample;
			if( d<=splatDepth && sample.weight>0 )
			{
				Point3D< Real > p = sample.data.p / sample.weight;
				Real w = sample.weight / samplesPerNode;
				_addWeightContribution( density , node , p , densityKey , w );
			}
		}
		return sample;
	};
	SetDensity( _spaceRoot );
#else // !FAST_SET_UP
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const TreeOctNode* node = samples[i].node;
		const ProjectiveData< OrientedPoint3D< Real > , Real >& sample = samples[i].sample;
		if( sample.weight>0 )
		{
			Point3D< Real > p = sample.data.p / sample.weight;
			Real w = sample.weight / samplesPerNode;
			for( TreeOctNode* _node=(TreeOctNode*)node ; _node ; _node=_node->parent ) if( _localDepth( _node )<=splatDepth ) _addWeightContribution( density , _node , p , densityKey , w );
		}
	}
#endif // FAST_SET_UP

	memoryUsage();
	return _density;
}
template< class Real >
template< int NormalDegree , int DensityDegree >
SparseNodeData< Point3D< Real > , NormalDegree > Octree< Real >::setNormalField( const std::vector< PointSample >& samples , const DensityEstimator< DensityDegree >& density , Real& pointWeightSum , bool forceNeumann )
{
	LocalDepth maxDepth = _localMaxDepth( _tree );
	PointSupportKey< DensityDegree > densityKey;
	PointSupportKey< NormalDegree > normalKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , normalKey.set( _localToGlobal( maxDepth ) );

	Real weightSum = 0;
	pointWeightSum = 0;
	SparseNodeData< Point3D< Real > , NormalDegree > normalField;
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const ProjectiveData< OrientedPoint3D< Real > , Real >& sample = samples[i].sample;
		if( sample.weight>0 )
		{
			Point3D< Real > p = sample.data.p / sample.weight , n = sample.data.n;
			weightSum += sample.weight;
			if( !_InBounds(p) ){ fprintf( stderr , "[WARNING] Octree:setNormalField: Point sample is out of bounds\n" ) ; continue; }
			pointWeightSum += _splatPointData< true >( density , p , n , normalField , densityKey , normalKey , 0 , maxDepth , 3 );
		}
	}
	pointWeightSum /= weightSum;
	memoryUsage();

	return normalField;
}
template< class Real >
template< int DataDegree , bool CreateNodes , int DensityDegree , class Data >
SparseNodeData< ProjectiveData< Data , Real > , DataDegree > Octree< Real >::setDataField( const std::vector< PointSample >& samples , std::vector< ProjectiveData< Data , Real > >& sampleData , const DensityEstimator< DensityDegree >* density )
{
	LocalDepth maxDepth = _localMaxDepth( _tree );
	PointSupportKey< DensityDegree > densityKey;
	PointSupportKey< DataDegree > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< ProjectiveData< Data , Real > , DataDegree > dataField;
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const ProjectiveData< OrientedPoint3D< Real > , Real >& sample = samples[i].sample;
		const ProjectiveData< Data , Real >& data = sampleData[i];
		Point3D< Real > p = sample.weight==0 ? sample.data.p : sample.data.p / sample.weight;
		if( !_InBounds(p) ){ fprintf( stderr , "[WARNING] Point is out of bounds: %f %f %f <- %f %f %f [%f]\n" , p[0] , p[1] , p[2] , sample.data.p[0] , sample.data.p[1] , sample.data.p[2] , sample.weight ) ; continue; }
		_multiSplatPointData< CreateNodes >( density , (TreeOctNode*)samples[i].node , p , data , dataField , densityKey , dataKey , 2 );
	}
	memoryUsage();
	return dataField;
}
template< class Real >
template< int MaxDegree , int FEMDegree , BoundaryType FEMBType , class HasDataFunctor >
void Octree< Real >::inalizeForBroodedMultigrid( LocalDepth fullDepth , const HasDataFunctor& F , std::vector< int >* map )
{
	if( FEMDegree>MaxDegree ) fprintf( stderr , "[ERROR] MaxDegree must be at least as large as the FEM degree: %d <= %d\n" , FEMDegree , MaxDegree );
	while( _localInset( 0 ) + BSplineEvaluationData< MaxDegree , BOUNDARY_FREE >::Begin( 0 )<0 || _localInset( 0 ) + BSplineEvaluationData< MaxDegree , BOUNDARY_FREE >::End( 0 )>(1<<_depthOffset) )
	{
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//          +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//          | | | | |    | | | | | | | | |
		// +-+-+    +-+-+-+-+    +-+-+-+-+-+-+-+-+
		// |*| |    | | | | |    | | | | | | | | |
		// +-o-+ -> +-+-o-+-+ -> +-+-+-+-o-+-+-+-+
		// | | |    | | |*| |    | | | | |*| | | |
		// +-+-+    +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//          | | | | |    | | | | | | | | |
		//          +-+-+-+-+    +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+
		//                       | | | | | | | | |
		//                       +-+-+-+-+-+-+-+-+

		TreeOctNode* newSpaceRootParent = TreeOctNode::NewBrood( _NodeInitializer );
		TreeOctNode* oldSpaceRootParent = _spaceRoot->parent;
		int corner = _depthOffset<=1 ? Cube::CORNERS-1 : 0;
		newSpaceRootParent[corner].children = _spaceRoot;
		oldSpaceRootParent->children = newSpaceRootParent;
		for( int c=0 ; c<Cube::CORNERS ; c++ ) _spaceRoot[c].parent = newSpaceRootParent + corner , newSpaceRootParent[c].parent = oldSpaceRootParent;
		_depthOffset++;
	}
	int d=0 , off[] = { 0 , 0 , 0 };
	TreeOctNode::ResetDepthAndOffset( _tree , d , off );
	_maxDepth = _localMaxDepth( _tree );

	// Make the low-resolution part of the tree be complete
	_fullDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( _maxDepth , fullDepth ) );
	_setFullDepth< MaxDegree , BOUNDARY_FREE >( _fullDepth );
	// Clear all the flags and make everything that is not low-res a ghost node
	for( TreeOctNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) node->nodeData.flags = 0 , SetGhostFlag( node , _localDepth( node )>_fullDepth );

	// Set the ghost nodes for the high-res part of the tree
	_clipTree( F );

	const int OverlapRadius = -BSplineOverlapSizes< MaxDegree , MaxDegree >::OverlapStart;
	typename TreeOctNode::NeighborKey< OverlapRadius , OverlapRadius > neighborKey;
	neighborKey.set( _localToGlobal( _maxDepth-1 ) );

	for( LocalDepth d=_maxDepth-1 ; d>=0 ; d-- )
		for( TreeOctNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( _localDepth( node )==d && IsActiveNode( node->children ) )
		{
			neighborKey.template getNeighbors< true >( node , _NodeInitializer );
			for( int i=0 ; i<neighborKey.Width ; i++ ) for( int j=0 ; j<neighborKey.Width ; j++ ) for( int k=0 ; k<neighborKey.Width ; k++ ) SetGhostFlag( neighborKey.neighbors[ _localToGlobal(d) ].neighbors[i][j][k] , false );
		}

	_sNodes.set( *_tree , map );
	_setValidityFlags< FEMDegree , FEMBType >();
	for( TreeOctNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( !IsActiveNode( node ) ) node->nodeData.nodeIndex = -1;
	memoryUsage();
}


template< class Real >
template< int FEMDegree , BoundaryType BType >
void Octree< Real >::_setValidityFlags( void )
{
	for( int i=0 ; i<_sNodes.size() ; i++ )
	{
		const unsigned char MASK = ~( TreeNodeData::SPACE_FLAG | TreeNodeData::FEM_FLAG );
		_sNodes.treeNodes[i]->nodeData.flags &= MASK;
		if( isValidSpaceNode( _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= TreeNodeData::SPACE_FLAG;
		if( isValidFEMNode< FEMDegree , BType >( _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= TreeNodeData::FEM_FLAG;
	}
}

// Trim off the branches of the tree (finer than _fullDepth) that don't contain data
template< class Real >
template< class HasDataFunctor >
void Octree< Real >::_clipTree( const HasDataFunctor& f )
{
	// Because we are doing things in a brooded fashion, if any of the children has data then the whole brood is active
	for( TreeOctNode* temp=_tree->nextNode() ; temp ; temp=_tree->nextNode(temp) ) if( temp->children && _localDepth( temp )>=_fullDepth )
	{
		bool hasData = false;
		for( int c=0 ; c<Cube::CORNERS && !hasData ; c++ ) hasData |= f( temp->children + c );
		for( int c=0 ; c<Cube::CORNERS ; c++ ) SetGhostFlag( temp->children+c , !hasData );
	}
}

template< class Real >
template< bool HasGradients >
bool Octree< Real >::_setInterpolationInfoFromChildren( TreeOctNode* node , SparseNodeData< PointData< Real , HasGradients > , 0 >& interpolationInfo ) const
{
	if( IsActiveNode( node->children ) )
	{
		bool hasChildData = false;
		PointData< Real , HasGradients > pData;
#if POINT_DATA_RES
		Point3D< Real > center;
		Real width;
		_centerAndWidth( node , center , width );
		for( int c=0 ; c<Cube::CORNERS ; c++ )
			if( _setInterpolationInfoFromChildren( node->children + c , interpolationInfo ) )
			{
				const PointData< Real , HasGradients >& _pData = interpolationInfo[ node->children + c ];
				for( int cc=0 ; cc<PointData< Real , HasGradients >::SAMPLES ; cc++ )
				{
					int x[3];
					PointData< Real , HasGradients >::SetIndices( _pData[cc].position / _pData[cc].weight , center , width , x );
					pData[ x[0] + x[1]*PointData< Real , HasGradients >::RES + x[2]*PointData< Real , HasGradients >::RES*PointData< Real , HasGradients >::RES ] += _pData[cc];
				}
				hasChildData = true;
			}
#else // !POINT_DATA_RES
		for( int c=0 ; c<Cube::CORNERS ; c++ )
			if( _setInterpolationInfoFromChildren( node->children + c , interpolationInfo ) )
			{
				pData += interpolationInfo[ node->children + c ];
				hasChildData = true;
			}
#endif // POINT_DATA_RES
		if( hasChildData && IsActiveNode( node ) ) interpolationInfo[ node ] += pData;
		return hasChildData;
	}
	else return interpolationInfo( node )!=NULL;
}
template< class Real >
template< bool HasGradients >
SparseNodeData< PointData< Real , HasGradients > , 0 > Octree< Real >::_densifyInterpolationInfo( const std::vector< PointSample >& samples , Real pointValue , int adaptiveExponent ) const
{
	SparseNodeData< PointData< Real , HasGradients > , 0 > iInfo;
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const TreeOctNode* node = samples[i].node;
		const ProjectiveData< OrientedPoint3D< Real > , Real >& pData = samples[i].sample;
		while( !IsActiveNode( node ) ) node = node->parent;
		if( pData.weight )
		{
#if POINT_DATA_RES
			Point3D< Real > center;
			Real width;
			_centerAndWidth( node , center , width );
			_PointDataAccumulator_< Real , HasGradients >::_AddToPointData_( iInfo[node] , pData.data.p , pointValue * pData.weight , pData.data.n , center , width , pData.weight );
#else // !POINT_DATA_RES
			_PointDataAccumulator_< Real , HasGradients >::_AddToPointData_( iInfo[node] , pData.data.p , pointValue * pData.weight , pData.data.n , pData.weight );
#endif // POINT_DATA_RES
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot, iInfo );
#pragma omp parallel for
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
#if POINT_DATA_RES
		for( int c=0 ; c<PointData< Real , HasGradients >::SAMPLES ; c++ )
		{
			Real w = iInfo[i][c].weight;
			iInfo[i][c] /= w ; iInfo[i][c].weight = w;
		}
#else // !POINT_DATA_RES
	{
		Real w = iInfo[i].weight;
		iInfo[i] /= w ; iInfo[i].weight = w;
	}
#endif // POINT_DATA_RES
	LocalDepth maxDepth = _localMaxDepth( _tree );

	// Set the average position and scale the weights
	for( const TreeOctNode* node=_tree->nextNode() ; node ; node=_tree->nextNode(node) ) if( IsActiveNode( node ) )
	{
		PointData< Real , HasGradients >* pData = iInfo( node );
		if( pData )
		{
			int e = _localDepth( node ) * adaptiveExponent - ( maxDepth ) * (adaptiveExponent-1);
#if POINT_DATA_RES
			for( int c=0 ; c<PointData< Real , HasGradients >::SAMPLES ; c++ ) if( (*pData)[c].weight )
			{
				if( e<0 ) (*pData)[c].weight /= Real( 1<<(-e) );
				else      (*pData)[c].weight *= Real( 1<<  e  );
			}
#else // !POINT_DATA_RES
			if( e<0 ) pData->weight /= Real( 1<<(-e) );
			else      pData->weight *= Real( 1<<  e  );
#endif // POINT_DATA_RES
		}
	}
	return iInfo;
}
////////////////
// VertexData //
////////////////
long long VertexData::CenterIndex( const TreeOctNode* node , int maxDepth )
{
	int idx[DIMENSION];
	return CenterIndex(node,maxDepth,idx);
}
long long VertexData::CenterIndex(const TreeOctNode* node,int maxDepth,int idx[DIMENSION])
{
	int d , o[3];
	node->depthAndOffset( d , o );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , o[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CenterIndex( int depth , const int offSet[DIMENSION] , int maxDepth , int idx[DIMENSION] )
{
	for(int i=0;i<DIMENSION;i++) idx[i]=BinaryNode::CornerIndex( maxDepth+1 , depth+1 , offSet[i]<<1 , 1 );
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::CornerIndex(const TreeOctNode* node,int cIndex,int maxDepth)
{
	int idx[DIMENSION];
	return CornerIndex(node,cIndex,maxDepth,idx);
}
long long VertexData::CornerIndex( const TreeOctNode* node , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	int d , o[3];
	node->depthAndOffset( d , o );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d , o[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndex( int depth , const int offSet[DIMENSION] , int cIndex , int maxDepth , int idx[DIMENSION] )
{
	int x[DIMENSION];
	Cube::FactorCornerIndex( cIndex , x[0] , x[1] , x[2] );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , depth , offSet[i] , x[i] );
	return CornerIndexKey( idx );
}
long long VertexData::CornerIndexKey( const int idx[DIMENSION] )
{
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth){
	int idx[DIMENSION];
	return FaceIndex(node,fIndex,maxDepth,idx);
}
long long VertexData::FaceIndex(const TreeOctNode* node,int fIndex,int maxDepth,int idx[DIMENSION])
{
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	int d,o[3];
	node->depthAndOffset(d,o);
	for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);}
	idx[dir]=BinaryNode::CornerIndex(maxDepth+1,d,o[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth ){ int idx[DIMENSION] ; return EdgeIndex( node , eIndex , maxDepth , idx ); }
long long VertexData::EdgeIndex( const TreeOctNode* node , int eIndex , int maxDepth , int idx[DIMENSION] )
{
	int o , i1 , i2;
	int d , off[3];
	node->depthAndOffset( d ,off );
	Cube::FactorEdgeIndex( eIndex , o , i1 , i2 );
	for( int i=0 ; i<DIMENSION ; i++ ) idx[i] = BinaryNode::CornerIndex( maxDepth+1 , d+1 , off[i]<<1 , 1 );
	switch(o)
	{
		case 0:
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 1:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[2] = BinaryNode::CornerIndex( maxDepth+1 , d , off[2] , i2 );
			break;
		case 2:
			idx[0] = BinaryNode::CornerIndex( maxDepth+1 , d , off[0] , i1 );
			idx[1] = BinaryNode::CornerIndex( maxDepth+1 , d , off[1] , i2 );
			break;
	};
	return (long long)(idx[0]) | (long long)(idx[1])<<VERTEX_COORDINATE_SHIFT | (long long)(idx[2])<<(2*VERTEX_COORDINATE_SHIFT);
}

//#include "MultiGridOctreeData.inl"
#include "MultiGridOctreeData.SortedTreeNodes.inl"
#include "MultiGridOctreeData.WeightedSamples.inl"
#include "MultiGridOctreeData.System.inl"
#include "MultiGridOctreeData.IsoSurface.inl"
//#include "MultiGridOctreeData.Evaluation.inl"
#endif // MULTI_GRID_OCTREE_DATA_INCLUDED
