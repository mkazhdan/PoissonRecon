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

#include <functional>
#include <cmath>
#include <climits>
#include "MyMiscellany.h"

/////////////////////
// FEMTreeNodeData //
/////////////////////
FEMTreeNodeData::FEMTreeNodeData( void ){ flags = 0 , isGeometrySupported = false; }
FEMTreeNodeData::~FEMTreeNodeData( void ) { }


/////////////
// FEMTree //
/////////////
template< unsigned int Dim , class Real >
double FEMTree< Dim , Real >::MemoryUsage( void )
{
	double mem = double( MemoryInfo::Usage() ) / (1<<20);
	_MaxMemoryUsage = std::max< double >( mem , _MaxMemoryUsage );
	_LocalMemoryUsage = std::max< double >( mem , _LocalMemoryUsage );
	return mem;
}

template< unsigned int Dim , class Real > FEMTree< Dim , Real >::FEMTree( size_t blockSize ) : _nodeInitializer( *this )
{
	if( blockSize )
	{
		nodeAllocators.resize( std::thread::hardware_concurrency() );
		for( size_t i=0 ; i<nodeAllocators.size() ; i++ )
		{
			nodeAllocators[i] = new Allocator< FEMTreeNode >();
			nodeAllocators[i]->set( blockSize );
		}
	}
	_nodeCount = 0;
	_tree = FEMTreeNode::NewBrood( nodeAllocators.size() ? nodeAllocators[0] : NULL , _nodeInitializer );
	_tree->template initChildren< false >( nodeAllocators.size() ? nodeAllocators[0] : NULL , _nodeInitializer ) , _spaceRoot = _tree->children;
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] _spaceRoot is the root of the tree until finalization" )
#endif // SHOW_WARNINGS
	_spaceRoot->parent = NULL;
	int offset[Dim];
	for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
	RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
	_depthOffset = 0;
	memset( _femSigs1 , -1 , sizeof( _femSigs1 ) );
	memset( _femSigs2 , -1 , sizeof( _femSigs2 ) );
}
template< unsigned int Dim , class Real >
FEMTree< Dim , Real >::FEMTree( FILE* fp , XForm< Real , Dim+1 > &xForm , size_t blockSize ) : _nodeInitializer( *this )
{
	if( blockSize )
	{
		nodeAllocators.resize( std::thread::hardware_concurrency() );
		for( size_t i=0 ; i<nodeAllocators.size() ; i++ )
		{
			nodeAllocators[i] = new Allocator< FEMTreeNode >();
			nodeAllocators[i]->set( blockSize );
		}
	}
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
	if( fp )
	{
		if( fread( xForm.coords , sizeof( Real ) , (Dim+1)*(Dim+1) , fp )!=(Dim+1)*(Dim+1) ) ERROR_OUT( "Failed to read transform" );
		if( fread( &_depthOffset , sizeof( int ) , 1 , fp )!=1 ) ERROR_OUT( "Failed to read depth offset" );
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _nodeInitializer );
		_tree->read( fp , nodeAllocator , _nodeInitializer );
		_maxDepth = _tree->maxDepth() - _depthOffset;

		_spaceRoot = _tree->children;

		if( _depthOffset>1 )
		{
			_spaceRoot = _tree->children + (1<<Dim)-1;
			for( int d=1 ; d<_depthOffset ; d++ )
				if( !_spaceRoot->children ) ERROR_OUT( "Expected children" );
				else _spaceRoot = _spaceRoot->children;
		}
		_sNodes.set( *_tree , NULL );
	}
	else
	{
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _nodeInitializer );
		_tree->template initChildren< false >( nodeAllocator , _nodeInitializer ) , _spaceRoot = _tree->children;
		int offset[Dim];
		for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
		RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
		_depthOffset = 0;
	}
}
template< unsigned int Dim , class Real > void FEMTree< Dim , Real >::write( FILE* fp , XForm< Real , Dim+1 > xForm ) const
{
	fwrite( xForm.coords , sizeof( Real ) , (Dim+1)*(Dim+1) , fp );
	fwrite( &_depthOffset , sizeof( int ) , 1 , fp );
	_tree->write( fp );
}

template< unsigned int Dim , class Real >
const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p ) const
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	while( node->children )
	{
		int cIndex = FEMTreeNode::ChildIndex( center , p );
		node = node->children + cIndex;
		width /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) center[d] += width/2;
			else                  center[d] -= width/2;
	}
	return node;
}
template< unsigned int Dim , class Real >
template< bool ThreadSafe >
RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* FEMTree< Dim , Real >::_leaf( Allocator< FEMTreeNode > *nodeAllocator , Point< Real , Dim > p , LocalDepth maxDepth )
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	LocalDepth d = _localDepth( node );
	while( ( d<0 && node->children ) || ( d>=0 && d<maxDepth ) )
	{
		if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , _nodeInitializer );
		int cIndex = FEMTreeNode::ChildIndex( center , p );
		node = node->children + cIndex;
		d++;
		width /= 2;
		for( int d=0 ; d<Dim ; d++ )
			if( (cIndex>>d) & 1 ) center[d] += width/2;
			else                  center[d] -= width/2;
	}
	return node;
}

template< unsigned int Dim , class Real > bool FEMTree< Dim , Real >::_InBounds( Point< Real , Dim > p ){ for( int d=0 ; d<Dim ; d++ ) if( p[d]<0 || p[d]>1 ) return false ; return true; }
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSignatures >
bool FEMTree< Dim , Real >::isValidFEMNode( UIntPack< FEMSignatures ... > , const FEMTreeNode* node ) const
{
	if( GetGhostFlag< Dim >( node ) ) return false;
	LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	return FEMIntegrator::IsValidFEMNode( UIntPack< FEMSignatures ... >() , d , off );
}
template< unsigned int Dim , class Real >
bool FEMTree< Dim , Real >::isValidSpaceNode( const FEMTreeNode* node ) const
{
	if( !node ) return false;
	LocalDepth d ; LocalOffset off ; _localDepthAndOffset( node , d , off );
	if( d<0 ) return false;
	int res = 1<<d;
	for( int dd=0 ; dd<Dim ; dd++ ) if( off[dd]<0 || off[dd]>=res ) return false;
	return true;
}

template< unsigned int Dim , class Real >
template< bool ThreadSafe , unsigned int ... Degrees >
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , Allocator< FEMTreeNode > *nodeAllocator , FEMTreeNode* node , LocalDepth depth )
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	bool refine = d<depth && ( d<0 || !FEMIntegrator::IsOutOfBounds( UIntPack< FEMDegreeAndBType< Degrees , BOUNDARY_FREE >::Signature ... >() , d , off ) );
	if( refine )
	{
		if( !node->children ) node->template initChildren< ThreadSafe >( nodeAllocator , _nodeInitializer );
		for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth< ThreadSafe >( UIntPack< Degrees ... >() , nodeAllocator , node->children+c , depth );
	}
}
template< unsigned int Dim , class Real >
template< bool ThreadSafe , unsigned int ... Degrees >
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , Allocator< FEMTreeNode > *nodeAllocator , LocalDepth depth )
{
	if( !_tree->children ) _tree->template initChildren< ThreadSafe >( nodeAllocator , _nodeInitializer );
	for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth< ThreadSafe >( UIntPack< Degrees ... >() , nodeAllocator , _tree->children+c , depth );
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
typename FEMTree< Dim , Real >::LocalDepth FEMTree< Dim , Real >::_getFullDepth( UIntPack< Degrees ... > , const FEMTreeNode* node ) const
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	bool refine = d<0 || !FEMIntegrator::IsOutOfBounds( UIntPack< FEMDegreeAndBType< Degrees , BOUNDARY_FREE >::Signature ... >() , d , off );

	if( refine )
	{
		if( !node->children ) return d;
		else
		{
			LocalDepth depth = INT_MAX;
			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				LocalDepth d = _getFullDepth( UIntPack< Degrees ... >() , node->children+c );
				if( d<depth ) depth = d;
			}
			return depth;
		}
	}
	else return INT_MAX;
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
typename FEMTree< Dim , Real >::LocalDepth FEMTree< Dim , Real >::getFullDepth( UIntPack< Degrees ... > ) const
{
	if( !_tree->children ) return -1;
	LocalDepth depth = INT_MAX;
	for( int c=0 ; c<(1<<Dim) ; c++ )
	{
		LocalDepth d = _getFullDepth( UIntPack< Degrees ... >() , _tree->children+c );
		if( d<depth ) depth = d;
	}
	return depth;
}

template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , bool CreateNodes , typename ProcessingNodeFunctor , typename ... DenseOrSparseNodeData , typename InitializeFunctor >
void FEMTree< Dim , Real >::processNeighbors( ProcessingNodeFunctor processingNode , std::tuple< DenseOrSparseNodeData *... > data , InitializeFunctor initialize )
{
	std::vector< FEMTreeNode * > nodes;
	nodes.reserve( _spaceRoot->nodes() );
	for( FEMTreeNode *node=_spaceRoot->nextNode() ; node ; node=_spaceRoot->nextNode(node) ) if( processingNode(node) ) nodes.push_back( node );
	processNeighbors< LeftRadius , RightRadius , CreateNodes >( &nodes[0] , nodes.size() , data , initialize );
}

template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , bool CreateNodes , typename ... DenseOrSparseNodeData , typename InitializeFunctor >
void FEMTree< Dim , Real >::processNeighbors( FEMTreeNode **nodes , size_t nodeCount, std::tuple< DenseOrSparseNodeData *... > data , InitializeFunctor initialize )
{
	int maxDepth = 0;
	for( size_t i=0 ; i<nodeCount ; i++ ) maxDepth = std::max( maxDepth , nodes[i]->depth() );
	std::vector< node_index_type > map( _nodeCount );
	for( node_index_type i=0 ; i<_nodeCount ; i++ ) map[i] = i;
	typedef typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > NeighborKey;

	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
	//	typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > neighborKey;
	NeighborKey neighborKey;
	neighborKey.set( maxDepth );
	for( size_t i=0 ; i<nodeCount ; i++ )
	{
		auto neighbors = neighborKey.template getNeighbors< CreateNodes , false >( nodes[i] , nodeAllocator , _nodeInitializer );
		for( unsigned int j=0 ; j<neighbors.neighbors.Size ; j++ ) if( neighbors.neighbors.data[j] ) initialize( neighbors.neighbors.data[j] );
	}

	_reorderDenseOrSparseNodeData< 0 >( GetPointer( map ) , _nodeCount , data );
}

template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , typename IsProcessingNodeFunctor , typename ProcessingKernel >
void FEMTree< Dim , Real >::processNeighboringLeaves( IsProcessingNodeFunctor isProcessingNode , ProcessingKernel kernel , bool processSubTree )
{
	std::vector< FEMTreeNode * > nodes;
	nodes.reserve( _spaceRoot->nodes() );
	for( FEMTreeNode *node=_spaceRoot->nextNode() ; node ; node=_spaceRoot->nextNode(node) ) if( isProcessingNode(node) ) nodes.push_back( node );
	processNeighboringLeaves< LeftRadius , RightRadius >( &nodes[0] , nodes.size() , kernel , processSubTree );
}

template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , typename ProcessingKernel >
void FEMTree< Dim , Real >::processNeighboringLeaves( FEMTreeNode **nodes , size_t nodeCount  , ProcessingKernel kernel , bool processSubTree )
{
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] may process the same leaf multiple times, if it is coarser" )
#endif // SHOW_WARNINGS
	typedef typename RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > NeighborKey;
	// Suppose that we have a node at index I and we want the leaf nodes supported on the (possibly virtual) node K away
	// Case 1: The K-th neighbor exists
	// ---> Iterate over the leaf nodes of the sub-tree rooted at the K-th neighbor
	// Case 2: The K-th neighbor does not exist
	// ---> The index of the K-th neighbor is I+K
	// ---> The index of the parent is floor( I/2 )
	// ---> The index of the K-th neighbors parent is floor( (I+K)/2 )
	// ---> The parent of the K-th neighbor is the [ floor( (I+k)/2 ) - floor( I/2 ) ]-th neighbor of the parent

	std::function< void ( FEMTreeNode * ) > ProcessSubTree = [&]( FEMTreeNode *node )
	{
		if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) ProcessSubTree( node->children+c );
		kernel( node );
	};

	unsigned int maxDepth=0;
	for( size_t i=0 ; i<nodeCount ; i++ ) maxDepth = std::max( maxDepth , (unsigned int)nodes[i]->depth() );

	std::vector< NeighborKey > neighborKeys( ThreadPool::NumThreads() );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( maxDepth );

	ThreadPool::Parallel_for( 0 , nodeCount , [&]( unsigned int t , size_t  i )
	{
		typedef StaticWindow< FEMTreeNode * , IsotropicUIntPack< Dim , LeftRadius+RightRadius+1 > > NeighborLeafNodes;
		NeighborLeafNodes neighborLeafNodes;
		neighborKeys[t].setLeafNeighbors( nodes[i] , neighborLeafNodes );
		for( int i=0 ; i<NeighborLeafNodes::Size ; i++ ) if( neighborLeafNodes.data[i] )
			if( processSubTree ) ProcessSubTree( neighborLeafNodes.data[i] );
			else kernel( neighborLeafNodes.data[i] );
	} );
}

template< unsigned int Dim , class Real >
template< unsigned int CoDim , unsigned int DensityDegree >
typename FEMTree< Dim , Real >::template DensityEstimator< DensityDegree >* FEMTree< Dim , Real >::setDensityEstimator( const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode )
{
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	splatDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( splatDepth , maxDepth ) );
	DensityEstimator< DensityDegree >* _density = new DensityEstimator< DensityDegree >( splatDepth , CoDim , samplesPerNode );
	DensityEstimator< DensityDegree >& density = *_density;
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	densityKey.set( _localToGlobal( splatDepth ) );

	std::vector< node_index_type > sampleMap( nodeCount() , -1 );
	ThreadPool::Parallel_for( 0 , samples.size() , [&]( unsigned int , size_t i ){ if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = (node_index_type)i; } );
	std::function< ProjectiveData< Point< Real , Dim > , Real > ( FEMTreeNode* ) > SetDensity = [&] ( FEMTreeNode* node )
	{
		ProjectiveData< Point< Real , Dim > , Real > sample;
		LocalDepth d = _localDepth( node );
		node_index_type idx = node->nodeData.nodeIndex;
		if( node->children )
			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				ProjectiveData< Point< Real , Dim > , Real > s = SetDensity( node->children + c );
				if( d<=splatDepth && s.weight>0 )
				{
					Point< Real , Dim > p = s.data / s.weight;
					_addWeightContribution< true , CoDim >( nodeAllocator , density , node , p , densityKey , s.weight );
				}
				sample += s;
			}
		else if( idx<(node_index_type)sampleMap.size() && sampleMap[idx]!=-1 )
		{
			sample = samples[ sampleMap[ idx ] ].sample;
			if( d<=splatDepth && sample.weight>0 )
			{
				Point< Real , Dim > p = sample.data / sample.weight;
				_addWeightContribution< true , CoDim >( nodeAllocator , density , node , p , densityKey , sample.weight );
			}
		}
		return sample;
	};
	SetDensity( _spaceRoot );

	MemoryUsage();
	return _density;
}

template< unsigned int Dim , class Real >
template< unsigned int ... DataSigs , unsigned int DensityDegree , class InData , class OutData >
SparseNodeData< OutData , UIntPack< DataSigs ... > > FEMTree< Dim , Real >::setInterpolatedDataField( UIntPack< DataSigs ... > , const std::vector< PointSample >& samples , const std::vector< InData >& data , const DensityEstimator< DensityDegree >* density , LocalDepth minDepth , LocalDepth maxDepth , Real minDepthCutoff , Real& pointWeightSum , std::function< bool ( InData , OutData& ) > ConversionFunction , std::function< Real ( InData ) > BiasFunction )
{
	std::function< bool ( InData , OutData & , Real & ) > ConversionAndBiasFunction = [&]( InData in , OutData &out , Real &bias )
	{
		if( ConversionFunction( in , out ) )
		{
			bias = BiasFunction( in );
			return true;
		}
		else return false;
	};
	return setInterpolatedDataField( UIntPack< DataSigs ... >() , samples , data , density , minDepth , maxDepth , minDepthCutoff , pointWeightSum , ConversionAndBiasFunction );
}
template< unsigned int Dim , class Real >
template< unsigned int ... DataSigs , unsigned int DensityDegree , class InData , class OutData >
SparseNodeData< OutData , UIntPack< DataSigs ... > > FEMTree< Dim , Real >::setInterpolatedDataField( UIntPack< DataSigs ... > , const std::vector< PointSample >& samples , const std::vector< InData >& data , const DensityEstimator< DensityDegree >* density , LocalDepth minDepth , LocalDepth maxDepth , Real minDepthCutoff , Real& pointWeightSum , std::function< bool ( InData , OutData & , Real & ) > ConversionAndBiasFunction )
{
	typedef PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > DensityKey;
	typedef UIntPack< FEMSignature< DataSigs >::Degree ... > DataDegrees;
	typedef PointSupportKey< UIntPack< FEMSignature< DataSigs >::Degree ... > > DataKey;
	std::vector< DensityKey > densityKeys( ThreadPool::NumThreads() );
	std::vector<    DataKey >    dataKeys( ThreadPool::NumThreads() );
	bool oneKey = DensityDegree==DataDegrees::Min() && DensityDegree==DataDegrees::Max();
	for( size_t i=0 ; i<densityKeys.size() ; i++ ) densityKeys[i].set( _localToGlobal( maxDepth ) );
	if( !oneKey ) for( size_t i=0 ; i<dataKeys.size() ; i++ ) dataKeys[i].set( _localToGlobal( maxDepth ) );
	Real weightSum = 0;
	pointWeightSum = 0;
	SparseNodeData< OutData , UIntPack< DataSigs ... > > dataField;
	Real _pointWeightSum = 0;
	ThreadPool::Parallel_for( 0 , samples.size() , [&]( unsigned int thread , size_t i )
	{
		DensityKey& densityKey = densityKeys[ thread ];
		DataKey& dataKey = dataKeys[ thread ];
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		if( sample.weight>0 )
		{
			Point< Real , Dim > p = sample.data / sample.weight;
			InData in = data[i] / sample.weight;
			OutData out;

			Real depthBias;
			if( !_InBounds(p) ) WARN( "Point sample is out of bounds" );
			else if( ConversionAndBiasFunction( in , out , depthBias ) )
			{
				AddAtomic( weightSum , sample.weight );
				out *= sample.weight;
				Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[ thread ] : NULL;
#if defined( __GNUC__ ) && __GNUC__ < 5
#ifdef SHOW_WARNINGS
#warning "you've got me gcc version<5"
#endif // SHOW_WARNINGS
				if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , true , DensityDegree , OutData >( nodeAllocator , *density , minDepthCutoff , p , out , dataField , densityKey , oneKey ? *( (DataKey*)&densityKey ) : dataKey , minDepth , maxDepth , Dim , depthBias ) * sample.weight );
#else // !__GNUC__ || __GNUC__ >=5
				if( density ) AddAtomic( _pointWeightSum , _splatPointData< true , true , DensityDegree , OutData , DataSigs ... >( nodeAllocator , *density , minDepthCutoff , p , out , dataField , densityKey , oneKey ? *( (DataKey*)&densityKey ) : dataKey , minDepth , maxDepth , Dim , depthBias ) * sample.weight );
#endif // __GNUC__ && __GNUC__ < 5
				else
				{
					Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#ifdef SHOW_WARNINGS
#warning "you've got me gcc version<5"
#endif // SHOW_WARNINGS
						_splatPointData< true , true , OutData >( nodeAllocator , _leaf< true >( nodeAllocator , p , maxDepth ) , p , out / (Real)pow( width , Dim ) , dataField , oneKey ? *( (DataKey*)&densityKey ) : dataKey );
#else // !__GNUC__ || __GNUC__ >=5
					_splatPointData< true , true , OutData , DataSigs ... >( nodeAllocator , _leaf< true >( nodeAllocator , p , maxDepth ) , p , out / (Real)pow( width , Dim ) , dataField , oneKey ? *( (DataKey*)&densityKey ) : dataKey );
#endif // __GNUC__ || __GNUC__ < 4
					AddAtomic( _pointWeightSum , sample.weight );
				}
			}
		}
	}
	);
	pointWeightSum = _pointWeightSum / weightSum;
	MemoryUsage();
	return dataField;
}

template< unsigned int Dim , class Real >
template< unsigned int DataSig , bool CreateNodes , unsigned int DensityDegree , class Data >
SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > FEMTree< Dim , Real >::setExtrapolatedDataField( const std::vector< PointSample >& samples , std::vector< Data >& sampleData , const DensityEstimator< DensityDegree >* density , bool nearest )
{
	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	PointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > dataField;
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		const Data& data = sampleData[i];
		Point< Real , Dim > p = sample.weight==0 ? sample.data : sample.data / sample.weight;
		if( !_InBounds(p) )
		{
			WARN( "Point is out of bounds" );
			continue;
		}
		if( nearest ) _nearestMultiSplatPointData< DensityDegree >( density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , 2 );
		else          _multiSplatPointData< CreateNodes , false , DensityDegree >( nodeAllocator , density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , dataKey , 2 );
	}
	MemoryUsage();
	return dataField;
}

template< unsigned int Dim , class Real >
template< unsigned int MaxDegree >
void FEMTree< Dim , Real >::_supportApproximateProlongation( void )
{
	// Refine the tree so that if an active element exists @{d}, all supporting elements exist @{d-1}
	const int OverlapRadius = -BSplineOverlapSizes< MaxDegree , MaxDegree >::OverlapStart;
	typedef typename FEMTreeNode::template NeighborKey< IsotropicUIntPack< Dim , OverlapRadius > , IsotropicUIntPack< Dim , OverlapRadius > > NeighborKey;

	std::vector< NeighborKey > neighborKeys( ThreadPool::NumThreads() );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth-1 ) );

	for( LocalDepth d=_maxDepth-1 ; d>_fullDepth ; d-- )
	{
		// Compute the set of nodes at depth d that have (non-ghost) children at depth d+1.
		std::vector< FEMTreeNode* > nodes;
		{
			auto NodeTerminationLambda = [&]( const FEMTreeNode *node ){ return _localDepth( node )==d; };
			for( FEMTreeNode* node=_tree->nextNode( NodeTerminationLambda , NULL ) ; node ; node=_tree->nextNode( NodeTerminationLambda , node ) ) if( _localDepth( node )==d && IsActiveNode( node->children ) ) nodes.push_back( node );
		}

		// Make sure that all finite elements whose support overlaps the support of the finite elements indexed by those nodes are in the tree.
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] This may be overkill as we only need to check if the support overlaps the support of the children" )
#endif // SHOW_WARNINGS
		ThreadPool::Parallel_for( 0 , nodes.size() , [&]( unsigned int thread , size_t i )
		{
			NeighborKey& neighborKey = neighborKeys[ thread ];
			FEMTreeNode *node = nodes[i];

			// Create the neighbors if they are not already in the tree
			neighborKey.template getNeighbors< true , true >( node , nodeAllocators.size() ? nodeAllocators[ thread ] : NULL , _nodeInitializer );

			// Mark the neighbors as active
			Pointer( FEMTreeNode* ) nodes = neighborKey.neighbors[ _localToGlobal(d) ].neighbors().data;
			unsigned int size = neighborKey.neighbors[ _localToGlobal(d) ].neighbors.Size;
			for( unsigned int i=0 ; i<size ; i++ ) SetGhostFlag( nodes[i] , false );
		}
		);
	}
}

template< unsigned int Dim , typename Real >
template< unsigned int SystemDegree >
void FEMTree< Dim , Real >::_markNonBaseDirichletElements( void )
{
	const int LeftSupportRadius = -BSplineSupportSizes< SystemDegree >::SupportStart;
	const int RightSupportRadius = BSplineSupportSizes< SystemDegree >::SupportEnd;
	typedef typename FEMTreeNode::template NeighborKey< IsotropicUIntPack< Dim , LeftSupportRadius > , IsotropicUIntPack< Dim , RightSupportRadius > > SupportKey;
	typedef StaticWindow< FEMTreeNode * , IsotropicUIntPack< Dim , LeftSupportRadius + RightSupportRadius + 1 > > NeighborLeaves;

	std::vector< NeighborLeaves > neighborLeaves( ThreadPool::NumThreads() );
	std::vector< SupportKey > supportKeys( ThreadPool::NumThreads() );
	for( int i=0 ; i<supportKeys.size() ; i++ ) supportKeys[i].set( _localToGlobal( _maxDepth ) );

	// Get the list of nodes @{_baseDepth)
	std::vector< FEMTreeNode * > baseNodes;
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( _localDepth( node )==_baseDepth ) baseNodes.push_back( node );

	// Process the sub-tree rooted at the node:
	// -- For each non-ghost node in the sub-tree check if the finite element associated to the node is supported on a Dirichlet node.
	//    If it is, mark the node as a Dirichlet element node.
	std::function< void ( FEMTreeNode * , SupportKey & , NeighborLeaves & ) > ProcessSubTree = [&]( FEMTreeNode *node , SupportKey &supportKey , NeighborLeaves &neighborLeaves )
	{
		if( !node->nodeData.getGhostFlag() )
		{
			supportKey.setLeafNeighbors( node , neighborLeaves );
			bool hasDirichletNeighbor = false;
			for( int i=0 ; i<NeighborLeaves::Size ; i++ ) if( neighborLeaves.data[i] && neighborLeaves.data[i]->nodeData.getDirichletNodeFlag() ) hasDirichletNeighbor = true;
			node->nodeData.setDirichletElementFlag( hasDirichletNeighbor );

			if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) ProcessSubTree( node->children+c , supportKey , neighborLeaves );
		}
	};

	ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t i ){ ProcessSubTree( baseNodes[i] , supportKeys[t] , neighborLeaves[t] ); } );
}

template< unsigned int Dim , typename Real >
template< unsigned int SystemDegree >
void FEMTree< Dim , Real >::_markBaseDirichletElements( void )
{
	const int LeftSupportRadius = BSplineSupportSizes< SystemDegree >::SupportEnd;
	const int RightSupportRadius = -BSplineSupportSizes< SystemDegree >::SupportStart;

	typedef typename FEMTreeNode::template NeighborKey< IsotropicUIntPack< Dim , LeftSupportRadius > , IsotropicUIntPack< Dim , RightSupportRadius > > SupportKey;
	std::vector< SupportKey > supportKeys( ThreadPool::NumThreads() );
	for( int i=0 ; i<supportKeys.size() ; i++ ) supportKeys[i].set( _localToGlobal( _baseDepth ) );

	std::vector< FEMTreeNode* > nodes;
	auto TerminationLambda = [&]( const FEMTreeNode *node ){ return _localDepth(node)==_baseDepth; };
	for( FEMTreeNode* node=_tree->nextNode( TerminationLambda , NULL ) ; node ; node=_tree->nextNode( TerminationLambda , node ) ) if( _localDepth( node )==_baseDepth && node->nodeData.getDirichletNodeFlag() ) nodes.push_back( node );

	ThreadPool::Parallel_for( 0 , nodes.size() , [&]( unsigned int thread , size_t i )
	{
		SupportKey &supportKey = supportKeys[ thread ];
		FEMTreeNode *node = nodes[i];
		supportKey.getNeighbors( node );
		for( LocalDepth d=0 ; d<=_baseDepth ; d++ )
		{
			Pointer( FEMTreeNode* ) _nodes = supportKey.neighbors[ _localToGlobal(d) ].neighbors().data;
			unsigned int size = supportKey.neighbors[ _localToGlobal(d) ].neighbors.Size;
			for( unsigned int i=0 ; i<size ; i++ ) if( _nodes[i] ) SetGhostFlag( _nodes[i] , false ) , _nodes[i]->nodeData.setDirichletElementFlag( true );
		}
	} );
}

template< unsigned int Dim , class Real >
template< unsigned int MaxDegree , unsigned int SystemDegree , typename HasDataFunctor , typename IsDirichletLeafFunctor , typename ... InterpolationInfos , typename ... DenseOrSparseNodeData >
void FEMTree< Dim , Real >::finalizeForMultigrid( LocalDepth baseDepth , LocalDepth fullDepth , const HasDataFunctor hasData , const IsDirichletLeafFunctor isDirichletLeaf , std::tuple< InterpolationInfos *... > interpolationInfos , std::tuple< DenseOrSparseNodeData *... > data )
{
	std::function< void ( FEMTreeNode * ) > pushFullDirichletFlag = [&]( FEMTreeNode *node )
	{
		if( node->children )
		{
			if( node->nodeData.getDirichletNodeFlag() )
			{
				for( int c=0 ; c<(1<<Dim) ; c++ ) node->children[c].nodeData.setDirichletNodeFlag( true );
				node->nodeData.setDirichletNodeFlag( false );
			}
			for( int c=0 ; c<(1<<Dim) ; c++ ) pushFullDirichletFlag( node->children + c );
		}
	};

	std::function< bool ( FEMTreeNode * ) > pullPartialDirichletFlag = [&]( FEMTreeNode *node )
	{
		if( node->children )
		{
			bool childDirichlet = false;
			for( int c=0 ; c<(1<<Dim) ; c++ ) childDirichlet |= pullPartialDirichletFlag( node->children + c );
			node->nodeData.setDirichletNodeFlag( childDirichlet );
		}
		return node->nodeData.getDirichletNodeFlag();
	};

	std::function< void ( FEMTreeNode * , node_index_type , bool ) > pushDirichletFlag = [&]( FEMTreeNode *node , node_index_type newNodeIndex , bool isDirichletNode )
	{
		if( node->nodeData.nodeIndex>=newNodeIndex ) node->nodeData.setDirichletNodeFlag( isDirichletNode );
		isDirichletNode = node->nodeData.getDirichletNodeFlag();
		if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) pushDirichletFlag( node->children+c , newNodeIndex , isDirichletNode );
	};

	if( baseDepth>fullDepth ) ERROR_OUT( "Base depth cannot exceed full depth: " , baseDepth , " <= " , fullDepth );
	_baseDepth = baseDepth;
	_spaceRoot->parent = _tree;

	Allocator< FEMTreeNode > *nodeAllocator = nodeAllocators.size() ? nodeAllocators[0] : NULL;
	_depthOffset = 1;
	while( _localInset( 0 ) + BSplineEvaluationData< FEMDegreeAndBType< MaxDegree >::Signature >::Begin( 0 )<0 || _localInset( 0 ) + BSplineEvaluationData< FEMDegreeAndBType< MaxDegree >::Signature >::End( 0 )>(1<<_depthOffset) )
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

		FEMTreeNode* newSpaceRootParent = FEMTreeNode::NewBrood( nodeAllocator , _nodeInitializer );
		FEMTreeNode* oldSpaceRootParent = _spaceRoot->parent;
		int corner = _depthOffset<=1 ? (1<<Dim)-1 : 0;
		newSpaceRootParent[corner].children = _spaceRoot;
		oldSpaceRootParent->children = newSpaceRootParent;
		for( int c=0 ; c<(1<<Dim) ; c++ ) _spaceRoot[c].parent = newSpaceRootParent + corner , newSpaceRootParent[c].parent = oldSpaceRootParent;
		_depthOffset++;
	}
	int d=0 , off[Dim];
	for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
	FEMTreeNode::ResetDepthAndOffset( _tree , d , off );
	_maxDepth = _spaceRoot->maxDepth();
	_fullDepth = std::min( fullDepth , _maxDepth );

	// Mark leaf nodes that are Dirichlet constraints so they do not get clipped out.
	// Need to do this before introducing new nodes into the tree  (since isDirichletLeaf depends on the structure at input).
	for( FEMTreeNode *leaf=_spaceRoot->nextLeaf() ; leaf ; leaf=_spaceRoot->nextLeaf(leaf) ) leaf->nodeData.setDirichletNodeFlag( isDirichletLeaf( leaf ) );

	// Make the low-resolution part of the tree be complete
	_setFullDepth< false >( IsotropicUIntPack< Dim , MaxDegree >() , nodeAllocator , _fullDepth );

	// Mark new leaf nodes
	pushFullDirichletFlag( _spaceRoot );

	// Pull the Dirichlet designator from the leaves so that nodes are now marked if they contain (possibly partial Dirichlet constraints
	pullPartialDirichletFlag( _spaceRoot );

	// Use the node Dirichlet designators to set the coarser finite element Dirichlet designators
	_markBaseDirichletElements< SystemDegree >();

	// Clear all the flags and make everything that is not low-res a ghost node
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) node->nodeData.flags &= ( FEMTreeNodeData::SCRATCH_FLAG | FEMTreeNodeData::DIRICHLET_NODE_FLAG | FEMTreeNodeData::DIRICHLET_ELEMENT_FLAG ) , SetGhostFlag( node , _localDepth( node )>fullDepth );

	// Clip off nodes that not have data and do not contain geometry or Dirichlet constraints below the exactDepth
	_clipTree( [&]( const FEMTreeNode *node ){ return hasData(node) || ( ( node->nodeData.getDirichletNodeFlag() || node->nodeData.getDirichletElementFlag() ) && _localDepth(node)<=_baseDepth ); } , fullDepth );

	// It is possible for the tree to have become shallower after clipping
	_maxDepth = _tree->maxDepth() - _depthOffset;

	node_index_type oldNodeCount = _nodeCount;

	// Refine the node so that finite elements @{depth-1} whose support overlaps finite elements @{depth} are in the tree
	_supportApproximateProlongation< MaxDegree >();

	// Mark new leaf nodes
	pushDirichletFlag( _spaceRoot , oldNodeCount , _spaceRoot->nodeData.getDirichletNodeFlag() );
	_markNonBaseDirichletElements< SystemDegree >();

	std::vector< node_index_type > map;
	_sNodes.set( *_tree , &map );
	_setSpaceValidityFlags();
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( !IsActiveNode( node ) ) node->nodeData.nodeIndex = -1;
	_reorderDenseOrSparseNodeData< 0 >( GetPointer( map ) , _sNodes.size() , data );
	_reorderInterpolationInfo< 0 >( GetPointer( map ) , _sNodes.size() , interpolationInfos );
	MemoryUsage();
}

template< unsigned int Dim , class Real >
template< class ... DenseOrSparseNodeData >
void FEMTree< Dim , Real >::resetIndices( std::tuple< DenseOrSparseNodeData *... > data )
{
	std::vector< node_index_type > map;
	_sNodes.set( *_tree , &map );
	_setSpaceValidityFlags();
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( !IsActiveNode< Dim >( node ) ) node->nodeData.nodeIndex = -1;
	_reorderDenseOrSparseNodeData< 0 >( &map[0] , _sNodes.size() , data );
	MemoryUsage();
}

template< unsigned int Dim , class Real >
void FEMTree< Dim , Real >::_setSpaceValidityFlags( void ) const
{
	ThreadPool::Parallel_for( 0 , _sNodes.size() , [&]( unsigned int , size_t i )
	{
		const unsigned char MASK = ~( FEMTreeNodeData::SPACE_FLAG );
		_sNodes.treeNodes[i]->nodeData.flags &= MASK;
		if( isValidSpaceNode( _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::SPACE_FLAG;
	}
	);
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs1 >
void FEMTree< Dim , Real >::_setFEM1ValidityFlags( UIntPack< FEMSigs1 ... > ) const
{
	bool needToReset;
	unsigned int femSigs1[] = { FEMSigs1 ... };
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock( m );
		needToReset = memcmp( femSigs1 , _femSigs1 , sizeof( _femSigs1 ) )!=0;
		if( needToReset ) memcpy( _femSigs1 , femSigs1 , sizeof( _femSigs1 ) );
	}
	if( needToReset )
		for( node_index_type i=0 ; i<(node_index_type)_sNodes.size() ; i++ )
		{
			const unsigned char MASK = ~( FEMTreeNodeData::FEM_FLAG_1 );
			_sNodes.treeNodes[i]->nodeData.flags &= MASK;
			if( isValidFEMNode( UIntPack< FEMSigs1 ... >() , _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::FEM_FLAG_1;
		}

}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs2 >
void FEMTree< Dim , Real >::_setFEM2ValidityFlags( UIntPack< FEMSigs2 ... > ) const
{
	bool needToReset;
	unsigned int femSigs2[] = { FEMSigs2 ... };
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
		needToReset = memcmp( femSigs2 , _femSigs2 , sizeof( _femSigs2 ) )!=0;
		if( needToReset ) memcpy( _femSigs2 , femSigs2 , sizeof( _femSigs2 ) );
	}
	if( needToReset )
		for( node_index_type i=0 ; i<(node_index_type)_sNodes.size() ; i++ )
		{
			const unsigned char MASK = ~( FEMTreeNodeData::FEM_FLAG_2 );
			_sNodes.treeNodes[i]->nodeData.flags &= MASK;
			if( isValidFEMNode( UIntPack< FEMSigs2 ... >() , _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::FEM_FLAG_2;
		}
}

template< unsigned int Dim , class Real >
template< class HasDataFunctor >
void FEMTree< Dim , Real >::_clipTree( const HasDataFunctor& f , LocalDepth fullDepth )
{
	std::vector< FEMTreeNode * > regularNodes;
	auto NodeTerminationLambda = [&]( const FEMTreeNode *node ){ return _localDepth( node )==fullDepth; };
	for( FEMTreeNode* temp=_tree->nextNode( NodeTerminationLambda , NULL ) ; temp ; temp=_tree->nextNode( NodeTerminationLambda , temp ) ) if( _localDepth( temp )==fullDepth ) regularNodes.push_back( temp );

	// Get the data status of each node
	std::vector< char > nodeHasData( nodeCount() , 0 );
	ThreadPool::Parallel_for( 0 , regularNodes.size() , [&]( unsigned int , size_t i )
	{
		for( FEMTreeNode* node=regularNodes[i]->nextNode() ; node ; node=regularNodes[i]->nextNode(node) ) nodeHasData[node->nodeData.nodeIndex] = f( node ) ? 1 : 0;
	} );

	// Pull the data status from the leaves
	std::function< char ( const FEMTreeNode * ) > PullHasDataFromChildren = [&]( const FEMTreeNode *node )
	{
		char hasData = nodeHasData[node->nodeData.nodeIndex];
		if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) hasData |= PullHasDataFromChildren( node->children+c );
		nodeHasData[node->nodeData.nodeIndex] = hasData;
		return hasData;
	};

	ThreadPool::Parallel_for( 0 , regularNodes.size() , [&]( unsigned int , size_t i ){ PullHasDataFromChildren( regularNodes[i] ); } );

	// Mark all children of a node as ghost if none of them have data
	ThreadPool::Parallel_for( 0 , regularNodes.size() , [&]( unsigned int , size_t i )
	{
		for( FEMTreeNode* node=regularNodes[i]->nextNode() ; node ; node=regularNodes[i]->nextNode(node) ) if( node->children )
		{
			char childHasData = 0;
			for( int c=0 ; c<(1<<Dim) ; c++ ) childHasData |= nodeHasData[node->children[c].nodeData.nodeIndex];
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetGhostFlag< Dim >( node->children+c , !childHasData );
		}
	} );
}

template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::_ExactPointAndDataInterpolationInfo< T , Data , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
	ThreadPool::Parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode* node , node_index_type &start )
	{
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

	node_index_type start = 0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointAndDataInfo< Dim , Real , Data , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.pointInfo.position = pData.data;
			_pData.pointInfo.weight = pData.weight;
			_pData.pointInfo.dualValues = _constraintDual( pData.data/pData.weight , sampleData[i]/pData.weight ) * pData.weight;
			_pData.data = sampleData[i];
		}
	}

	ThreadPool::Parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i  )
	{
		Real w = _iData[i].pointInfo.weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].pointInfo.weight = w;
		else            _iData[i].pointInfo.weight = w * ( 1<<tree._maxDepth );
		_iData[i].pointInfo.dualValues *= _iData[i].pointInfo.weight;
	}
	);
}

template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< T , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
	ThreadPool::Parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode* node , node_index_type &start )
	{
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

	node_index_type start=0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointInfo< Dim , Real , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.position = pData.data;
			_pData.dualValues = _constraintDual( pData.data/pData.weight ) * pData.weight;
			_pData.weight = pData.weight;
		}
	}

	ThreadPool::Parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i )
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
	);
}
template< unsigned int Dim , class Real >
template< unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< double , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
	ThreadPool::Parallel_for( 0 , tree.nodesSize() , [&]( unsigned int , size_t i ){ _sampleSpan[i] = std::pair< node_index_type , node_index_type >( 0 , 0 ); } );
	for( node_index_type i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , node_index_type & ) > SetRange = [&] ( FEMTreeNode *node , node_index_type &start )
	{
		std::pair< node_index_type , node_index_type >& span = _sampleSpan[ node->nodeData.nodeIndex ];
		if( tree._isValidSpaceNode( node->children ) )
		{
			for( int c=0 ; c<(1<<Dim) ; c++ ) SetRange( node->children + c , start );
			span.first  = _sampleSpan[ node->children[0           ].nodeData.nodeIndex ].first;
			span.second = _sampleSpan[ node->children[ (1<<Dim)-1 ].nodeData.nodeIndex ].second;
		}
		else
		{
			span.second = start + span.second - span.first;
			span.first = start;
			start += span.second - span.first;
		}
	};

	node_index_type start = 0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( node_index_type i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) )
		{
			const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
			DualPointInfo< Dim , Real , T , PointD >& _pData = _iData[ _sampleSpan[ leaf->nodeData.nodeIndex ].second++ ];
			_pData.position = pData.data;
			_pData.dualValues = _constraintDual( pData.data/pData.weight ) * pData.weight;
			_pData.weight = pData.weight;
		}
	}

	ThreadPool::Parallel_for( 0 , _iData.size() , [&]( unsigned int , size_t i )
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
	);
}
template< unsigned int Dim , class Real >
template< typename T >
bool FEMTree< Dim , Real >::_setInterpolationInfoFromChildren( FEMTreeNode* node , SparseNodeData< T , IsotropicUIntPack< Dim , FEMTrivialSignature > >& interpolationInfo ) const
{
	if( IsActiveNode< Dim >( node->children ) )
	{
		bool hasChildData = false;
		T t = {};
		for( int c=0 ; c<(1<<Dim) ; c++ )
			if( _setInterpolationInfoFromChildren( node->children + c , interpolationInfo ) )
			{
				t += interpolationInfo[ node->children + c ];
				hasChildData = true;
			}
		if( hasChildData && IsActiveNode< Dim >( node ) ) interpolationInfo[ node ] += t;
		return hasChildData;
	}
	else return interpolationInfo( node )!=NULL;
}
template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual >
SparseNodeData< DualPointInfo< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstraintDual constraintDual , int adaptiveExponent ) const
{
	SparseNodeData< DualPointInfo< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointInfo< Dim , Real , T , PointD >& _pData = iInfo[node];
			_pData.position += pData.data;
			_pData.weight += pData.weight;
			_pData.dualValues += constraintDual( pData.data/pData.weight ) * pData.weight;
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

	ThreadPool::Parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
	{
		Real w = iInfo[i].weight;
		iInfo[i] /= w ; iInfo[i].weight = w;
	}
	);
	LocalDepth maxDepth = _spaceRoot->maxDepth();

	// Set the average position and scale the weights
	for( const FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode(node) ) if( IsActiveNode< Dim >( node ) )
	{
		DualPointInfo< Dim , Real , T , PointD >* pData = iInfo( node );
		if( pData )
		{
			int e = _localDepth( node ) * adaptiveExponent - ( maxDepth ) * (adaptiveExponent-1);
			if( e<0 ) pData->weight /= Real( 1<<(-e) );
			else      pData->weight *= Real( 1<<  e  );
			pData->dualValues *= pData->weight;
		}
	}
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual >
SparseNodeData< DualPointAndDataInfo< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , int adaptiveExponent ) const
{
	SparseNodeData< DualPointAndDataInfo< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
	for( node_index_type i=0 ; i<(node_index_type)samples.size() ; i++ )
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointAndDataInfo< Dim , Real , Data , T , PointD >& _pData = iInfo[node];
			_pData.pointInfo.position += pData.data;
			_pData.pointInfo.dualValues += constraintDual( pData.data/pData.weight , sampleData[i]/pData.weight ) * pData.weight;
			_pData.pointInfo.weight += pData.weight;
			_pData.data += sampleData[i];
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

	ThreadPool::Parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
	{
		Real w = iInfo[i].pointInfo.weight;
		iInfo[i] /= w ; iInfo[i].pointInfo.weight = w;
	}
	);
	LocalDepth maxDepth = _spaceRoot->maxDepth();

	// Set the average position and scale the weights
	for( const FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode(node) ) if( IsActiveNode< Dim >( node ) )
	{
		DualPointAndDataInfo< Dim , Real , Data , T , PointD >* pData = iInfo( node );
		if( pData )
		{
			int e = _localDepth( node ) * adaptiveExponent - ( maxDepth ) * (adaptiveExponent-1);
			if( e<0 ) pData->pointInfo.weight /= Real( 1<<(-e) );
			else      pData->pointInfo.weight *= Real( 1<<  e  );
			pData->pointInfo.dualValues *= pData->pointInfo.weight;
		}
	}
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual >
SparseNodeData< DualPointInfoBrood< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstraintDual constraintDual , bool noRescale ) const
{
	SparseNodeData< DualPointInfoBrood< Dim , Real , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
	for( node_index_type i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointInfoBrood< Dim , Real , T , PointD >& _pData = iInfo[node];
			Point< Real , Dim > p = pData.data/pData.weight;
			int cIdx = _childIndex( node , p );
			_pData[cIdx].position += pData.data;
			_pData[cIdx].weight += pData.weight;
			_pData[cIdx].dualValues += constraintDual( p ) * pData.weight;
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

	ThreadPool::Parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
	{
		iInfo[i].finalize();
		for( size_t c=0 ; c<iInfo[i].size() ; c++ )
		{
			iInfo[i][c].position /= iInfo[i][c].weight;
			if( !noRescale )
			{
				iInfo[i][c].weight     *= ( 1<<_maxDepth );
				iInfo[i][c].dualValues *= ( 1<<_maxDepth );
			}
		}
	}
	);
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual >
SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , bool noRescale ) const
{
	SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
	for( node_index_type i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* node = samples[i].node;
		const ProjectiveData< Point< Real , Dim > , Real >& pData = samples[i].sample;
		while( !IsActiveNode< Dim >( node ) ) node = node->parent;
		if( pData.weight )
		{
			DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD >& _pData = iInfo[node];
			Point< Real , Dim > p = pData.data/pData.weight;
			int cIdx = _childIndex( node , p );
			_pData[cIdx].pointInfo.position += pData.data;
			_pData[cIdx].pointInfo.dualValues += constraintDual( p , sampleData[i]/pData.weight ) * pData.weight;
			_pData[cIdx].pointInfo.weight += pData.weight;
			_pData[cIdx].data += sampleData[i];
		}
	}

	// Set the interior values
	_setInterpolationInfoFromChildren( _spaceRoot , iInfo );

	ThreadPool::Parallel_for( 0 , iInfo.size() , [&]( unsigned int , size_t i )
	{
		iInfo[i].finalize();
		for( size_t c=0 ; c<iInfo[i].size() ; c++ )
		{
			iInfo[i][c].pointInfo.position /= iInfo[i][c].pointInfo.weight;
			iInfo[i][c].data /= iInfo[i][c].pointInfo.weight;
			if( !noRescale )
			{
				iInfo[i][c].pointInfo.weight     *= ( 1<<_maxDepth );
				iInfo[i][c].pointInfo.dualValues *= ( 1<<_maxDepth );
				iInfo[i][c].data                 *= ( 1<<_maxDepth );
			}
		}
	}
	);
	return iInfo;
}



template< unsigned int Dim , class Real >
std::vector< node_index_type > FEMTree< Dim , Real >::merge( FEMTree* tree )
{
	std::vector< node_index_type > map;
	if( _depthOffset!=tree->_depthOffset ) ERROR_OUT( "depthOffsets don't match: %d != %d" , _depthOffset , tree->_depthOffset );

	// Compute the next available index
	node_index_type nextIndex = 0;
	for( const FEMTreeNode* node=_tree->nextNode() ; node!=NULL ; node=_tree->nextNode( node ) ) nextIndex = std::max< node_index_type >( nextIndex , node->nodeData.nodeIndex+1 );

	// Set the size of the map
	{
		node_index_type mapSize = 0;
		for( const FEMTreeNode* node=tree->_tree->nextNode() ; node!=NULL ; node=tree->_tree->nextNode( node ) ) mapSize = std::max< node_index_type >( mapSize , node->nodeData.nodeIndex+1 );
		map.resize( mapSize );
	}

	std::function< void ( FEMTreeNode* , FEMTreeNode* , std::vector< node_index_type > & , node_index_type & ) > MergeNodes = [&]( FEMTreeNode* node1 , FEMTreeNode* node2 , std::vector< node_index_type > &map , node_index_type &nextIndex )
	{
		if( node1 && node2 )
		{
			if( node2->nodeData.nodeIndex>=0 )
			{
				if( node1->nodeData.nodeIndex<0 ) node1->nodeData.nodeIndex = nextIndex++;
				map[ node2->nodeData.nodeIndex ] = node1->nodeData.nodeIndex;
			}
			if( node1->children && node2->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( node1->children+c , node2->children+c , map , nextIndex );
			else if( node2->children )
			{
				for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( NULL , node2->children+c , map , nextIndex );
				node1->children = node2->children;
				node2->children = NULL;
				for( int c=0 ; c<(1<<Dim) ; c++ ) node1->children[c].parent = node1;
			}
		}
		else if( node2 )
		{
			if( node2->nodeData.nodeIndex>=0 ){ map[ node2->nodeData.nodeIndex ] = nextIndex ; node2->nodeData.nodeIndex = nextIndex++; }
			if( node2->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) MergeNodes( NULL , node2->children+c , map , nextIndex );
		}
	};

	MergeNodes( _tree , tree->_tree , map , nextIndex );
	return map;
}

