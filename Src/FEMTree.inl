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
FEMTreeNodeData::FEMTreeNodeData( void ){ flags = 0; }
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

template< unsigned int Dim , class Real > FEMTree< Dim , Real >::FEMTree( int blockSize )
{
	if( blockSize>0 )
	{
		nodeAllocator = new Allocator< FEMTreeNode >();
		nodeAllocator->set( blockSize );
	}
	else nodeAllocator = NULL;
	_nodeCount = 0;
	_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
	_tree->initChildren( nodeAllocator , _NodeInitializer( *this ) ) , _spaceRoot = _tree->children;
	int offset[Dim];
	for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
	RegularTreeNode< Dim , FEMTreeNodeData >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
	_depthOffset = 0;
	memset( _femSigs1 , -1 , sizeof( _femSigs1 ) );
	memset( _femSigs2 , -1 , sizeof( _femSigs2 ) );
	memset( _refinableSigs , -1 , sizeof( _refinableSigs ) );
}
template< unsigned int Dim , class Real >
FEMTree< Dim , Real >::FEMTree( FILE* fp , int blockSize )
{
	if( blockSize>0 )
	{
		nodeAllocator = new Allocator< FEMTreeNode >();
		nodeAllocator->set( blockSize );
	}
	else nodeAllocator = NULL;
	if( fp )
	{
		if( fread( &_depthOffset , sizeof( int ) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] FEMTree::FEMTree: failed to read depth offset\n" ) , exit( 0 );
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
		_tree->read( fp , nodeAllocator , _NodeInitializer( *this ) );
		_maxDepth = _tree->maxDepth() - _depthOffset;

		_spaceRoot = _tree->children;

		if( _depthOffset>1 )
		{
			_spaceRoot = _tree->children + (1<<Dim)-1;
			for( int d=1 ; d<_depthOffset ; d++ )
				if( !_spaceRoot->children ) fprintf( stderr , "[ERROR] FEMTree::FEMTree expected children\n" ) , exit( 0 );
				else _spaceRoot = _spaceRoot->children;
		}
		_sNodes.set( *_tree , NULL );
	}
	else
	{
		_tree = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
		_tree->initChildren( nodeAllocator , _NodeInitializer( *this ) ) , _spaceRoot = _tree->children;
		int offset[Dim];
		for( int d=0 ; d<Dim ; d++ ) offset[d] = 0;
		RegularTreeNode< Dim , FEMTreeNodeData >::ResetDepthAndOffset( _spaceRoot , 0 , offset );
		_depthOffset = 0;
	}
}
template< unsigned int Dim , class Real > void FEMTree< Dim , Real >::write( FILE* fp ) const
{
	fwrite( &_depthOffset , sizeof( int ) , 1 , fp );
	_tree->write( fp );
}

template< unsigned int Dim , class Real >
const RegularTreeNode< Dim , FEMTreeNodeData >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p ) const
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
RegularTreeNode< Dim , FEMTreeNodeData >* FEMTree< Dim , Real >::leaf( Point< Real , Dim > p , LocalDepth maxDepth )
{
	if( !_InBounds( p ) ) return NULL;
	Point< Real , Dim > center;
	for( int d=0 ; d<Dim ; d++ ) center[d] = (Real)0.5;
	Real width = Real(1.0);
	FEMTreeNode* node = _spaceRoot;
	LocalDepth d = _localDepth( node );
	while( ( d<0 && node->children ) || ( d>=0 && d<maxDepth ) )
	{
		if( !node->children ) node->initChildren( nodeAllocator , _NodeInitializer( *this ) );
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
template< unsigned int ... Degrees >
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , FEMTreeNode* node , LocalDepth depth )
{
	LocalDepth d ; LocalOffset off;
	_localDepthAndOffset( node , d , off );
	bool refine = d<depth && ( d<0 || !FEMIntegrator::IsOutOfBounds( UIntPack< FEMDegreeAndBType< Degrees , BOUNDARY_FREE >::Signature ... >() , d , off ) );
	if( refine )
	{
		if( !node->children ) node->initChildren( nodeAllocator , _NodeInitializer( *this ) );
		for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , node->children+c , depth );
	}
}
template< unsigned int Dim , class Real >
template< unsigned int ... Degrees >
void FEMTree< Dim , Real >::_setFullDepth( UIntPack< Degrees ... > , LocalDepth depth )
{
	if( !_tree->children ) _tree->initChildren( nodeAllocator , _NodeInitializer( *this ) );
	for( int c=0 ; c<(1<<Dim) ; c++ ) _setFullDepth( UIntPack< Degrees ... >() , _tree->children+c , depth );
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
template< unsigned int LeftRadius , unsigned int RightRadius , class ... DenseOrSparseNodeData > 
void FEMTree< Dim , Real >::thicken( FEMTreeNode** nodes , size_t nodeCount, DenseOrSparseNodeData* ... data )
{
	std::vector< int > map( _nodeCount );
	for( int i=0 ; i<_nodeCount ; i++ ) map[i] = i;
	{
		int d=0 , off[Dim];
		for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
		FEMTreeNode::ResetDepthAndOffset( _tree , d , off );
	}
	typename RegularTreeNode< Dim , FEMTreeNodeData >::template NeighborKey< IsotropicUIntPack< Dim , LeftRadius > , IsotropicUIntPack< Dim , RightRadius > > neighborKey;
	neighborKey.set( _tree->maxDepth() );
	for( int i=0 ; i<nodeCount ; i++ ) neighborKey.template getNeighbors< true >( nodes[i] , nodeAllocator , _NodeInitializer( *this ) );
	{
		int d=0 , off[Dim];
		for( int d=0 ; d<Dim ; d++ ) off[d] = 0;
		FEMTreeNode::ResetDepthAndOffset( _spaceRoot , d , off );
	}

	_reorderDenseOrSparseNodeData( &map[0] , _nodeCount , data ... );
}
template< unsigned int Dim , class Real >
template< unsigned int LeftRadius , unsigned int RightRadius , class IsThickenNode , class ... DenseOrSparseNodeData > 
void FEMTree< Dim , Real >::thicken( IsThickenNode F , DenseOrSparseNodeData* ... data )
{
	std::vector< FEMTreeNode* > nodes;
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( IsActiveNode( node ) && F( node ) ) nodes.push_back( node );
	thicken< LeftRadius , RightRadius >( &nodes[0] , nodes.size() , data ... );
}

template< unsigned int Dim , class Real >
template< unsigned int DensityDegree >
typename FEMTree< Dim , Real >::template DensityEstimator< DensityDegree >* FEMTree< Dim , Real >::setDensityEstimator( const std::vector< PointSample >& samples , LocalDepth splatDepth , Real samplesPerNode , int coDimension )
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	splatDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( splatDepth , maxDepth ) );
	DensityEstimator< DensityDegree >* _density = new DensityEstimator< DensityDegree >( splatDepth , coDimension );
	DensityEstimator< DensityDegree >& density = *_density;
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	densityKey.set( _localToGlobal( splatDepth ) );

	std::vector< int > sampleMap( nodeCount() , -1 );
#pragma omp parallel for
	for( int i=0 ; i<samples.size() ; i++ ) if( samples[i].sample.weight>0 ) sampleMap[ samples[i].node->nodeData.nodeIndex ] = i;
	std::function< ProjectiveData< Point< Real , Dim > , Real > ( FEMTreeNode* ) > SetDensity = [&] ( FEMTreeNode* node )
	{
		ProjectiveData< Point< Real , Dim > , Real > sample;
		LocalDepth d = _localDepth( node );
		int idx = node->nodeData.nodeIndex;
		if( node->children )
			for( int c=0 ; c<(1<<Dim) ; c++ )
			{
				ProjectiveData< Point< Real , Dim > , Real > s = SetDensity( node->children + c );
				if( d<=splatDepth && s.weight>0 )
				{
					Point< Real , Dim > p = s.data / s.weight;
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
				Point< Real , Dim > p = sample.data / sample.weight;
				Real w = sample.weight / samplesPerNode;
				_addWeightContribution( density , node , p , densityKey , w );
			}
		}
		return sample;
	};
	SetDensity( _spaceRoot );

	MemoryUsage();
	return _density;
}
template< unsigned int Dim , class Real >
template< unsigned int ... NormalSigs , unsigned int DensityDegree , class Data >
SparseNodeData< Point< Real , Dim > , UIntPack< NormalSigs ... > > FEMTree< Dim , Real >::setNormalField( UIntPack< NormalSigs ... > , const std::vector< PointSample >& samples , const std::vector< Data >& normalData , const DensityEstimator< DensityDegree >* density , Real& pointWeightSum , std::function< Real ( Real ) > BiasFunction )
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	typedef PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > DensityKey;
	typedef UIntPack< FEMSignature< NormalSigs >::Degree ... > NormalDegrees;
	typedef PointSupportKey< UIntPack< FEMSignature< NormalSigs >::Degree ... > > NormalKey;
	std::vector< DensityKey > densityKeys( omp_get_max_threads() );
	std::vector<  NormalKey >  normalKeys( omp_get_max_threads() );
	bool oneKey = DensityDegree==NormalDegrees::Min() && DensityDegree==NormalDegrees::Max();
	for( int i=0 ; i<densityKeys.size() ; i++ ) densityKeys[i].set( _localToGlobal( maxDepth ) );
	if( !oneKey ) for( int i=0 ; i<normalKeys.size() ; i++ ) normalKeys[i].set( _localToGlobal( maxDepth ) );

	Real weightSum = 0;
	pointWeightSum = 0;
	SparseNodeData< Point< Real , Dim > , UIntPack< NormalSigs ... > > normalField;
	Real _pointWeightSum = 0;
#pragma omp parallel for reduction( + : weightSum , _pointWeightSum )
	for( int i=0 ; i<samples.size() ; i++ )
	{
		DensityKey& densityKey = densityKeys[ omp_get_thread_num() ];
		NormalKey& normalKey = normalKeys[ omp_get_thread_num() ];
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		if( sample.weight>0 )
		{
			Point< Real , Dim > p = sample.data / sample.weight , n = std::get< 0 >( normalData[i].data ).data;
			Real l = (Real)Length( n );
			// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
			if( !l ) continue;
			Real confidence = l / sample.weight;
			n *= sample.weight / l;
			Real depthBias = BiasFunction( confidence );
			weightSum += sample.weight;
			if( !_InBounds(p) )
			{
				fprintf( stderr , "[WARNING] FEMTree:setNormalField: Point sample is out of bounds:" );
				for( int d=0 ; d<Dim ; d++ ) fprintf( stderr , " %g" , p[d] );
				fprintf( stderr , "\n" );
				continue;
			}
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#else // !__GNUC__ || __GNUC__ >=5
			if( density ) _pointWeightSum += _splatPointData< true , DensityDegree , Point< Real , Dim > , NormalSigs ... >( *density , p , n , normalField , densityKey , oneKey ? *( (NormalKey*)&densityKey ) : normalKey , 0 , maxDepth , Dim , depthBias ) * sample.weight;
#endif // __GNUC__ || __GNUC__ < 4
			else
			{
				Real width = (Real)( 1.0 / ( 1<<maxDepth ) );
#if defined( __GNUC__ ) && __GNUC__ < 5
#warning "you've got me gcc version<5"
				_splatPointData< true , Point< Real , Dim > >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#else // !__GNUC__ || __GNUC__ >=5
				_splatPointData< true , Point< Real , Dim > , NormalSigs ... >( leaf( p , maxDepth ) , p , n / (Real)pow( width , Dim ) , normalField , oneKey ? *( (NormalKey*)&densityKey ) : normalKey );
#endif // __GNUC__ || __GNUC__ < 4
				_pointWeightSum += sample.weight;
			}
		}
	}
	pointWeightSum = _pointWeightSum / weightSum;
	MemoryUsage();
	return normalField;
}
template< unsigned int Dim , class Real >
template< unsigned int DataSig , bool CreateNodes , unsigned int DensityDegree , class Data >
SparseNodeData< Data , IsotropicUIntPack< Dim , DataSig > > FEMTree< Dim , Real >::setSingleDepthDataField( const std::vector< PointSample >& samples , const std::vector< Data >& sampleData , const DensityEstimator< DensityDegree >* density )
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	PointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< Data , IsotropicUIntPack< Dim , DataSig > > dataField;
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		const Data& data = sampleData[i];
		Point< Real , Dim > p = sample.weight==0 ? sample.data : sample.data / sample.weight;
		if( !_InBounds(p) )
		{
			fprintf( stderr , "[WARNING] Point is out of bounds:" );
			for( int d=0 ; d<Dim ; d++ ) fprintf( stderr , " %f" , p[d] );
			fprintf( stderr , " <-" );
			for( int d=0 ; d<Dim ; d++ ) fprintf( stderr , " %f" , sample.data[d] );
			fprintf( stderr , " [%f]\n" , sample.weight );
			continue;
		}
		if( density ) _splatPointData< CreateNodes , DensityDegree , DataSig >( *density             , p , data * sample.weight , dataField , densityKey , dataKey , 0 , maxDepth , Dim );
		else          _splatPointData< CreateNodes ,                 DataSig >( leaf( p , maxDepth ) , p , data * sample.weight , dataField , dataKey );
	}
	MemoryUsage();
	return dataField;
}
template< unsigned int Dim , class Real >
template< unsigned int DataSig , bool CreateNodes , unsigned int DensityDegree , class Data >
SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > FEMTree< Dim , Real >::setDataField( const std::vector< PointSample >& samples , std::vector< Data >& sampleData , const DensityEstimator< DensityDegree >* density , bool nearest )
{
	LocalDepth maxDepth = _spaceRoot->maxDepth();
	PointSupportKey< IsotropicUIntPack< Dim , DensityDegree > > densityKey;
	PointSupportKey< IsotropicUIntPack< Dim , FEMSignature< DataSig >::Degree > > dataKey;
	densityKey.set( _localToGlobal( maxDepth ) ) , dataKey.set( _localToGlobal( maxDepth ) );

	SparseNodeData< ProjectiveData< Data , Real > , IsotropicUIntPack< Dim , DataSig > > dataField;
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const ProjectiveData< Point< Real , Dim > , Real >& sample = samples[i].sample;
		const Data& data = sampleData[i];
		Point< Real , Dim > p = sample.weight==0 ? sample.data : sample.data / sample.weight;
		if( !_InBounds(p) )
		{
			fprintf( stderr , "[WARNING] Point is out of bounds:" );
			for( int d=0 ; d<Dim ; d++ ) fprintf( stderr , " %f" , p[d] );
			fprintf( stderr , " <-" );
			for( int d=0 ; d<Dim ; d++ ) fprintf( stderr , " %f" , sample.data[d] );
			fprintf( stderr , " [%f]\n" , sample.weight );
			continue;
		}
		if( nearest ) _nearestMultiSplatPointData< DensityDegree >( density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , 2 );
		else          _multiSplatPointData< CreateNodes , DensityDegree >( density , (FEMTreeNode*)samples[i].node , p , ProjectiveData< Data , Real >( data , sample.weight ) , dataField , densityKey , dataKey , 2 );
	}
	MemoryUsage();
	return dataField;
}
template< unsigned int Dim , class Real >
template< unsigned int MaxDegree , class HasDataFunctor , class ... DenseOrSparseNodeData >
void FEMTree< Dim , Real >::finalizeForMultigrid( LocalDepth fullDepth , const HasDataFunctor F , DenseOrSparseNodeData* ... data )
{
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

		FEMTreeNode* newSpaceRootParent = FEMTreeNode::NewBrood( nodeAllocator , _NodeInitializer( *this ) );
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
	// Make the low-resolution part of the tree be complete
	fullDepth = std::max< LocalDepth >( 0 , std::min< LocalDepth >( _maxDepth , fullDepth ) );
	_setFullDepth( IsotropicUIntPack< Dim , MaxDegree >() , fullDepth );
	// Clear all the flags and make everything that is not low-res a ghost node
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) node->nodeData.flags = 0 , SetGhostFlag< Dim >( node , _localDepth( node )>fullDepth );

	// Set the ghost nodes for the high-res part of the tree
	_clipTree( F , fullDepth );

	const int OverlapRadius = -BSplineOverlapSizes< MaxDegree , MaxDegree >::OverlapStart;
	int maxDepth = _tree->maxDepth( );
	typedef typename FEMTreeNode::template NeighborKey< IsotropicUIntPack< Dim , OverlapRadius > , IsotropicUIntPack< Dim , OverlapRadius > > NeighborKey;

	std::vector< NeighborKey > neighborKeys( omp_get_max_threads() );
	for( int i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth-1 ) );

	for( LocalDepth d=_maxDepth-1 ; d>=0 ; d-- )
	{
		std::vector< FEMTreeNode* > nodes;
		for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( _localDepth( node )==d && IsActiveNode< Dim >( node->children ) ) nodes.push_back( node );
#pragma omp parallel for
		for( int i=0 ; i<nodes.size() ; i++ )
		{
			NeighborKey& neighborKey = neighborKeys[ omp_get_thread_num() ];
			FEMTreeNode* node = nodes[i];
			neighborKey.template getNeighbors< true >( node , nodeAllocator , _NodeInitializer( *this ) );
			Pointer( FEMTreeNode* ) nodes = neighborKey.neighbors[ _localToGlobal(d) ].neighbors().data;
			unsigned int size = neighborKey.neighbors[ _localToGlobal(d) ].neighbors.Size;
			for( unsigned int i=0 ; i<size ; i++ ) SetGhostFlag< Dim >( nodes[i] , false );
		}
	}
	std::vector< int > map;
	_sNodes.set( *_tree , &map );
	_setSpaceValidityFlags();
	for( FEMTreeNode* node=_tree->nextNode() ; node ; node=_tree->nextNode( node ) ) if( !IsActiveNode< Dim >( node ) ) node->nodeData.nodeIndex = -1;
	_reorderDenseOrSparseNodeData( &map[0] , _sNodes.size() , data ... );
	MemoryUsage();
}

template< unsigned int Dim , class Real >
void FEMTree< Dim , Real >::_setSpaceValidityFlags( void ) const
{
	for( int i=0 ; i<_sNodes.size() ; i++ )
	{
		const unsigned char MASK = ~( FEMTreeNodeData::SPACE_FLAG );
		_sNodes.treeNodes[i]->nodeData.flags &= MASK;
		if( isValidSpaceNode( _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::SPACE_FLAG;
	}
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs1 >
void FEMTree< Dim , Real >::_setFEM1ValidityFlags( UIntPack< FEMSigs1 ... > ) const
{
	bool needToReset;
	unsigned int femSigs1[] = { FEMSigs1 ... };
#pragma omp critical (set_fem_1_validity_flags)
	{
		needToReset = memcmp( femSigs1 , _femSigs1 , sizeof( _femSigs1 ) )!=0;
		if( needToReset ) memcpy( _femSigs1 , femSigs1 , sizeof( _femSigs1 ) );
	}
	if( needToReset )
		for( int i=0 ; i<_sNodes.size() ; i++ )
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
#pragma omp critical (set_fem_2_validity_flags)
	{
		needToReset = memcmp( femSigs2 , _femSigs2 , sizeof( _femSigs2 ) )!=0;
		if( needToReset ) memcpy( _femSigs2 , femSigs2 , sizeof( _femSigs2 ) );
	}
	if( needToReset )
		for( int i=0 ; i<_sNodes.size() ; i++ )
		{
			const unsigned char MASK = ~( FEMTreeNodeData::FEM_FLAG_2 );
			_sNodes.treeNodes[i]->nodeData.flags &= MASK;
			if( isValidFEMNode( UIntPack< FEMSigs2 ... >() , _sNodes.treeNodes[i] ) ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::FEM_FLAG_2;
		}
}
template< unsigned int Dim , class Real >
template< unsigned int ... FEMSigs >
void FEMTree< Dim , Real >::_setRefinabilityFlags( UIntPack< FEMSigs ... > ) const
{
	bool needToReset;
	unsigned int refinableSigs[] = { FEMSigs ... };
#pragma omp critical (set_refinability_flags)
	{
		needToReset = memcmp( refinableSigs , _refinableSigs , sizeof( _refinableSigs ) )!=0;
		if( needToReset ) memcpy( _refinableSigs , refinableSigs , sizeof( _refinableSigs ) );
	}
	if( needToReset )
	{
		typedef typename FEMTreeNode::template ConstNeighborKey< UIntPack< ( - BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ) ... > , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleEnd ... > > UpSampleKey;
		typedef typename FEMTreeNode::template ConstNeighbors< UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... > > UpSampleNeighbors;
		static const int UpSampleStart[] = { BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleStart ... };
		std::vector< UpSampleKey > neighborKeys( omp_get_max_threads() );
		for( size_t i=0 ; i<neighborKeys.size() ; i++ ) neighborKeys[i].set( _localToGlobal( _maxDepth ) );

		for( int d=0 ; d<_maxDepth ; d++ )
#pragma omp parallel for
			for( int i=_sNodesBegin(d) ; i<_sNodesEnd(d) ; i++ )
			{
				UpSampleKey& neighborKey = neighborKeys[ omp_get_thread_num() ];

				// Clear the refinability flag
				const unsigned char MASK = ~( FEMTreeNodeData::REFINABLE_FLAG );
				_sNodes.treeNodes[i]->nodeData.flags &= MASK;

				LocalDepth d ; LocalOffset pOff;
				_localDepthAndOffset( _sNodes.treeNodes[i] , d , pOff );

				// Get the supporting child neighbors
				neighborKey.getNeighbors( _sNodes.treeNodes[i] );
				UpSampleNeighbors neighbors;
				neighborKey.getChildNeighbors( 0 , _localToGlobal( d ) , neighbors );

				// Check if the child neighbors exist (i.e. that the children nodes are not ghost-nodes if they correspond to valid coefficients)
				bool refinable = true;
				LocalOffset cOff;
				WindowLoop< Dim >::Run
				(
					IsotropicUIntPack< Dim , 0 >() , UIntPack< BSplineSupportSizes< FEMSignature< FEMSigs >::Degree >::UpSampleSize ... >() ,
					[&]( int d , int i ){ cOff[d] = pOff[d]*2 + UpSampleStart[d] + i; } ,
					[&]( const FEMTreeNode* node ){ if( GetGhostFlag< Dim >( node ) && FEMIntegrator::IsValidFEMNode( UIntPack< FEMSigs ... >() , d+1 , cOff ) ) refinable = false; } ,
					neighbors.neighbors()
				);
				if( refinable ) _sNodes.treeNodes[i]->nodeData.flags |= FEMTreeNodeData::REFINABLE_FLAG;
			}
	}
}
template< unsigned int Dim , class Real >
template< class HasDataFunctor >
void FEMTree< Dim , Real >::_clipTree( const HasDataFunctor& f , LocalDepth fullDepth )
{
	for( FEMTreeNode* temp=_tree->nextNode() ; temp ; temp=_tree->nextNode(temp) ) if( temp->children && _localDepth( temp )>=fullDepth )
	{
		bool hasData = false;
		for( int c=0 ; c<(1<<Dim) && !hasData ; c++ ) hasData |= f( temp->children + c );
		for( int c=0 ; c<(1<<Dim) ; c++ ) SetGhostFlag< Dim >( temp->children+c , !hasData );
	}
}

template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::_ExactPointAndDataInterpolationInfo< T , Data , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , int& ) > SetRange = [&] ( FEMTreeNode* node , int& start )
	{
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
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

	int start = 0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)_iData.size() ; i++ )
	{
		Real w = _iData[i].pointInfo.weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].pointInfo.weight = w;
		else            _iData[i].pointInfo.weight = w * ( 1<<tree._maxDepth );
		_iData[i].pointInfo.dualValues *= _iData[i].pointInfo.weight;
	}
}
template< unsigned int Dim , class Real >
template< typename T , unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< T , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , int& ) > SetRange = [&] ( FEMTreeNode* node , int& start )
	{
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
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

	int start = 0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)_iData.size() ; i++ )
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
}
template< unsigned int Dim , class Real >
template< unsigned int PointD , typename ConstraintDual , typename SystemDual >
void FEMTree< Dim , Real >::ExactPointInterpolationInfo< double , PointD , ConstraintDual , SystemDual >::_init( const class FEMTree< Dim , Real >& tree , const std::vector< PointSample >& samples , bool noRescale )
{
	_sampleSpan.resize( tree.nodesSize() );
#pragma omp parallel for
	for( int i=0 ; i<tree.nodesSize() ; i++ ) _sampleSpan[i] = std::pair< int , int >( 0 , 0 );
	for( int i=0 ; i<samples.size() ; i++ )
	{
		const FEMTreeNode* leaf = samples[i].node;
		while( leaf && !tree._isValidSpaceNode( leaf ) ) leaf = leaf->parent;
		if( leaf && tree._isValidSpaceNode( leaf ) ) _sampleSpan[ leaf->nodeData.nodeIndex ].second++;
	}
	_iData.resize( samples.size() );

	std::function< void ( FEMTreeNode* , int& ) > SetRange = [&] ( FEMTreeNode* node , int& start )
	{
		std::pair< int , int >& span = _sampleSpan[ node->nodeData.nodeIndex ];
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

	int start = 0;
	SetRange( tree._spaceRoot , start );
	for( FEMTreeNode* node=tree._spaceRoot->nextNode() ; node ; node=tree._spaceRoot->nextNode(node) )
		if( tree._isValidSpaceNode( node ) && !tree._isValidSpaceNode( node->children ) ) _sampleSpan[ node->nodeData.nodeIndex ].second = _sampleSpan[ node->nodeData.nodeIndex ].first;

	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)_iData.size() ; i++ )
	{
		Real w = _iData[i].weight;
		_iData[i] /= w;
		if( noRescale ) _iData[i].weight = w;
		else            _iData[i].weight = w * ( 1<<tree._maxDepth );
		_iData[i].dualValues *= _iData[i].weight;
	}
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
	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
	{
		Real w = iInfo[i].weight;
		iInfo[i] /= w ; iInfo[i].weight = w;
	}
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
	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
	{
		Real w = iInfo[i].pointInfo.weight;
		iInfo[i] /= w ; iInfo[i].pointInfo.weight = w;
	}
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
	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
	{
		iInfo[i].finalize();
		for( int c=0 ; c<(int)iInfo[i].size() ; c++ )
		{
			iInfo[i][c].position /= iInfo[i][c].weight;
			if( !noRescale )
			{
				iInfo[i][c].weight     *= ( 1<<_maxDepth );
				iInfo[i][c].dualValues *= ( 1<<_maxDepth );
			}
		}
	}
	return iInfo;
}
template< unsigned int Dim , class Real >
template< typename T , typename Data , unsigned int PointD , typename ConstraintDual >
SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > FEMTree< Dim , Real >::_densifyChildInterpolationInfoAndSetDualConstraints( const std::vector< PointSample >& samples , ConstPointer( Data ) sampleData , ConstraintDual constraintDual , bool noRescale ) const
{
	SparseNodeData< DualPointAndDataInfoBrood< Dim , Real , Data , T , PointD > , IsotropicUIntPack< Dim , FEMTrivialSignature > > iInfo;
	for( int i=0 ; i<samples.size() ; i++ )
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

#pragma omp parallel for
	for( int i=0 ; i<(int)iInfo.size() ; i++ )
	{
		iInfo[i].finalize();
		for( int c=0 ; c<(int)iInfo[i].size() ; c++ )
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
	return iInfo;
}



template< unsigned int Dim , class Real >
std::vector< int > FEMTree< Dim , Real >::merge( FEMTree* tree )
{
	std::vector< int > map;
	if( _depthOffset!=tree->_depthOffset ) fprintf( stderr , "[ERROR] FEMTree::merge: depthOffsets don't match: %d != %d\n" , _depthOffset , tree->_depthOffset ) , exit( 0 );

	// Compute the next available index
	int nextIndex = 0;
	for( const FEMTreeNode* node=_tree->nextNode() ; node!=NULL ; node=_tree->nextNode( node ) ) nextIndex = std::max< int >( nextIndex , node->nodeData.nodeIndex+1 );

	// Set the size of the map
	{
		int mapSize = 0;
		for( const FEMTreeNode* node=tree->_tree->nextNode() ; node!=NULL ; node=tree->_tree->nextNode( node ) ) mapSize = std::max< int >( mapSize , node->nodeData.nodeIndex+1 );
		map.resize( mapSize );
	}

	std::function< void ( FEMTreeNode* , FEMTreeNode* , std::vector< int >& , int& ) > MergeNodes = [&]( FEMTreeNode* node1 , FEMTreeNode* node2 , std::vector< int >& map , int& nextIndex )
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

