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

#ifndef REGULAR_TREE_NODE_INCLUDED
#define REGULAR_TREE_NODE_INCLUDED

#include "Allocator.h"
#include "BinaryNode.h"
#include "Window.h"
#include <functional>

#ifdef USE_DEEP_TREE_NODES
template< unsigned int Dim , class NodeData , class DepthAndOffsetType=unsigned int >
#else // !USE_DEEP_TREE_NODES
template< unsigned int Dim , class NodeData , class DepthAndOffsetType=unsigned short >
#endif // USE_DEEP_TREE_NODES
struct RegularTreeNode
{
private:
	DepthAndOffsetType _depth , _offset[Dim];
public:

	RegularTreeNode* parent;
	RegularTreeNode* children;
	NodeData nodeData;

	RegularTreeNode( std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );
	static RegularTreeNode* NewBrood( Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );
	static void ResetDepthAndOffset( RegularTreeNode* root , int d , int off[Dim] );
	int initChildren( Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );

	void cleanChildren( Allocator< RegularTreeNode >* nodeAllocator );
	~RegularTreeNode( void );

	// The merge functor takes two objects of type NodeData and returns an object of type NodeData
	// [NOTE] We are assuming that the merge functor is symmetric, f(a,b) = f(b,a), and implicity satisfies f(a) = a
	template< class MergeFunctor >
	void merge( RegularTreeNode* node , MergeFunctor& f );

	void depthAndOffset( int& depth , int offset[Dim] ) const; 
	void centerIndex( int index[Dim] ) const;
	int depth( void ) const;
	template< class Real > void centerAndWidth( Point< Real , Dim >& center , Real& width ) const;
	template< class Real > void startAndWidth( Point< Real , Dim >& start , Real& width ) const;
	template< class Real > bool isInside( Point< Real , Dim > p ) const;

	size_t leaves( void ) const;
	size_t maxDepthLeaves( int maxDepth ) const;
	size_t nodes( void ) const;
	int maxDepth( void ) const;

	const RegularTreeNode* root( void ) const;

	const RegularTreeNode* nextLeaf( const RegularTreeNode* currentLeaf=NULL ) const;
	RegularTreeNode* nextLeaf( RegularTreeNode* currentLeaf=NULL );
	const RegularTreeNode* nextNode( const RegularTreeNode* currentNode=NULL ) const;
	RegularTreeNode* nextNode( RegularTreeNode* currentNode=NULL );
	const RegularTreeNode* nextBranch( const RegularTreeNode* current ) const;
	RegularTreeNode* nextBranch( RegularTreeNode* current );
	const RegularTreeNode* prevBranch( const RegularTreeNode* current ) const;
	RegularTreeNode* prevBranch( RegularTreeNode* current );

	void setFullDepth( int maxDepth , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );

	void printLeaves( void ) const;
	void printRange( void ) const;

	template< class Real > static int ChildIndex( const Point< Real , Dim >& center , const Point< Real , Dim > &p );

	bool write( const char* fileName ) const;
	bool write( FILE* fp ) const;
	bool read( const char* fileName , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );
	bool read( FILE* fp , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );

	template< typename Pack > struct Neighbors{};
	template< unsigned int ... Widths >
	struct Neighbors< UIntPack< Widths ... > >
	{
		typedef StaticWindow< RegularTreeNode* , UIntPack< Widths ... > > Window;
		Window neighbors;
		Neighbors( void );
		void clear( void );
	};
	template< typename Pack > struct ConstNeighbors{};
	template< unsigned int ... Widths >
	struct ConstNeighbors< UIntPack< Widths ... > >
	{
		typedef StaticWindow< const RegularTreeNode* , UIntPack< Widths ... > > Window;
		Window neighbors;
		ConstNeighbors( void );
		void clear( void );
	};

	template< typename LeftPack , typename RightPack > struct NeighborKey{};
	template< unsigned int ... LeftRadii , unsigned int ... RightRadii >
	struct NeighborKey< UIntPack< LeftRadii ... > , UIntPack< RightRadii ... > >
	{
	protected:
		static_assert( sizeof...(LeftRadii)==sizeof...(RightRadii) , "[ERROR] Left and right radii dimensions don't match" );
		static const unsigned int CenterIndex = WindowIndex< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > , UIntPack< LeftRadii ... > >::Index;
		int _depth;

		template< bool CreateNodes , unsigned int ... _PLeftRadii , unsigned int ... _PRightRadii , unsigned int ... _CLeftRadii , unsigned int ... _CRightRadii >
		static unsigned int _NeighborsLoop( UIntPack< _PLeftRadii ... > , UIntPack< _PRightRadii ... > , UIntPack< _CLeftRadii ... > , UIntPack< _CRightRadii ... > , ConstWindowSlice< RegularTreeNode* , UIntPack< ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< RegularTreeNode* , UIntPack< ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int cIdx , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer );
		template< bool CreateNodes , unsigned int ... _PLeftRadii , unsigned int ... _PRightRadii , unsigned int ... _CLeftRadii , unsigned int ... _CRightRadii >
		static unsigned int _NeighborsLoop( UIntPack< _PLeftRadii ... > , UIntPack< _PRightRadii ... > , UIntPack< _CLeftRadii ... > , UIntPack< _CRightRadii ... > ,      WindowSlice< RegularTreeNode* , UIntPack< ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< RegularTreeNode* , UIntPack< ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int cIdx , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer );

		template< bool CreateNodes , typename PLeft , typename PRight , typename CLeft , typename CRight > struct _Run{};

		template< bool CreateNodes , unsigned int _PLeftRadius , unsigned int ... _PLeftRadii , unsigned int _PRightRadius , unsigned int ... _PRightRadii , unsigned int _CLeftRadius , unsigned int ... _CLeftRadii , unsigned int _CRightRadius , unsigned int ... _CRightRadii >
		struct _Run< CreateNodes , UIntPack< _PLeftRadius , _PLeftRadii ... > , UIntPack< _PRightRadius , _PRightRadii ... > , UIntPack< _CLeftRadius , _CLeftRadii ... > , UIntPack< _CRightRadius , _CRightRadii ... > >
		{
			static unsigned int Run( ConstWindowSlice< RegularTreeNode* , UIntPack< _PLeftRadius+_PRightRadius+1 , ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< RegularTreeNode* , UIntPack< _CLeftRadius+_CRightRadius+1 , ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int* c , int cornerIndex , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer );
		};
		template< bool CreateNodes , unsigned int _PLeftRadius , unsigned int _PRightRadius , unsigned int _CLeftRadius , unsigned int _CRightRadius >
		struct _Run< CreateNodes , UIntPack< _PLeftRadius > , UIntPack< _PRightRadius > , UIntPack< _CLeftRadius > , UIntPack< _CRightRadius > >
		{
			static unsigned int Run( ConstWindowSlice< RegularTreeNode* , UIntPack< _PLeftRadius+_PRightRadius+1 > > pNeighbors , WindowSlice< RegularTreeNode* , UIntPack< _CLeftRadius+_CRightRadius+1 > > cNeighbors , int* c , int cornerIndex , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer );
		};
	public:
		typedef Neighbors< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > > NeighborType;
		NeighborType* neighbors;


		NeighborKey( void );
		NeighborKey( const NeighborKey& key );
		~NeighborKey( void );
		int depth( void ) const { return _depth; }

		void set( int depth );

		template< bool CreateNodes >
		typename RegularTreeNode< Dim , NodeData , DepthAndOffsetType >::template Neighbors< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > >& getNeighbors( RegularTreeNode* node , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );

		NeighborType& getNeighbors( const RegularTreeNode* node ) { return getNeighbors< false >( (RegularTreeNode*)node , NULL , std::function< void ( RegularTreeNode& ) >() ); }

		template< bool CreateNodes , unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > ,       RegularTreeNode* node , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );
		template< unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > , const RegularTreeNode* node , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors ){ return getNeighbors< false >( UIntPack< _LeftRadii ... >() , UIntPack< _RightRadii ... >() , (RegularTreeNode*)node , NULL , std::function< void ( RegularTreeNode& ) >() ); }
		template< bool CreateNodes , unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > ,       RegularTreeNode* node , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& pNeighbors , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() );
		template< unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > , const RegularTreeNode* node , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& pNeighbors , Neighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors ){ return getNeighbors< false >( UIntPack< _LeftRadii ... >() , UIntPack< _RightRadii ... >() , (RegularTreeNode*)node , NULL , std::function< void ( RegularTreeNode& ) >() ); }

		template< bool CreateNodes >
		unsigned int getChildNeighbors( int cIdx , int d , NeighborType& childNeighbors , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() ) const;
		unsigned int getChildNeighbors( int cIdx , int d , NeighborType& childNeighbors ) const { return getChildNeighbors< false >( cIdx , d , childNeighbors , NULL , std::function< void ( RegularTreeNode& ) >() ); }

		template< bool CreateNodes , class Real >
		unsigned int getChildNeighbors( Point< Real , Dim > p , int d , NeighborType& childNeighbors , Allocator< RegularTreeNode >* nodeAllocator , std::function< void ( RegularTreeNode& ) > Initializer=std::function< void ( RegularTreeNode& ) >() ) const;
		template< class Real >
		unsigned int getChildNeighbors( Point< Real , Dim > p , int d , NeighborType& childNeighbors ) const { return getChildNeighbors< false , Real >( p , d , childNeighbors , NULL , std::function< void ( RegularTreeNode& ) >() ); }
	};

	template< typename LeftPack , typename RightPack > struct ConstNeighborKey{};

	template< unsigned int ... LeftRadii , unsigned int ... RightRadii >
	struct ConstNeighborKey< UIntPack< LeftRadii ... > , UIntPack< RightRadii ... > >
	{
	protected:
		static_assert( sizeof...(LeftRadii)==sizeof...(RightRadii) , "[ERROR] Left and right radii dimensions don't match" );
		static const unsigned int CenterIndex = WindowIndex< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > , UIntPack< LeftRadii ... > >::Index;
		int _depth;

		template< unsigned int ... _PLeftRadii , unsigned int ... _PRightRadii , unsigned int ... _CLeftRadii , unsigned int ... _CRightRadii >
		static unsigned int _NeighborsLoop( UIntPack< _PLeftRadii ... > , UIntPack< _PRightRadii ... > , UIntPack< _CLeftRadii ... > , UIntPack< _CRightRadii ... > , ConstWindowSlice< const RegularTreeNode* , UIntPack< ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< const RegularTreeNode* , UIntPack< ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int cIdx );
		template< unsigned int ... _PLeftRadii , unsigned int ... _PRightRadii , unsigned int ... _CLeftRadii , unsigned int ... _CRightRadii >
		static unsigned int _NeighborsLoop( UIntPack< _PLeftRadii ... > , UIntPack< _PRightRadii ... > , UIntPack< _CLeftRadii ... > , UIntPack< _CRightRadii ... > , WindowSlice< const RegularTreeNode* , UIntPack< ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< const RegularTreeNode* , UIntPack< ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int cIdx );

		template< typename PLeft , typename PRight , typename CLeft , typename CRight > struct _Run{};

		template< unsigned int _PLeftRadius , unsigned int ... _PLeftRadii , unsigned int _PRightRadius , unsigned int ... _PRightRadii , unsigned int _CLeftRadius , unsigned int ... _CLeftRadii , unsigned int _CRightRadius , unsigned int ... _CRightRadii >
		struct _Run< UIntPack< _PLeftRadius , _PLeftRadii ... > , UIntPack< _PRightRadius , _PRightRadii ... > , UIntPack< _CLeftRadius , _CLeftRadii ... > , UIntPack< _CRightRadius , _CRightRadii ... > >
		{
			static unsigned int Run( ConstWindowSlice< const RegularTreeNode* , UIntPack< _PLeftRadius + _PRightRadius + 1 , ( _PLeftRadii+_PRightRadii+1 ) ... > > pNeighbors , WindowSlice< const RegularTreeNode* , UIntPack< _CLeftRadius + _CRightRadius + 1 , ( _CLeftRadii+_CRightRadii+1 ) ... > > cNeighbors , int* c , int cornerIndex );
		};
		template< unsigned int _PLeftRadius , unsigned int _PRightRadius , unsigned int _CLeftRadius , unsigned int _CRightRadius >
		struct _Run< UIntPack< _PLeftRadius > , UIntPack< _PRightRadius > , UIntPack< _CLeftRadius > , UIntPack< _CRightRadius > >
		{
			static unsigned int Run( ConstWindowSlice< const RegularTreeNode* , UIntPack< _PLeftRadius+_PRightRadius+1 > > pNeighbors , WindowSlice< const RegularTreeNode* , UIntPack< _CLeftRadius+_CRightRadius+1 > > cNeighbors , int* c , int cornerIndex );
		};

	public:

		typedef ConstNeighbors< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > > NeighborType;
		NeighborType* neighbors;

		ConstNeighborKey( void );
		ConstNeighborKey( const ConstNeighborKey& key );
		~ConstNeighborKey( void );
#if 1
		ConstNeighborKey& operator = ( const ConstNeighborKey& key );
#endif

		int depth( void ) const { return _depth; }
		void set( int depth );

		typename RegularTreeNode< Dim , NodeData , DepthAndOffsetType >::template ConstNeighbors< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > >& getNeighbors( const RegularTreeNode* node );
		template< unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > , const RegularTreeNode* node , ConstNeighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors );
		template< unsigned int ... _LeftRadii , unsigned int ... _RightRadii >
		void getNeighbors( UIntPack< _LeftRadii ... > , UIntPack< _RightRadii ... > , const RegularTreeNode* node , ConstNeighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& pNeighbors , ConstNeighbors< UIntPack< ( _LeftRadii + _RightRadii + 1 ) ... > >& neighbors );
		unsigned int getChildNeighbors( int cIdx , int d , NeighborType& childNeighbors ) const;
		template< class Real >
		unsigned int getChildNeighbors( Point< Real , Dim > p , int d , ConstNeighbors< UIntPack< ( LeftRadii + RightRadii + 1 ) ... > >& childNeighbors ) const;
	};

	int width( int maxDepth ) const;
};

#include "RegularTree.inl"

#endif // REGULAR_TREE_NODE_INCLUDED
