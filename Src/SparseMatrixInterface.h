#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define FORCE_TWO_BYTE_ALIGNMENT 1
#include "Array.h"

#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(push)
#pragma pack(2)
#endif // FORCE_TWO_BYTE_ALIGNMENT
template< class T , class IndexType=int >
struct MatrixEntry
{
	MatrixEntry( void )             { N =-1 , Value = 0; }
	MatrixEntry( IndexType i )      { N = i , Value = 0; }
	MatrixEntry( IndexType n , T v ){ N = n , Value = v; }
	IndexType N;
	T Value;
};
#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(pop)
#endif // FORCE_TWO_BYTE_ALIGNMENT


enum
{
	MULTIPLY_ADD = 1 ,
	MULTIPLY_NEGATE = 2
};

template< class T , class const_iterator > class SparseMatrixInterface
{
public:
	virtual const_iterator begin( size_t row ) const = 0;
	virtual const_iterator end  ( size_t row ) const = 0;
	virtual size_t rows   ( void )             const = 0;
	virtual size_t rowSize( size_t idx )       const = 0;

	size_t entries( void ) const;

	double squareNorm( void ) const;
	double squareASymmetricNorm( void ) const;
	double squareASymmetricNorm( int& idx1 , int& idx2 ) const;

	template< class T2 > void multiply      (           ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const;
	template< class T2 > void multiply      (                Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiply      (         ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }
	template< class T2 > void multiplyScaled( T scale ,      Pointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag=0 ) const { multiplyScaled( scale , ( ConstPointer(T2) )( In ) , Out , multiplyFlag ); }

	void setDiagonal( Pointer( T ) diagonal ) const;
	void setDiagonalR( Pointer( T ) diagonal ) const;
	template< class T2 > void jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const;
	template< class T2 > void gsIteration(                                                              ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
	template< class T2 > void gsIteration( const              std::vector< int >  & multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x ,                bool dReciprocal ) const;
	template< class T2 > void gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const;
};

// Assuming that the SPDOperator class defines:
//		auto SPDOperator::()( ConstPointer( T ) , Pointer( T ) ) const
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot );

#include "SparseMatrixInterface.inl"
#endif // SPARSE_MATRIX_INTERFACE_INCLUDED
