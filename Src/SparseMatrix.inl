/* -*- C++ -*-
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

#include <float.h>
#include <complex>
#include <unordered_map>

///////////////////
//  SparseMatrix //
///////////////////
///////////////////////////////////////
// SparseMatrix Methods and Memebers //
///////////////////////////////////////

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( void )
{
	rowSizes = NullPointer( size_t );
	rowNum = 0;
	_entries = NullPointer( Pointer( MatrixEntry< T , IndexType > ) );
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( size_t rowNum )
{
	this->rowNum = 0;
	rowSizes = NullPointer( size_t );
	_entries= NullPointer( Pointer( MatrixEntry< T , IndexType > ) );
	resize( rowNum );
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix& M )
{
	rowSizes = NullPointer( size_t );
	rowNum = 0;
	_entries = NullPointer( Pointer( MatrixEntry< T , IndexType > ) );
	resize( M.rowNum );
	for( int i=0 ; i<rowNum ; i++ )
	{
		setRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = M._entries[i][j];
	}
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( SparseMatrix&& M )
{
	rowSizes = NullPointer( size_t );
	rowNum = 0;
	_entries = NullPointer( Pointer( MatrixEntry< T , IndexType > ) );

	Swap( *this , M );
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M )
{
	rowSizes = NullPointer( size_t );
	rowNum = 0;
	_entries = NULL;
	resize( M.rowNum );
	for( int i=0 ; i<rowNum ; i++ )
	{
		setRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = MatrixEntry< T , IndexType >( M._entries[i][j].N , T( M._entries[i][j].Value ) );
	}
}

template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::copy( const SparseMatrix< T2 , IndexType2 >& M  )
{
	resize( M.rowNum );
	for ( int i=0 ; i<rowNum ; i++)
	{
		setRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ )
		{
			int idx = M._entries[i][j].N;
			_entries[i][j] = MatrixEntry< T , IndexType >( idx , T( M[i][j].Value ) );
		}
	}
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = ( SparseMatrix< T , IndexType >&& M )
{
	Swap( *this , M );
	return *this;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = ( const SparseMatrix< T , IndexType >& M )
{
	resize( M.rowNum );
	for( int i=0 ; i<rowNum ; i++ )
	{
		setRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j]=M._entries[i][j];
	}
	return *this;
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = (const SparseMatrix< T2 , IndexType2 >& M)
{
	resize( M.rowNum );
	for( int i=0 ; i<rowNum ; i++ )
	{
		setRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = MatrixEntry< T , IndexType >( M._entries[i][j].N , T( M._entries[i][j].Value ) );
	}
	return *this;
}

template< class T , class IndexType >
template< class T2 >
void SparseMatrix< T , IndexType >::operator() ( const T2* in , T2* out ) const { Interface::multiply( in , out ); }


template< class T , class IndexType > SparseMatrix< T , IndexType >::~SparseMatrix( void ) { resize( 0 ); }

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::resize( size_t r )
{
	if( rowNum>0 )
	{
		for( int i=0 ; i<rowNum ; i++ ) FreePointer( _entries[i] );
		FreePointer( _entries );
		FreePointer( rowSizes );
	}
	rowNum = r;
	if( r )
	{
		rowSizes = AllocPointer< size_t >( r ) , memset( rowSizes , 0 , sizeof(size_t)*r );
		_entries = AllocPointer< Pointer( MatrixEntry< T , IndexType > ) >( r );
		for( int i=0 ; i<r ; i++ ) _entries[i] = NullPointer( MatrixEntry< T , IndexType > );
	}
}

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::setRowSize( size_t row , size_t count )
{
	if( row>=0 && row<rowNum )
	{
		FreePointer( _entries[row] );
		if( count>0 )
		{
			_entries[ row ] = AllocPointer< MatrixEntry< T , IndexType > >( count );
			memset( _entries[ row ] , 0 , sizeof( MatrixEntry< T , IndexType > )*count );
		}
		rowSizes[row] = count;
	}
	else fprintf( stderr , "[ERROR] SparseMatrix::setRowSize: Row is out of bounds: 0 <= %d < %d\n" , (int)row , (int)rowNum ) , exit( 0 );
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::resetRowSize( size_t row , size_t count )
{
	if( row>=0 && row<rowNum )
	{
		size_t oldCount = rowSizes[row];
		_entries[row] = ReAllocPointer< MatrixEntry< T, IndexType > >( _entries[row] , count );
		if( count>oldCount ) memset( _entries[row]+oldCount , 0 , sizeof( MatrixEntry< T , IndexType > ) * ( count - oldCount ) );
		rowSizes[row] = count;
	}
	else fprintf( stderr , "[ERROR] SparseMatrix::ResetRowSize: Row is out of bounds: 0 <= %d < %d\n" , (int)row , (int)rowNum ) , exit( 0 );
}

template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::Identity( size_t dim )
{
	SparseMatrix I;
	I.resize( dim );
	for( int i=0 ; i<dim ; i++ ) I.setRowSize( i , 1 ) , I[i][0] = MatrixEntry< T , IndexType >( (IndexType)i , (T)1 );
	return I;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator *= ( T s )
{
#pragma omp parallel for
	for( int i=0 ; i<rowNum ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j].Value *= s;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator /= ( T s ){ return (*this) * ( (T)1./s ); }
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator *= ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) * B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator += ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) + B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator -= ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) - B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator * ( T s ) const
{
	SparseMatrix out = (*this);
	return out *= s;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator / ( T s ) const { return (*this) * ( (T)1. / s ); }
template< class T , class IndexType >
Pointer( T ) SparseMatrix< T , IndexType >::operator * ( const Pointer( T ) in ) const
{
	Pointer( T ) out = AllocPointer< T >( rowNum );
	MultiplyParallel( in , out , omp_get_num_procs() , 0 );
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator * ( const SparseMatrix< T , IndexType >& B ) const
{
	SparseMatrix out;
	const SparseMatrix& A = *this;
	size_t aCols = 0 , aRows = A.rowNum;
	size_t bCols = 0 , bRows = B.rowNum;
	for( int i=0 ; i<A.rowNum ; i++ ) for( int j=0 ; j<A.rowSizes[i] ; j++ ) if( aCols<=A[i][j].N ) aCols = A[i][j].N+1;
	for( int i=0 ; i<B.rowNum ; i++ ) for( int j=0 ; j<B.rowSizes[i] ; j++ ) if( bCols<=B[i][j].N ) bCols = B[i][j].N+1;
	if( bRows<aCols ) fprintf( stderr , "[Error] SparseMatrix::operator *: Matrix sizes do not support multiplication %lld x %lld * %lld x %lld\n" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols ) , exit( 0 );

	out.resize( (int)aRows );
#pragma omp parallel for
	for( int i=0 ; i<aRows ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		for( int j=0 ; j<A.rowSizes[i] ; j++ )
		{
			IndexType idx1 = A[i][j].N;
			T AValue = A[i][j].Value;
			for( int k=0 ; k<B.rowSizes[idx1] ; k++ )
			{
				IndexType idx2 = B[idx1][k].N;
				T BValue = B[idx1][k].Value;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = AValue * BValue;
				else iter->second += AValue * BValue;
			}
		}
		out.setRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator + ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	size_t rowNum = std::max< size_t >( A.rowNum , B.rowNum );
	SparseMatrix out;

	out.resize( rowNum );
#pragma omp parallel for
	for( int i=0 ; i<rowNum ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		if( i<A.rowNum )
			for( int j=0 ; j<A.rowSizes[i] ; j++ )
			{
				IndexType idx = A[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = A[i][j].Value;
				else iter->second += A[i][j].Value;
			}
		if( i<B.rowNum )
			for( int j=0 ; j<B.rowSizes[i] ; j++ )
			{
				IndexType idx = B[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = B[i][j].Value;
				else iter->second += B[i][j].Value;
			}
		out.setRowSize( i , row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator - ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	size_t rowNum = std::max< size_t >( A.rowNum , B.rowNum );
	SparseMatrix out;

	out.resize( rowNum );
#pragma omp parallel for
	for( int i=0 ; i<rowNum ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		if( i<A.rowNum )
			for( int j=0 ; j<A.rowSizes[i] ; j++ )
			{
				IndexType idx = A[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = A[i][j].Value;
				else iter->second += A[i][j].Value;
			}
		if( i<B.rowNum )
			for( int j=0 ; j<B.rowSizes[i] ; j++ )
			{
				IndexType idx = B[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = -B[i][j].Value;
				else iter->second -= B[i][j].Value;
			}
		out.setRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::transpose( T (*TransposeFunction)( const T& ) ) const
{
	SparseMatrix A;
	const SparseMatrix& At = *this;
	size_t aRows = 0 , aCols = At.rowNum;
	for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
#pragma omp parallel for
	for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
#pragma omp atomic
		A.rowSizes[ At[i][j].N ]++;
#pragma omp parallel for
	for( int i=0 ; i<A.rowNum ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.setRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction ) for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
	{
		int ii = At[i][j].N;
		A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
	}
	else for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
	{
		int ii = At[i][j].N;
		A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
	}
	return A;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::transpose( size_t aRows , T (*TransposeFunction)( const T& ) ) const
{
	SparseMatrix A;
	const SparseMatrix& At = *this;
	size_t _aRows = 0 , aCols = At.rowNum;
	for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aCols<=At[i][j].N ) _aRows = At[i][j].N+1;
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] SparseMatrix::transpose: prescribed output dimension too low: %d < %zu\n" , (int)aRows , _aRows );
		return false;
	}

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
#pragma omp parallel for
	for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
#pragma omp atomic
		A.rowSizes[ At[i][j].N ]++;
#pragma omp parallel for
	for( int i=0 ; i<A.rowNum ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.setRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
		}
	else
		for( int i=0 ; i<At.rowNum ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
		}
	return A;
}

template< class T , class IndexType >
template< class A_const_iterator , class B_const_iterator >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::Multiply( const SparseMatrixInterface< T , A_const_iterator >& A , const SparseMatrixInterface< T , B_const_iterator >& B )
{
	SparseMatrix M;
	size_t aCols = 0 , aRows = A.rows();
	size_t bCols = 0 , bRows = B.rows();
	for( int i=0 ; i<A.rows() ; i++ ) for( A_const_iterator iter=A.begin(i) ; iter!=A.end(i) ; iter++ ) if( aCols<=iter->N ) aCols = iter->N+1;
	for( int i=0 ; i<B.rows() ; i++ ) for( B_const_iterator iter=B.begin(i) ; iter!=B.end(i) ; iter++ ) if( bCols<=iter->N ) bCols = iter->N+1;
	if( bRows<aCols )
	{
		fprintf( stderr , "[ERROR] Multiply: Matrix sizes do not support multiplication %lld x %lld * %lld x %lld\n" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols );
		exit( 0 );
	}

	M.resize( (int)aRows );
#pragma omp parallel for
	for( int i=0 ; i<aRows ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		for( A_const_iterator iterA=A.begin(i) ; iterA!=A.end(i) ; iterA++ )
		{
			IndexType idx1 = iterA->N;
			T AValue = iterA->Value;
			for( B_const_iterator iterB=B.begin(idx1) ; iterB!=B.end(idx1) ; iterB++ )
			{
				IndexType idx2 = iterB->N;
				T BValue = iterB->Value;
				T temp = BValue * AValue; // temp = A( i , idx1 ) * B( idx1 , idx2 )
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = temp;
				else iter->second += temp;
			}
		}
		M.setRowSize( i , (int)row.size() );
		M.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ )
			M[i][ M.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return M;
}
template< class T , class IndexType >
template< class const_iterator >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::Transpose( const SparseMatrixInterface< T , const_iterator >& At , T (*TransposeFunction)( const T& ) )
{
	SparseMatrix< T , IndexType > A;
	size_t aRows = 0 , aCols = At.rows();
	for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aRows<=iter->N ) aRows = iter->N+1;

	A.resize( aRows );
	for( size_t i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( size_t i=0 ; i<A.rows ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.setRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			size_t ii = (size_t)iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( (IndexType)i , TransposeFunction( iter->Value ) );
		}
	else
		for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			size_t ii = (size_t)iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( (IndexType)i , iter->Value );
		}
	return A;
}
template< class T , class IndexType >
template< class const_iterator >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::Transpose( const SparseMatrixInterface< T , const_iterator >& At , size_t outRows , T (*TransposeFunction)( const T& ) )
{
	SparseMatrix< T , IndexType > A;
	size_t _aRows = 0 , aCols = At.rows() , aRows = outRows;
	for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aCols<=iter->N ) _aRows = iter->N+1;
	if( _aRows>aRows ) fprintf( stderr , "[ERROR] Transpose: prescribed output dimension too low: %d < %zu\n" , aRows , _aRows ) , exit( 0 );

	A.resize( aRows );
	for( size_t i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( size_t i=0 ; i<A.rows ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.setRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			size_t ii = (size_t)iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( (IndexType)i , TransposeFunction( iter->Value ) );
		}
	else
		for( size_t i=0 ; i<At.rows() ; i++ ) for( const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			size_t ii = (size_t)iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( (IndexType)i , iter->Value );
		}
	return true;
}
