
template< class T , class const_iterator > size_t SparseMatrixInterface< T , const_iterator >::entries( void ) const
{
	size_t entries = 0;
	for( size_t i=0 ; i<rows() ; i++ ) entries += rowSize( i );
	return entries;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) n += iter->Value * iter->Value;
	}
	return n;

}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter1 = begin( i ) ; iter1!=e ; iter1++ )
		{
			int j = iter1->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			n += (iter1->Value-value) * (iter1->Value-value);
		}
	}
	return n;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::squareASymmetricNorm( int& idx1 , int& idx2 ) const
{
	double n=0;
	double max=0;
	for( size_t i=0 ; i<rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ )
		{
			int j = iter->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			double temp = (iter->Value-value) * (iter->Value-value);
			n += temp;
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for
	for( int i=0 ; i<rows() ; i++ )
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += (T2)( _in[ iter->N ] * iter->Value );
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::multiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for
	for( int i=0 ; i<rows() ; i++ )
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		temp *= scale;
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}

template< class T , class const_iterator >
void SparseMatrixInterface< T , const_iterator >::setDiagonal( Pointer( T ) diagonal ) const
{
#pragma omp parallel for
	for( int i=0 ; i<rows() ; i++ )
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
	}
}

template< class T , class const_iterator >
void SparseMatrixInterface< T , const_iterator >::setDiagonalR( Pointer( T ) diagonal ) const
{
#pragma omp parallel for
	for( int i=0 ; i<rows() ; i++ )
	{
		diagonal[i] = (T)0;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
		if( diagonal[i] ) diagonal[i] = (T)( 1./diagonal[i] );
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::jacobiIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , ConstPointer( T2 ) in , Pointer( T2 ) out , bool dReciprocal ) const
{
	multiply( in , out );
	if( dReciprocal )
#pragma omp parallel for
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) * diagonal[i];
	else
#pragma omp parallel for
		for( int i=0 ; i<rows() ; i++ ) out[i] = in[i] + ( b[i] - out[i] ) / diagonal[i];
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
{
	if( dReciprocal )
	{
#define ITERATE( j )                                                                                \
	{                                                                                               \
		T2 _b = b[j];                                                                               \
		const_iterator e = end( j );                                                                \
		for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
		x[j] += _b * diagonal[j];                                                                   \
	}
		if( forward ) for( int j=0 ; j<int( rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( j )                                                                                \
	{                                                                                               \
		T2 _b = b[j];                                                                               \
		const_iterator e = end( j );                                                                \
		for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
		x[j] += _b / diagonal[j];                                                                   \
	}

		if( forward ) for( int j=0 ; j<int( rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#undef ITERATE
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< int >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool dReciprocal ) const
{
	if( dReciprocal )
#pragma omp parallel for
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
		{
			int jj = multiColorIndices[j];
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b * diagonal[jj];
		}
	else
#pragma omp parallel for
		for( int j=0 ; j<(int)multiColorIndices.size() ; j++ )
		{
			int jj = multiColorIndices[j];
			T2 _b = b[jj];
			const_iterator e = end( jj );
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;
			x[jj] += _b / diagonal[jj];
		}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::gsIteration( const std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward , bool dReciprocal ) const
{
#ifdef _WIN32
#define SetOMPParallel __pragma( omp parallel for )
#else // !_WIN32
#define SetOMPParallel _Pragma( "omp parallel for" )
#endif // _WIN32

	if( dReciprocal )
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
		{                                                                                                \
			int jj = indices[k];                                                                         \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b * diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
		{                                                                                                \
			int jj = indices[k];                                                                         \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b / diagonal[jj];                                                                  \
		}                                                                                                \
	}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
#undef SetOMPParallel
}
template< class SPDFunctor , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot )
{
	eps *= eps;
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );

	Real delta_new = 0 , delta_0;
	M( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] );

	delta_0 = delta_new;
	if( delta_new<=eps )
	{
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		M( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#pragma omp parallel for reduction( + : dDotQ )
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#pragma omp parallel for
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
			M( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
		}
		else
#pragma omp parallel for reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );

		Real beta = delta_new / delta_old;
#pragma omp parallel for
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
	}
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	return ii;
}
template< class SPDFunctor , class Preconditioner , class T , typename Real , class TDotTFunctor > int SolveCG( const SPDFunctor& M , const Preconditioner& P , int dim , ConstPointer( T ) b , int iters , Pointer( T ) x , double eps , TDotTFunctor Dot  )
{
	eps *= eps;
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );
	Pointer( T ) Pb = AllocPointer< T >( dim );
	Pointer( T ) temp = AllocPointer< T >( dim );

	auto PM = [&] ( ConstPointer(T) x , Pointer(T) y )
	{
		M( x , temp );
		P( ( ConstPointer(T) )temp , y );
	};

	Real delta_new = 0 , delta_0;
	P( b , Pb );
	PM( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] );

	delta_0 = delta_new;
	if( delta_new<=eps )
	{
		FreePointer( Pb );
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		FreePointer( temp );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		PM( ( ConstPointer( T ) )d , q );
		Real dDotQ = 0;
#pragma omp parallel for reduction( + : dDotQ )
		for( int i=0 ; i<dim ; i++ ) dDotQ += Dot( d[i] , q[i] );
		if( !dDotQ ) break;

		Real alpha = delta_new / dDotQ;
		Real delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
#pragma omp parallel for
			for( int i=0 ; i<dim ; i++ ) x[i] += (T)( d[i] * alpha );
			PM( ( ConstPointer( T ) )x , r );
#pragma omp parallel for reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] = Pb[i] - r[i] , delta_new += Dot( r[i] , r[i] ) , x[i] += (T)( d[i] * alpha );
		}
		else
#pragma omp parallel for reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] -=(T)( q[i] * alpha ) , delta_new += Dot( r[i] , r[i] ) ,  x[i] += (T)( d[i] * alpha );

		Real beta = delta_new / delta_old;
#pragma omp parallel for
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + (T)( d[i] * beta );
	}
	FreePointer( Pb );
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	FreePointer( temp );
	return ii;
}
