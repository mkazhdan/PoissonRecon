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

namespace VertexFactory
{
	int ToPlyType( TypeOnDisk typeOnDisk )
	{
		switch( typeOnDisk )
		{
			case TypeOnDisk::CHAR:   return PLY::Type<          char >();
			case TypeOnDisk::UCHAR:  return PLY::Type< unsigned char >();
			case TypeOnDisk::INT:    return PLY::Type<           int >();
			case TypeOnDisk::UINT:   return PLY::Type< unsigned  int >();
			case TypeOnDisk::FLOAT:  return PLY::Type<         float >();
			case TypeOnDisk::DOUBLE: return PLY::Type<        double >();
			default: ERROR_OUT( "Unrecognized type: " , typeOnDisk );
		}
		return -1;
	}

	TypeOnDisk FromPlyType( int plyType )
	{
		switch( plyType )
		{
			case PLY_INT:    return TypeOnDisk::INT;
			case PLY_UINT:   return TypeOnDisk::UINT;
			case PLY_CHAR:   return TypeOnDisk::CHAR;
			case PLY_UCHAR:  return TypeOnDisk::UCHAR;
			case PLY_FLOAT:  return TypeOnDisk::FLOAT;
			case PLY_DOUBLE: return TypeOnDisk::DOUBLE;
			default: ERROR_OUT( "Unrecognized type: " , plyType );
		}
		return TypeOnDisk::UNKNOWN;
	}

	template<> TypeOnDisk GetTypeOnDisk<          char >( void ){ return TypeOnDisk::CHAR;  }
	template<> TypeOnDisk GetTypeOnDisk< unsigned char >( void ){ return TypeOnDisk::UCHAR; }
	template<> TypeOnDisk GetTypeOnDisk<           int >( void ){ return TypeOnDisk::INT;  }
	template<> TypeOnDisk GetTypeOnDisk< unsigned  int >( void ){ return TypeOnDisk::UINT; }
	template<> TypeOnDisk GetTypeOnDisk<         float >( void ){ return TypeOnDisk::FLOAT;  }
	template<> TypeOnDisk GetTypeOnDisk<        double >( void ){ return TypeOnDisk::DOUBLE; }

	template< typename Real >
	template< typename Type >
	bool VertexIO< Real >::_ReadBinary( FILE *fp , Real &r )
	{
		Type t;
		if( fread( &t , sizeof(Type) , 1 , fp )!=1 ) return false;
		r = (Real)t;
		return true;
	}

	template< typename Real >
	template< typename Type >
	void VertexIO< Real >::_WriteBinary( FILE *fp , Real r ){ Type t = (Type)r ; fwrite(  &t , sizeof(Type) , 1 , fp ); }

	template< typename Real >
	bool VertexIO< Real >::ReadASCII( FILE *fp , TypeOnDisk , Real &s )
	{
		double d;
		if( fscanf( fp , " %lf"  , &d )!=1 ) return false;
		s = (Real)d;
		return true;
	}

	template< typename Real >
	bool VertexIO< Real >::ReadBinary( FILE *fp , TypeOnDisk typeOnDisk , Real &s )
	{
		if( TypeOnDisk()==typeOnDisk ) return fread( &s , sizeof(Real) , 1 , fp )==1;
		switch( typeOnDisk )
		{
			case TypeOnDisk::CHAR:   return _ReadBinary<          char >( fp , s );
			case TypeOnDisk::UCHAR:  return _ReadBinary< unsigned char >( fp , s );
			case TypeOnDisk::INT:    return _ReadBinary<           int >( fp , s );
			case TypeOnDisk::UINT:   return _ReadBinary< unsigned  int >( fp , s );
			case TypeOnDisk::FLOAT:  return _ReadBinary<         float >( fp , s );
			case TypeOnDisk::DOUBLE: return _ReadBinary<        double >( fp , s );
			default: ERROR_OUT( "Unrecognized type: " , typeOnDisk );
		}
		return true;
	}
	template< typename Real >
	void VertexIO< Real >::WriteASCII( FILE *fp , TypeOnDisk typeOnDisk , const Real &s )
	{
		switch( typeOnDisk )
		{
			case TypeOnDisk::CHAR:   fprintf( fp , " %d" , (         char)s ) ; break;
			case TypeOnDisk::UCHAR:  fprintf( fp , " %u" , (unsigned char)s ) ; break;
			case TypeOnDisk::INT:    fprintf( fp , " %d" , (          int)s ) ; break;
			case TypeOnDisk::UINT:   fprintf( fp , " %u" , (unsigned  int)s ) ; break;
			case TypeOnDisk::FLOAT:  fprintf( fp , " %f" , (        float)s ) ; break;
			case TypeOnDisk::DOUBLE: fprintf( fp , " %f" , (       double)s ) ; break;
			default: ERROR_OUT( "Unrecongized type: " , typeOnDisk );
		}
	}

	template< typename Real >
	void VertexIO< Real >::WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , const Real &s )
	{
		if( TypeOnDisk()==typeOnDisk ) fwrite( &s , sizeof(Real) , 1 , fp );
		switch( typeOnDisk )
		{
			case TypeOnDisk::CHAR:   _WriteBinary<          char >( fp , s ) ; break;
			case TypeOnDisk::UCHAR:  _WriteBinary< unsigned char >( fp , s ) ; break;
			case TypeOnDisk::INT:    _WriteBinary<           int >( fp , s ) ; break;
			case TypeOnDisk::UINT:   _WriteBinary< unsigned  int >( fp , s ) ; break;
			case TypeOnDisk::FLOAT:  _WriteBinary<         float >( fp , s ) ; break;
			case TypeOnDisk::DOUBLE: _WriteBinary<        double >( fp , s ) ; break;
			default: ERROR_OUT( "Unrecongized type: " , typeOnDisk );
		}
	}

	template< typename Real >
	bool VertexIO< Real >::ReadASCII( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , Real *s )
	{
		for( size_t i=0 ; i<sz ; i++ ) if( !ReadASCII( fp , typeOnDisk , s[i] ) ) return false;
		return true;
	}

	template< typename Real >
	bool VertexIO< Real >::ReadBinary( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , Real *s )
	{
		if( TypeOnDisk()==typeOnDisk ) return fread( &s , sizeof(Real) , sz , fp )==sz;
		else for( size_t i=0 ; i<sz ; i++ ) if( !ReadBinary( fp , typeOnDisk , s[i] ) ) return false;
		return true;
	}
	template< typename Real >
	void VertexIO< Real >::WriteASCII( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s )
	{
		for( size_t i=0 ; i<sz ; i++ ) WriteASCII( fp , typeOnDisk , s[i] );
	}

	template< typename Real >
	void VertexIO< Real >::WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s )
	{
		if( TypeOnDisk()==typeOnDisk ) fwrite( &s , sizeof(Real) , sz , fp );
		else for( size_t i=0 ; i<sz ; i++ ) WriteBinary( fp , typeOnDisk , s[i] );
	}

	/////////////////////
	// PositionFactory //
	/////////////////////
	template< typename Real , unsigned int Dim >
	PlyProperty PositionFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty PositionFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real , unsigned int Dim >
	PlyProperty PositionFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty PositionFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}

	template<> const std::string PositionFactory<  float , 1 >::_PlyNames[] = { "x" };
	template<> const std::string PositionFactory< double , 1 >::_PlyNames[] = { "x" };
	template<> const std::string PositionFactory<  float , 2 >::_PlyNames[] = { "x" , "y" };
	template<> const std::string PositionFactory< double , 2 >::_PlyNames[] = { "x" , "y" };
	template<> const std::string PositionFactory<  float , 3 >::_PlyNames[] = { "x" , "y" , "z" };
	template<> const std::string PositionFactory< double , 3 >::_PlyNames[] = { "x" , "y" , "z" };

	///////////////////
	// NormalFactory //
	///////////////////
	template< typename Real , unsigned int Dim >
	PlyProperty NormalFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty NormalFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real , unsigned int Dim >
	PlyProperty NormalFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty NormalFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}

	template<> const std::string NormalFactory<  float , 1 >::_PlyNames[] = { "nx" };
	template<> const std::string NormalFactory< double , 1 >::_PlyNames[] = { "nx" };
	template<> const std::string NormalFactory<  float , 2 >::_PlyNames[] = { "nx" , "ny" };
	template<> const std::string NormalFactory< double , 2 >::_PlyNames[] = { "nx" , "ny" };
	template<> const std::string NormalFactory<  float , 3 >::_PlyNames[] = { "nx" , "ny" , "nz" };
	template<> const std::string NormalFactory< double , 3 >::_PlyNames[] = { "nx" , "ny" , "nz" };

	////////////////////
	// TextureFactory //
	////////////////////
	template< typename Real , unsigned int Dim >
	PlyProperty TextureFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty TextureFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real , unsigned int Dim >
	PlyProperty TextureFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}
	template< typename Real , unsigned int Dim >
	PlyProperty TextureFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}

	template<> const std::string TextureFactory<  float , 1 >::_PlyNames[] = { "u" };
	template<> const std::string TextureFactory< double , 1 >::_PlyNames[] = { "u" };
	template<> const std::string TextureFactory<  float , 2 >::_PlyNames[] = { "u" , "v" };
	template<> const std::string TextureFactory< double , 2 >::_PlyNames[] = { "u" , "v" };
	template<> const std::string TextureFactory<  float , 3 >::_PlyNames[] = { "u" , "v" , "w" };
	template<> const std::string TextureFactory< double , 3 >::_PlyNames[] = { "u" , "v" , "w" };

	/////////////////////
	// RGBColorFactory //
	/////////////////////
	template< typename Real >
	PlyProperty RGBColorFactory< Real >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real >
	PlyProperty RGBColorFactory< Real >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real >
	PlyProperty RGBColorFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx%3] ) );
	}
	template< typename Real >
	PlyProperty RGBColorFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}

	template< typename Real > const std::string RGBColorFactory< Real >::_PlyNames[] = { "red" , "green" , "blue" , "r" , "g" , "b" };

	//////////////////////
	// RGBAColorFactory //
	//////////////////////
	template< typename Real >
	PlyProperty RGBAColorFactory< Real >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real >
	PlyProperty RGBAColorFactory< Real >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real >
	PlyProperty RGBAColorFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx%3] ) );
	}
	template< typename Real >
	PlyProperty RGBAColorFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( VertexType , coords[idx] ) );
	}

	template< typename Real > const std::string RGBAColorFactory< Real >::_PlyNames[] = { "red" , "green" , "blue" , "alpha" , "r" , "g" , "b" , "a" };

	//////////////////
	// ValueFactory //
	//////////////////
	template< typename Real >
	PlyProperty ValueFactory< Real >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>= plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real >
	PlyProperty ValueFactory< Real >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real >
	PlyProperty ValueFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , 0 );
	}
	template< typename Real >
	PlyProperty ValueFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , 0 );
	}

	template< typename Real > const std::string ValueFactory< Real >::_PlyNames[] = { "value" };

	////////////////////
	// DynamicFactory //
	////////////////////
	template< typename Real >
	DynamicFactory< Real >::DynamicFactory( const std::vector< std::pair< std::string , TypeOnDisk > > &namesAndTypesOnDisk ) : _namesAndTypesOnDisk(namesAndTypesOnDisk)
	{
	}
	template< typename Real >
	DynamicFactory< Real >::DynamicFactory( const std::vector< PlyProperty > &plyProperties )
	{
		_namesAndTypesOnDisk.resize( plyProperties.size() );
		for( int i=0 ; i<plyProperties.size() ; i++ ) _namesAndTypesOnDisk[i] = std::pair< std::string , TypeOnDisk >( plyProperties[i].name , FromPlyType( plyProperties[i].external_type ) );
	}

	template< typename Real >
	bool DynamicFactory< Real >::readASCII( FILE *fp , VertexType &dt ) const
	{
		for( unsigned int i=0 ; i<dt.dim() ; i++ ) if( !VertexIO< Real >::ReadASCII( fp , _namesAndTypesOnDisk[i].second , dt[i] ) ) return false;
		return true;
	}
	template< typename Real >
	bool DynamicFactory< Real >::readBinary( FILE *fp , VertexType &dt ) const
	{
		for( unsigned int i=0 ; i<dt.dim() ; i++ ) if( !VertexIO< Real >::ReadBinary( fp , _namesAndTypesOnDisk[i].second , dt[i] ) ) return false;
		return true;
	}
	template< typename Real >
	void DynamicFactory< Real >::writeASCII( FILE *fp , const VertexType &dt ) const
	{
		for( unsigned int i=0 ; i<dt.dim() ; i++ ) VertexIO< Real >::WriteASCII( fp , _namesAndTypesOnDisk[i].second , dt[i] );
	}
	template< typename Real >
	void DynamicFactory< Real >::writeBinary( FILE *fp , const VertexType &dt ) const
	{
		for( unsigned int i=0 ; i<dt.dim() ; i++ ) VertexIO< Real >::WriteBinary( fp , _namesAndTypesOnDisk[i].second , dt[i] );
	}

	template< typename Real >
	PlyProperty DynamicFactory< Real >::plyReadProperty( unsigned int idx ) const
	{
		if( idx>=plyReadNum() ) ERROR_OUT( "read property out of bounds" );
		return PlyProperty( _namesAndTypesOnDisk[idx].first , ToPlyType( _namesAndTypesOnDisk[idx].second ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}
	template< typename Real >
	PlyProperty DynamicFactory< Real >::plyWriteProperty( unsigned int idx ) const
	{
		if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" );
		return PlyProperty( _namesAndTypesOnDisk[idx].first , ToPlyType( _namesAndTypesOnDisk[idx].second ) , PLY::Type< Real >() , sizeof(Real)*idx );
	}

	template< typename Real >
	bool DynamicFactory< Real >::operator == ( const DynamicFactory< Real > &factory ) const
	{
		if( size()!=factory.size() ) return false;
		for( int i=0 ; i<size() ; i++ ) if( _namesAndTypesOnDisk[i].first!=factory._namesAndTypesOnDisk[i].first || _namesAndTypesOnDisk[i].second!=factory._namesAndTypesOnDisk[i].second ) return false;
		return true;
	}

	/////////////
	// Factory //
	/////////////
	template< typename Real , typename ... Factories >
	template< unsigned int I >
	typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type Factory< Real , Factories ... >::_plyReadProperty( unsigned int idx , size_t offset ) const
	{
		if( idx<this->template get<I>().plyReadNum() )
		{
			PlyProperty prop = this->template get<I>().plyReadProperty(idx);
			prop.offset += (int)offset;
			return prop;
		}
		else return _plyReadProperty<I+1>( idx - this->template get<I>().plyReadNum() , offset + this->template get<I>().bufferSize() );
	}

	template< typename Real , typename ... Factories >
	template< unsigned int I >
	typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type Factory< Real , Factories ... >::_plyWriteProperty( unsigned int idx , size_t offset ) const
	{
		if( idx<this->template get<I>().plyWriteNum() )
		{
			PlyProperty prop = this->template get<I>().plyWriteProperty(idx);
			prop.offset += (int)offset;
			return prop;
		}
		else return _plyWriteProperty<I+1>( idx - this->template get<I>().plyWriteNum() , offset + this->template get<I>().bufferSize() );
	}

	template< typename Real , typename ... Factories >
	template< unsigned int I >
	typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type Factory< Real , Factories ... >::_plyStaticReadProperty( unsigned int idx ) const
	{
		if( idx<this->template get<I>().plyReadNum() )
		{
			VertexType v;
			PlyProperty prop = this->template get<I>().plyStaticReadProperty( idx );
			prop.offset += (int)( (size_t)&v.template get<I>() - (size_t)&v );
			return prop;
		}
		else return _plyStaticReadProperty<I+1>( idx - this->template get<I>().plyReadNum() );
	}

	template< typename Real , typename ... Factories >
	template< unsigned int I >
	typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type Factory< Real , Factories ... >::_plyStaticWriteProperty( unsigned int idx ) const
	{
		if( idx<this->template get<I>().plyWriteNum() )
		{
			VertexType v;
			PlyProperty prop = this->template get<I>().plyStaticWriteProperty( idx );
			prop.offset += (int)( (size_t)&v.template get<I>() - (size_t)&v );
			return prop;
		}
		else return _plyStaticWriteProperty<I+1>( idx - this->template get<I>().plyWriteNum() ) ;
	}
}
