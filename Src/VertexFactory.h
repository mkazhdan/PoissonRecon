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
#ifndef VERTEX_FACTORY_INCLUDED
#define VERTEX_FACTORY_INCLUDED

#include "Ply.h"

namespace VertexFactory
{
	// Mimicking the scalar types in Ply.inl
	// Assuming Real is something that is castable to/from a double
	enum TypeOnDisk
	{
		CHAR ,
		UCHAR ,
		INT ,
		UINT ,
		FLOAT ,
		DOUBLE ,
		UNKNOWN
	};

	int ToPlyType( TypeOnDisk typeOnDisk );
	TypeOnDisk FromPlyType( int plyType );
	template< typename Real > TypeOnDisk GetTypeOnDisk( void );

	template< typename Real >
	struct VertexIO
	{
		static bool ReadASCII  ( FILE *fp , TypeOnDisk typeOnDisk ,       Real &s );
		static bool ReadBinary ( FILE *fp , TypeOnDisk typeOnDisk ,       Real &s );
		static void WriteASCII ( FILE *fp , TypeOnDisk typeOnDisk , const Real &s );
		static void WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , const Real &s );

		static bool ReadASCII  ( FILE *fp , TypeOnDisk typeOnDisk , size_t sz ,       Real *s );
		static bool ReadBinary ( FILE *fp , TypeOnDisk typeOnDisk , size_t sz ,       Real *s );
		static void WriteASCII ( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s );
		static void WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s );

	protected:
		template< typename Type > static bool  _ReadBinary( FILE *fp , Real &s );
		template< typename Type > static void _WriteBinary( FILE *fp , Real  s );
	};

	template< typename _VertexType , typename FactoryType >
	struct _Factory
	{
		typedef _VertexType VertexType;
		virtual VertexType operator()( void ) const = 0;

		// Reading/writing methods for PLY format
		virtual unsigned int  plyReadNum( void ) const = 0;
		virtual unsigned int plyWriteNum( void ) const = 0;
		virtual bool plyValidReadProperties( const bool *flags ) const = 0;

		virtual PlyProperty  plyReadProperty( unsigned int idx ) const = 0;
		virtual PlyProperty plyWriteProperty( unsigned int idx ) const = 0;

		// Reading/writing methods for ASCII/Binary files
		virtual bool   readASCII( FILE *fp ,       VertexType &dt ) const = 0;
		virtual bool  readBinary( FILE *fp ,       VertexType &dt ) const = 0;
		virtual void  writeASCII( FILE *fp , const VertexType &dt ) const = 0;
		virtual void writeBinary( FILE *fp , const VertexType &dt ) const = 0;

		virtual size_t bufferSize( void ) const = 0;
		virtual void toBuffer( const VertexType &dt , char *buffer ) const = 0;
		virtual void fromBuffer( const char *buffer , VertexType &dt ) const = 0;

		virtual bool isStaticallyAllocated( void ) const = 0;
		virtual PlyProperty  plyStaticReadProperty( unsigned int idx ) const = 0;
		virtual PlyProperty plyStaticWriteProperty( unsigned int idx ) const = 0;

		virtual bool operator == ( const FactoryType &factory ) const = 0;
		bool operator != ( const FactoryType &factory ) const { return !( (*this)==factory ); }
	};

	// An empty factory
	template< typename Real >
	struct EmptyFactory : _Factory< EmptyVectorType< Real > , EmptyFactory< Real > >
	{
		typedef typename _Factory< EmptyVectorType< Real > , EmptyFactory< Real > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename X > Transform( X ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return 0; }
		unsigned int plyWriteNum( void ) const { return 0; }
		bool plyValidReadProperties( const bool *flags ) const { return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const { if( idx>= plyReadNum() ) ERROR_OUT(  "read property out of bounds" ) ; return PlyProperty(); }
		PlyProperty plyWriteProperty( unsigned int idx ) const { if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" ) ; return PlyProperty(); }
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return true; }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return true; }
		void  writeASCII( FILE *fp , const VertexType &dt ) const {}
		void writeBinary( FILE *fp , const VertexType &dt ) const {};

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const { if( idx>= plyReadNum() ) ERROR_OUT(  "read property out of bounds" ) ; return PlyProperty(); }
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const { if( idx>=plyWriteNum() ) ERROR_OUT( "write property out of bounds" ) ; return PlyProperty(); }

		size_t bufferSize( void ) const { return 0; }
		void toBuffer( const VertexType &dt , char *buffer ) const {}
		void fromBuffer( const char* buffer , VertexType &dt ) const {}

		bool operator == ( const EmptyFactory &factory ) const { return true; }

	};

	// The position factory
	template< typename Real , unsigned int Dim >
	struct PositionFactory : _Factory< Point< Real , Dim > , PositionFactory< Real , Dim > >
	{
		typedef typename _Factory< Point< Real , Dim > , PositionFactory< Real , Dim > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){ _xForm = XForm< Real , Dim+1 >::Identity(); }
			Transform( XForm< Real , Dim > xForm ){ _xForm = XForm< Real , Dim+1 >::Identity() ; for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _xForm(i,j) = xForm(i,j); }
			Transform( XForm< Real , Dim+1 > xForm ) : _xForm(xForm) {}
			VertexType operator()( VertexType dt ) const { return _xForm * dt; }
		protected:
			XForm< Real , Dim+1 > _xForm;
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return Dim; }
		unsigned int plyWriteNum( void ) const { return Dim; }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<Dim ; d++ ) if( !flags[d] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , Dim , &dt[0] ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , Dim , &dt[0] ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		PositionFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< Real >() ) : _typeOnDisk( typeOnDisk ) {}

		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const PositionFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }
	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The vertex factory
	template< typename Real , unsigned int Dim >
	struct NormalFactory : public _Factory< Point< Real , Dim > , NormalFactory< Real , Dim > >
	{
		typedef typename _Factory< Point< Real , Dim > , NormalFactory< Real , Dim > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){ _xForm = XForm< Real , Dim >::Identity(); }
			Transform( XForm< Real , Dim > xForm ){ _xForm = xForm.inverse().transpose(); }
			Transform( XForm< Real , Dim+1 > xForm )
			{
				for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _xForm(i,j) = xForm(i,j);
				_xForm = _xForm.inverse().transpose();
			}
			VertexType operator()( VertexType dt ) const { return _xForm * dt; }
		protected:
			XForm< Real , Dim > _xForm;
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return Dim; }
		unsigned int plyWriteNum( void ) const { return Dim; }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<Dim ; d++ ) if( !flags[d] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , Dim , &dt[0] ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , Dim , &dt[0] ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		NormalFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< Real >() ) : _typeOnDisk( typeOnDisk ) {}

		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const NormalFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }
	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The texture factory
	template< typename Real , unsigned int Dim >
	struct TextureFactory : public _Factory< Point< Real , Dim > , TextureFactory< Real , Dim > >
	{
		typedef typename _Factory< Point< Real , Dim > , TextureFactory< Real , Dim > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename XForm > Transform( XForm xForm ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return Dim; }
		unsigned int plyWriteNum( void ) const { return Dim; }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<Dim ; d++ ) if( !flags[d] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , Dim , &dt[0] ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , Dim , &dt[0] ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , Dim , &dt[0] ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		TextureFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< Real >() ) : _typeOnDisk( typeOnDisk ) {}
		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const TextureFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }

	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The rgb color factory
	template< typename Real >
	struct RGBColorFactory : public _Factory< Point< Real , 3 > , RGBColorFactory< Real > >
	{
		typedef typename _Factory< Point< Real , 3 > , RGBColorFactory< Real > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename XForm > Transform( XForm xForm ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return 6; }
		unsigned int plyWriteNum( void ) const { return 3; }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<3 ; d++ ) if( !flags[d] && !flags[d+3] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , 3 , &dt[0] ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , 3 , &dt[0] ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , 3 , &dt[0] ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , 3 , &dt[0] ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		RGBColorFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< unsigned char >() ) : _typeOnDisk( typeOnDisk ) {}
		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const RGBColorFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }
	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The rgba color factory
	template< typename Real >
	struct RGBAColorFactory : public _Factory< Point< Real , 4 > , RGBAColorFactory< Real > >
	{
		typedef typename _Factory< Point< Real , 4 > , RGBAColorFactory< Real > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename XForm > Transform( XForm xForm ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return VertexType(); }

		unsigned int  plyReadNum( void ) const { return 8; }
		unsigned int plyWriteNum( void ) const { return 4; }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<4 ; d++ ) if( !flags[d] && !flags[d+4] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , 4 , &dt[0] ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , 4 , &dt[0] ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , 4 , &dt[0] ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , 4 , &dt[0] ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		RGBAColorFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< unsigned char >() ) : _typeOnDisk( typeOnDisk ) {}
		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const RGBAColorFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }
	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The value factory
	template< typename Real >
	struct ValueFactory : public _Factory< Real , ValueFactory< Real > >
	{
		typedef typename _Factory< Real , ValueFactory< Real > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename XForm > Transform( XForm xForm ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return (Real)0; }

		unsigned int  plyReadNum( void ) const { return 1; }
		unsigned int plyWriteNum( void ) const { return 1; }
		bool plyValidReadProperties( const bool *flags ) const { if( !flags[0] ) return false ; return true ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >:: ReadASCII( fp , _typeOnDisk , 1 , &dt ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return VertexIO< Real >::ReadBinary( fp , _typeOnDisk , 1 , &dt ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const { VertexIO< Real >:: WriteASCII( fp , _typeOnDisk , 1 , &dt ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { VertexIO< Real >::WriteBinary( fp , _typeOnDisk , 1 , &dt ); }

		bool isStaticallyAllocated( void ) const{ return true; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const;
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const;

		ValueFactory( TypeOnDisk typeOnDisk=GetTypeOnDisk< Real >() ) : _typeOnDisk( typeOnDisk ) {}
		size_t bufferSize( void ) const { return sizeof( VertexType ); }
		void toBuffer( const VertexType &dt , char *buffer ) const { *( (VertexType*)buffer ) = dt; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { dt = *( (VertexType*)buffer ); }

		bool operator == ( const ValueFactory &factory ) const { return _typeOnDisk==factory._typeOnDisk; }
	protected:
		TypeOnDisk _typeOnDisk;
		static const std::string _PlyNames[];
	};

	// The named array factory
	template< typename Real >
	struct DynamicFactory : public _Factory< Point< Real > , DynamicFactory< Real > >
	{
		typedef typename _Factory< Point< Real > , DynamicFactory< Real > >::VertexType VertexType;

		struct Transform
		{
			Transform( void ){}
			template< typename XForm > Transform( XForm xForm ){}
			VertexType operator()( VertexType dt ) const { return dt; }
		};

		VertexType operator()( void ) const { return Point< Real >( _namesAndTypesOnDisk.size() ); }
		unsigned int  plyReadNum( void ) const { return (unsigned int )_namesAndTypesOnDisk.size(); }
		unsigned int plyWriteNum( void ) const { return (unsigned int )_namesAndTypesOnDisk.size(); }
		bool plyValidReadProperties( const bool *flags ) const { for( int d=0 ; d<_namesAndTypesOnDisk.size() ; d++ ) if( !flags[d] ) return false ; return true ; }
		PlyProperty plyReadProperty( unsigned int idx ) const;
		PlyProperty plyWriteProperty( unsigned int idx ) const;

		bool   readASCII( FILE *fp ,       VertexType &dt ) const;
		bool  readBinary( FILE *fp ,       VertexType &dt ) const;
		void  writeASCII( FILE *fp , const VertexType &dt ) const;
		void writeBinary( FILE *fp , const VertexType &dt ) const;

		bool isStaticallyAllocated( void ) const{ return false; }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const { ERROR_OUT( "does not support static allocation" ) ; return PlyProperty(); }
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const { ERROR_OUT( "does not support static allocation" ) ; return PlyProperty(); }

		size_t size( void ) const { return _namesAndTypesOnDisk.size(); }

		DynamicFactory( const std::vector< std::pair< std::string , TypeOnDisk > > &namesAndTypesOnDisk );
		DynamicFactory( const std::vector< PlyProperty > &plyProperties );
		size_t bufferSize( void ) const { return sizeof(Real) * _namesAndTypesOnDisk.size(); }
		void toBuffer( const VertexType &dt , char *buffer ) const { for( size_t i=0 ; i<_namesAndTypesOnDisk.size() ; i++ ) ( (Real*)buffer )[i] = dt[i]; }
		void fromBuffer( const char *buffer , VertexType &dt ) const { for( size_t i=0 ; i<_namesAndTypesOnDisk.size() ; i++ ) dt[i] = ( (const Real*)buffer )[i]; }

		bool operator == ( const DynamicFactory &factory ) const;
	protected:
		std::vector< std::pair< std::string , TypeOnDisk > > _namesAndTypesOnDisk;
	};

	// A factory for building the union of multiple components
	template< typename Real , typename ... Factories >
	struct Factory :  public _Factory< VectorTypeUnion< Real , typename Factories::VertexType ... > , Factory< Real , Factories ... > >
	{
	protected:
		typedef std::tuple< Factories ... > _FactoryTuple;
	public:
		typedef typename _Factory< VectorTypeUnion< Real , typename Factories::VertexType ... > , Factory< Real , Factories ... > >::VertexType VertexType;
		template< unsigned int I > using FactoryType = typename std::tuple_element< I , _FactoryTuple >::type;
		template< unsigned int I >       FactoryType< I >& get( void )       { return std::get< I >( _factoryTuple ); }
		template< unsigned int I > const FactoryType< I >& get( void ) const { return std::get< I >( _factoryTuple ); }

		struct Transform
		{
		protected:
			typedef std::tuple< typename Factories::Transform ... > _TransformTuple;
		public:
			Transform( void ){}
			template< typename XForm >
			Transform( XForm xForm ){ _set<0>( xForm ); }
			VertexType operator()( VertexType in ) const { VertexType out ; _xForm< 0 >( in , out ) ; return out; }

			template< unsigned int I > using TransformType = typename std::tuple_element< I , _TransformTuple >::type;
			template< unsigned int I >       TransformType< I >& get( void )       { return std::get< I >( _transformTuple ); }
			template< unsigned int I > const TransformType< I >& get( void ) const { return std::get< I >( _transformTuple ); }
		protected:
			_TransformTuple _transformTuple;
			template< unsigned int I , typename XForm > typename std::enable_if< I!=sizeof...(Factories) >::type _set( XForm xForm ){ this->template get<I>() = xForm ; _set< I+1 >( xForm ); }
			template< unsigned int I , typename XForm > typename std::enable_if< I==sizeof...(Factories) >::type _set( XForm xForm ){}
			template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) >::type _xForm( const VertexType &in , VertexType &out ) const { out.template get<I>() = this->template get<I>()( in.template get<I>() ) ; _xForm< I+1 >( in , out ); }
			template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) >::type _xForm( const VertexType &in , VertexType &out ) const {}
		};

		Factory( void ){}
		Factory( Factories ... factories ) : _factoryTuple( factories ... ){ }

		VertexType operator()( void ) const{ return _VertexType( std::make_index_sequence< sizeof...(Factories) >() ); }

		unsigned int  plyReadNum( void ) const { return  _plyReadNum<0>(); }
		unsigned int plyWriteNum( void ) const { return _plyWriteNum<0>(); }
		bool plyValidReadProperties( const bool *flags ) const { return _plyValidReadProperties<0>( flags ) ; }
		template< unsigned int I > bool plyValidReadProperties( const bool* flags ) const { return get< I >().plyValidReadProperties( flags + _readOffset< I >() ) ; }
		PlyProperty  plyReadProperty( unsigned int idx ) const { return  _plyReadProperty<0>( idx , 0 ); }
		PlyProperty plyWriteProperty( unsigned int idx ) const { return _plyWriteProperty<0>( idx , 0 ); }
		bool   readASCII( FILE *fp ,       VertexType &dt ) const { return  _readASCII<0>( fp , dt ); }
		bool  readBinary( FILE *fp ,       VertexType &dt ) const { return _readBinary<0>( fp , dt ); }
		void  writeASCII( FILE *fp , const VertexType &dt ) const {  _writeASCII<0>( fp , dt ); }
		void writeBinary( FILE *fp , const VertexType &dt ) const { _writeBinary<0>( fp , dt ); }

		bool isStaticallyAllocated( void ) const { return _isStaticallyAllocated<0>(); }
		PlyProperty  plyStaticReadProperty( unsigned int idx ) const { return  _plyStaticReadProperty<0>( idx ); }
		PlyProperty plyStaticWriteProperty( unsigned int idx ) const { return _plyStaticWriteProperty<0>( idx ); }

		size_t bufferSize( void ) const { return _bufferSize<0>(); }
		void toBuffer( const VertexType &dt , char *buffer ) const { _toBuffer<0>( dt , buffer ); }
		void fromBuffer( const char *buffer , VertexType &dt ) const { _fromBuffer<0>( buffer , dt ); }

		bool operator == ( const Factory &factory ) const { return _equal<0>( factory ); }
	protected:
		_FactoryTuple _factoryTuple;
		template< size_t ... Is > VertexType _VertexType( std::index_sequence< Is ... > ) const{ return VertexType( get<Is>().operator()() ... ); }

		template< typename Factory1 , typename Factory2 >
		static typename std::enable_if< !std::is_same< Factory1 , Factory1 >::value , bool >::type _EqualFactories( const Factory1 &f1 , const Factory1 &f2 ){ return false; }
		template< typename Factory1 , typename Factory2 >
		static typename std::enable_if<  std::is_same< Factory1 , Factory1 >::value , bool >::type _EqualFactories( const Factory1 &f1 , const Factory1 &f2 ){ return f1==f2; }


		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , unsigned int >::type _plyReadNum( void ) const { return get<I>().plyReadNum() + _plyReadNum< I+1 >(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , unsigned int >::type _plyReadNum( void ) const { return 0; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , unsigned int >::type _plyWriteNum( void ) const { return get<I>().plyWriteNum() + _plyWriteNum< I+1 >(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , unsigned int >::type _plyWriteNum( void ) const { return 0; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , bool >::type _plyValidReadProperties( const bool *flags ) const { return get< I >().plyValidReadProperties( flags ) && _plyValidReadProperties< I+1 >( flags + get< I >().plyReadNum() ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , bool >::type _plyValidReadProperties( const bool *flags ) const { return true; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type _plyReadProperty( unsigned int idx , size_t offset ) const;
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type _plyWriteProperty( unsigned int idx , size_t offset ) const;
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , PlyProperty >::type _plyReadProperty( unsigned int idx , size_t offset ) const { ERROR_OUT( "read property out of bounds" ) ; return PlyProperty(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , PlyProperty >::type _plyWriteProperty( unsigned int idx , size_t offset ) const { ERROR_OUT( "write property out of bounds" ) ; return PlyProperty(); }
		template< unsigned int I > typename std::enable_if< I==0 , unsigned int >::type _readOffset( void ) const { return 0; }
		template< unsigned int I > typename std::enable_if< I!=0 , unsigned int >::type _readOffset( void ) const { return _readOffset< I-1 >() + get< I-1 >().plyReadNum(); }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , bool >::type _readASCII( FILE *fp , VertexType &dt ) const { return this->template get<I>().readASCII( fp , dt.template get<I>() ) && _readASCII< I+1 >( fp , dt ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , bool >::type _readASCII( FILE *fp , VertexType &dt ) const { return true; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , bool >::type _readBinary( FILE *fp , VertexType &dt ) const { return this->template get<I>().readBinary( fp , dt.template get<I>() ) && _readBinary< I+1 >( fp , dt ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , bool >::type _readBinary( FILE *fp , VertexType &dt ) const { return true; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) >::type  _writeASCII( FILE *fp , const VertexType &dt ) const { this->template get<I>().writeASCII( fp , dt.template get<I>() ) ; _writeASCII< I+1 >( fp , dt ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) >::type  _writeASCII( FILE *fp , const VertexType &dt ) const {}
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) >::type _writeBinary( FILE *fp , const VertexType &dt ) const { this->template get<I>().writeBinary( fp , dt.template get<I>() ) ; _writeBinary< I+1 >( fp , dt ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) >::type _writeBinary( FILE *fp , const VertexType &dt ) const {}

		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , bool >::type _isStaticallyAllocated( void ) const { return this->template get< I >().isStaticallyAllocated() && _isStaticallyAllocated< I+1 >(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , bool >::type _isStaticallyAllocated( void ) const { return true; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type _plyStaticReadProperty ( unsigned int idx ) const;
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , PlyProperty >::type _plyStaticWriteProperty( unsigned int idx ) const;
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , PlyProperty >::type _plyStaticReadProperty ( unsigned int idx ) const { ERROR_OUT(  "read property out of bounds" ) ; return PlyProperty(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , PlyProperty >::type _plyStaticWriteProperty( unsigned int idx ) const { ERROR_OUT( "write property out of bounds" ) ; return PlyProperty(); }

		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , size_t >::type _bufferSize( void ) const { return this->template get<I>().bufferSize() + _bufferSize< I+1 >(); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , size_t >::type _bufferSize( void ) const { return 0; }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) >::type _toBuffer( const VertexType &dt , char *buffer ) const { this->template get<I>().toBuffer( dt.template get<I>() , buffer ) ; _toBuffer< I+1 >( dt , buffer + this->template get<I>().bufferSize() ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) >::type _toBuffer( const VertexType &dt , char *buffer ) const {}
		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) >::type _fromBuffer( const char * buffer , VertexType &dt ) const {this->template  get<I>().fromBuffer( buffer , dt.template get<I>() ) ; _fromBuffer< I+1 >( buffer + this->template get<I>().bufferSize() , dt ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) >::type _fromBuffer( const char * buffer , VertexType &dt ) const {}

		template< unsigned int I > typename std::enable_if< I!=sizeof...(Factories) , bool >::type _equal( const Factory &factory ) const { return _EqualFactories< FactoryType<I> , Factory >( this->template get< I >() , factory.template get<I>() ) && _equal< I+1 >( factory ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(Factories) , bool >::type _equal( const Factory &factory ) const { return true; }
	};
}

#include "VertexFactory.inl"

#endif // VERTEX_FACTORY_INCLUDED