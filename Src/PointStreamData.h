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
#ifndef POINT_STREAM_DATA_INCLUDED
#define POINT_STREAM_DATA_INCLUDED

#include <algorithm>
#include "Ply.h"

template< class Real > using Color = Point< Real , 3 >;
template< class Real > void SetColorValues( const Color< Real >& color , unsigned char c[3] ){ for( int i=0 ; i<3 ; i++ ) c[i] = (unsigned char)std::max< int >( 0 , std::min< int >( 255 , (int)( color[i]+0.5 ) ) ); }
template< class Real > void SetColorValues( const Color< Real >& color , RGBColor& c ){ for( int i=0 ; i<3 ; i++ ) c[i] = (unsigned char)std::max< int >( 0 , std::min< int >( 255 , (int)( color[i]+0.5 ) ) ); }

template< class Real , int Dim >
struct Normal
{
	Point< Real , Dim > normal;
	Normal( void ){;}
	Normal( Point< Real , Dim > n ) : normal( n ) {;}
	Normal( const Normal& n ) : normal( n.normal ) {;}

	Normal& operator += ( const Normal& p ){ normal += p.normal ; return *this; }
	Normal& operator -= ( const Normal& p ){ normal -= p.normal ; return *this; }
	Normal& operator *= ( Real s ){ normal *= s ; return *this; }
	Normal& operator /= ( Real s ){ normal /= s ; return *this; }
	Normal  operator +  ( const Normal& p ) const { return Normal( normal+p.normal ); }
	Normal  operator -  ( const Normal& p ) const { return Normal( normal-p.normal ); }
	Normal  operator *  ( Real s ) const { return Normal( normal*s ); }
	Normal  operator /  ( Real s ) const { return Normal( normal/s ); }
};
template< class Real , int Dim >
struct NormalAndColor
{
	Point< Real , Dim > normal;
	Color< Real > color;
	NormalAndColor( void ){;}
	NormalAndColor( Point< Real , Dim > n , Color< Real > c ) : normal(n) , color(c) {;}

	NormalAndColor& operator += ( const NormalAndColor& p ){ color += p.color , normal += p.normal ; return *this; }
	NormalAndColor& operator -= ( const NormalAndColor& p ){ color -= p.color , normal -= p.normal ; return *this; }
	NormalAndColor& operator *= ( Real s ){ color *= s , normal *= s ; return *this; }
	NormalAndColor& operator /= ( Real s ){ color /= s , normal /= s ; return *this; }
	NormalAndColor  operator +  ( const NormalAndColor& p ) const { return NormalAndColor( normal+p.normal , color+p.color ); }
	NormalAndColor  operator -  ( const NormalAndColor& p ) const { return NormalAndColor( normal-p.normal , color-p.color ); }
	NormalAndColor  operator *  ( Real s ) const { return NormalAndColor( normal*s , color*s ); }
	NormalAndColor  operator /  ( Real s ) const { return NormalAndColor( normal/s , color/s ); }
};

template< class Real , int Dim >
struct NormalInfo
{
	typedef Normal< Real , Dim > Type;
	static Type ReadASCII( FILE* fp )
	{
		float n[3];
		if( fscanf( fp , " %f %f %f " , &n[0] , &n[1] , &n[2] )!=3 ) fprintf( stderr , "[ERROR] Failed to read normal\n" ) , exit( 0 );
		return Type( Point< Real , Dim >( n[0] , n[1] , n[2] ) );
	};
	static bool ReadBinary( FILE* fp , Point< Real , Dim >& p , Type& d )
	{
		struct TypeOnDisk { Point< float , Dim > point , normal; };
		TypeOnDisk t;
		if( fread( &t , sizeof(TypeOnDisk) , 1 , fp )!=1 ) return false;
		p = Point< Real , Dim >(t.point);
		d.normal = Point< Real , Dim >(t.normal);
		return true;
	}
	static void WriteASCII( FILE* fp , const Type& d )
	{
		fprintf( fp , " %f %f %f " , d.normal[0] , d.normal[1] , d.normal[2] );
	};

	static void WriteBinary( FILE* fp , const Point< Real , Dim >& p , const Type& d )
	{
		struct TypeOnDisk { Point< float , Dim > point , normal; };
		TypeOnDisk t;
		t.point = Point< float , Dim >( p ) , t.normal = Point< float , Dim >( d.normal );
		fwrite( &t , sizeof(TypeOnDisk) , 1 , fp );
	}

	static bool ValidPlyProperties( const bool* props )
	{
		if( !( props[0] && props[1] && props[2] ) ) printf( "%s %s %s\n" , props[0] ? "true" : "false" , props[1] ? "true" : "false" , props[2] ? "true" : "false" );
		return ( props[0] && props[1] && props[2] );
	}
	const static PlyProperty PlyProperties[];
	const static int PlyPropertyNum = 3;

	template< class Vertex >
	struct VertexSetter
	{
		static void SetValue( Vertex& v , Real w ){ ; }
		static void SetData ( Vertex& v , const Type& nc ){ ; }
	};
	template< bool HasNormal , bool HasValue , bool HasColor >
	struct VertexSetter< FullPlyVertex< float , Dim , HasNormal , HasValue , HasColor > >
	{
		typedef FullPlyVertex< float , Dim , HasNormal , HasValue , HasColor > Vertex;
		static void SetValue( Vertex& v , Real w ){ if( HasValue ) v.value() = (float)w; }
		static void SetData ( Vertex& v , const Type& nc ){ if( HasNormal ) v.normal() = Point< Real , Dim >( nc.normal ); }
	};

	static Real ProcessDataWithConfidence( const Point< Real , Dim >& p , Type& data )
	{
		Real l = (Real)Length( data.normal );
		if( !l || l!=l ) return (Real)-1.;
		return l;
	}
	static Real ProcessData( const Point< Real , Dim >& p , Type& data )
	{
		Real l = (Real)Length( data.normal );
		if( !l || l!=l ) return (Real)-1.;
		data.normal /= l;
		return (Real)1.;
	}

	struct Transform
	{
		Transform( const XForm< Real , Dim+1 >& xForm ) : _pointXForm( xForm )
		{
			for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _normalXForm(i,j) = _pointXForm(i,j);
			_normalXForm = _normalXForm.transpose().inverse();
			_normalXForm /= (Real)pow( fabs( _normalXForm.determinant() ) , 1./Dim );
		}
		void operator()( Point< Real , Dim >& p , Type& n ) const { p = _pointXForm * p , n.normal = _normalXForm * n.normal; }
	protected:
		XForm< Real , Dim+1 > _pointXForm;
		XForm< Real , Dim > _normalXForm;
	};
};
template<>
const PlyProperty NormalInfo< float , 3 >::PlyProperties[] =
{
	{ "nx" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "ny" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "nz" , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty NormalInfo< double , 3 >::PlyProperties[] =
{
	{ "nx" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "ny" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "nz" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};

template< class Real , int Dim >
struct NormalAndColorInfo
{
	typedef NormalAndColor< Real , Dim > Type;

	static Type ReadASCII( FILE* fp )
	{
		float n[3];
		unsigned char c[3];
		if( fscanf( fp , " %f %f %f %c %c %c " , &n[0] , &n[1] , &n[2] , &c[0] , &c[1] , &c[2] )!=6 ) fprintf( stderr , "[ERROR] Failed to read normal and color\n" ) , exit( 0 );
		return Type( Point< Real , Dim >( n[0] , n[1] , n[2] ) , Color< Real >( (Real)c[0] , (Real)c[1] , (Real)c[2] ) );
	};

	static bool ReadBinary( FILE* fp , Point< Real , Dim >& p , Type& d )
	{
		struct TypeOnDisk
		{
			Point< float , Dim > point;
			Point< float , Dim > normal;
			unsigned char color[3];
		};
		TypeOnDisk t;
		if( fread( &t , sizeof( TypeOnDisk ) , 1 , fp )!=1 ) return false;
		p = Point< Real , Dim >(t.point);
		d.normal = Point< Real , Dim >(t.normal);
		for( int c=0 ; c<3 ; c++ ) d.color[c] = (Real)t.color[c];
		return true;
	}
	static void WriteASCII( FILE* fp , const Type& d )
	{
		unsigned char c[3];
		SetColorValues( d.color , c );
		fprintf( fp , " %f %f %f  %d %d %d " , d.normal[0] , d.normal[1] , d.normal[2] , c[0] , c[1] , c[2] );
	};
	static void WriteBinary( FILE* fp , const Point< Real , Dim >& p , const Type& d )
	{
		struct TypeOnDisk
		{
			Point< float , Dim > point;
			Point< float , Dim > normal;
			unsigned char color[3];
		};
		TypeOnDisk t;
		t.point = Point< float , Dim >( p ) , t.normal = Point< float , Dim >( d.normal );
		SetColorValues( d.color , t.color );
		fwrite( &t , sizeof(TypeOnDisk) , 1 , fp );
	}

	static bool ValidPlyProperties( const bool* props ){ return ( props[0] && props[1] && props[2] ) && ( ( props[3] || props[6] ) && ( props[4] || props[7] ) && ( props[5] || props[8] ) ); }
	const static PlyProperty PlyProperties[];
	const static int PlyPropertyNum = 9;

	template< class Vertex >
	struct VertexSetter
	{
		static void SetValue( Vertex& v , Real w ){}
		static void SetData ( Vertex& v , const Type& nc ){}
	};
	template< bool HasNormal , bool HasValue , bool HasColor >
	struct VertexSetter< FullPlyVertex< float , Dim , HasNormal , HasValue , HasColor > >
	{
		typedef FullPlyVertex< float , Dim , HasNormal , HasValue , HasColor > Vertex;
		static void SetValue( Vertex& v , Real w ){ if( HasValue ) v.value() = (float)w; }
		static void SetData ( Vertex& v , const Type& nc )
		{
			if( HasNormal ) v.normal() = Point< Real , Dim >( nc.normal );
			if( HasColor ) SetColorValues( nc.color , v.color() );
		}
	};

	static Real ProcessDataWithConfidence( const Point< Real , Dim >& p , Type& data )
	{
		Real l = (Real)Length( data.normal );
		if( !l || l!=l ) return (Real)-1.;
		return l;
	}
	static Real ProcessData( const Point< Real , Dim >& p , Type& data )
	{
		Real l = (Real)Length( data.normal );
		if( !l || l!=l ) return (Real)-1.;
		data.normal /= l;
		return (Real)1.;
	}

	struct Transform
	{
		Transform( const XForm< Real , Dim+1 >& xForm ) : _pointXForm( xForm )
		{
			for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _normalXForm(i,j) = _pointXForm(i,j);
			_normalXForm = _normalXForm.transpose().inverse();
			_normalXForm /= (Real)pow( fabs( _normalXForm.determinant() ) , 1./Dim );
		}
		void operator()( Point< Real , Dim >& p , Type& nc ) const { p = _pointXForm * p , nc.normal = _normalXForm * nc.normal; }
	protected:
		XForm< Real , Dim+1 > _pointXForm;
		XForm< Real , Dim > _normalXForm;
	};
};
template<>
const PlyProperty NormalAndColorInfo< float , 3 >::PlyProperties[] =
{
	{ "nx"    , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "ny"    , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "nz"    , PLY_FLOAT , PLY_FLOAT , int( offsetof( Type , normal.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "r"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "g"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "b"     , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "red"   , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[0] ) ) , 0 , 0 , 0 , 0 } , 
	{ "green" , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "blue"  , PLY_UCHAR , PLY_FLOAT , int( offsetof( Type , color .coords[2] ) ) , 0 , 0 , 0 , 0 }
};
template<>
const PlyProperty NormalAndColorInfo< double , 3 >::PlyProperties[] =
{
	{ "nx"    , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "ny"    , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "nz"    , PLY_FLOAT , PLY_DOUBLE , int( offsetof( Type , normal.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "r"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "g"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "b"     , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[2] ) ) , 0 , 0 , 0 , 0 } ,
	{ "red"   , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[0] ) ) , 0 , 0 , 0 , 0 } , 
	{ "green" , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "blue"  , PLY_UCHAR , PLY_DOUBLE , int( offsetof( Type , color .coords[2] ) ) , 0 , 0 , 0 , 0 }
};
#endif // POINT_STREAM_DATA_INCLUDED