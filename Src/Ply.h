/*

 Header for PLY polygon files.
 
  - Greg Turk, March 1994
  
   A PLY file contains a single polygonal _object_.
   
	An object is composed of lists of _elements_.  Typical elements are
	vertices, faces, edges and materials.
	
	 Each type of element for a given object has one or more _properties_
	 associated with the element type.  For instance, a vertex element may
	 have as properties three floating-point values x,y,z and three unsigned
	 chars for red, green and blue.
	 
	  ---------------------------------------------------------------
	  
	   Copyright (c) 1994 The Board of Trustees of The Leland Stanford
	   Junior University.  All rights reserved.   
	   
		Permission to use, copy, modify and distribute this software and its   
		documentation for any purpose is hereby granted without fee, provided   
		that the above copyright notice and this permission notice appear in   
		all copies of this software and that you do not sell the software.   
		
		 THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
		 EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
		 WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   
		 
*/


#ifndef __PLY_H__
#define __PLY_H__

#define NEW_PLY_CODE

#ifndef WIN32
#define _strdup strdup
#endif

#ifdef __cplusplus
extern "C" {
#endif
	
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
    
#define PLY_ASCII         1      /* ascii PLY file */
#define PLY_BINARY_BE     2      /* binary PLY file, big endian */
#define PLY_BINARY_LE     3      /* binary PLY file, little endian */
#define PLY_BINARY_NATIVE 4      /* binary PLY file, same endianness as current architecture */
    
#define PLY_OKAY    0           /* ply routine worked okay */
#define PLY_ERROR  -1           /* error in ply routine */
	
	/* scalar data types supported by PLY format */
	
#define PLY_START_TYPE 0
#define PLY_CHAR       1
#define PLY_SHORT      2
#define PLY_INT        3
#define PLY_UCHAR      4
#define PLY_USHORT     5
#define PLY_UINT       6
#define PLY_FLOAT      7
#define PLY_DOUBLE     8
#define PLY_INT_8      9
#define PLY_UINT_8     10
#define PLY_INT_16     11
#define PLY_UINT_16    12
#define PLY_INT_32     13
#define PLY_UINT_32    14
#define PLY_FLOAT_32   15
#define PLY_FLOAT_64   16
	
#define PLY_END_TYPE   17
	
#define  PLY_SCALAR  0
#define  PLY_LIST    1
	
#define PLY_STRIP_COMMENT_HEADER 0

typedef struct PlyProperty {    /* description of a property */
	char *name;                           /* property name */
	int external_type;                    /* file's data type */
	int internal_type;                    /* program's data type */
	int offset;                           /* offset bytes of prop in a struct */
	
	int is_list;                          /* 1 = list, 0 = scalar */
	int count_external;                   /* file's count type */
	int count_internal;                   /* program's count type */
	int count_offset;                     /* offset byte for list count */
}
PlyProperty;

typedef struct PlyElement {     /* description of an element */
	char *name;                   /* element name */
	int num;                      /* number of elements in this object */
	int size;                     /* size of element (bytes) or -1 if variable */
	int nprops;                   /* number of properties for this element */
	PlyProperty **props;          /* list of properties in the file */
	char *store_prop;             /* flags: property wanted by user? */
	int other_offset;             /* offset to un-asked-for props, or -1 if none*/
	int other_size;               /* size of other_props structure */
} PlyElement;

typedef struct PlyOtherProp {   /* describes other properties in an element */
	char *name;                   /* element name */
	int size;                     /* size of other_props */
	int nprops;                   /* number of properties in other_props */
	PlyProperty **props;          /* list of properties in other_props */
} PlyOtherProp;

typedef struct OtherData { /* for storing other_props for an other element */
	void *other_props;
} OtherData;

typedef struct OtherElem {     /* data for one "other" element */
	char *elem_name;             /* names of other elements */
	int elem_count;              /* count of instances of each element */
	OtherData **other_data;      /* actual property data for the elements */
	PlyOtherProp *other_props;   /* description of the property data */
} OtherElem;

typedef struct PlyOtherElems {  /* "other" elements, not interpreted by user */
	int num_elems;                /* number of other elements */
	OtherElem *other_list;        /* list of data for other elements */
} PlyOtherElems;

typedef struct PlyFile {        /* description of PLY file */
	FILE *fp;                     /* file pointer */
	int file_type;                /* ascii or binary */
	float version;                /* version number of file */
	int nelems;                   /* number of elements of object */
	PlyElement **elems;           /* list of elements */
	int num_comments;             /* number of comments */
	char **comments;              /* list of comments */
	int num_obj_info;             /* number of items of object information */
	char **obj_info;              /* list of object info items */
	PlyElement *which_elem;       /* which element we're currently writing */
	PlyOtherElems *other_elems;   /* "other" elements from a PLY file */
} PlyFile;
	
	/* memory allocation */
extern char *my_alloc();
#define myalloc(mem_size) my_alloc((mem_size), __LINE__, __FILE__)

#ifndef ALLOCN
#define REALLOCN(PTR,TYPE,OLD_N,NEW_N)							\
{										\
	if ((OLD_N) == 0)                                           		\
{   ALLOCN((PTR),TYPE,(NEW_N));}                            		\
	else									\
{								    		\
	(PTR) = (TYPE *)realloc((PTR),(NEW_N)*sizeof(TYPE));			\
	if (((PTR) == NULL) && ((NEW_N) != 0))					\
{									\
	fprintf(stderr, "Memory reallocation failed on line %d in %s\n", 	\
	__LINE__, __FILE__);                             		\
	fprintf(stderr, "  tried to reallocate %d->%d\n",       		\
	(OLD_N), (NEW_N));                              		\
	exit(-1);								\
}									\
	if ((NEW_N)>(OLD_N))							\
	memset((char *)(PTR)+(OLD_N)*sizeof(TYPE), 0,			\
	((NEW_N)-(OLD_N))*sizeof(TYPE));				\
}										\
}

#define  ALLOCN(PTR,TYPE,N) 					\
{ (PTR) = (TYPE *) calloc(((unsigned)(N)),sizeof(TYPE));\
	if ((PTR) == NULL) {    				\
	fprintf(stderr, "Memory allocation failed on line %d in %s\n", \
	__LINE__, __FILE__);                           \
	exit(-1);                                             \
	}							\
}


#define FREE(PTR)  { free((PTR)); (PTR) = NULL; }
#endif


/*** delcaration of routines ***/

extern PlyFile *ply_write(FILE *, int, const char **, int);
extern PlyFile *ply_open_for_writing( const char *, int, const char **, int, float *);
extern void ply_describe_element(PlyFile *, char *, int, int, PlyProperty *);
extern void ply_describe_property(PlyFile *, const char *, PlyProperty *);
extern void ply_element_count(PlyFile *, const char *, int);
extern void ply_header_complete(PlyFile *);
extern void ply_put_element_setup(PlyFile *, const char *);
extern void ply_put_element(PlyFile *, void *);
extern void ply_put_comment(PlyFile *, char *);
extern void ply_put_obj_info(PlyFile *, char *);
extern PlyFile *ply_read(FILE *, int *, char ***);
extern PlyFile *ply_open_for_reading( const char *, int *, char ***, int *, float *);
extern PlyProperty **ply_get_element_description(PlyFile *, char *, int*, int*);
extern void ply_get_element_setup( PlyFile *, char *, int, PlyProperty *);
extern int ply_get_property(PlyFile *, char *, PlyProperty *);
extern PlyOtherProp *ply_get_other_properties(PlyFile *, char *, int);
extern void ply_get_element(PlyFile *, void *);
extern char **ply_get_comments(PlyFile *, int *);
extern char **ply_get_obj_info(PlyFile *, int *);
extern void ply_close(PlyFile *);
extern void ply_get_info(PlyFile *, float *, int *);
extern PlyOtherElems *ply_get_other_element (PlyFile *, char *, int);
extern void ply_describe_other_elements ( PlyFile *, PlyOtherElems *);
extern void ply_put_other_elements (PlyFile *);
extern void ply_free_other_elements (PlyOtherElems *);
extern void ply_describe_other_properties(PlyFile *, PlyOtherProp *, int);

extern int equal_strings(const char *, const char *);

#ifdef __cplusplus
}
#endif
#include "Geometry.h"
#include <vector>

template< class Real > int PLYType( void );
template<> inline int PLYType< int           >( void ){ return PLY_INT   ; }
template<> inline int PLYType<          char >( void ){ return PLY_CHAR  ; }
template<> inline int PLYType< unsigned char >( void ){ return PLY_UCHAR ; }
template<> inline int PLYType<        float  >( void ){ return PLY_FLOAT ; }
template<> inline int PLYType<        double >( void ){ return PLY_DOUBLE; }
template< class Real > inline int PLYType( void ){ fprintf( stderr , "[ERROR] Unrecognized type\n" ) , exit( 0 ); }

typedef struct PlyFace
{
	unsigned char nr_vertices;
	int *vertices;
	int segment;
} PlyFace;
static PlyProperty face_props[] =
{
	{ _strdup( "vertex_indices" ) , PLY_INT , PLY_INT , offsetof( PlyFace , vertices ) , 1 , PLY_UCHAR, PLY_UCHAR , offsetof(PlyFace,nr_vertices) },
};


///////////////
// PlyVertex //
///////////////
#ifdef NEW_PLY_CODE
static char*       PlyPositionNames[] = { "x" , "y" , "z" ,  "w" };
static char*         PlyNormalNames[] = { "nx" , "ny" , "nz" , "nw" };
static char*          PlyValueNames[] = { "value" };
static char*          PlyColorNames[] = { "red" , "green" , "blue" , "alpha" };
static char* PlyAlternateColorNames[] = { "r" , "g" , "b" , "a" };
#else // !NEW_PLY_CODE
static char*       PlyPositionNames[] = { _strdup(  "x" ) , _strdup(  "y" ) , _strdup(  "z" ) , _strdup(  "w" ) };
static char*         PlyNormalNames[] = { _strdup( "nx" ) , _strdup( "ny" ) , _strdup( "nz" ) , _strdup( "nw" ) };
static char*          PlyValueNames[] = { _strdup( "value" ) };
static char*          PlyColorNames[] = { _strdup( "red"  ), _strdup( "green" ) , _strdup( "blue" ) , _strdup( "alpha" ) };
static char* PlyAlternateColorNames[] = { _strdup( "r" ) , _strdup( "g" ) , _strdup( "b" ) , _strdup( "a" ) };
#endif // NEW_PLY_CODE
inline PlyProperty MakePlyProperty( char* name , int internalType , int externalType , int offset )
{
	PlyProperty p;
	p.name = name;
	p.external_type = externalType;
	p.internal_type = internalType;
	p.offset = offset;
	p.is_list = p.count_external = p.count_internal = p.count_offset = 0;
	return p;
}

// The "Wrapper" class indicates the class to cast to/from in order to support linear operations.

struct RGBColor
{
	unsigned char c[3];
	unsigned char& operator[]( int idx )       { return c[idx]; }
	unsigned char  operator[]( int idx ) const { return c[idx]; }
	RGBColor( void ){ c[0] = c[1] = c[2] = 0; }
	RGBColor( const RGBColor& rgb ){ memcpy( c , rgb.c , sizeof(unsigned char) * 3 ); }
	RGBColor& operator = ( const RGBColor& rgb ){ memcpy( c , rgb.c , sizeof(unsigned char) * 3 ) ; return *this; }
};
///////////////////
// FullPlyVertex //
///////////////////
template< class _Real , int Dim , bool HasNormal , bool HasValue , bool HasColor >
class FullPlyVertex
{
public:
	typedef _Real Real;
protected:
	static PlyProperty _Properties[];
	static const int _NormalSize = HasNormal ? sizeof( Point< Real , Dim > ) : 0;
	static const int  _ValueSize = HasValue  ? sizeof( Real )                : 0;
	static const int  _ColorSize = HasColor  ? sizeof( RGBColor )            : 0;
	static const int _Size = _NormalSize + _ValueSize + _ColorSize;

	static const int _NormalOffset =  0;
	static const int  _ValueOffset = _NormalOffset + _NormalSize;
	static const int  _ColorOffset =  _ValueOffset +  _ValueSize;

	char _vertexData[ _Size==0 ? 1 :_Size ];
public:
	struct _FullPlyVertex
	{
		static const int  PointCount =               Dim      ;
		static const int NormalCount = ( HasNormal ? Dim : 0 );
		static const int  ValueCount = ( HasValue  ?   1 : 0 );
		static const int  ColorCount = ( HasColor  ?   3 : 0 );
		static const int Count = NormalCount + ValueCount + ColorCount;

		static const int NormalOffset = 0;
		static const int  ValueOffset = NormalOffset + NormalCount;
		static const int  ColorOffset =  ValueOffset +  ValueCount;

		Point< Real , Count==0 ? 1 : Count > vertexData;
		Point< Real , Dim > point;
		Point< Real , Dim >& normal( void ){ return *( ( Point< Real , Dim >*)(&vertexData[0] + NormalOffset) ); }
		Real               & value ( void ){ return *( ( Real               *)(&vertexData[0] +  ValueOffset) ); }
		Point< Real , 3   >& color ( void ){ return *( ( Point< Real , 3   >*)(&vertexData[0] +  ColorOffset) ); }
		const Point< Real , Dim >& normal( void ) const { return *( ( Point< Real , Dim >*)(&vertexData[0] + NormalOffset) ); }
		const Real               & value ( void ) const { return *( ( Real               *)(&vertexData[0] +  ValueOffset) ); }
		const Point< Real , 3   >& color ( void ) const { return *( ( Point< Real , 3   >*)(&vertexData[0] +  ColorOffset) ); }

		_FullPlyVertex( void ) {;}
		_FullPlyVertex( Point< Real , Dim > p , Point< Real , Count > vData ){ point = p , vertexData = vData; }
		_FullPlyVertex( FullPlyVertex p )
		{
			point = p.point;
			if( HasNormal ) for( int i=0 ; i<Dim ; i++ ) vertexData[NormalOffset+i] = p.normal()[i];
			if( HasValue  )                              vertexData[ ValueOffset  ] = p.value();
			if( HasColor  ) for( int i=0 ; i<3 ; i++ )   vertexData[ ColorOffset+i] = (Real)p.color()[i];
		}
		operator FullPlyVertex()
		{
			FullPlyVertex p;
			p.point = point;
			if( HasNormal ) for( int i=0 ; i<Dim ; i++ ) p.normal()[i] = vertexData[NormalOffset+i];
			if( HasValue  )                              p.value ()    = vertexData[ ValueOffset  ];
			if( HasColor  ) for( int i=0 ; i<3 ; i++ )   p.color ()[i] = (unsigned char)std::max< int >( 0 , std::min< int >( 255 , (int)( vertexData[ColorOffset+i]+0.5 ) ) );
			return p;
		}
		_FullPlyVertex& operator += ( const _FullPlyVertex& p ){ point += p.point , vertexData += p.vertexData ; return *this; }
		_FullPlyVertex& operator -= ( const _FullPlyVertex& p ){ point -= p.point , vertexData -= p.vertexData ; return *this; }
		_FullPlyVertex& operator *= ( Real s )                 { point *= s , vertexData *= s            ; return *this; }
		_FullPlyVertex& operator /= ( Real s )                 { point /= s , vertexData /= s            ; return *this; }
		_FullPlyVertex  operator +  ( const _FullPlyVertex& p ) const { return _FullPlyVertex( point + p.point , vertexData + p.vertexData ); }
		_FullPlyVertex  operator -  ( const _FullPlyVertex& p ) const { return _FullPlyVertex( point - p.point , vertexData - p.vertexData ); }
		_FullPlyVertex  operator *  ( Real s )                  const { return _FullPlyVertex( point * s , vertexData * s            ); }
		_FullPlyVertex  operator /  ( Real s )                  const { return _FullPlyVertex( point / s , vertexData / s            ); }
	};

	typedef _FullPlyVertex Wrapper;

	const static int  ReadComponents = Dim + ( HasNormal ? Dim : 0 ) + ( HasValue ? 1 : 0 ) + ( HasColor ? 6 : 0 );
	const static int WriteComponents = Dim + ( HasNormal ? Dim : 0 ) + ( HasValue ? 1 : 0 ) + ( HasColor ? 3 : 0 );

	static PlyProperty* Properties( void );

	Point< Real , Dim > point;
	Point< Real , Dim >& normal( void ){ return *( ( Point< Real , Dim >* )( _vertexData + _NormalOffset ) ); }
	Real&                 value( void ){ return *( ( Real* )               ( _vertexData +  _ValueOffset ) ); }
	RGBColor&             color( void ){ return *( ( RGBColor* )           ( _vertexData +  _ColorOffset ) ); }
	FullPlyVertex( void ){ memset( _vertexData , 0 , _Size ); }
	FullPlyVertex( Point< Real , Dim > p )  : FullPlyVertex() { point = p; }
	FullPlyVertex( const FullPlyVertex& v ) : FullPlyVertex() { point = v.point , memcpy( _vertexData , v._vertexData , _Size ); }
	FullPlyVertex& operator = ( const FullPlyVertex& v )      { point = v.point , memcpy( _vertexData , v._vertexData , _Size ) ; return *this; }
};
template< class Real , int Dim , bool HasNormal , bool HasValue , bool HasColor , class _Real >
FullPlyVertex< Real , Dim , HasNormal , HasValue , HasColor > operator * ( XForm< _Real , Dim+1 > xForm , FullPlyVertex< Real , Dim , HasNormal , HasValue , HasColor > p )
{
	FullPlyVertex< Real , Dim , HasNormal , HasValue , HasColor > _p = p;
	_p.point = xForm * p.point;
	if( HasNormal ) _p.normal() = xForm.inverse().transpose() * p.normal();
	return _p;
}

template< class Real , int Dim , bool HasNormal , bool HasValue , bool HasColor > PlyProperty FullPlyVertex< Real , Dim , HasNormal , HasValue , HasColor >::_Properties[ FullPlyVertex::ReadComponents ];
template< class Real , int Dim , bool HasNormal , bool HasValue , bool HasColor >
PlyProperty* FullPlyVertex< Real , Dim , HasNormal , HasValue , HasColor >::Properties( void )
{
	int idx = 0;

	// Primary values (for writing)
	int vertexDataOffset = (size_t)(&( ( (FullPlyVertex*)0 )->_vertexData ));
	int      pointOffset = (size_t)(&( ( (FullPlyVertex*)0 )->point ));
	for( int d=0 ; d<Dim ; d++ )                 _Properties[idx++] = MakePlyProperty( PlyPositionNames[d] , PLYType< Real >()          , PLYType< Real >()          ,      pointOffset                 + sizeof( Real )*d );
	if( HasNormal ) for( int d=0 ; d<Dim ; d++ ) _Properties[idx++] = MakePlyProperty(   PlyNormalNames[d] , PLYType< Real >()          , PLYType< Real >()          , vertexDataOffset + _NormalOffset + sizeof( Real )*d );
	if( HasValue )                               _Properties[idx++] = MakePlyProperty(    PlyValueNames[0] , PLYType< Real >()          , PLYType< Real >()          , vertexDataOffset +  _ValueOffset );
	if( HasColor ) for( int c=0 ; c<3 ; c++ )    _Properties[idx++] = MakePlyProperty(    PlyColorNames[c] , PLYType< unsigned char >() , PLYType< unsigned char >() , vertexDataOffset +  _ColorOffset + sizeof( unsigned char )*c );

	// Alternative values (for reading or writing)
	if( HasColor ) for( int c=0 ; c<3 ; c++ ) _Properties[idx++] = MakePlyProperty( PlyAlternateColorNames[c] , PLYType< unsigned char >() , PLYType< unsigned char >() , vertexDataOffset + _ColorOffset + sizeof( unsigned char )*c );

	return _Properties;
}

template< class Real , int Dim > using PlyVertex              = FullPlyVertex< Real , Dim , false , false , false >;
template< class Real , int Dim > using PlyOrientedVertex      = FullPlyVertex< Real , Dim , true  , false , false >;
template< class Real , int Dim > using PlyValueVertex         = FullPlyVertex< Real , Dim , false , true  , false >;
template< class Real , int Dim > using PlyColorVertex         = FullPlyVertex< Real , Dim , false , false , true  >;
template< class Real , int Dim > using PlyOrientedColorVertex = FullPlyVertex< Real , Dim , true  , false , true  >;
template< class Real , int Dim > using PlyColorAndValueVertex = FullPlyVertex< Real , Dim , false , true  , true  >;

template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >*  mesh , int file_type , const Point< float , Dim >& translate , float scale , char** comments=NULL , int commentNum=0 , XForm< Real , Dim+1 > xForm=XForm< Real , Dim+1 >::Identity() );

template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >*  mesh , int file_type , char** comments=NULL , int commentNum=0 , XForm< Real , Dim+1 > xForm=XForm< Real , Dim+1 >::Identity() );

inline bool PlyReadHeader( const char* fileName , PlyProperty* properties , int propertyNum , bool* readFlags , int& file_type )
{
	int nr_elems;
	char **elist;
	float version;
	PlyFile* ply;
	char* elem_name;
	int num_elems;
	int nr_props;
	PlyProperty** plist;

	ply = ply_open_for_reading( fileName , &nr_elems , &elist , &file_type , &version );
	if( !ply ) return false;

	for( int i=0 ; i<nr_elems ; i++ )
	{
		elem_name = elist[i];
		plist = ply_get_element_description( ply , elem_name , &num_elems , &nr_props );
		if( !plist )
		{
			for( int i=0 ; i<nr_elems ; i++ )
			{
				free( ply->elems[i]->name );
				free( ply->elems[i]->store_prop );
				for( int j=0 ; j<ply->elems[i]->nprops ; j++ )
				{
					free( ply->elems[i]->props[j]->name );
					free( ply->elems[i]->props[j] );
				}
				free( ply->elems[i]->props );
			}
			for( int i=0 ; i<nr_elems ; i++ ) free( ply->elems[i] );
			free( ply->elems );
			for( int i=0 ; i<ply->num_comments ; i++ ) free( ply->comments[i] );
			free( ply->comments );
			for( int i=0 ; i<ply->num_obj_info ; i++ ) free( ply->obj_info[i] );
			free( ply->obj_info );
			ply_free_other_elements( ply->other_elems );
			
			for( int i=0 ; i<nr_elems ; i++ ) free( elist[i] );
			free( elist );
			ply_close( ply );
			return 0;
		}		
		if( equal_strings( "vertex" , elem_name ) )
			for( int i=0 ; i<propertyNum ; i++ )
				if( readFlags ) readFlags[i] = ply_get_property( ply , elem_name , &properties[i] )!=0;

		for( int j=0 ; j<nr_props ; j++ )
		{
			free( plist[j]->name );
			free( plist[j] );
		}
		free( plist );
	}  // for each type of element
	
	for( int i=0 ; i<nr_elems ; i++ )
	{
		free( ply->elems[i]->name );
		free( ply->elems[i]->store_prop );
		for( int j=0 ; j<ply->elems[i]->nprops ; j++ )
		{
			free( ply->elems[i]->props[j]->name );
			free( ply->elems[i]->props[j] );
		}
		if( ply->elems[i]->props && ply->elems[i]->nprops ) free(ply->elems[i]->props);
	}
	for( int i=0 ; i<nr_elems ; i++ ) free(ply->elems[i]);
	free( ply->elems) ;
	for( int i=0 ; i<ply->num_comments ; i++ ) free( ply->comments[i] );
	free( ply->comments );
	for( int i=0 ; i<ply->num_obj_info ; i++ ) free( ply->obj_info[i] );
	free( ply->obj_info );
	ply_free_other_elements(ply->other_elems);
	
	
	for( int i=0 ; i<nr_elems ; i++ ) free( elist[i] );
	free( elist );
	ply_close( ply );
	return true;
}
inline bool PlyReadHeader( const char* fileName , PlyProperty* properties , int propertyNum , bool* readFlags )
{
	int file_type;
	return PlyReadHeader( fileName , properties , propertyNum , readFlags , file_type );
}


template<class Vertex>
int PlyReadPolygons( const char* fileName,
					std::vector<Vertex>& vertices,std::vector<std::vector<int> >& polygons,
					PlyProperty* properties,int propertyNum,
					int& file_type,
					char*** comments=NULL,int* commentNum=NULL , bool* readFlags=NULL );

template<class Vertex>
int PlyWritePolygons( const char* fileName,
					 const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons ,
					 PlyProperty* properties,int propertyNum,
					 int file_type,
					 char** comments=NULL,const int& commentNum=0);

template<class Vertex>
int PlyWritePolygons( const char* fileName,
					 const std::vector<Vertex>& vertices , const std::vector< std::vector< int > >& polygons,
					 PlyProperty* properties,int propertyNum,
					 int file_type,
					 char** comments,const int& commentNum)
{
	int nr_vertices=int(vertices.size());
	int nr_faces=int(polygons.size());
	float version;
	const char *elem_names[] = { "vertex" , "face" };
	PlyFile *ply = ply_open_for_writing( fileName , 2 , elem_names , file_type , &version );
	if (!ply){return 0;}
	
	//
	// describe vertex and face properties
	//
	ply_element_count(ply, "vertex", nr_vertices);
	for(int i=0;i<propertyNum;i++)
		ply_describe_property(ply, "vertex", &properties[i]);
	
	ply_element_count(ply, "face", nr_faces);
	ply_describe_property(ply, "face", &face_props[0]);
	
	// Write in the comments
	if(comments && commentNum)
		for(int i=0;i<commentNum;i++)
			ply_put_comment(ply,comments[i]);

	ply_header_complete(ply);
	
	// write vertices
	ply_put_element_setup(ply, "vertex");
	for (int i=0; i < int(vertices.size()); i++)
		ply_put_element(ply, (void *) &vertices[i]);

	// write faces
	PlyFace ply_face;
	int maxFaceVerts=3;
	ply_face.nr_vertices = 3;
	ply_face.vertices = new int[3];

	ply_put_element_setup(ply, "face");
	for (int i=0; i < nr_faces; i++)
	{
		if(int(polygons[i].size())>maxFaceVerts)
		{
			delete[] ply_face.vertices;
			maxFaceVerts=int(polygons[i].size());
			ply_face.vertices=new int[maxFaceVerts];
		}
		ply_face.nr_vertices=int(polygons[i].size());
		for(int j=0;j<ply_face.nr_vertices;j++)
			ply_face.vertices[j]=polygons[i][j];
		ply_put_element(ply, (void *) &ply_face);
	}

	delete[] ply_face.vertices;
	ply_close(ply);
	return 1;
}
template<class Vertex>
int PlyReadPolygons( const char* fileName,
					std::vector<Vertex>& vertices , std::vector<std::vector<int> >& polygons ,
					 PlyProperty* properties , int propertyNum ,
					int& file_type ,
					char*** comments , int* commentNum , bool* readFlags )
{
	int nr_elems;
	char **elist;
	float version;
	int i,j,k;
	PlyFile* ply;
	char* elem_name;
	int num_elems;
	int nr_props;
	PlyProperty** plist;
	PlyFace ply_face;

	ply = ply_open_for_reading( fileName , &nr_elems , &elist , &file_type , &version );
	if(!ply) return 0;

	if( comments )
	{
		(*comments)=new char*[*commentNum+ply->num_comments];
		for( int i=0 ; i<ply->num_comments ; i++ ) (*comments)[i] = _strdup(ply->comments[i]);
		*commentNum = ply->num_comments;
	}

	for (i=0; i < nr_elems; i++) {
		elem_name = elist[i];
		plist = ply_get_element_description(ply, elem_name, &num_elems, &nr_props);
		if(!plist)
		{
			for(i=0;i<nr_elems;i++){
				free(ply->elems[i]->name);
				free(ply->elems[i]->store_prop);
				for(j=0;j<ply->elems[i]->nprops;j++){
					free(ply->elems[i]->props[j]->name);
					free(ply->elems[i]->props[j]);
				}
				free(ply->elems[i]->props);
			}
			for(i=0;i<nr_elems;i++){free(ply->elems[i]);}
			free(ply->elems);
			for(i=0;i<ply->num_comments;i++) free(ply->comments[i]);
			free(ply->comments);
			for(i=0;i<ply->num_obj_info;i++) free(ply->obj_info[i]);
			free(ply->obj_info);
			ply_free_other_elements (ply->other_elems);
			
			for(i=0;i<nr_elems;i++){free(elist[i]);}
			free(elist);
			ply_close(ply);
			return 0;
		}		
		if (equal_strings("vertex", elem_name))
		{
			for( int i=0 ; i<propertyNum ; i++)
			{
				int hasProperty = ply_get_property(ply,elem_name,&properties[i]);
				if( readFlags ) readFlags[i] = (hasProperty!=0);
			}
			vertices.resize(num_elems);
			for (j=0; j < num_elems; j++)	ply_get_element (ply, (void *) &vertices[j]);
		}
		else if (equal_strings("face", elem_name))
		{
			ply_get_property (ply, elem_name, &face_props[0]);
			polygons.resize(num_elems);
			for (j=0; j < num_elems; j++)
			{
				ply_get_element (ply, (void *) &ply_face);
				polygons[j].resize(ply_face.nr_vertices);
				for(k=0;k<ply_face.nr_vertices;k++)	polygons[j][k]=ply_face.vertices[k];
				delete[] ply_face.vertices;
			}  // for, read faces
		}  // if face
		else{ply_get_other_element (ply, elem_name, num_elems);}

		for(j=0;j<nr_props;j++){
			free(plist[j]->name);
			free(plist[j]);
		}
		free(plist);
	}  // for each type of element
	
	for(i=0;i<nr_elems;i++){
		free(ply->elems[i]->name);
		free(ply->elems[i]->store_prop);
		for(j=0;j<ply->elems[i]->nprops;j++){
			free(ply->elems[i]->props[j]->name);
			free(ply->elems[i]->props[j]);
		}
		if(ply->elems[i]->props && ply->elems[i]->nprops){free(ply->elems[i]->props);}
	}
	for(i=0;i<nr_elems;i++){free(ply->elems[i]);}
	free(ply->elems);
	for(i=0;i<ply->num_comments;i++){free(ply->comments[i]);}
	free(ply->comments);
	for(i=0;i<ply->num_obj_info;i++){free(ply->obj_info[i]);}
	free(ply->obj_info);
	ply_free_other_elements (ply->other_elems);
	
	
	for(i=0;i<nr_elems;i++){free(elist[i]);}
	free(elist);
	ply_close(ply);
	return 1;
}

template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >* mesh , int file_type , const Point< float , Dim >& translate , float scale , char** comments , int commentNum , XForm< Real , Dim+1 > xForm )
{
	int i;
	int nr_vertices=int(mesh->outOfCorePointCount()+mesh->inCorePoints.size());
	int nr_faces=mesh->polygonCount();
	float version;
	const char *elem_names[] = { "vertex" , "face" };
	PlyFile *ply = ply_open_for_writing( fileName , 2 , elem_names , file_type , &version );
	if( !ply ) return 0;

	mesh->resetIterator();
	
	//
	// describe vertex and face properties
	//
	ply_element_count( ply , "vertex" , nr_vertices );
	for( int i=0 ; i<Vertex::Components ; i++ ) ply_describe_property( ply , "vertex" , &Vertex::Properties[i] );
	
	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		Vertex vertex = xForm * ( mesh->inCorePoints[i] * scale + translate );
		ply_put_element(ply, (void *) &vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		Vertex vertex;
		mesh->nextOutOfCorePoint( vertex );
		vertex = xForm * ( vertex * scale + translate );
		ply_put_element(ply, (void *) &vertex);		
	}  // for, write vertices
	
	// write faces
	std::vector< CoredVertexIndex > polygon;
	ply_put_element_setup( ply , "face" );
	for( i=0 ; i<nr_faces ; i++ )
	{
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		mesh->nextPolygon( polygon );
		ply_face.nr_vertices = int( polygon.size() );
		ply_face.vertices = new int[ polygon.size() ];
		for( int i=0 ; i<int(polygon.size()) ; i++ )
			if( polygon[i].inCore ) ply_face.vertices[i] = polygon[i].idx;
			else                    ply_face.vertices[i] = polygon[i].idx + int( mesh->inCorePoints.size() );
		ply_put_element( ply, (void *) &ply_face );
		delete[] ply_face.vertices;
	}  // for, write faces
	
	ply_close( ply );
	return 1;
}
template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >* mesh , int file_type , char** comments , int commentNum , XForm< Real , Dim+1 > xForm )
{
	int i;
	int nr_vertices=int(mesh->outOfCorePointCount()+mesh->inCorePoints.size());
	int nr_faces=mesh->polygonCount();
	float version;
	const char *elem_names[] = { "vertex" , "face" };
	PlyFile *ply = ply_open_for_writing( fileName , 2 , elem_names , file_type , &version );
	if( !ply ) return 0;

	mesh->resetIterator();
	
	//
	// describe vertex and face properties
	//
	ply_element_count( ply , "vertex" , nr_vertices );
	for( int i=0 ; i<Vertex::WriteComponents ; i++ ) ply_describe_property( ply , "vertex" , &Vertex::Properties()[i] );
	
	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		Vertex vertex = xForm * mesh->inCorePoints[i];
		ply_put_element(ply, (void *) &vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		Vertex vertex;
		mesh->nextOutOfCorePoint( vertex );
		vertex = xForm * ( vertex );
		ply_put_element(ply, (void *) &vertex);		
	}  // for, write vertices
	
	// write faces
	std::vector< CoredVertexIndex > polygon;
	ply_put_element_setup( ply , "face" );
	for( i=0 ; i<nr_faces ; i++ )
	{
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		mesh->nextPolygon( polygon );
		ply_face.nr_vertices = int( polygon.size() );
		ply_face.vertices = new int[ polygon.size() ];
		for( int i=0 ; i<int(polygon.size()) ; i++ )
			if( polygon[i].inCore ) ply_face.vertices[i] = polygon[i].idx;
			else                    ply_face.vertices[i] = polygon[i].idx + int( mesh->inCorePoints.size() );
		ply_put_element( ply, (void *) &ply_face );
		delete[] ply_face.vertices;
	}  // for, write faces
	
	ply_close( ply );
	return 1;
}
inline int PlyDefaultFileType(void){return PLY_ASCII;}

#endif /* !__PLY_H__ */
