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
#include <vector>
#include "PlyFile.h"
#include "Geometry.h"

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


struct RGBColor
{
	unsigned char c[3];
	unsigned char& operator[]( int idx )       { return c[idx]; }
	unsigned char  operator[]( int idx ) const { return c[idx]; }
	RGBColor( void ){ c[0] = c[1] = c[2] = 0; }
	RGBColor( const RGBColor& rgb ){ memcpy( c , rgb.c , sizeof(unsigned char) * 3 ); }
	RGBColor& operator = ( const RGBColor& rgb ){ memcpy( c , rgb.c , sizeof(unsigned char) * 3 ) ; return *this; }
};
///////////////
// PlyVertex //
///////////////
template< typename _Real , int Dim , typename _RealOnDisk=float >
class PlyVertex
{
public:
	typedef _Real Real;

	PlyVertex& operator += ( const PlyVertex& p ){ point += p.point ; return *this; }
	PlyVertex& operator -= ( const PlyVertex& p ){ point -= p.point ; return *this; }
	PlyVertex& operator *= ( Real s )            { point *= s ; return *this; }
	PlyVertex& operator /= ( Real s )            { point /= s ; return *this; }
	PlyVertex  operator +  ( const PlyVertex& p ) const { return PlyVertex( point + p.point ); }
	PlyVertex  operator -  ( const PlyVertex& p ) const { return PlyVertex( point - p.point ); }
	PlyVertex  operator *  ( Real s )             const { return PlyVertex( point * s ); }
	PlyVertex  operator /  ( Real s )             const { return PlyVertex( point / s ); }

	const static int PlyReadNum = Dim;
	const static int PlyWriteNum = Dim;

	static const PlyProperty* PlyReadProperties( void ){ return _PlyProperties; }
	static const PlyProperty* PlyWriteProperties( void ){ return _PlyProperties; }

	Point< Real , Dim > point;
	PlyVertex( void ) {}
	PlyVertex( Point< Real , Dim > p ) : point(p) { }

	struct Transform
	{
		Transform( void ){}
		Transform( const XForm< Real , Dim+1 >& xForm ) : _pointXForm(xForm) { }
		PlyVertex operator() ( const PlyVertex& p ) const
		{
			PlyVertex _p;
			_p.point = _pointXForm * p.point;
			return _p;
		}
	protected:
		XForm< Real , Dim+1 > _pointXForm;
	};

protected:
	static const PlyProperty _PlyProperties[];
};

template<>
const PlyProperty PlyVertex< float , 2 , float >::_PlyProperties[] =
{
	{ "x" , PLY_FLOAT , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_FLOAT , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< double , 2 , float >::_PlyProperties[] =
{
	{ "x" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< float , 2 , double >::_PlyProperties[] =
{
	{ "x" , PLY_DOUBLE , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_DOUBLE , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< double , 2 , double >::_PlyProperties[] =
{
	{ "x" , PLY_DOUBLE , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_DOUBLE , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
};

template<>
const PlyProperty PlyVertex< float , 3 , float >::_PlyProperties[] =
{
	{ "x" , PLY_FLOAT , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_FLOAT , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "z" , PLY_FLOAT , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< double , 3 , float >::_PlyProperties[] =
{
	{ "x" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "z" , PLY_FLOAT , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< float , 3 , double >::_PlyProperties[] =
{
	{ "x" , PLY_DOUBLE , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_DOUBLE , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "z" , PLY_DOUBLE , PLY_FLOAT , int( offsetof( PlyVertex , point.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};
template<>
const PlyProperty PlyVertex< double , 3 , double >::_PlyProperties[] =
{
	{ "x" , PLY_DOUBLE , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[0] ) ) , 0 , 0 , 0 , 0 } ,
	{ "y" , PLY_DOUBLE , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[1] ) ) , 0 , 0 , 0 , 0 } ,
	{ "z" , PLY_DOUBLE , PLY_DOUBLE , int( offsetof( PlyVertex , point.coords[2] ) ) , 0 , 0 , 0 , 0 } ,
};

///////////////////////
// PlyVertexWithData //
///////////////////////
template< typename _Real , int Dim , typename Data , typename _RealOnDisk=float >
class PlyVertexWithData
{
public:
	typedef _Real Real;

	PlyVertexWithData& operator += ( const PlyVertexWithData& p ){ point += p.point , data += p.data ; return *this; }
	PlyVertexWithData& operator -= ( const PlyVertexWithData& p ){ point -= p.point , data -= p.data ; return *this; }
	PlyVertexWithData& operator *= ( Real s )                    { point *= s , data *= s ; return *this; }
	PlyVertexWithData& operator /= ( Real s )                    { point /= s , data /= s ; return *this; }
	PlyVertexWithData  operator +  ( const PlyVertexWithData& p ) const { return PlyVertexWithData( point + p.point , data + p.data ); }
	PlyVertexWithData  operator -  ( const PlyVertexWithData& p ) const { return PlyVertexWithData( point - p.point , data - p.data ); }
	PlyVertexWithData  operator *  ( Real s )                     const { return PlyVertexWithData( point * s , data * s ); }
	PlyVertexWithData  operator /  ( Real s )                     const { return PlyVertexWithData( point / s , data / s ); }

	const static int PlyReadNum = Data::PlyReadNum + Dim;
	const static int PlyWriteNum = Data::PlyWriteNum + Dim;

	static const PlyProperty* PlyReadProperties( void ){ _SetReadProperties() ; return _PlyReadProperties; }
	static const PlyProperty* PlyWriteProperties( void ){ _SetWriteProperties() ; return _PlyWriteProperties; }


	Point< Real , Dim > point;
	Data data;
	PlyVertexWithData( void ) {}
	PlyVertexWithData( Point< Real , Dim > p , Data d ) : point(p) , data(d) { }

	struct Transform
	{
		Transform( void ){}
		Transform( const XForm< Real , Dim+1 >& xForm ) : _pointXForm(xForm) , _dataXForm(xForm) { }
		PlyVertexWithData operator() ( const PlyVertexWithData& p ) const
		{
			PlyVertexWithData _p;
			_p.point = _pointXForm * p.point;
			_p.data = _dataXForm( p.data );
			return _p;
		}
	protected:
		XForm< Real , Dim+1 > _pointXForm;
		typename Data::Transform _dataXForm;
	};

protected:
	static void _SetReadProperties( void );
	static void _SetWriteProperties( void );
	static PlyProperty _PlyReadProperties[];
	static PlyProperty _PlyWriteProperties[];
};
template< typename Real , int Dim , typename Data , typename RealOnDisk > PlyProperty PlyVertexWithData< Real , Dim , Data , RealOnDisk >::_PlyReadProperties[ PlyReadNum ];
template< typename Real , int Dim , typename Data , typename RealOnDisk > PlyProperty PlyVertexWithData< Real , Dim , Data , RealOnDisk >::_PlyWriteProperties[ PlyWriteNum ];
template< typename Real , int Dim , typename Data , typename RealOnDisk >
void PlyVertexWithData< Real , Dim , Data , RealOnDisk >::_SetReadProperties( void )
{
	{
		const PlyProperty * ReadProps = PlyVertex< Real , Dim , RealOnDisk >::PlyReadProperties();
		for( int d=0 ; d<PlyVertex< Real , Dim , RealOnDisk >::PlyReadNum ; d++ ) _PlyReadProperties[d] = ReadProps[d];
	}
	{
		const PlyProperty * ReadProps = Data::PlyReadProperties();
		for( int d=0 ; d<Data::PlyReadNum ; d++ )
		{
			_PlyReadProperties[d+PlyVertex< Real , Dim , RealOnDisk >::PlyReadNum ] = ReadProps[d];
			_PlyReadProperties[d+PlyVertex< Real , Dim , RealOnDisk >::PlyReadNum ].offset += (int)offsetof( PlyVertexWithData , data );
		}
	}
}
template< typename Real , int Dim , typename Data , typename RealOnDisk >
void PlyVertexWithData< Real , Dim , Data , RealOnDisk >::_SetWriteProperties( void )
{
	{
		const PlyProperty * WriteProps = PlyVertex< Real , Dim , RealOnDisk >::PlyWriteProperties();
		for( int d=0 ; d<PlyVertex< Real , Dim , RealOnDisk >::PlyWriteNum ; d++ ) _PlyWriteProperties[d] = WriteProps[d];
	}
	{
		const PlyProperty * WriteProps = Data::PlyWriteProperties();
		for( int d=0 ; d<Data::PlyWriteNum ; d++ )
		{
			_PlyWriteProperties[d+PlyVertex< Real , Dim , RealOnDisk >::PlyWriteNum ] = WriteProps[d];
			_PlyWriteProperties[d+PlyVertex< Real , Dim , RealOnDisk >::PlyWriteNum ].offset += (int)offsetof( PlyVertexWithData , data );
		}
	}
}

template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >*  mesh , int file_type , const Point< float , Dim >& translate , float scale , char** comments=NULL , int commentNum=0 , XForm< Real , Dim+1 > xForm=XForm< Real , Dim+1 >::Identity() );

template< class Vertex , class Real , int Dim >
int PlyWritePolygons( const char* fileName , CoredMeshData< Vertex >*  mesh , int file_type , char** comments=NULL , int commentNum=0 , XForm< Real , Dim+1 > xForm=XForm< Real , Dim+1 >::Identity() );

inline bool PlyReadHeader( char* fileName , const PlyProperty* properties , int propertyNum , bool* readFlags , int& file_type )
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
inline bool PlyReadHeader( char* fileName , const PlyProperty* properties , int propertyNum , bool* readFlags )
{
	int file_type;
	return PlyReadHeader( fileName , properties , propertyNum , readFlags , file_type );
}


template<class Vertex>
int PlyReadPolygons( const char* fileName,
					std::vector<Vertex>& vertices,std::vector<std::vector<int> >& polygons,
					const PlyProperty* properties,int propertyNum,
					int& file_type,
					char*** comments=NULL,int* commentNum=NULL , bool* readFlags=NULL );

template<class Vertex>
int PlyWritePolygons( const char* fileName,
					 const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons ,
					 const PlyProperty* properties,int propertyNum,
					 int file_type,
					 char** comments=NULL,const int& commentNum=0);

template<class Vertex>
int PlyWritePolygons( const char* fileName,
					 const std::vector<Vertex>& vertices , const std::vector< std::vector< int > >& polygons,
					 const PlyProperty* properties,int propertyNum,
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
					const PlyProperty* properties , int propertyNum ,
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
	typename Vertex::Transform _xForm( xForm );
	const PlyProperty* PlyWriteProperties = Vertex::PlyWriteProperties();
	for( int i=0 ; i<Vertex::PlyWriteNum ; i++ ) ply_describe_property( ply , "vertex" , &PlyWriteProperties[i] );

	ply_element_count( ply , "face" , nr_faces );
	ply_describe_property( ply , "face" , &face_props[0] );
	
	// Write in the comments
	for( i=0 ; i<commentNum ; i++ ) ply_put_comment( ply , comments[i] );

	ply_header_complete( ply );
	
	// write vertices
	ply_put_element_setup( ply , "vertex" );
	for( i=0 ; i<int( mesh->inCorePoints.size() ) ; i++ )
	{
		Vertex vertex = _xForm( mesh->inCorePoints[i] );
		ply_put_element(ply, (void *) &vertex);
	}
	for( i=0; i<mesh->outOfCorePointCount() ; i++ )
	{
		Vertex vertex;
		mesh->nextOutOfCorePoint( vertex );
		vertex = _xForm( vertex );
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
