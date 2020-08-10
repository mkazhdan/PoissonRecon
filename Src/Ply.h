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

#ifndef PLY_INCLUDED
#define PLY_INCLUDED

#include <vector>
#include <string>
#include <functional>
#include "PlyFile.h"
#include "Geometry.h"
#include "CoredMesh.h"
#include "MyMiscellany.h"

namespace PLY
{
	// Converts from C-type to PLY type
	template< class Scalar > int Type( void );

	// Converts from C-type to PLY name
	template< typename Integer > struct Traits{ static const std::string name; };

	// A structure representing a face
	template< typename Index >
	struct Face
	{
		unsigned int nr_vertices;
		Index *vertices;

		static PlyProperty Properties[];
	};

	int DefaultFileType( void );

	// PLY read header functionality

	// Get the properties (and return the file type)
	int ReadVertexHeader( std::string fileName , std::vector< PlyProperty > &properties );

	// Test which properties are represented by elements of the vertex factory (and return the file type)
	template< typename VertexFactory >
	int ReadVertexHeader( std::string fileName , const VertexFactory &vFactory , bool *readFlags );

	// Test which properties are represented by elements of the vertex factory and add the others to the property list (and return the file type)
	template< typename VertexFactory >
	int ReadVertexHeader( std::string fileName , const VertexFactory &vFactory , bool *readFlags , std::vector< PlyProperty > &unprocessedProperties );

	// PLY write mesh functionality
	template< typename VertexFactory , typename Index , class Real , int Dim , typename OutputIndex=int >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , CoredMeshData< typename VertexFactory::VertexType , Index > *mesh , int file_type , const std::vector< std::string >& comments , std::function< typename VertexFactory::VertexType ( typename VertexFactory::VertexType ) > xForm = []( typename VertexFactory::VertexType v ){ return v; } );

	template< typename VertexFactory , typename Index >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::vector< Index > > &polygons , int file_type , const std::vector< std::string > &comments );

	// PLY read mesh functionality
	template< typename VertexFactory , typename Index >
	void ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::vector< Index > >& polygons , int &file_type , std::vector< std::string > &comments , bool* readFlags=NULL );
}
#include "Ply.inl"
#endif // PLY_INCLUDED
