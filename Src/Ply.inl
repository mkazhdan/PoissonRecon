/*
Copyright (c) 2020, Michael Kazhdan
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

namespace PLY
{
	inline int DefaultFileType( void ){ return PLY_ASCII; }

	template<> inline int Type< int           >( void ){ return PLY_INT   ; }
	template<> inline int Type<          char >( void ){ return PLY_CHAR  ; }
	template<> inline int Type< unsigned char >( void ){ return PLY_UCHAR ; }
	template<> inline int Type<        float  >( void ){ return PLY_FLOAT ; }
	template<> inline int Type<        double >( void ){ return PLY_DOUBLE; }
	template< class Scalar > inline int Type( void )
	{
		ERROR_OUT( "Unrecognized scalar type" );
		return -1;
	}

	template<> const std::string Traits<          int >::name="int";
	template<> const std::string Traits< unsigned int >::name="unsigned int";
	template<> const std::string Traits<          long >::name="long";
	template<> const std::string Traits< unsigned long >::name="unsigned long";
	template<> const std::string Traits<          long long >::name="long long";
	template<> const std::string Traits< unsigned long long >::name="unsigned long long";

	template<>
	PlyProperty Face<          int       >::Properties[] = { PlyProperty( "vertex_indices" , PLY_INT       , PLY_INT       , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face< unsigned int       >::Properties[] = { PlyProperty( "vertex_indices" , PLY_UINT      , PLY_UINT      , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face<          long long >::Properties[] = { PlyProperty( "vertex_indices" , PLY_LONGLONG  , PLY_LONGLONG  , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };
	template<>
	PlyProperty Face< unsigned long long >::Properties[] = { PlyProperty( "vertex_indices" , PLY_ULONGLONG , PLY_ULONGLONG , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };

	// Read
	template< typename VertexFactory >
	inline int ReadVertexHeader( std::string fileName , const VertexFactory &vFactory , bool *readFlags )
	{
		int fileType;
		std::vector< std::string > elist;
		float version;

		PlyFile *ply = PlyFile::Read( fileName , elist , fileType , version );
		if( !ply ) THROW( "could not create read ply file: " , fileName );

		for( int i=0 ; i<(int)elist.size() ; i++ ) if( elist[i]=="vertex" ) for( unsigned int j=0 ; j<vFactory.plyReadNum() ; j++ )
		{
			PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(j) : vFactory.plyReadProperty(j);
			readFlags[j] = ( ply->get_property( elist[i] , &prop )!=0 );
		}

		delete ply;
		return fileType;
	}

	template< typename VertexFactory >
	inline int ReadVertexHeader( std::string fileName , const VertexFactory &vFactory , bool *readFlags , std::vector< PlyProperty > &unprocessedProperties )
	{
		int fileType;
		std::vector< std::string > elist;
		float version;

		PlyFile *ply = PlyFile::Read( fileName, elist, fileType, version );
		if( !ply ) ERROR_OUT( "Failed to open ply file for reading: " , fileName );

		size_t numElems;
		std::vector< PlyProperty > plist = ply->get_element_description( "vertex" , numElems );
		if( !plist.size() ) ERROR_OUT( "Failed to get element description: vertex" );
		for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++ ) readFlags[i] = false;

		for( int i=0 ; i<plist.size() ; i++ )
		{
			bool found = false;
			for( unsigned int j=0 ; j<vFactory.plyReadNum() ; j++ )
			{
				PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(j) : vFactory.plyReadProperty(j);
				if( prop.name==plist[i].name ) found = readFlags[j] = true;
			}
			if( !found ) unprocessedProperties.push_back( plist[i] );
		}
		delete ply;
		return fileType;
	}

	inline int ReadVertexHeader( std::string fileName , std::vector< PlyProperty > &properties )
	{
		int fileType;
		std::vector< std::string > elist;
		float version;

		PlyFile *ply = PlyFile::Read( fileName, elist, fileType, version );
		if( !ply ) ERROR_OUT( "Failed to open ply file for reading: " , fileName );

		size_t numElems;
		std::vector< PlyProperty > plist = ply->get_element_description( "vertex" , numElems );
		for( int i=0 ; i<plist.size() ; i++ ) properties.push_back( plist[i] );
		delete ply;
		return fileType;
	}


	template< typename VertexFactory , typename Index >
	void ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::vector< Index > > &polygons , int &file_type , std::vector< std::string > &comments , bool *readFlags )
	{
		std::vector< std::string > elist = { std::string( "vertex" ) , std::string( "face" ) };
		float version;

		PlyFile *ply = PlyFile::Read( fileName , elist , file_type , version );
		if( !ply ) ERROR_OUT( "Could not create ply file for reading: " , fileName );

		comments.reserve( comments.size() + ply->comments.size() );
		for( int i=0 ; i<ply->comments.size() ; i++ ) comments.push_back( ply->comments[i] );

		for( int i=0 ; i<elist.size() ; i++ )
		{
			std::string &elem_name = elist[i];
			size_t num_elems;
			std::vector< PlyProperty > plist = ply->get_element_description( elem_name , num_elems );
			if( !plist.size() ) ERROR_OUT( "Could not read element properties: " , elem_name );
			if( elem_name=="vertex" )
			{
				for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
				{
					PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
					int hasProperty = ply->get_property( elem_name , &prop );
					if( readFlags ) readFlags[i] = (hasProperty!=0);
				}
				vertices.resize( num_elems , vFactory() );
				char *buffer = new char[ vFactory.bufferSize() ];
				for( size_t j=0 ; j<num_elems ; j++ )
				{
					if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
					else
					{
						ply->get_element( (void *)buffer );
						vFactory.fromBuffer( buffer , vertices[j] );
					}
				}
				delete[] buffer;
			}
			else if( elem_name=="face" )
			{
				ply->get_property( elem_name , Face< Index >::Properties );
				polygons.resize( num_elems );
				for( size_t j=0 ; j<num_elems ; j++ )
				{
					Face< Index > ply_face;
					ply->get_element( (void *)&ply_face );
					polygons[j].resize( ply_face.nr_vertices );
					for( unsigned int k=0 ; k<ply_face.nr_vertices ; k++ ) polygons[j][k] = ply_face.vertices[k];
					free( ply_face.vertices );
				}  // for, read faces
			}  // if face
			else ply->get_other_element( elem_name , num_elems );
		}  // for each type of element

		delete ply;
	}

	template< typename VertexFactory , typename Index >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::vector< Index > > &polygons , int file_type , const std::vector< std::string > &comments )
	{
		size_t nr_vertices = vertices.size();
		size_t nr_faces = polygons.size();
		float version;
		std::vector< std::string > elem_names = { std::string( "vertex" ) , std::string( "face" ) };
		PlyFile *ply = PlyFile::Write( fileName , elem_names , file_type , version );
		if( !ply ) ERROR_OUT( "Could not create ply file for writing: " , fileName );

		//
		// describe vertex and face properties
		//
		ply->element_count( "vertex", nr_vertices );
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
		{
			PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
			ply->describe_property( "vertex" , &prop );
		}
		ply->element_count( "face" , nr_faces );
		ply->describe_property( "face" , Face< Index >::Properties );

		// Write in the comments
		for( int i=0 ; i<comments.size() ; i++ ) ply->put_comment( comments[i] );
		ply->header_complete();

		// write vertices
		ply->put_element_setup( elem_names[0] );

		char *buffer = new char[ vFactory.bufferSize() ];
		for( size_t j=0 ; j<vertices.size() ; j++ )
		{
			if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
			else
			{
				vFactory.toBuffer( vertices[j] , buffer );
				ply->put_element( (void *)buffer );
			}
		}
		delete[] buffer;

		// write faces
		Face< Index > ply_face;
		int maxFaceVerts=3;
		ply_face.nr_vertices = 3;
		ply_face.vertices = new Index[3];

		ply->put_element_setup( elem_names[1] );
		for( size_t i=0 ; i<nr_faces ; i++ )
		{
			if( (int)polygons[i].size()>maxFaceVerts )
			{
				delete[] ply_face.vertices;
				maxFaceVerts = (int)polygons[i].size();
				ply_face.vertices=new Index[ maxFaceVerts ];
			}
			ply_face.nr_vertices = (int)polygons[i].size();
			for( size_t j=0 ; j<ply_face.nr_vertices ; j++ ) ply_face.vertices[j] = polygons[i][j];
			ply->put_element( (void *)&ply_face );
		}

		delete[] ply_face.vertices;
		delete ply;
	}

	template< typename VertexFactory , typename Index , class Real , int Dim , typename OutputIndex >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , CoredMeshData< typename VertexFactory::VertexType , Index >* mesh , int file_type , const std::vector< std::string > &comments , std::function< typename VertexFactory::VertexType ( typename VertexFactory::VertexType ) > xForm )
	{
		if( mesh->outOfCoreVertexNum()+mesh->inCoreVertices.size()>(size_t)std::numeric_limits< OutputIndex >::max() )
		{
			if( std::is_same< Index , OutputIndex >::value ) ERROR_OUT( "more vertices than can be represented using " , Traits< Index >::name );
			WARN( "more vertices than can be represented using " , Traits< OutputIndex >::name , " using " , Traits< Index >::name , " instead" );
			return WritePolygons< VertexFactory , Index , Real , Dim , Index >( fileName , vFactory , mesh , file_type , comments , xForm );
		}
		size_t nr_vertices = mesh->outOfCoreVertexNum()+mesh->inCoreVertices.size();
		size_t nr_faces = mesh->polygonNum();
		float version;
		std::vector< std::string > elem_names = { std::string( "vertex" ) , std::string( "face" ) };
		PlyFile *ply = PlyFile::Write( fileName , elem_names , file_type , version );
		if( !ply ) ERROR_OUT( "Could not create ply file for writing: " , fileName );

		mesh->resetIterator();

		//
		// describe vertex and face properties
		//
		ply->element_count( "vertex" , nr_vertices );
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
		{
			PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
			ply->describe_property( "vertex" , &prop );
		}
		ply->element_count( "face" , nr_faces );
		ply->describe_property( "face" , Face< OutputIndex >::Properties );

		// Write in the comments
		for( int i=0 ; i<comments.size() ; i++ ) ply->put_comment( comments[i] );
		ply->header_complete();

		// write vertices
		ply->put_element_setup( "vertex" );
		if( vFactory.isStaticallyAllocated() )
		{
			for( size_t i=0 ; i<mesh->inCoreVertices.size() ; i++ )
			{
				typename VertexFactory::VertexType vertex = xForm( mesh->inCoreVertices[i] );
				ply->put_element( (void *)&vertex );
			}
			for( size_t i=0; i<mesh->outOfCoreVertexNum() ; i++ )
			{
				typename VertexFactory::VertexType vertex = vFactory();
				mesh->nextOutOfCoreVertex( vertex );
				vertex = xForm( vertex );
				ply->put_element( (void *)&vertex );
			}
		}
		else
		{
			char *buffer = new char[ vFactory.bufferSize() ];
			for( size_t i=0 ; i<mesh->inCoreVertices.size() ; i++ )
			{
				typename VertexFactory::VertexType vertex = xForm( mesh->inCoreVertices[i] );
				vFactory.toBuffer( vertex , buffer );
				ply->put_element( (void *)buffer );
			}
			for( size_t i=0; i<mesh->outOfCoreVertexNum() ; i++ )
			{
				typename VertexFactory::VertexType vertex = vFactory();
				mesh->nextOutOfCoreVertex( vertex );
				vFactory.toBuffer( xForm( vertex ) , buffer );
				ply->put_element( (void *)buffer );
			}
			delete[] buffer;
		}

	   // write faces
		std::vector< CoredVertexIndex< Index > > polygon;
		ply->put_element_setup( "face" );
		for( size_t i=0 ; i<nr_faces ; i++ )
		{
			//
			// create and fill a struct that the ply code can handle
			//
			Face< OutputIndex > ply_face;
			mesh->nextPolygon( polygon );
			ply_face.nr_vertices = int( polygon.size() );
			ply_face.vertices = new OutputIndex[ polygon.size() ];
			for( int j=0 ; j<int(polygon.size()) ; j++ )
				if( polygon[j].inCore ) ply_face.vertices[j] = (OutputIndex)polygon[j].idx;
				else                    ply_face.vertices[j] = (OutputIndex)( polygon[j].idx + mesh->inCoreVertices.size() );
			ply->put_element( (void *)&ply_face );
			delete[] ply_face.vertices;
		}  // for, write faces

		delete ply;
	}
}
