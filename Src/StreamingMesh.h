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

#ifndef STREAMING_MESH_INCLUDED
#define STREAMING_MESH_INCLUDED

#include "Geometry.h"
#include "MyMiscellany.h"

template< typename Vertex , typename Index >
class StreamingVertices
{
public:
	virtual void resetIterator( void ) = 0;
	virtual Index addVertex( const Vertex &v ) = 0;
	virtual bool nextVertex( Vertex &v ) = 0;
	virtual size_t vertexNum( void ) const = 0;
};

template< class Vertex , typename Index >
class StreamingCurve : public StreamingVertices< Vertex , Index >
{
public:
	virtual void resetIterator( void ) = 0;
	virtual void addEdge_s( unsigned int thread , Index v1 , Index v2 ) = 0;
	virtual bool nextEdge( std::pair< Index , Index > &e ) = 0;
	virtual size_t edgeNum( void ) const = 0;
};

template< class Vertex , typename Index >
class StreamingMesh : public StreamingVertices< Vertex , Index >
{
public:
	virtual ~StreamingMesh( void ){}
	virtual void resetIterator( void ) = 0;
	virtual void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices ) = 0;
	virtual bool nextPolygon( std::vector< Index >& vertices ) = 0;
	virtual size_t polygonNum( void ) const = 0;
};

template< typename Vertex , typename Index >
class VectorStreamingVertices : public StreamingVertices< Vertex , Index >
{
	std::vector< Vertex > _vertices;
	Index _vertexIndex;
public:
	VectorStreamingVertices( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v );

	bool nextVertex( Vertex &v );

	size_t vertexNum( void ) const;

	Vertex &vertex( unsigned int idx ){ return _vertices[idx]; }
	const Vertex &vertex( unsigned int idx ) const { return _vertices[idx]; }
};

template< class Vertex , typename Index >
class VectorStreamingCurve : public StreamingCurve< Vertex , Index > , public VectorStreamingVertices< Vertex , Index >
{
	std::vector< std::vector< std::pair< Index , Index > > > _edges;
	unsigned int _threadIndex;
	Index _edgeIndex;
public:
	VectorStreamingCurve( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v ){ return VectorStreamingVertices< Vertex , Index >::addVertex( v ); }
	void addEdge_s( unsigned int thread , Index v1 , Index v2 );

	bool nextVertex( Vertex &v ){ return VectorStreamingVertices< Vertex , Index >::nextVertex(v); }
	bool nextEdge( std::pair< Index , Index > &e );

	size_t vertexNum( void ) const { return VectorStreamingVertices< Vertex , Index >::vertexNum(); }
	size_t edgeNum( void ) const;
};

template< class Vertex , typename Index >
class VectorStreamingMesh : public StreamingMesh< Vertex , Index > , public VectorStreamingVertices< Vertex , Index >
{
	std::vector< std::vector< std::vector< Index > > > _polygons;
	unsigned int _threadIndex;
	Index _polygonIndex;
public:
	VectorStreamingMesh( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v ){ return VectorStreamingVertices< Vertex , Index >::addVertex( v ); }
	void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices );

	bool nextVertex( Vertex &v ){ return VectorStreamingVertices< Vertex , Index >::nextVertex(v); }
	bool nextPolygon( std::vector< Index > &vertices );

	size_t vertexNum( void ) const { return VectorStreamingVertices< Vertex , Index >::vertexNum(); }
	size_t polygonNum( void ) const;
};

class BufferedReadWriteFile
{
	bool tempFile;
	FILE* _fp;
	Pointer( char ) _buffer;
	char _fileName[2048];
	size_t _bufferIndex , _bufferSize;
public:
	BufferedReadWriteFile( const char* fileName=NULL , const char* fileHeader="" , unsigned int bufferSize=(1<<20) )
	{
		_bufferIndex = 0;
		_bufferSize = bufferSize;
		if( fileName ) strcpy( _fileName , fileName ) , tempFile = false , _fp = fopen( _fileName , "w+b" );
		else
		{
			if( fileHeader && strlen(fileHeader) ) sprintf( _fileName , "%sXXXXXX" , fileHeader );
			else strcpy( _fileName , "XXXXXX" );
#ifdef _WIN32
			_mktemp( _fileName );
			_fp = fopen( _fileName , "w+b" );
#else // !_WIN32
			_fp = fdopen( mkstemp( _fileName ) , "w+b" );
#endif // _WIN32
			tempFile = true;
		}
		if( !_fp ) ERROR_OUT( "Failed to open file: " , _fileName );
		_buffer = AllocPointer< char >( _bufferSize );
	}
	~BufferedReadWriteFile( void )
	{
		FreePointer( _buffer );
		fclose( _fp );
		if( tempFile ) remove( _fileName );
	}
	bool write( ConstPointer( char ) data , size_t size )
	{
		if( !size ) return true;
		ConstPointer( char ) _data = data;
		size_t sz = _bufferSize - _bufferIndex;
		while( sz<=size )
		{
			memcpy( _buffer+_bufferIndex , _data , sz );
			fwrite( _buffer , 1 , _bufferSize , _fp );
			_data += sz;
			size -= sz;
			_bufferIndex = 0;
			sz = _bufferSize;
		}
		if( size )
		{
			memcpy( _buffer+_bufferIndex , _data , size );
			_bufferIndex += size;
		}
		return true;
	}
	bool read( Pointer( char ) data , size_t size )
	{
		if( !size ) return true;
		Pointer( char ) _data = data;
		size_t sz = _bufferSize - _bufferIndex;
		while( sz<=size )
		{
			if( size && !_bufferSize ) return false;
			memcpy( _data , _buffer+_bufferIndex , sz );
			_bufferSize = fread( _buffer , 1 , _bufferSize , _fp );
			_data += sz;
			size -= sz;
			_bufferIndex = 0;
			if( !size ) return true;
			sz = _bufferSize;
		}
		if( size )
		{
			if( !_bufferSize ) return false;
			memcpy( _data , _buffer+_bufferIndex , size );
			_bufferIndex += size;
		}
		return true;
	}
	void reset( void )
	{
		if( _bufferIndex ) fwrite( _buffer , 1 , _bufferIndex , _fp );
		_bufferIndex = 0;
		fseek( _fp , 0 , SEEK_SET );
		_bufferIndex = 0;
		_bufferSize = fread( _buffer , 1 , _bufferSize , _fp );
	}
};

template< typename VertexFactory , typename Index >
class FileStreamingVertices : public StreamingVertices< typename VertexFactory::VertexType , Index >
{
	BufferedReadWriteFile *_vertexFile;
	Index _vertexNum;
	const VertexFactory &_factory;
	Pointer( char ) _buffer;
public:
	using Vertex = typename VertexFactory::VertexType;
	FileStreamingVertices( const VertexFactory &vFactory , const char* fileHeader="" );
	~FileStreamingVertices( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v );

	bool nextVertex( Vertex &v );

	size_t vertexNum( void ) const;
};


template< typename VertexFactory , typename Index >
class FileStreamingCurve : public StreamingCurve< typename VertexFactory::VertexType , Index > , public FileStreamingVertices< VertexFactory , Index >
{
	std::vector< BufferedReadWriteFile* > _edgeFiles;
	unsigned int _threadIndex;
public:
	using Vertex = typename VertexFactory::VertexType;
	FileStreamingCurve( const char* fileHeader="" );
	~FileStreamingCurve( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v ){ return FileStreamingVertices< VertexFactory , Index >::addVertex( v ); }
	void addEdge_s( unsigned int thread , Index v1 , Index v2 );

	bool nextVertex( Vertex &v ){ return FileStreamingVertices< VertexFactory , Index >::nextVertex(v); }
	bool nextEdge( std::pair< Index , Index > &e );

	size_t vertexNum( void ) const { return FileStreamingVertices< VertexFactory , Index >::vertexNum(); }
	size_t edgeNum( void ) const;
};

template< typename VertexFactory , typename Index >
class FileStreamingMesh : public StreamingMesh< typename VertexFactory::VertexType , Index > , public FileStreamingVertices< VertexFactory , Index >
{
	std::vector< Index > _polygonNum;
	std::vector< BufferedReadWriteFile * > _polygonFiles;
	unsigned int _threadIndex;
public:
	using Vertex = typename VertexFactory::VertexType;
	FileStreamingMesh( const VertexFactory &vFactory , const char* fileHeader="" );
	~FileStreamingMesh( void );

	void resetIterator( void );

	Index addVertex( const Vertex &v ){ return FileStreamingVertices< VertexFactory , Index >::addVertex( v ); }
	void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices );

	bool nextVertex( Vertex &v ) { return FileStreamingVertices< VertexFactory , Index >::nextVertex(v); }
	bool nextPolygon( std::vector< Index >& vertices );

	size_t vertexNum( void ) const { return FileStreamingVertices< VertexFactory , Index >::vertexNum(); }
	size_t polygonNum( void ) const;
};

#include "StreamingMesh.inl"

#endif // STREAMING_MESH_INCLUDED
