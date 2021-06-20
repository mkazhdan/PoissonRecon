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

#ifndef CORED_MESH_INCLUDED
#define CORED_MESH_INCLUDED

#include "Geometry.h"
#include "MyMiscellany.h"

template< typename Index >
class CoredPointIndex
{
public:
	Index index;
	char inCore;

	bool operator == (const CoredPointIndex& cpi) const {return (index==cpi.index) && (inCore==cpi.inCore);};
	bool operator != (const CoredPointIndex& cpi) const {return (index!=cpi.index) || (inCore!=cpi.inCore);};
};

template< typename Index > struct CoredEdgeIndex{ CoredPointIndex< Index > idx[2]; };

template< typename Index >
struct CoredVertexIndex
{
	Index idx;
	bool inCore;
};

template< class Vertex , typename Index >
class CoredCurveData
{
public:
	std::vector< Vertex > inCoreVertices;
	virtual void resetIterator( void ) = 0;

	virtual Index addOutOfCoreVertex  (                       const Vertex& v ) = 0;
	virtual Index addOutOfCoreVertex_s( unsigned int thread , const Vertex& v ) = 0;
	virtual void addEdge_s( unsigned int thread , CoredVertexIndex< Index > v1 , CoredVertexIndex< Index > v2 ) = 0;
	virtual void addEdge_s( unsigned int thread , Index v1 , Index v2 ) = 0;

	virtual Index nextOutOfCoreVertex( Vertex &v )=0;
	virtual Index nextEdge( CoredVertexIndex< Index >& v1 , CoredVertexIndex< Index >& v2 ) = 0;

	virtual size_t outOfCoreVertexNum(void)=0;
	virtual size_t edgeCount( void ) = 0;
};

template< class Vertex , typename Index >
class CoredMeshData
{
public:
	virtual ~CoredMeshData( void ){}
	std::vector< Vertex > inCoreVertices;
	virtual void resetIterator( void ) = 0;

	virtual Index addOutOfCoreVertex  (                       const Vertex &v ) = 0;
	virtual Index addOutOfCoreVertex_s( unsigned int thread , const Vertex &v ) = 0;
	virtual void addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices ) = 0;
	virtual void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices ) = 0;

	virtual Index nextOutOfCoreVertex( Vertex &v )=0;
	virtual Index nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices ) = 0;

	virtual size_t outOfCoreVertexNum( void )=0;
	virtual size_t polygonNum( void ) = 0;
};

template< class Vertex , typename Index >
class CoredVectorCurveData : public CoredCurveData< Vertex , Index >
{
	std::vector< Vertex > oocPoints;
	std::vector< std::pair< Index , Index > > edges;
	unsigned int threadIndex;
	Index edgeIndex;
	Index oocPointIndex;
public:
	CoredVectorCurveData( void );

	void resetIterator( void );

	Index addOutOfCoreVertex  (                       const Vertex &v );
	Index addOutOfCoreVertex_s( unsigned int thread , const Vertex &v );
	void addEdge_s( unsigned int thread , CoredVertexIndex< Index > v1 , CoredVertexIndex< Index > v2 );
	void addEdge_s( unsigned int thread , Index v1 , Index v2 );

	Index nextOutOfCoreVertex( Vertex &v );
	Index nextEdge( CoredVertexIndex< Index > &v1 , CoredVertexIndex< Index > &v2 );

	size_t outOfCoreVertexNum( void );
	size_t edgeCount( void );
};
template< class Vertex , typename Index >
class CoredVectorMeshData : public CoredMeshData< Vertex , Index >
{
	std::vector< Vertex > oocPoints;
	std::vector< std::vector< std::vector< Index > > > polygons;
	unsigned int threadIndex;
	Index polygonIndex;
	Index oocPointIndex;
public:
	CoredVectorMeshData( void );

	void resetIterator( void );

	Index addOutOfCoreVertex  (                       const Vertex &v );
	Index addOutOfCoreVertex_s( unsigned int thread , const Vertex &v );
	void addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices );
	void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices );

	Index nextOutOfCoreVertex( Vertex &v );
	Index nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices );

	size_t outOfCoreVertexNum( void );
	size_t polygonNum( void );
};
class BufferedReadWriteFile
{
	bool tempFile;
	FILE* _fp;
	char *_buffer , _fileName[1024];
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
		_buffer = (char*) malloc( _bufferSize );
	}
	~BufferedReadWriteFile( void )
	{
		free( _buffer );
		fclose( _fp );
		if( tempFile ) remove( _fileName );
	}
	bool write( const void* data , size_t size )
	{
		if( !size ) return true;
		const char* _data = (char*) data;
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
	bool read( void* data , size_t size )
	{
		if( !size ) return true;
		char *_data = (char*) data;
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
template< class Vertex , typename Index >
class CoredFileCurveData : public CoredCurveData< Vertex , Index >
{
	BufferedReadWriteFile *oocPointFile;
	Index oocPoints;
	std::vector< BufferedReadWriteFile* > edgeFiles;
	unsigned int threadIndex;
public:
	CoredFileCurveData( const char* fileHeader="" );
	~CoredFileCurveData( void );

	void resetIterator( void );

	Index addOutOfCoreVertex  (                       const Vertex &v );
	Index addOutOfCoreVertex_s( unsigned int thread , const Vertex &v );
	void addEdge_s( unsigned int thread , CoredVertexIndex< Index > v1 , CoredVertexIndex< Index > v2 );
	void addEdge_s( unsigned int thread , Index v1 , Index v2 );

	Index nextOutOfCoreVertex( Vertex &v );
	Index nextEdge( CoredVertexIndex< Index > &v1 , CoredVertexIndex< Index > &v2 );

	size_t outOfCoreVertexNum( void );
	size_t edgeCount( void );
};

template< typename Index , typename VertexFactory >
class CoredFileMeshData : public CoredMeshData< typename VertexFactory::VertexType , Index >
{
	typedef typename VertexFactory::VertexType VertexType;
	BufferedReadWriteFile *_oocVertexFile;
	Index _oocVertexNum;
	std::vector< Index > _polygonNum;
	std::vector< BufferedReadWriteFile * > _polygonFiles;
	unsigned int _threadIndex;
	const VertexFactory &_factory;
	char *_buffer;
public:
	CoredFileMeshData( const VertexFactory &vFactory , const char* fileHeader="" );
	~CoredFileMeshData( void );

	void resetIterator( void );

	Index addOutOfCoreVertex  (                       const VertexType& v );
	Index addOutOfCoreVertex_s( unsigned int thread , const VertexType& v );
	void addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices );
	void addPolygon_s( unsigned int thread , const std::vector< Index >& vertices );

	Index nextOutOfCoreVertex( VertexType &v );
	Index nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices );

	size_t outOfCoreVertexNum( void );
	size_t polygonNum( void );
};

#include "CoredMesh.inl"

#endif // CORED_MESH
