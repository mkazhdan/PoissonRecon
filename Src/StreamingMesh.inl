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

/////////////////////////////
// VectorStreamingVertices //
/////////////////////////////
template< typename Vertex , typename Index >
VectorStreamingVertices< Vertex , Index >::VectorStreamingVertices( void ) { _vertexIndex = 0; }

template< typename Vertex , typename Index >
void VectorStreamingVertices< Vertex , Index >::resetIterator ( void ) { _vertexIndex = 0; }

template< typename Vertex , typename Index >
Index VectorStreamingVertices< Vertex , Index >::addVertex( const Vertex& p )
{
	_vertices.push_back(p);
	return ( Index )_vertices.size()-1;
}

template< typename Vertex , typename Index >
bool VectorStreamingVertices< Vertex , Index >::nextVertex( Vertex &p )
{
	if( _vertexIndex<(Index)_vertices.size() )
	{
		p = _vertices[ _vertexIndex++ ];
		return true;
	}
	else return false;
}

template< typename Vertex , typename Index >
size_t VectorStreamingVertices< Vertex , Index >::vertexNum( void ) const { return _vertices.size(); }

//////////////////////////
// VectorStreamingCurve //
//////////////////////////
template< class Vertex , typename Index >
VectorStreamingCurve< Vertex , Index >::VectorStreamingCurve( void ) { _threadIndex = _edgeIndex = 0 ; _edges.resize( std::thread::hardware_concurrency() ); }

template< class Vertex , typename Index >
void VectorStreamingCurve< Vertex , Index >::resetIterator ( void ) { _threadIndex = _edgeIndex = 0 ; VectorStreamingVertices< Vertex , Index >::resetIterator(); }


template< class Vertex , typename Index >
void VectorStreamingCurve< Vertex , Index >::addEdge_s( unsigned int thread , Index v1 , Index v2 )
{
	_edges[ thread ].push_back( std::make_pair( v1 , v2 ) );
}

template< class Vertex , typename Index >
bool VectorStreamingCurve< Vertex , Index >::nextEdge( std::pair< Index , Index > &e )
{
	while( true )
	{
		if( _threadIndex<(int)_edges.size() )
		{
			if( _edgeIndex<(Index)( _edges[_threadIndex].size() ) )
			{
				e = _edges[_threadIndex][ _edgeIndex++ ];
				return true;
			}
			else _threadIndex++ , _edgeIndex = 0;
		}
		else return false;
	}
}

template< class Vertex , typename Index >
size_t VectorStreamingCurve< Vertex , Index >::edgeNum( void ) const
{
	size_t count = 0;
	for( size_t i=0 ; i<_edges.size() ; i++ ) count += _edges[i].size();
	return count;
}

/////////////////////////
// VectorStreamingMesh //
/////////////////////////
template< class Vertex , typename Index >
VectorStreamingMesh< Vertex , Index >::VectorStreamingMesh( void ) { _threadIndex = 0 ; _polygonIndex = 0 ; _polygons.resize( std::thread::hardware_concurrency() ); }

template< class Vertex , typename Index >
void VectorStreamingMesh< Vertex , Index >::resetIterator ( void ){ _threadIndex = 0 ; _polygonIndex = 0 ; VectorStreamingVertices< Vertex , Index >::resetIterator(); }

template< class Vertex , typename Index >
void VectorStreamingMesh< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< Index >& polygon )
{
	_polygons[ thread ].push_back( polygon );
}

template< class Vertex , typename Index >
bool VectorStreamingMesh< Vertex , Index >::nextPolygon( std::vector< Index >& vertices )
{
	while( true )
	{
		if( _threadIndex<(int)_polygons.size() )
		{
			if( _polygonIndex<(Index)( _polygons[_threadIndex].size() ) )
			{
				vertices = _polygons[ _threadIndex ][ _polygonIndex++ ];
				return true;
			}
			else _threadIndex++ , _polygonIndex = 0;
		}
		else return false;
	}
}

template< class Vertex , typename Index >
size_t VectorStreamingMesh< Vertex , Index >::polygonNum( void ) const
{
	size_t count = 0;
	for( size_t i=0 ; i<_polygons.size() ; i++ ) count += _polygons[i].size();
	return count;
}

///////////////////////
// FileStreamingMesh //
///////////////////////
template< typename VertexFactory , typename Index >
FileStreamingVertices< VertexFactory , Index >::FileStreamingVertices( const VertexFactory &factory , const char* fileHeader ) : _factory(factory)
{
	_vertexNum = 0;

	char _fileHeader[1024];
	sprintf( _fileHeader , "%s_points_" , fileHeader );
	_vertexFile = new BufferedReadWriteFile( NULL , _fileHeader );
	_buffer = NewPointer< char >( factory.bufferSize() );
}

template< typename VertexFactory , typename Index >
FileStreamingVertices< VertexFactory , Index >::~FileStreamingVertices( void )
{
	delete _vertexFile;
	DeletePointer( _buffer );
}

template< typename VertexFactory , typename Index >
void FileStreamingVertices< VertexFactory , Index >::resetIterator( void )
{
	_vertexFile->reset();
#if 1
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Why are we not resetting _vertexNum to zero?" )
#endif // SHOW_WARNINGS
#else
	_vertexNum = 0;
#endif
}

template< typename VertexFactory , typename Index >
Index FileStreamingVertices< VertexFactory , Index >::addVertex( const Vertex &v )
{
	_factory.toBuffer( v , _buffer );
	_vertexFile->write( _buffer , _factory.bufferSize() );
	return _vertexNum++;
}

template< typename VertexFactory , typename Index >
bool FileStreamingVertices< VertexFactory , Index >::nextVertex( Vertex &v )
{
	if( _vertexFile->read( _buffer , _factory.bufferSize() ) )
	{
		_factory.fromBuffer( _buffer , v );
		return true;
	}
	else return false;
}

template< typename VertexFactory , typename Index >
size_t FileStreamingVertices< VertexFactory , Index >::vertexNum( void ) const { return _vertexNum; }

///////////////////////
// FileStreamingMesh //
///////////////////////
template< typename VertexFactory , typename Index >
FileStreamingMesh< VertexFactory , Index >::FileStreamingMesh( const VertexFactory &factory , const char* fileHeader ) : FileStreamingVertices< VertexFactory , Index >( factory , fileHeader )
{
	_threadIndex = 0;
	_polygonNum.resize( std::thread::hardware_concurrency() );
	for( unsigned int i=0 ; i<_polygonNum.size() ; i++ ) _polygonNum[i] = 0;

	char _fileHeader[1024];
	_polygonFiles.resize( std::thread::hardware_concurrency() );
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ )
	{
		sprintf( _fileHeader , "%s_polygons_t%d_" , fileHeader , i );
		_polygonFiles[i] = new BufferedReadWriteFile( NULL , _fileHeader );
	}
}

template< typename VertexFactory , typename Index >
FileStreamingMesh< VertexFactory , Index >::~FileStreamingMesh( void )
{
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ ) delete _polygonFiles[i];
}

template< typename VertexFactory , typename Index >
void FileStreamingMesh< VertexFactory , Index >::resetIterator ( void )
{
	_threadIndex = 0;
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ ) _polygonFiles[i]->reset();
	FileStreamingVertices< VertexFactory , Index >::resetIterator();
}

template< typename VertexFactory , typename Index >
void FileStreamingMesh< VertexFactory , Index >::addPolygon_s( unsigned int thread , const std::vector< Index >& vertices )
{
	unsigned int vSize = (unsigned int)vertices.size();
	_polygonFiles[thread]->write( ( ConstPointer( char ) )GetPointer( vSize ), sizeof(unsigned int) );
	_polygonFiles[thread]->write( ( ConstPointer( char ) )GetPointer( vertices ) , sizeof(Index) * vSize );
	_polygonNum[thread]++;
}

template< typename VertexFactory , typename Index >
bool FileStreamingMesh< VertexFactory , Index >::nextPolygon( std::vector< Index > &vertices )
{
	while( true )
	{
		if( _threadIndex<(unsigned int)_polygonFiles.size() )
		{
			unsigned int pSize;
			if( _polygonFiles[ _threadIndex ]->read( ( Pointer( char ) )GetPointer( pSize ) , sizeof(unsigned int) ) )
			{
				vertices.resize( pSize );
				if( _polygonFiles[ _threadIndex ]->read( ( Pointer( char ) )GetPointer( vertices ) , sizeof(Index)*pSize ) ) return true;
				ERROR_OUT( "Failed to read polygon from file" );
			}
			else _threadIndex++;
		}
		else return false;
	}
}

template< typename VertexFactory , typename Index >
size_t FileStreamingMesh< VertexFactory , Index >::polygonNum( void ) const
{
	size_t count = 0;
	for( size_t i=0 ; i<_polygonNum.size() ; i++ ) count += _polygonNum[i];
	return count;
}
