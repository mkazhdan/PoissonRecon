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


/////////////////////////
// CoredVectorMeshData //
/////////////////////////
template< class Vertex , typename Index >
CoredVectorMeshData< Vertex , Index >::CoredVectorMeshData( void ) { oocPointIndex = polygonIndex = threadIndex = 0 ; polygons.resize( std::thread::hardware_concurrency() ); }
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::resetIterator ( void ) { oocPointIndex = polygonIndex = threadIndex = 0; }
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::addOutOfCoreVertex( const Vertex& p )
{
	oocPoints.push_back(p);
	return ( Index )oocPoints.size()-1;
}
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::addOutOfCoreVertex_s( unsigned int thread , const Vertex& p )
{
	size_t sz;
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
		sz = oocPoints.size();
		oocPoints.push_back(p);
	}
	return (Index)sz;
}
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< Index >& polygon )
{
	polygons[ thread ].push_back( polygon );
}
template< class Vertex , typename Index >
void CoredVectorMeshData< Vertex , Index >::addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices )
{
	std::vector< Index > polygon( vertices.size() );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
	return addPolygon_s( thread , polygon );
}
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::nextOutOfCoreVertex( Vertex &p )
{
	if( oocPointIndex<(Index)oocPoints.size() )
	{
		p=oocPoints[oocPointIndex++];
		return 1;
	}
	else return 0;
}
template< class Vertex , typename Index >
Index CoredVectorMeshData< Vertex , Index >::nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices )
{
	while( true )
	{
		if( threadIndex<(int)polygons.size() )
		{
			if( polygonIndex<(Index)( polygons[threadIndex].size() ) )
			{
				std::vector< Index >& polygon = polygons[threadIndex][ polygonIndex++ ];
				vertices.resize( polygon.size() );
				for( int i=0 ; i<int(polygon.size()) ; i++ )
					if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
					else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
				return 1;
			}
			else threadIndex++ , polygonIndex = 0;
		}
		else return 0;
	}
}
template< class Vertex , typename Index >
size_t CoredVectorMeshData< Vertex , Index >::outOfCoreVertexNum( void ){ return oocPoints.size(); }
template< class Vertex , typename Index >
size_t CoredVectorMeshData< Vertex , Index >::polygonNum( void )
{
	size_t count = 0;
	for( size_t i=0 ; i<polygons.size() ; i++ ) count += polygons[i].size();
	return count;
}

///////////////////////
// CoredFileMeshData //
///////////////////////
template< typename Index , typename VertexFactory >
CoredFileMeshData< Index , VertexFactory >::CoredFileMeshData( const VertexFactory &factory , const char* fileHeader ) : _factory(factory)
{
	_threadIndex = 0;
	_oocVertexNum = 0;
	_polygonNum.resize( std::thread::hardware_concurrency() );
	for( unsigned int i=0 ; i<_polygonNum.size() ; i++ ) _polygonNum[i] = 0;

	char _fileHeader[1024];
	sprintf( _fileHeader , "%s_points_" , fileHeader );
	_oocVertexFile = new BufferedReadWriteFile( NULL , _fileHeader );
	_polygonFiles.resize( std::thread::hardware_concurrency() );
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ )
	{
		sprintf( _fileHeader , "%s_polygons_t%d_" , fileHeader , i );
		_polygonFiles[i] = new BufferedReadWriteFile( NULL , _fileHeader );
	}

	_buffer = new char[ factory.bufferSize() ];
}

template< typename Index , typename VertexFactory >
CoredFileMeshData< Index , VertexFactory >::~CoredFileMeshData( void )
{
	delete _oocVertexFile;
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ ) delete _polygonFiles[i];
	delete[] _buffer;
}

template< typename Index , typename VertexFactory >
void CoredFileMeshData< Index , VertexFactory >::resetIterator ( void )
{
	_oocVertexFile->reset();
	_threadIndex = 0;
	for( unsigned int i=0 ; i<_polygonFiles.size() ; i++ ) _polygonFiles[i]->reset();
}

template< typename Index , typename VertexFactory >
Index CoredFileMeshData< Index , VertexFactory >::addOutOfCoreVertex( const VertexType& v )
{
	_factory.toBuffer( v , _buffer );
	_oocVertexFile->write( _buffer , _factory.bufferSize() );
	return _oocVertexNum++;
}

template< typename Index , typename VertexFactory >
Index CoredFileMeshData< Index , VertexFactory >::addOutOfCoreVertex_s( unsigned int thread , const VertexType& v )
{
	Index sz;
	{
		static std::mutex m;
		std::lock_guard< std::mutex > lock(m);
		sz = _oocVertexNum;
		_factory.toBuffer( v , _buffer );
		_oocVertexFile->write( _buffer , _factory.bufferSize() );
		_oocVertexNum++;
	}
	return sz;
}

template< typename Index , typename VertexFactory >
void CoredFileMeshData< Index , VertexFactory >::addPolygon_s( unsigned int thread , const std::vector< Index >& vertices )
{
	unsigned int vSize = (unsigned int)vertices.size();
	_polygonFiles[thread]->write( &vSize , sizeof(unsigned int) );
	_polygonFiles[thread]->write( &vertices[0] , sizeof(Index) * vSize );
	_polygonNum[thread]++;
}

template< typename Index , typename VertexFactory >
void CoredFileMeshData< Index , VertexFactory >::addPolygon_s( unsigned int thread , const std::vector< CoredVertexIndex< Index > >& vertices )
{
	std::vector< Index > polygon( vertices.size() );
	for( unsigned int i=0 ; i<(unsigned int)vertices.size() ; i++ ) 
		if( vertices[i].inCore ) polygon[i] =  vertices[i].idx;
		else                     polygon[i] = -vertices[i].idx-1;
	return addPolygon_s( thread , polygon );
}

template< typename Index , typename VertexFactory >
Index CoredFileMeshData< Index , VertexFactory >::nextOutOfCoreVertex( VertexType& v )
{
	if( _oocVertexFile->read( _buffer , _factory.bufferSize() ) )
	{
		_factory.fromBuffer( _buffer , v );
		return 1;
	}
	else return 0;
}

template< typename Index , typename VertexFactory >
Index CoredFileMeshData< Index , VertexFactory >::nextPolygon( std::vector< CoredVertexIndex< Index > >& vertices )
{
	while( true )
	{
		if( _threadIndex<(unsigned int)_polygonFiles.size() )
		{
			unsigned int pSize;
			if( _polygonFiles[ _threadIndex ]->read( &pSize , sizeof(unsigned int) ) )
			{
				std::vector< Index > polygon( pSize );
				if( _polygonFiles[ _threadIndex ]->read( &polygon[0] , sizeof(Index)*pSize ) )
				{
					vertices.resize( pSize );
					for( unsigned int i=0 ; i<(unsigned int)polygon.size() ; i++ )
						if( polygon[i]<0 ) vertices[i].idx = -polygon[i]-1 , vertices[i].inCore = false;
						else               vertices[i].idx =  polygon[i]   , vertices[i].inCore = true;
					return 1;
				}
				ERROR_OUT( "Failed to read polygon from file" );
			}
			else _threadIndex++;
		}
		else return 0;
	}
}

template< typename Index , typename VertexFactory >
size_t CoredFileMeshData< Index , VertexFactory >::outOfCoreVertexNum( void ){ return _oocVertexNum; }

template< typename Index , typename VertexFactory >
size_t CoredFileMeshData< Index , VertexFactory >::polygonNum( void )
{
	size_t count = 0;
	for( size_t i=0 ; i<_polygonNum.size() ; i++ ) count += _polygonNum[i];
	return count;
}