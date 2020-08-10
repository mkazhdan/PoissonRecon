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

///////////////////////////
// MemoryInputDataStream //
///////////////////////////
template< typename Data >
MemoryInputDataStream< Data >::MemoryInputDataStream( size_t size , const Data *data ) : _data(data) , _size(size) , _current(0) {}
template< typename Data >
MemoryInputDataStream< Data >::~MemoryInputDataStream( void ){ ; }
template< typename Data >
void MemoryInputDataStream< Data >::reset( void ) { _current=0; }
template< typename Data >
bool MemoryInputDataStream< Data >::next( Data &d )
{
	if( _current>=_size ) return false;
	d = _data[_current++];
	return true;
}

////////////////////////////
// MemoryOutputDataStream //
////////////////////////////
template< typename Data >
MemoryOutputDataStream< Data >::MemoryOutputDataStream( size_t size , Data *data ) : _data(data) , _size(size) , _current(0) {}
template< typename Data >
MemoryOutputDataStream< Data >::~MemoryOutputDataStream( void ){ ; }
template< typename Data >
void MemoryOutputDataStream< Data >::next( const Data &d )
{
	_data[_current++] = d;
}

//////////////////////////
// ASCIIInputDataStream //
//////////////////////////
template< typename Factory >
ASCIIInputDataStream< Factory >::ASCIIInputDataStream( const char* fileName , const Factory &factory ) : _factory( factory )
{
	_fp = fopen( fileName , "r" );
	if( !_fp ) ERROR_OUT( "Failed to open file for reading: %s" , fileName );
}
template< typename Factory >
ASCIIInputDataStream< Factory >::~ASCIIInputDataStream( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< typename Factory >
void ASCIIInputDataStream< Factory >::reset( void ) { fseek( _fp , SEEK_SET , 0 ); }
template< typename Factory >
bool ASCIIInputDataStream< Factory >::next( Data &d )
{
	d = _factory();
	return _factory.readASCII( _fp , d );
}

///////////////////////////
// ASCIIOutputDataStream //
///////////////////////////
template< typename Factory >
ASCIIOutputDataStream< Factory >::ASCIIOutputDataStream( const char* fileName , const Factory &factory ) : _factory( factory )
{
	_fp = fopen( fileName , "w" );
	if( !_fp ) ERROR_OUT( "Failed to open file for writing: %s" , fileName );
}
template< typename Factory >
ASCIIOutputDataStream< Factory >::~ASCIIOutputDataStream( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< typename Factory >
void ASCIIOutputDataStream< Factory >::next( const Data &d )
{
	_factory.writeASCII( _fp , d );
}

///////////////////////////
// BinaryInputDataStream //
///////////////////////////
template< typename Factory >
BinaryInputDataStream< Factory >::BinaryInputDataStream( const char* fileName , const Factory &factory ) : _factory(factory)
{
	_fp = fopen( fileName , "rb" );
	if( !_fp ) ERROR_OUT( "Failed to open file for reading: %s" , fileName );
}
template< typename Factory >
void BinaryInputDataStream< Factory >::reset( void ) { fseek( _fp , SEEK_SET , 0 ); }
template< typename Factory >
bool BinaryInputDataStream< Factory >::next( Data &d )
{
	d = _factory();
	return _factory.readBinary( _fp , d );
}

////////////////////////////
// BinaryOutputDataStream //
////////////////////////////
template< typename Factory >
BinaryOutputDataStream< Factory >::BinaryOutputDataStream( const char* fileName , const Factory &factory ) : _factory(factory)
{
	_fp = fopen( fileName , "wb" );
	if( !_fp ) ERROR_OUT( "Failed to open file for writing: %s" , fileName );
}
template< typename Factory >
void BinaryOutputDataStream< Factory >::next( const Data &d ){ return _factory.writeBinary( _fp , d ); }

////////////////////////
// PLYInputDataStream //
////////////////////////
template< typename Factory >
PLYInputDataStream< Factory >::PLYInputDataStream( const char* fileName , const Factory &factory ) : _factory(factory)
{
	_fileName = new char[ strlen( fileName )+1 ];
	strcpy( _fileName , fileName );
	_ply = NULL;
	if( factory.bufferSize() ) _buffer = new char[ factory.bufferSize() ];
	else _buffer = NULL;
	reset();
}
template< typename Factory >
void PLYInputDataStream< Factory >::reset( void )
{
	int fileType;
	float version;
	std::vector< PlyProperty > plist;
	if( _ply ) _free();
	_ply = PlyFile::Read( _fileName, _elist, fileType, version );
	if( !_ply ) ERROR_OUT( "Failed to open ply file for reading: %s" , _fileName );

	bool foundData = false;
	for( int i=0 ; i<_elist.size() ; i++ )
	{
		size_t num_elems;
		std::string &elem_name = _elist[i];
		plist = _ply->get_element_description( elem_name , num_elems );
		if( !plist.size() ) ERROR_OUT( "Failed to get element description: %s" , elem_name );

		if( elem_name=="vertex" )
		{
			foundData = true;
			_pCount = num_elems , _pIdx = 0;

			bool* properties = new bool[ _factory.plyReadNum() ];
			for( unsigned int i=0 ; i<_factory.plyReadNum() ; i++ )
			{
				PlyProperty prop = _factory.isStaticallyAllocated() ? _factory.plyStaticReadProperty(i) : _factory.plyReadProperty(i);
				if( !_ply->get_property( elem_name , &prop ) ) properties[i] = false;
				else                                           properties[i] = true;
			}
			bool valid = _factory.plyValidReadProperties( properties );
			delete[] properties;
			if( !valid ) ERROR_OUT( "Failed to validate properties in file" );
		}
	}
	if( !foundData ) ERROR_OUT( "Could not find data in ply file" );
}
template< typename Factory >
void PLYInputDataStream< Factory >::_free( void ){ delete _ply; }

template< typename Factory >
PLYInputDataStream< Factory >::~PLYInputDataStream( void )
{
	_free();
	if( _fileName ) delete[] _fileName , _fileName = NULL;
	if( _buffer ) delete[] _buffer , _buffer = NULL;
}
template< typename Factory >
bool PLYInputDataStream< Factory >::next( Data &d )
{
	if( _pIdx<_pCount )
	{
		d = _factory();
		if( _factory.isStaticallyAllocated() ) _ply->get_element( (void *)&d );
		else
		{
			_ply->get_element( (void *)_buffer );
			_factory.fromBuffer( _buffer , d );
		}
		_pIdx++;
		return true;
	}
	else return false;
}

/////////////////////////
// PLYOutputDataStream //
/////////////////////////
template< typename Factory >
PLYOutputDataStream< Factory >::PLYOutputDataStream( const char* fileName , const Factory &factory , size_t count , int fileType ) : _factory(factory)
{
	float version;
	std::vector< std::string > elem_names = { std::string( "vertex" ) };
	_ply = PlyFile::Write( fileName , elem_names , fileType , version );
	if( !_ply ) ERROR_OUT( "Failed to open ply file for writing: " , fileName );

	_pIdx = 0;
	_pCount = count;
	_ply->element_count( "vertex" , _pCount );
	for( unsigned int i=0 ; i<_factory.plyWriteNum() ; i++ )
	{
		PlyProperty prop = _factory.isStaticallyAllocated() ? _factory.plyStaticWriteProperty(i) : _factory.plyWriteProperty(i);
		_ply->describe_property( "vertex" , &prop );
	}
	_ply->header_complete();
	_ply->put_element_setup( "vertex" );
	if( _factory.bufferSize() ) _buffer = new char[ _factory.bufferSize() ];
	else                        _buffer = NULL;
}
template< typename Factory >
PLYOutputDataStream< Factory >::~PLYOutputDataStream( void )
{
	if( _pIdx!=_pCount ) ERROR_OUT( "Streamed points not equal to total count: " , _pIdx , " != " , _pCount );
	delete _ply;
	if( _buffer ) delete[] _buffer , _buffer = NULL;
}
template< typename Factory >
void PLYOutputDataStream< Factory >::next( const Data &d )
{
	if( _pIdx==_pCount ) ERROR_OUT( "Trying to add more points than total: " , _pIdx , " < " , _pCount );
	if( _factory.isStaticallyAllocated() ) _ply->put_element( (void *)&d );
	else
	{
		_factory.toBuffer( d , _buffer );
		_ply->put_element( (void *)_buffer );
	}
	_pIdx++;
}

