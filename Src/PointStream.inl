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

////////////////////////////
// MemoryInputPointStream //
////////////////////////////
template< class Real , int Dim >
MemoryInputPointStream< Real , Dim >::MemoryInputPointStream( size_t pointCount , const Point< Real , Dim >* points ){ _points = points , _pointCount = pointCount , _current = 0; }
template< class Real , int Dim >
MemoryInputPointStream< Real , Dim >::~MemoryInputPointStream( void ){ ; }
template< class Real , int Dim >
void MemoryInputPointStream< Real , Dim >::reset( void ) { _current=0; }
template< class Real , int Dim >
bool MemoryInputPointStream< Real , Dim >::nextPoint( Point< Real , Dim >& p )
{
	if( _current>=_pointCount ) return false;
	p = _points[_current];
	_current++;
	return true;
}

///////////////////////////
// ASCIIInputPointStream //
///////////////////////////
template< class Real , int Dim >
ASCIIInputPointStream< Real , Dim >::ASCIIInputPointStream( const char* fileName )
{
	_fp = fopen( fileName , "r" );
	if( !_fp ) fprintf( stderr , "Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
}
template< class Real , int Dim >
ASCIIInputPointStream< Real , Dim >::~ASCIIInputPointStream( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< class Real , int Dim >
void ASCIIInputPointStream< Real , Dim >::reset( void ) { fseek( _fp , SEEK_SET , 0 ); }
template< class Real , int Dim >
bool ASCIIInputPointStream< Real , Dim >::nextPoint( Point< Real , Dim >& p )
{
	float c;
	for( int d=0 ; d<Dim ; d++ )
		if( fscanf( _fp , " %f " , &c )!=1 ) return false;
		else p[d] = (Real)c;
	return true;
}

////////////////////////////
// ASCIIOutputPointStream //
////////////////////////////
template< class Real , int Dim >
ASCIIOutputPointStream< Real , Dim >::ASCIIOutputPointStream( const char* fileName )
{
	_fp = fopen( fileName , "w" );
	if( !_fp ) fprintf( stderr , "Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
}
template< class Real , int Dim >
ASCIIOutputPointStream< Real , Dim >::~ASCIIOutputPointStream( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< class Real , int Dim >
void ASCIIOutputPointStream< Real , Dim >::nextPoint( const Point< Real , Dim >& p )
{
	for( int d=0 ; d<Dim ; d++ ) fprintf( _fp , " %f" , (float)p[d] ); 
	fprintf( _fp , "\n" );
}

////////////////////////////
// BinaryInputPointStream //
////////////////////////////
template< class Real , int Dim >
BinaryInputPointStream< Real , Dim >::BinaryInputPointStream( const char* fileName , bool (*readPoint)( FILE* , Point< Real , Dim >& ) )
{
	_readPoint = readPoint!=NULL ? readPoint : _DefaultReadPoint;
	_fp = fopen( fileName , "rb" );
	if( !_fp ) fprintf( stderr , "Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
}

/////////////////////////////
// BinaryOutputPointStream //
/////////////////////////////
template< class Real , int Dim >
BinaryOutputPointStream< Real , Dim >::BinaryOutputPointStream( const char* fileName , void (*writePoint)( FILE* , const Point< Real , Dim >& ) )
{
	_writePoint = writePoint!=NULL ? writePoint : _DefaultWritePoint;
	_fp = fopen( fileName , "wb" );
	if( !_fp ) fprintf( stderr , "Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
}

/////////////////////////
// PLYInputPointStream //
/////////////////////////
template< class Real , int Dim >
PLYInputPointStream< Real , Dim >::PLYInputPointStream( const char* fileName )
{
	_fileName = new char[ strlen( fileName )+1 ];
	strcpy( _fileName , fileName );
	_ply = NULL;
	reset();
}
template< class Real , int Dim >
void PLYInputPointStream< Real , Dim >::reset( void )
{
	int fileType;
	float version;
	PlyProperty** plist;
	if( _ply ) _free();
	_ply = ply_open_for_reading( _fileName, &_nr_elems, &_elist, &fileType, &version );
	if( !_ply )
	{
		fprintf( stderr, "[ERROR] Failed to open ply file for reading: %s\n" , _fileName );
		exit( 0 );
	}
	bool foundVertices = false;
	for( int i=0 ; i<_nr_elems ; i++ )
	{
		int num_elems;
		int nr_props;
		char* elem_name = _elist[i];
		plist = ply_get_element_description( _ply , elem_name , &num_elems , &nr_props );
		if( !plist )
		{
			fprintf( stderr , "[ERROR] Failed to get element description: %s\n" , elem_name );
			exit( 0 );
		}	

		if( equal_strings( "vertex" , elem_name ) )
		{
			foundVertices = true;
			_pCount = num_elems , _pIdx = 0;
			for( int i=0 ; i<PlyVertex< Real , Dim >::ReadComponents ; i++ ) 
				if( !ply_get_property( _ply , elem_name , &(PlyVertex< Real , Dim >::Properties()[i]) ) )
				{
					fprintf( stderr , "[ERROR] Failed to find property in ply file: %s\n" , PlyVertex< Real , Dim >::Properties()[i].name );
					exit( 0 );
				}
		}
		for( int j=0 ; j<nr_props ; j++ )
		{
			free( plist[j]->name );
			free( plist[j] );
		}
		free( plist );
		if( foundVertices ) break;
	}
	if( !foundVertices )
	{
		fprintf( stderr , "[ERROR] Could not find vertices in ply file\n" );
		exit( 0 );
	}
}
template< class Real , int Dim >
void PLYInputPointStream< Real , Dim >::_free( void )
{
	if( _ply ) ply_close( _ply ) , _ply = NULL;
	if( _elist )
	{
		for( int i=0 ; i<_nr_elems ; i++ ) free( _elist[i] );
		free( _elist );
	}
}
template< class Real , int Dim >
PLYInputPointStream< Real , Dim >::~PLYInputPointStream( void )
{
	_free();
	if( _fileName ) delete[] _fileName , _fileName = NULL;
}
template< class Real , int Dim >
bool PLYInputPointStream< Real , Dim >::nextPoint( Point< Real , Dim >& p )
{
	if( _pIdx<_pCount )
	{
		PlyVertex< Real , Dim > v;
		ply_get_element( _ply, (void *)&v );
		p = v.point;
		_pIdx++;
		return true;
	}
	else return false;
}

////////////////////////////////////
// MemoryInputPointStreamWithData //
////////////////////////////////////
template< class Real , int Dim , class Data >
MemoryInputPointStreamWithData< Real , Dim , Data >::MemoryInputPointStreamWithData( size_t pointCount , const std::pair< Point< Real , Dim > , Data >* points ){ _points = points , _pointCount = pointCount , _current = 0; }
template< class Real , int Dim , class Data >
MemoryInputPointStreamWithData< Real , Dim , Data >::~MemoryInputPointStreamWithData( void ){ ; }
template< class Real , int Dim , class Data >
void MemoryInputPointStreamWithData< Real , Dim , Data >::reset( void ) { _current=0; }
template< class Real , int Dim , class Data >
bool MemoryInputPointStreamWithData< Real , Dim , Data >::nextPoint( Point< Real , Dim >& p , Data& d )
{
	if( _current>=_pointCount ) return false;
	p = _points[_current].first;
	d = _points[_current].second;
	_current++;
	return true;
}

///////////////////////////////////
// ASCIIInputPointStreamWithData //
///////////////////////////////////
template< class Real , int Dim , class Data >
ASCIIInputPointStreamWithData< Real , Dim , Data >::ASCIIInputPointStreamWithData( const char* fileName , Data (*readData)( FILE* ) ) : _readData( readData )
{
	_fp = fopen( fileName , "r" );
	if( !_fp ) fprintf( stderr , "Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
}
template< class Real , int Dim , class Data >
ASCIIInputPointStreamWithData< Real , Dim , Data >::~ASCIIInputPointStreamWithData( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< class Real , int Dim , class Data >
void ASCIIInputPointStreamWithData< Real , Dim , Data >::reset( void ) { fseek( _fp , SEEK_SET , 0 ); }
template< class Real , int Dim , class Data >
bool ASCIIInputPointStreamWithData< Real , Dim , Data >::nextPoint( Point< Real , Dim >& p , Data& d )
{
	float c;
	for( int dd=0 ; dd<Dim ; dd++ ) 
		if( fscanf( _fp , " %f " , &c )!=1 ) return false;
		else p[dd] = c;
	d = _readData( _fp );
	return true;
}

////////////////////////////////////
// ASCIIOutputPointStreamWithData //
////////////////////////////////////
template< class Real , int Dim , class Data >
ASCIIOutputPointStreamWithData< Real , Dim , Data >::ASCIIOutputPointStreamWithData( const char* fileName , void (*writeData)( FILE* , const Data& ) ) : _writeData( writeData )
{
	_fp = fopen( fileName , "w" );
	if( !_fp ) fprintf( stderr , "Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
}
template< class Real , int Dim , class Data >
ASCIIOutputPointStreamWithData< Real , Dim , Data >::~ASCIIOutputPointStreamWithData( void )
{
	fclose( _fp );
	_fp = NULL;
}
template< class Real , int Dim , class Data >
void ASCIIOutputPointStreamWithData< Real , Dim , Data >::nextPoint( const Point< Real , Dim >& p , const Data& d )
{
	for( int d=0 ; d<Dim ; d++ )  fprintf( _fp , " %f" , (float)p[d] );
	fprintf( _fp , " " );
	_writeData( _fp , d );
	fprintf( _fp , "\n" );
}

////////////////////////////////////
// BinaryInputPointStreamWithData //
////////////////////////////////////
template< class Real , int Dim , class Data >
BinaryInputPointStreamWithData< Real , Dim , Data >::BinaryInputPointStreamWithData( const char* fileName , bool (*readPointAndData)( FILE* , Point< Real , Dim >& , Data& ) )
{
	_readPointAndData = readPointAndData!=NULL ? readPointAndData : _DefaultReadPointAndData;
	_fp = fopen( fileName , "rb" );
	if( !_fp ) fprintf( stderr , "Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
}

/////////////////////////////////////
// BinaryOutputPointStreamWithData //
/////////////////////////////////////
template< class Real , int Dim , class Data >
BinaryOutputPointStreamWithData< Real , Dim , Data >::BinaryOutputPointStreamWithData( const char* fileName , void (*writePointAndData)( FILE* , const Point< Real , Dim >& , const Data& ) )
{
	_writePointAndData = writePointAndData!=NULL ? writePointAndData : _DefaultWritePointAndData;
	_fp = fopen( fileName , "wb" );
	if( !_fp ) fprintf( stderr , "Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
}

/////////////////////////////////
// PLYInputPointStreamWithData //
/////////////////////////////////
template< class Real , int Dim , class Data >
PLYInputPointStreamWithData< Real , Dim , Data >::PLYInputPointStreamWithData( const char* fileName , const PlyProperty* dataProperties , int dataPropertiesCount , bool (*validationFunction)( const bool* ) ) : _dataPropertiesCount( dataPropertiesCount ) , _validationFunction( validationFunction )
{
	_dataProperties = new PlyProperty[ _dataPropertiesCount ];
	memcpy( _dataProperties , dataProperties , sizeof(PlyProperty) * _dataPropertiesCount );
	for( int i=0 ; i<_dataPropertiesCount ; i++ ) _dataProperties[i].offset += sizeof( PlyVertex< Real , Dim > );
	_fileName = new char[ strlen( fileName )+1 ];
	strcpy( _fileName , fileName );
	_ply = NULL;
	reset();
}
template< class Real , int Dim , class Data >
void PLYInputPointStreamWithData< Real , Dim , Data >::reset( void )
{
	int fileType;
	float version;
	PlyProperty** plist;
	if( _ply ) _free();
	_ply = ply_open_for_reading( _fileName , &_nr_elems , &_elist , &fileType , &version );
	if( !_ply ) fprintf( stderr, "[ERROR] Failed to open ply file for reading: %s\n" , _fileName ) , exit( 0 );

	bool foundVertices = false;
	for( int i=0 ; i<_nr_elems ; i++ )
	{
		int num_elems;
		int nr_props;
		char* elem_name = _elist[i];
		plist = ply_get_element_description( _ply , elem_name , &num_elems , &nr_props );
		if( !plist ) fprintf( stderr , "[ERROR] Failed to get element description: %s\n" , elem_name ) , exit( 0 );

		if( equal_strings( "vertex" , elem_name ) )
		{
			foundVertices = true;
			_pCount = num_elems , _pIdx = 0;
			for( int i=0 ; i<PlyVertex< Real , Dim >::ReadComponents ; i++ ) 
				if( !ply_get_property( _ply , elem_name , &(PlyVertex< Real , Dim >::Properties()[i]) ) )
					fprintf( stderr , "[ERROR] Failed to find property in ply file: %s\n" , PlyVertex< Real , Dim >::Properties()[i].name ) , exit( 0 );
			if( _validationFunction )
			{
				bool* properties = new bool[_dataPropertiesCount];
				for( int i=0 ; i<_dataPropertiesCount ; i++ )
					if( !ply_get_property( _ply , elem_name , &(_dataProperties[i]) ) ) properties[i] = false;
					else                                                                properties[i] = true;
				bool valid = _validationFunction( properties );
				delete[] properties;
				if( !valid ) fprintf( stderr , "[ERROR] Failed to validate properties in file\n" ) , exit( 0 );
			}
			else
			{
				for( int i=0 ; i<_dataPropertiesCount ; i++ )
					if( !ply_get_property( _ply , elem_name , &(_dataProperties[i]) ) )
						fprintf( stderr , "[WARNING] Failed to find property in ply file: %s\n" , _dataProperties[i].name );
			}
		}
		for( int j=0 ; j<nr_props ; j++ )
		{
			free( plist[j]->name );
			free( plist[j] );
		}
		free( plist );
		if( foundVertices ) break;
	}
	if( !foundVertices )
	{
		fprintf( stderr , "[ERROR] Could not find vertices in ply file\n" );
		exit( 0 );
	}
}
template< class Real , int Dim , class Data >
void PLYInputPointStreamWithData< Real , Dim , Data >::_free( void )
{
	if( _ply ) ply_close( _ply ) , _ply = NULL;
	if( _elist )
	{
		for( int i=0 ; i<_nr_elems ; i++ ) free( _elist[i] );
		free( _elist );
	}
}
template< class Real , int Dim , class Data >
PLYInputPointStreamWithData< Real , Dim , Data >::~PLYInputPointStreamWithData( void )
{
	_free();
	if( _fileName ) delete[] _fileName , _fileName = NULL;
	if( _dataProperties ) delete[] _dataProperties , _dataProperties = NULL;
}
template< class Real , int Dim , class Data >
bool PLYInputPointStreamWithData< Real , Dim , Data >::nextPoint( Point< Real , Dim >& p , Data& d )
{
	if( _pIdx<_pCount )
	{
		_PlyVertexWithData v;
		ply_get_element( _ply, (void *)&v );
		p = v.point;
		d = v.data;
		_pIdx++;
		return true;
	}
	else return false;
}

//////////////////////////
// PLYOutputPointStream //
//////////////////////////
template< class Real , int Dim >
PLYOutputPointStream< Real , Dim >::PLYOutputPointStream( const char* fileName , size_t count , int fileType )
{
	static const char *elem_names[] = { "vertex" };
	float version;
	char* _fileName = new char[ strlen(fileName)+1];
	strcpy( _fileName , fileName );
	_ply = ply_open_for_writing( fileName , 1 , elem_names , fileType , &version );
	delete[] _fileName;
	if( !_ply )
	{
		fprintf( stderr, "[ERROR] Failed to open ply file for writing: %s\n" , fileName );
		exit( 0 );
	}

	_pIdx = 0;
	_pCount = count;
	ply_element_count( _ply, "vertex" , _pCount );
	for( int i=0 ; i<PlyVertex< Real , Dim >::WriteComponents ; i++ ) ply_describe_property( _ply , "vertex" , &PlyVertex< Real , Dim >::WriteProperties()[i] );

	ply_header_complete( _ply );
	ply_put_element_setup( _ply , "vertex" );
}
template< class Real , int Dim >
PLYOutputPointStream< Real , Dim >::~PLYOutputPointStream( void )
{
	if( _pIdx!=_pCount )
	{
		fprintf( stderr , "[ERROR] Streamed points not equal to total count: %d!=%d\n" , _pIdx , _pCount );
		exit( 0 );
	}
	ply_close( _ply );
	_ply = NULL;
}
template< class Real , int Dim >
void PLYOutputPointStream< Real , Dim >::nextPoint( const Point< Real , Dim >& p )
{
	if( _pIdx==_pCount )
	{
		fprintf( stderr , "[ERROR] Trying to add more points than total: %d<%d\n" , _pIdx , _pCount );
		exit( 0 );
	}
	PlyVertex< Real , Dim > op;
	op.point = p;
	ply_put_element( _ply , (void *) &op );
	_pIdx++;
}

//////////////////////////////////
// PLYOutputPointStreamWithData //
//////////////////////////////////
template< class Real , int Dim , class Data >
PLYOutputPointStreamWithData< Real , Dim , Data >::PLYOutputPointStreamWithData( const char* fileName , size_t count , int fileType , const PlyProperty* dataProperties , int dataPropertiesCount )
{
	static const char *elem_names[] = { "vertex" };
	float version;
	char* _fileName = new char[ strlen(fileName)+1];
	strcpy( _fileName , fileName );
	_ply = ply_open_for_writing( _fileName , 1 , &elem_names[0] , fileType , &version );
	delete[] _fileName;
	if( !_ply )
	{
		fprintf( stderr, "[ERROR] Failed to open ply file for writing: %s\n" , fileName );
		exit( 0 );
	}

	_pIdx = 0;
	_pCount = (int)count;
	ply_element_count( _ply, "vertex" , _pCount );
	for( int i=0 ; i<PlyVertex< Real , Dim >::WriteComponents ; i++ ) ply_describe_property( _ply , "vertex" , &PlyVertex< Real , Dim >::Properties()[i] );
	for( int i=0 ; i<dataPropertiesCount ; i++ )
	{
		PlyProperty prop = dataProperties[i];
		prop.offset += sizeof( PlyVertex< Real , Dim > );
		ply_describe_property( _ply , "vertex" , &prop );
	}

	ply_header_complete( _ply );
	ply_put_element_setup( _ply , "vertex" );
}
template< class Real , int Dim , class Data >
PLYOutputPointStreamWithData< Real , Dim , Data >::~PLYOutputPointStreamWithData( void )
{
	if( _pIdx!=_pCount )
	{
		fprintf( stderr , "[ERROR] Streamed points not equal to total count: %d!=%d\n" , _pIdx , _pCount );
		exit( 0 );
	}
	ply_close( _ply );
	_ply = NULL;
}
template< class Real , int Dim , class Data >
void PLYOutputPointStreamWithData< Real , Dim , Data >::nextPoint( const Point< Real , Dim >& p , const Data& d )
{
	if( _pIdx==_pCount )
	{
		fprintf( stderr , "[ERROR] Trying to add more points than total: %d<%d\n" , _pIdx , _pCount );
		exit( 0 );
	}
	_PlyVertexWithData op;
	op.point = p;
	op.data = d;
	ply_put_element( _ply , (void *) &op );
	_pIdx++;
}
