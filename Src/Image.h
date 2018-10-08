#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED

#define SUPPORT_TILES

#include <string.h>

struct ImageReader
{
	virtual unsigned int nextRow( unsigned char* row ) = 0;
	static unsigned char* Read( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
	{
		ImageReader* reader = Get( fileName );
		width = reader->width() , height = reader->height() , channels = reader->channels();
		unsigned char* pixels = new unsigned char[ width*height*channels ];
		for( unsigned int j=0 ; j<height ; j++ ) reader->nextRow( pixels + j*width*channels );
		delete reader;
		return pixels;
	}
	static unsigned char* ReadColor( const char* fileName , unsigned int& width , unsigned int& height )
	{
		unsigned int channels;
		ImageReader* reader = Get( fileName );
		width = reader->width() , height = reader->height();
		if( channels!=1 && channels!=3 ) fprintf( stderr , "[ERROR] ImageReader::ReadColor requres one- or three-channel input\n" ) , exit( 0 );
		unsigned char* pixels = new unsigned char[ width*height*3 ];
		unsigned char* pixelRow = new unsigned char[ width*channels];
		for( unsigned int j=0 ; j<height ; j++ )
		{
			reader->nextRow( pixelRow );
			if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(unsigned char)*width*3 );
			else if( channels==1 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i];
		}
		delete[] pixelRow;
		delete reader;
		return pixels;
	}

	static ImageReader* Get( const char* fileName );
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	virtual ~ImageReader( void ){ }
	unsigned int width( void ) const { return _width; }
	unsigned int height( void ) const { return _height; }
	unsigned int channels( void ) const { return _channels; }
protected:
	unsigned int _width , _height , _channels;
};
struct ImageWriterParams
{
	static const char* DefaultTileExtension;
	unsigned int quality;
#ifdef SUPPORT_TILES
	const char* tileExtension;
	unsigned int tileWidth , tileHeight;
	ImageWriterParams* tileParams;
	ImageWriterParams( void ) : quality( 100 ) , tileExtension( DefaultTileExtension ) , tileWidth( 4096 ) , tileHeight( 4096 ) , tileParams( NULL ){};
#else // !SUPPORT_TILES
	ImageWriterParams( void ) : quality( 100 ){};
#endif // SUPPORT_TILES
};
const char* ImageWriterParams::DefaultTileExtension = "jpg";
struct ImageWriter
{
	virtual unsigned int nextRow( const unsigned char* row ) = 0;
	virtual unsigned int nextRows( const unsigned char* rows , unsigned int rowNum ){ unsigned int row ; for( unsigned int r=0 ; r<rowNum ; r++ ) row = nextRow( rows + _width * _channels * r ) ; return row; }
	static bool Write( const char* fileName , const unsigned char* pixels , unsigned int width , unsigned int height , unsigned int channels , ImageWriterParams params=ImageWriterParams() )
	{
		ImageWriter* writer = Get( fileName , width , height , channels , params );
		if( !writer ) return false;
		for( unsigned int j=0 ; j<height ; j++ ) writer->nextRow( pixels + j*width*channels );
		delete writer;
		return true;
	}
	static ImageWriter* Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , ImageWriterParams params=ImageWriterParams() );
	virtual ~ImageWriter( void ){ }
	unsigned int width( void ) const { return _width; }
	unsigned int height( void ) const { return _height; }
	unsigned int channels( void ) const { return _channels; }
protected:
	unsigned int _width , _height , _channels;
};

#ifdef SUPPORT_TILES
struct TiledImageReader : public ImageReader
{
	unsigned int nextRow( unsigned char* row );
	TiledImageReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~TiledImageReader( void );
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
protected:
	ImageReader** _tileReaders;
	char** _tileNames;
	unsigned int _tileRows , _tileColumns , _currentPixelRow , _currentTileRow , *_tileWidths , *_tileHeights;
};
struct TiledImageWriter : public ImageWriter
{
	unsigned int nextRow( const unsigned char* row );
	TiledImageWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , ImageWriterParams params );
	~TiledImageWriter( void );
protected:
	ImageWriter** _tileWriters;
	char** _tileNames;
	unsigned int _tileWidth , _tileHeight , _tileRows , _tileColumns , _currentPixelRow;
	ImageWriterParams _params;
};
#endif // SUPPORT_TILES


// [WARNING] Need to include "png.h" before "jpeg.h" so that "setjmp.h" is not already included (?)
#include "PNG.h"
#include "JPEG.h"

struct FileNameParser
{
#if defined( _WIN32 ) || defined( _WIN64 )
	static const char Separator = (char)'\\';
#else // !_WIN
	static const char Separator = (char)'/';
#endif // _WIN
	static inline char* Extension  ( const char* fileName ){ return __Split( fileName , '.' , false , false ); }
	static inline char* Header     ( const char* fileName ){ return __Split( fileName , '.' , true , false ); }
	static inline char* Local      ( const char* fileName ){ return __Split( fileName , Separator , false , false ); }
	static inline char* Dir        ( const char* fileName ){ return __Split( fileName , Separator , true , false ); }
	static inline char* LocalHeader( const char* fileName )
	{
		char* localFileName = Local( fileName );
		if( !localFileName ) fprintf( stderr , "[ERROR] GetLocalFileHeader: Couldn't get local file name: %s\n" , fileName ) , exit( 0 );
		char* localFileHeader = Header( localFileName );
		delete[] localFileName;
		return localFileHeader;
	}

protected:
	static inline char* __Split( const char* fileName , char splitChar , bool front , bool first )
	{
		int position;
		char* out;
		if( first ){ for( position=0 ; position<strlen(fileName) ; position++ ) if( fileName[position]==splitChar ) break; }
		else       { for( position=(int)strlen(fileName)-1 ; position>=0 ; position-- ) if( fileName[position]==splitChar ) break; }

		if( front )
		{
			if( position==-1 ) out = NULL;
			else
			{
				out = new char[ strlen(fileName)+1 ];
				strcpy( out , fileName );
				out[ position ] = 0;
			}
		}
		else
		{
			if( position==strlen(fileName) ) out = NULL;
			else
			{
				out = new char[ strlen(fileName)-position ];
				strcpy( out , fileName+position+1 );
			}
		}
		return out;
	}
};

inline ImageReader* ImageReader::Get( const char* fileName )
{
	unsigned int width , height , channels;
	ImageReader* reader = NULL;
	char* ext = FileNameParser::Extension( fileName );
#ifdef WIN32
	if     ( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) reader = new       JPEGReader( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" )                              ) reader = new        PNGReader( fileName , width , height , channels );
	else if( !_stricmp( ext , "iGrid" )                            ) reader = new TiledImageReader( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) reader = new       JPEGReader( fileName , width , height , channels );
	else if( !strcasecmp( ext , "png" )                           ) reader = new        PNGReader( fileName , width , height , channels );
	else if( !strcasecmp( ext , "iGrid" )                         ) reader = new TiledImageReader( fileName , width , height , channels );
#endif // WIN32
	reader->_width = width;
	reader->_height = height;
	reader->_channels = channels;

	delete[] ext;
	return reader;
}
inline void ImageReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	char* ext = FileNameParser::Extension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" ) )                          PNGReader::GetInfo( fileName , width , height , channels );
	else if( !_stricmp( ext , "iGrid" ) )                 TiledImageReader::GetInfo( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels );
	else if( !strcasecmp( ext , "png" ) )                            PNGReader::GetInfo( fileName , width , height , channels );
	else if( !strcasecmp( ext , "iGrid" ) )                   TiledImageReader::GetInfo( fileName , width , height , channels );
#endif // WIN32
	delete[] ext;
}
inline ImageWriter* ImageWriter::Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , ImageWriterParams params )
{
	ImageWriter* writer = NULL;
	char* ext = FileNameParser::Extension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , params.quality );
	else if( !_stricmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , params.quality );
#ifdef SUPPORT_TILES
	else if( !_stricmp( ext , "iGrid" ) ) writer = new TiledImageWriter( fileName , width , height , channels , params );
#endif // SUPPORT_TILES
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , params.quality );
	else if( !strcasecmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , params.quality );
#ifdef SUPPORT_TILES
	else if( !strcasecmp( ext , "iGrid" ) ) writer = new TiledImageWriter( fileName , width , height , channels , params );
#endif // SUPPORT_TILES
#endif // WIN32
	else fprintf( stderr , "[ERROR] ImageWriter::Get: Unrecognized file extension: %s\n" , ext ) , exit( 0 );
	writer->_width = width;
	writer->_height = height;
	writer->_channels = channels;

	delete[] ext;
	return writer;
}

#ifdef SUPPORT_TILES
bool TiledImageReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	char* fileDir = FileNameParser::Dir( fileName );
	unsigned int *_tileHeights , *_tileWidths;
	unsigned int _tileRows , _tileColumns , _channels;
	FILE* fp = fopen( fileName , "r" );
	if( !fp ){ fprintf( stderr , "[WARNING] TiledImageReader::GetInfo: Couldn't open file for reading: %s\n" , fileName ) ; return false; }
	{
		char line[1024];
		if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Failed to read column line from: %s\n" , fileName ) , exit( 0 );
		line[strlen(line)-1] = 0;
		if( sscanf( line , "Columns: %d" , &_tileColumns )!=1 ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Failed to read column count from: %s (%s)\n" , fileName , line ) , exit( 0 );
		if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Failed to read row line from: %s\n" , fileName ) , exit( 0 );
		line[strlen(line)-1] = 0;
		if( sscanf( line , "Rows: %d" , &_tileRows )!=1 ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Failed to read row count from: %s (%s)\n" , fileName , line ) , exit( 0 );
		_tileHeights = new unsigned int[ _tileRows+1 ];
		_tileWidths  = new unsigned int[ _tileColumns+1 ];

		char tileName[1024];
		for( unsigned int r=0 ; r<_tileRows ; r++ ) for( unsigned int c=0 ; c<_tileColumns ; c++ )
		{
			if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Failed to read tile name from: %s\n" , fileName ) , exit( 0 );
			line[strlen(line)-1] = 0;
			if( fileDir ) sprintf( tileName , "%s%c%s" , fileDir , FileNameParser::Separator , line );
			else          sprintf( tileName , "%s" , line );

			unsigned int _w , _h , _c;
			ImageReader::GetInfo( tileName , _w , _h , _c );
			if( !r && !c ) _channels = _c;
			else if( _channels!=_c ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Number of color channels don't match: %d != %d\n" , _channels , _c ) , exit( 0 );
			if( !r ) _tileWidths[c+1] = _w;
			else if( _tileWidths[c+1]!=_w ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Images in the same column must have the same width: %d != %d\n" , _tileWidths[c+1] , _w ) , exit( 0 );
			if( !c ) _tileHeights[r+1] = _h;
			else if( _tileHeights[r+1]!=_h ) fprintf( stderr , "[ERROR] TiledImageReader::GetInfo: Images in the same row must have the same heights: %d != %d\n" , _tileHeights[r+1] , _h ) , exit( 0 );
		}
	}
	fclose( fp );
	if( fileDir ) delete[] fileDir;
	_tileWidths[0] = _tileHeights[0] = 0;
	for( unsigned int c=0 ; c<_tileColumns ; c++ ) _tileWidths[c+1] += _tileWidths[c];
	for( unsigned int r=0 ; r<_tileRows ; r++ ) _tileHeights[r+1] += _tileHeights[r];
	width = _tileWidths[_tileColumns] , height = _tileHeights[_tileRows] , channels = _channels;
	return true;
}

TiledImageReader::TiledImageReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	char* fileDir = FileNameParser::Dir( fileName );
	FILE* fp = fopen( fileName , "r" );
	if( !fp ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Couldn't open file for reading: %s\n" , fileName ) , exit( 0 );
	{
		char line[1024];
		if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Failed read column line from: %s\n" , fileName ) , exit( 0 );
		line[strlen(line)-1] = 0;
		if( sscanf( line , "Columns: %d" , &_tileColumns )!=1 ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Failed to read column count from: %s (%s)\n" , fileName , line ) , exit( 0 );
		if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Failed read row line from: %s\n" , fileName ) , exit( 0 );
		line[strlen(line)-1] = 0;
		if( sscanf( line , "Rows: %d" , &_tileRows )!=1 ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Failed to read row count from: %s (%s)\n" , fileName , line ) , exit( 0 );

		_tileReaders = new ImageReader*[ _tileColumns ];
		_tileHeights = new unsigned int[ _tileRows+1 ];
		_tileWidths  = new unsigned int[ _tileColumns+1 ];

		_tileNames = new char*[ _tileColumns * _tileRows ];
		char tileName[1024];
		for( unsigned int r=0 ; r<_tileRows ; r++ ) for( unsigned int c=0 ; c<_tileColumns ; c++ )
		{
			if( !fgets( line , 1024 , fp ) ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Failed read tile name from: %s\n" , fileName ) , exit( 0 );
			line[strlen(line)-1] = 0;
			if( fileDir ) sprintf( tileName , "%s%c%s" , fileDir , FileNameParser::Separator , line );
			else          sprintf( tileName , "%s" , line );
			_tileNames[r*_tileColumns+c] = new char[ strlen(tileName)+1 ];
			strcpy( _tileNames[r*_tileColumns+c] , tileName );
		}
	}
	fclose( fp );
	if( fileDir ) delete[] fileDir;
	for( unsigned int r=0 ; r<_tileRows ; r++ ) for( unsigned int c=0 ; c<_tileColumns ; c++ )
	{
		unsigned int _w , _h , _c;
		ImageReader::GetInfo( _tileNames[r*_tileColumns+c] , _w , _h , _c );
		if( !r && !c ) _channels = _c;
		else if( _channels!=_c ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Number of color channels don't match: %d != %d\n" , _channels , _c ) , exit( 0 );
		if( !r ) _tileWidths[c+1] = _w;
		else if( _tileWidths[c+1]!=_w ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Images in the same column must have the same width: %d != %d\n" , _tileWidths[c+1] , _w ) , exit( 0 );
		if( !c ) _tileHeights[r+1] = _h;
		else if( _tileHeights[r+1]!=_h ) fprintf( stderr , "[ERROR] TiledImageReader::TiledImageReader: Images in the same row must have the same heights: %d != %d\n" , _tileHeights[r+1] , _h ) , exit( 0 );
	}
	_tileWidths[0] = _tileHeights[0] = 0;
	for( unsigned int c=0 ; c<_tileColumns ; c++ ) _tileWidths[c+1] += _tileWidths[c];
	for( unsigned int r=0 ; r<_tileRows ; r++ ) _tileHeights[r+1] += _tileHeights[r];
	width = _width = _tileWidths[_tileColumns] , height = _height = _tileHeights[_tileRows] , channels = _channels;
	_currentPixelRow = _currentTileRow = 0;
}
TiledImageReader::~TiledImageReader( void )
{
	delete[] _tileReaders;
	for( unsigned int i=0 ; i<_tileColumns*_tileRows ; i++ ) delete[] _tileNames[i];
	delete[] _tileNames;
	delete[] _tileWidths;
	delete[] _tileHeights;
}
unsigned TiledImageReader::nextRow( unsigned char* row )
{
	// If it's the first row, set up the readers
	if( _currentPixelRow==_tileHeights[ _currentTileRow ] ) for( unsigned int c=0 ; c<_tileColumns ; c++ ) _tileReaders[c] = ImageReader::Get( _tileNames[ _currentTileRow * _tileColumns + c ] );

	// Read the row fragments
	for( unsigned int c=0 ; c<_tileColumns ; c++ ) _tileReaders[c]->nextRow( row + c * _tileWidths[c] * _channels );

	// If it's the last row of the tile, free up the readers
	if( _currentPixelRow==_tileHeights[_currentTileRow+1]-1 )
	{
		for( unsigned int c=0 ; c<_tileColumns ; c++ ) delete _tileReaders[c];
		_currentTileRow++;
	}

	return _currentPixelRow++;
}

TiledImageWriter::TiledImageWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , ImageWriterParams params )
{
	_width = width , _height = height , _channels = channels , _tileWidth = params.tileWidth , _tileHeight = params.tileHeight;
	_tileColumns = ( _width + ( _tileWidth-1 ) ) / _tileWidth , _tileRows = ( _height + ( _tileHeight-1 ) ) / _tileHeight;
	_tileWriters = new ImageWriter*[ _tileColumns ];
	_tileNames = new char*[ _tileColumns * _tileRows ];
	if( params.tileParams ) _params = *params.tileParams;

	char tileName[1024];
	char* tileHeader = FileNameParser::Header( fileName );
	for( unsigned int r=0 ; r<_tileRows ; r++ ) for( unsigned int c=0 ; c<_tileColumns ; c++ )
	{
		sprintf( tileName , "%s.%d.%d.%s" , tileHeader , c , r , params.tileExtension );
		_tileNames[r*_tileColumns+c] = new char[ strlen(tileName)+1 ];
		strcpy( _tileNames[r*_tileColumns+c] , tileName );
	}
	delete[] tileHeader;
	FILE* fp = fopen( fileName , "w" );
	if( !fp ) fprintf( stderr , "[ERROR] TiledImageWriter::TiledImageWriter: Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
	fprintf( fp , "Columns: %d\n" , _tileColumns );
	fprintf( fp , "Rows: %d\n" , _tileRows );
	for( unsigned int i=0 ; i<_tileRows*_tileColumns ; i++ )
	{
		char* localTileName = FileNameParser::Local( _tileNames[i] );
		fprintf( fp , "%s\n" , localTileName );
		delete[] localTileName;
	}
	fclose( fp );
	_currentPixelRow = 0;
}
TiledImageWriter::~TiledImageWriter( void )
{
	delete[] _tileWriters;
	for( unsigned int i=0 ; i<_tileColumns*_tileRows ; i++ ) delete[] _tileNames[i];
	delete[] _tileNames;
}
unsigned int TiledImageWriter::nextRow( const unsigned char* row )
{
	unsigned int r = _currentPixelRow / _tileHeight;
	if( ( _currentPixelRow % _tileHeight )==0 )
	{
		for( unsigned int c=0 ; c<_tileColumns ; c++ )
			_tileWriters[c] = ImageWriter::Get( _tileNames[ r * _tileColumns + c ] , std::min< unsigned int >( _tileWidth , _width - _tileWidth*c ) , std::min< unsigned int >( _tileHeight , _height - _tileHeight*r ) , _channels , _params );
	}
	for( int c=0 ; c<(int)_tileColumns ; c++ ) _tileWriters[c]->nextRow( row + c * _tileWidth * _channels );
	if( ( _currentPixelRow % _tileHeight )==( _tileHeight-1 ) || _currentPixelRow==(_height-1) ) for( unsigned int c=0 ; c<_tileColumns ; c++ ) delete _tileWriters[c];

	return _currentPixelRow++;
}
#endif // SUPPORT_TILES

#endif // IMAGE_INCLUDED
