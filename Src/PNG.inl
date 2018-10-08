#include <stdio.h>
#include <vector>
#ifdef _WIN32
#include "PNG/png.h"
#else // !_WIN32
#include <png.h>
#endif // _WIN32

inline PNGReader::PNGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;

	_png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !_png_ptr ) fprintf( stderr , "[ERROR] PNGReader: failed to create png pointer\n" ) , exit( 0 );
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) fprintf( stderr , "[ERROR] PNGReader: failed to create info pointer\n" ) , exit( 0 );

	_end_info = png_create_info_struct( _png_ptr );
	if( !_end_info ) fprintf( stderr , "[ERROR] PNGReader: failed to create end pointer\n" ) , exit( 0 );


	_fp = fopen( fileName , "rb" );
	if( !_fp ) fprintf( stderr , "[ERROR] PNGReader: Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
	png_init_io( _png_ptr , _fp );

	png_read_info( _png_ptr, _info_ptr );

	width = png_get_image_width( _png_ptr , _info_ptr );
	height = png_get_image_height( _png_ptr, _info_ptr );
	channels = png_get_channels( _png_ptr , _info_ptr );
	int bit_depth=png_get_bit_depth( _png_ptr , _info_ptr );
	int color_type = png_get_color_type( _png_ptr , _info_ptr );
	if( bit_depth==16 )
	{
		fprintf( stderr , "[WARNING] Converting 16-bit image to 8-bit image\n" );
		_scratchRow = new unsigned char[ channels*width*2 ];
	}
	else
	{
		if( bit_depth!=8 ) fprintf( stderr , "[ERROR] PNGReader: expected 8 bits per channel\n" ) , exit( 0 );
		_scratchRow = NULL;
	}
	if( color_type==PNG_COLOR_TYPE_PALETTE ) png_set_expand( _png_ptr ) , printf( "Expanding PNG color pallette\n" );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}
inline unsigned int PNGReader::nextRow( unsigned char* row )
{
	if( _scratchRow )
	{
		int width = png_get_image_width( _png_ptr , _info_ptr );
		int channels = png_get_channels( _png_ptr , _info_ptr );

		png_read_row( _png_ptr , (png_bytep)_scratchRow , NULL );
#pragma omp parallel for
		for( int i=0 ; i<width*channels ; i++ ) row[i] = _scratchRow[2*i];
	}
	else png_read_row( _png_ptr , (png_bytep)row , NULL );
	return _currentRow++;
}

PNGReader::~PNGReader( void )
{
	if( _scratchRow ) delete[] _scratchRow;
	_scratchRow = NULL;
	png_destroy_read_struct( &_png_ptr , &_info_ptr , &_end_info );
}

inline bool PNGReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_infop end_info ;
	FILE* fp;

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !png_ptr ) fprintf( stderr , "[ERROR] PNGReader: failed to create png pointer\n" ) , exit( 0 );
	info_ptr = png_create_info_struct( png_ptr );
	if( !info_ptr ) fprintf( stderr , "[ERROR] PNGReader: failed to create info pointer\n" ) , exit( 0 );
	end_info = png_create_info_struct( png_ptr );
	if( !end_info ) fprintf( stderr , "[ERROR] PNGReader: failed to create end pointer\n" ) , exit( 0 );

	fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "[ERROR] PNGReader: Failed to open file for reading: %s\n" , fileName ) , exit( 0 );
	png_init_io( png_ptr , fp );

	png_read_info( png_ptr, info_ptr );

	width = png_get_image_width( png_ptr , info_ptr );
	height = png_get_image_height( png_ptr, info_ptr );
	channels = png_get_channels( png_ptr , info_ptr );

	png_destroy_read_struct( &png_ptr , &info_ptr , &end_info );
	fclose( fp );
	return true;
}

#if 1
PNGWriter::PNGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
{
	_currentRow = 0;

	_png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
	if( !_png_ptr )	fprintf( stderr , "[ERROR] PNGWriter: Failed to create png write struct\n" ) , exit(0);
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) fprintf( stderr , "[ERROR] PNGWriter: Failed to create png info struct\n") , exit(0);

	_fp = fopen( fileName , "wb" );
	if( !_fp ) fprintf( stderr , "[ERROR] PNGWriter: Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
	png_init_io( _png_ptr , _fp );

	png_set_compression_level( _png_ptr , Z_BEST_SPEED );

	int pngColorType;
	switch( channels )
	{
		case 1: pngColorType = PNG_COLOR_TYPE_GRAY ; break;
		case 3: pngColorType = PNG_COLOR_TYPE_RGB  ; break;
		case 4: pngColorType = PNG_COLOR_TYPE_RGBA ; break;
		default: fprintf( stderr , "[ERROR] PNGWriter: Only 1, 3, or 4 channel PNGs are supported\n" ) , exit( 0 );
	};
	png_set_IHDR( _png_ptr , _info_ptr, width , height, 8 , pngColorType , PNG_INTERLACE_NONE , PNG_COMPRESSION_TYPE_DEFAULT , PNG_FILTER_TYPE_DEFAULT );
	png_write_info( _png_ptr , _info_ptr );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}
PNGWriter::~PNGWriter( void )
{
	png_write_end( _png_ptr , NULL );
	png_destroy_write_struct( &_png_ptr , &_info_ptr );
	fclose( _fp );
}
unsigned int PNGWriter::nextRow( const unsigned char* row )
{
	png_write_row( _png_ptr , (png_bytep)row );
	return _currentRow++;
}
unsigned int PNGWriter::nextRows( const unsigned char* rows , unsigned int rowNum )
{
	for( unsigned int r=0 ; r<rowNum ; r++ ) png_write_row( _png_ptr , (png_bytep)( rows + r * 3 * sizeof( unsigned char ) * _png_ptr->width ) );
	return _currentRow += rowNum;
}
#else

void PNGWriteColor( const char* fileName , const unsigned char* pixels , int width , int height )
{
	FILE* fp = fopen( fileName , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName ) , exit( 0 );
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!png_ptr)	return;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return;
	png_init_io(png_ptr, fp);
	// turn off compression or set another filter
	// png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
	png_set_IHDR(png_ptr, info_ptr, width , height ,
		8,PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	if (0) {                    // high-level write
		std::vector<unsigned char> matrix( width * height * 3 );
		std::vector<png_bytep> row_pointers( height );
		for(int y=0;y<height;y++)
		{
			row_pointers[y]=&matrix[y*width*3];
			unsigned char* buf=&matrix[y*width*3];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3 + z ];
		}
		png_set_rows(png_ptr, info_ptr, &row_pointers[0]);
		int png_transforms=0;
		png_write_png(png_ptr, info_ptr, png_transforms, NULL);
	} else {                    // low-level write
		png_write_info(png_ptr, info_ptr);
		// png_set_filler(png_ptr, 0, PNG_FILLER_AFTER);
		//  but no way to provide GRAY data with RGBA fill, so pack each row
		std::vector<unsigned char> buffer(width*3);
		for(int y=0;y<height;y++)
		{
			unsigned char* buf=&buffer[0];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3+z ];
			png_bytep row_pointer=&buffer[0];
			png_write_row(png_ptr, row_pointer);
		}
	}
}
#endif