#ifndef PNG_INCLUDED
#define PNG_INCLUDED

#include "PNG/png.h"

struct PNGReader : public ImageReader
{
	PNGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~PNGReader( void );
	unsigned int nextRow( unsigned char* row );
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
protected:
	png_structp _png_ptr;
	png_infop _info_ptr;
	png_infop _end_info ;
	FILE* _fp;
	unsigned char* _scratchRow;
	unsigned int _currentRow;
};

struct PNGWriter : public ImageWriter
{
	PNGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	~PNGWriter( void );
	unsigned int nextRow( const unsigned char* row );
	unsigned int nextRows( const unsigned char* rows , unsigned int rowNum );
protected:
	FILE* _fp;
	png_structp _png_ptr;
	png_infop _info_ptr;
	unsigned int _currentRow;
};

#include "PNG.inl"
#endif //PNG_INCLUDED
