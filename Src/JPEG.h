#ifndef JPEG_INCLUDED
#define JPEG_INCLUDED
#include "Image.h"

#include <setjmp.h>

#ifdef _WIN32
#include <windows.h>
#include "JPEG/jpeglib.h"
#include "JPEG/jerror.h"
#include "JPEG/jmorecfg.h"
#else // !_WIN32
#include <jpeglib.h>
#include <jerror.h>
#include <jmorecfg.h>
#endif // _WIN32

struct my_error_mgr
{
	struct jpeg_error_mgr pub;    // "public" fields
	jmp_buf setjmp_buffer;        // for return to caller
};
typedef struct my_error_mgr * my_error_ptr;

struct JPEGReader : public ImageReader
{
	JPEGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~JPEGReader( void );
	unsigned int nextRow( unsigned char* row );
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
protected:
	FILE* _fp;
	struct jpeg_decompress_struct _cInfo;
	struct my_error_mgr _jErr;
	unsigned int _currentRow;
};

struct JPEGWriter : public ImageWriter
{
	JPEGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	~JPEGWriter( void );
	unsigned int nextRow( const unsigned char* row );
	unsigned int nextRows( const unsigned char* rows , unsigned int rowNum );
protected:
	FILE* _fp;
	struct jpeg_compress_struct _cInfo;
	struct my_error_mgr _jErr;
	unsigned int _currentRow;
};

#include "JPEG.inl"
#endif //JPEG_INCLUDED
