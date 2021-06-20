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

#ifndef VERTEX_STREAM_INCLUDED
#define VERTEX_STREAM_INCLUDED

#include <functional>
#include "Ply.h"
#include "Geometry.h"

template< typename Data >
class InputDataStream
{
public:
	virtual ~InputDataStream( void ){}
	virtual void reset( void ) = 0;
	virtual bool next( Data &d ) = 0;
	virtual size_t next( Data *d , size_t count )
	{
		size_t c=0;
		for( size_t i=0 ; i<count ; i++ , c++ ) if( !next( d[i] ) ) break;
		return c;
	}
};

template< typename Data >
class OutputDataStream
{
public:
	virtual ~OutputDataStream( void ){}
	virtual void next( const Data &d ) = 0;
	virtual void next( const Data *d , size_t count ){ for( size_t i=0 ; i<count ; i++ ) next( d[i] ); }
};

template< typename Data >
class TransformedInputDataStream : public InputDataStream< Data >
{
	std::function< void ( Data& ) > _xForm;
	InputDataStream< Data >& _stream;
public:
	TransformedInputDataStream( std::function< void ( Data & ) > xForm , InputDataStream< Data > &stream ) : _xForm(xForm) , _stream(stream) {;}
	virtual void reset( void ){ _stream.reset(); }
	virtual bool next( Data &d )
	{
		bool ret = _stream.next( d );
		_xForm( d );
		return ret;
	}
};

template< typename Data >
class TransformedOutputDataStream : public OutputDataStream< Data >
{
	std::function< void ( Data & ) > _xForm;
	OutputDataStream< Data > &_stream;
public:
	TransformedOutputDataStream( std::function< void ( Data & ) > xForm , OutputDataStream< Data >& stream ) : _xForm(xForm) , _stream(stream) {;}
	virtual void reset( void ){ _stream.reset(); }
	virtual bool next( const Data& p )
	{
		Data _p = p;
		_xForm( _p );
		return _stream.next( _p );
	}
};

template< typename Data >
class MemoryInputDataStream : public InputDataStream< Data >
{
	const Data *_data;
	size_t _size;
	size_t _current;
public:
	MemoryInputDataStream( size_t size , const Data *data );
	~MemoryInputDataStream( void );
	void reset( void );
	bool next( Data &d );
};

template< typename Data >
class MemoryOutputDataStream : public OutputDataStream< Data >
{
	Data *_data;
	size_t _size;
	size_t _current;
public:
	MemoryOutputDataStream( size_t size , Data *data );
	~MemoryOutputDataStream( void );
	void next( const Data &d );
};

template< typename Factory >
class ASCIIInputDataStream : public InputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	FILE *_fp;
public:
	ASCIIInputDataStream( const char* fileName , const Factory &factory );
	~ASCIIInputDataStream( void );
	void reset( void );
	bool next( Data &d );
};

template< typename Factory >
class ASCIIOutputDataStream : public OutputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	FILE *_fp;
public:
	ASCIIOutputDataStream( const char* fileName , const Factory &factory );
	~ASCIIOutputDataStream( void );
	void next( const Data &d );
};

template< typename Factory >
class BinaryInputDataStream : public InputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	FILE* _fp;
public:
	BinaryInputDataStream( const char* filename , const Factory &factory );
	~BinaryInputDataStream( void ){ fclose( _fp ) , _fp=NULL; }
	void reset( void );
	bool next( Data &d );
};

template< typename Factory >
class BinaryOutputDataStream : public OutputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	FILE* _fp;
public:
	BinaryOutputDataStream( const char* filename , const Factory &factory );
	~BinaryOutputDataStream( void ){ fclose( _fp ) , _fp=NULL; }
	void reset( void ){ fseek( _fp , SEEK_SET , 0 ); }
	void next( const Data &d );
};

template< typename Factory >
class PLYInputDataStream : public InputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	char* _fileName;
	PlyFile *_ply;
	std::vector< std::string > _elist;
	char *_buffer;

	size_t _pCount , _pIdx;
	void _free( void );
public:
	PLYInputDataStream( const char* fileName , const Factory &factory );
	~PLYInputDataStream( void );
	void reset( void );
	bool next( Data &d );
};

template< typename Factory >
class PLYOutputDataStream : public OutputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Data;
	const Factory &_factory;
	PlyFile *_ply;
	size_t _pCount , _pIdx;
	char *_buffer;
public:
	PLYOutputDataStream( const char* fileName , const Factory &factory , size_t count , int fileType=PLY_BINARY_NATIVE );
	~PLYOutputDataStream( void );
	void next( const Data &d );
};

#include "VertexStream.inl"

#endif // VERTEX_STREAM_INCLUDED
