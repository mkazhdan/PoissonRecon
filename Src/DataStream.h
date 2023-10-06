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

#ifndef DATA_STREAM_INCLUDED
#define DATA_STREAM_INCLUDED

#include <mutex>
#include <vector>
#include <atomic>

// Pre-declare so we can make friends
template< typename Data > struct MultiInputDataStream;
template< typename Data > struct MultiOutputDataStream;

////////////////////////////////////////
// Abstract input/output data streams //
////////////////////////////////////////

// An input stream containing "Data" types
// Supporting:
// -- Resetting the stream to the start
// -- Trying to read the next element from the stream
template< typename Data >
struct InputDataStream
{
	friend struct MultiInputDataStream< Data >;

	virtual ~InputDataStream( void ){}

	// Reset to the start of the stream
	virtual void reset( void ) = 0;

	bool read( Data &d ){ return base_read(d); }
	bool read( unsigned int thread , Data &d ){ return base_read(thread,d); }

protected:
	std::mutex _insertionMutex;

	// Read in data in a single-threaded context
	virtual bool base_read( Data &d ) = 0;

	// Read in data in a multi-threaded context
	virtual bool base_read( unsigned int thread , Data &d )
	{
		std::lock_guard< std::mutex > lock( _insertionMutex );
		return base_read(d);
	}
};

// An output stream containing "Data" types
// Supporting:
// -- Writing the next element to the stream
// -- Writing the next element to the stream by a particular thread
template< typename Data >
struct OutputDataStream
{
	friend struct MultiOutputDataStream< Data >;

	OutputDataStream( void ) : _size(0) {}
	virtual ~OutputDataStream( void ){}

	// Returns the number of items written to the stream
	size_t size( void ) const { return _size; }

	void write( const Data &d ){ base_write(d) ; _size++; }
	void write( unsigned int thread , const Data &d ){ base_write( thread , d ) ; _size++; }

protected:
	std::mutex _insertionMutex;
	std::atomic< size_t > _size;

	// Write out data in a single-threaded context
	virtual void base_write( const Data &d ) = 0;
	// Write out data in a multi-threaded context
	virtual void base_write( unsigned int thread , const Data &d )
	{
		std::lock_guard< std::mutex > lock( _insertionMutex );
		return base_write(d);
	}

};

//////////////////////////////////////////
// Multi-streams for multi-threaded I/O //
//////////////////////////////////////////
template< typename Data >
struct MultiInputDataStream : public InputDataStream< Data >
{
	MultiInputDataStream( InputDataStream< Data > **streams , size_t N ) : _current(0) , _streams( streams , streams+N ) {}
	MultiInputDataStream( const std::vector< InputDataStream< Data > * > &streams ) : _current(0) , _streams( streams ) {}
	void reset( void ){ for( unsigned int i=0 ; i<_streams.size() ; i++ ) _streams[i]->reset(); }

protected:
	std::vector< InputDataStream< Data > * > _streams;
	unsigned int _current;

	bool base_read( unsigned int t , Data &d ){ return _streams[t]->base_read(d); }
	bool base_read( Data &d )
	{
		while( _current<_streams.size() )
		{
			if( _streams[_current]->read( d ) ) return true;
			else _current++;
		}
		return false;
	}
};

template< typename Data >
struct MultiOutputDataStream : public OutputDataStream< Data >
{
	MultiOutputDataStream( OutputDataStream< Data > **streams , size_t N ) : _current(0) , _streams( streams , streams+N ) {}
	MultiOutputDataStream( const std::vector< OutputDataStream< Data > * > &streams ) : _current(0) , _streams( streams ) {}

protected:
	std::vector< OutputDataStream< Data > * > _streams;
	unsigned int _current;

	void base_write( const Data &d ){ _streams[0]->base_write(d); }
	void base_write( unsigned int t , const Data &d ){ _streams[t]->base_write(d); }
};


#endif // DATA_STREAM_INCLUDED
