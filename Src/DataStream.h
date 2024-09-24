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

namespace PoissonRecon
{

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
		unsigned int numStreams( void ) const { return (unsigned int)_streams.size(); }

	protected:
		std::vector< InputDataStream< Data > * > _streams;
		unsigned int _current;

		MultiInputDataStream( void ) {}
		void _init( InputDataStream< Data > **streams , size_t N )
		{
			_streams.resize( N );
			for( unsigned int i=0 ; i<N ; i++ ) _streams[i] = streams[i];
		}
		void _init( const std::vector< InputDataStream< Data > * > &streams ){ _streams = streams; }
		bool base_read( unsigned int t , Data &d ){ return _streams[t]->base_read(d); }
		bool base_read(                  Data &d )
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
		MultiOutputDataStream( OutputDataStream< Data > **streams , size_t N ) : _streams( streams , streams+N ) {}
		MultiOutputDataStream( const std::vector< OutputDataStream< Data > * > &streams ) : _streams( streams ) {}
		unsigned int numStreams( void ) const { return (unsigned int)_streams.size(); }

	protected:
		std::vector< OutputDataStream< Data > * > _streams;

		MultiOutputDataStream( void ) {}
		void _init( OutputDataStream< Data > **streams , size_t N )
		{
			_streams.resize( N );
			for( unsigned int i=0 ; i<N ; i++ ) _streams[i] = streams[i];
		}
		void _init( const std::vector< OutputDataStream< Data > * > &streams ){ _streams = streams; }
		void base_write(                  const Data &d ){ _streams[0]->base_write(d); }
		void base_write( unsigned int t , const Data &d ){ _streams[t]->base_write(d); }
	};

	template< typename Index , typename Data >
	struct IndexedInputDataStream : public InputDataStream< Data >
	{
		IndexedInputDataStream( MultiInputDataStream< std::pair< Index , Data > > &multiStream ) : _multiStream( multiStream ) , _firstTime(true)
		{
			_nextValues.resize( multiStream.numStreams() );
		}
		void reset( void ){ _multiStream.reset() , _firstTime = true; }

	protected:
		struct _NextValue
		{
			std::pair< Index , Data > data;
			unsigned int streamIndex;
			bool validData;

			// Returns true if v1>v2 so that it's a min-heap
			static bool Compare( const _NextValue &v1 , const _NextValue &v2 )
			{
				if( !v2.validData ) return false;
				else if( !v1.validData && v2.validData ) return true;
				else return v1.data.first>v2.data.first;
			};
		};

		MultiInputDataStream< std::pair< Index , Data > > &_multiStream;
		std::vector< _NextValue > _nextValues;
		bool _firstTime;

		void _init( const Data &data )
		{
			for( unsigned int i=0 ; i<_nextValues.size() ; i++ ) _nextValues[i].data.second = data;
			for( unsigned int i=0 ; i<_nextValues.size() ; i++ )
			{
				_nextValues[i].validData = _multiStream.read( i , _nextValues[i].data );
				_nextValues[i].streamIndex = i;
			}
			std::make_heap( _nextValues.begin() , _nextValues.end() , _NextValue::Compare );
		}

		bool base_read( unsigned int t , Data &d )
		{
			ERROR_OUT( "Multi-threaded read not supported" );
			return false;
		}
		bool base_read( Data &d )
		{
			if( _firstTime ) _init(d) , _firstTime = false;
			std::pop_heap( _nextValues.begin() , _nextValues.end() , _NextValue::Compare );
			_NextValue &next = _nextValues.back();
			if( !next.validData ) return false;
			d = next.data.second;
			next.validData = _multiStream.read( next.streamIndex , next.data );

			std::push_heap( _nextValues.begin() , _nextValues.end() , _NextValue::Compare );
			return true;
		}
	};
}

#endif // DATA_STREAM_INCLUDED
