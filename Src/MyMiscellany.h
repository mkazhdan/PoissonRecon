/*
Copyright (c) 2017, Michael Kazhdan
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
#ifndef MY_MISCELLANY_INCLUDED
#define MY_MISCELLANY_INCLUDED

#include <iostream>
#include <sstream>
#include <filesystem>
#include <thread>
#include <string.h>
#include <sys/timeb.h>
#include <cstdio>
#include <ctime>
#include <chrono>
#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#include <Windows.h>
#include <Psapi.h>
#else // !_WIN32 && !_WIN64
#include <unistd.h>
#include <sys/time.h> 
#include <sys/resource.h> 
#endif // _WIN32 || _WIN64
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <vector>
#if defined(_WIN32) || defined( _WIN64 )
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "Array.h"

namespace PoissonRecon
{
	////////////////////////////
	// Formatted float output //
	////////////////////////////
	struct StreamFloatPrecision
	{
		StreamFloatPrecision( std::ostream &str , unsigned int precision , bool scientific=false ) : _str(str)
		{
			_defaultPrecision = (int)_str.precision();
			_str.precision( precision );
			if( scientific ) _str << std::scientific;
			else             _str << std::fixed;
		}
		~StreamFloatPrecision( void )
		{
			_str << std::defaultfloat;
			_str.precision( _defaultPrecision );
		}
	protected:
		int _defaultPrecision;
		std::ostream &_str;
	};

	////////////////
	// Time Stuff //
	////////////////
	inline double Time( void )
	{
#ifdef WIN32
		struct _timeb t;
		_ftime( &t );
		return double( t.time ) + double( t.millitm ) / 1000.0;
#else // WIN32
		struct timeval t;
		gettimeofday( &t , NULL );
		return t.tv_sec + double( t.tv_usec ) / 1000000;
#endif // WIN32
	}

	struct Timer
	{
		Timer( void ){ _startCPUClock = std::clock() , _startWallClock = std::chrono::system_clock::now(); }
		double cpuTime( void ) const{ return (std::clock() - _startCPUClock) / (double)CLOCKS_PER_SEC; };
		double wallTime( void ) const{ std::chrono::duration<double> diff = (std::chrono::system_clock::now() - _startWallClock) ; return diff.count(); }
		std::string operator()( bool showCpuTime , unsigned int precision=1 ) const
		{
			std::stringstream ss;
			StreamFloatPrecision sfp( ss , precision );
			ss << wallTime() << " (s)";
			if( showCpuTime ) ss << " / " << cpuTime() << " (s)";
			return ss.str();
		}
		friend std::ostream &operator << ( std::ostream &os , const Timer &timer ){ return os << timer(false); }
	protected:
		std::clock_t _startCPUClock;
		std::chrono::time_point< std::chrono::system_clock > _startWallClock;
	};

	///////////////
	// I/O Stuff //
	///////////////
#if defined( _WIN32 ) || defined( _WIN64 )
	const char FileSeparator = '\\';
#else // !_WIN
	const char FileSeparator = '/';
#endif // _WIN

#ifndef SetTempDirectory
#if defined( _WIN32 ) || defined( _WIN64 )
#define SetTempDirectory( tempDir , sz ) GetTempPath( (sz) , (tempDir) )
#else // !_WIN32 && !_WIN64
#define SetTempDirectory( tempDir , sz ) if( std::getenv( "TMPDIR" ) ) strcpy( tempDir , std::getenv( "TMPDIR" ) );
#endif // _WIN32 || _WIN64
#endif // !SetTempDirectory

#if defined( _WIN32 ) || defined( _WIN64 )
	inline void FSync( FILE *fp )
	{
		//	FlushFileBuffers( (HANDLE)_fileno( fp ) );
		_commit( _fileno( fp ) );
	}
#else // !_WIN32 && !_WIN64
	inline void FSync( FILE *fp )
	{
		fsync( fileno( fp ) );
	}
#endif // _WIN32 || _WIN64

	template< typename Value > bool SetAtomic( volatile Value *value , Value newValue , Value oldValue );
	template< typename Data > void AddAtomic( Data& a , Data b );

	////////////////////
	// MKThread Stuff //
	////////////////////
	struct ThreadPool
	{
		enum ParallelType
		{
#ifdef _OPENMP
			OPEN_MP ,
#endif // _OPENMP
#ifdef SANITIZED_PR
#else // !SANITIZED_PR
			THREAD_POOL ,
#endif // SANITIZED_PR
			ASYNC ,
			NONE
		};
		static const std::vector< std::string > ParallelNames;

		enum ScheduleType
		{
			STATIC ,
			DYNAMIC
		};
		static const std::vector< std::string > ScheduleNames;

		static size_t DefaultChunkSize;
		static ScheduleType DefaultSchedule;

		template< typename ... Functions >
		static void ParallelSections( const Functions & ... functions )
		{
			std::vector< std::future< void > > futures( sizeof...(Functions) );
			_ParallelSections( &futures[0] , functions ... );
			for( size_t t=0 ; t<futures.size() ; t++ ) futures[t].get();
		}

		static void Parallel_for( size_t begin , size_t end , const std::function< void ( unsigned int , size_t ) > &iterationFunction , ScheduleType schedule=DefaultSchedule , size_t chunkSize=DefaultChunkSize )
		{
			if( begin>=end ) return;
			size_t range = end - begin;
			size_t chunks = ( range + chunkSize - 1 ) / chunkSize;
			unsigned int threads = (unsigned int)NumThreads();
			std::atomic< size_t > index;
			index.store( 0 );


			if( range<chunkSize || _ParallelType==NONE || threads==1 )
			{
				for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
				return;
			}

			auto _ChunkFunction = [ &iterationFunction , begin , end , chunkSize ]( unsigned int thread , size_t chunk )
				{
					const size_t _begin = begin + chunkSize*chunk;
					const size_t _end = std::min< size_t >( end , _begin+chunkSize );
					for( size_t i=_begin ; i<_end ; i++ ) iterationFunction( thread , i );
				};
			auto _StaticThreadFunction = [ &_ChunkFunction , chunks , threads ]( unsigned int thread )
				{
					for( size_t chunk=thread ; chunk<chunks ; chunk+=threads ) _ChunkFunction( thread , chunk );
				};
			auto _DynamicThreadFunction = [ &_ChunkFunction , chunks , &index ]( unsigned int thread )
				{
					size_t chunk;
					while( ( chunk=index.fetch_add(1) )<chunks ) _ChunkFunction( thread , chunk );
				};

			if     ( schedule==STATIC  ) _ThreadFunction = _StaticThreadFunction;
			else if( schedule==DYNAMIC ) _ThreadFunction = _DynamicThreadFunction;

			if( false ){}
#ifdef _OPENMP
			else if( _ParallelType==OPEN_MP )
			{
				if( schedule==STATIC )
#pragma omp parallel for num_threads( threads ) schedule( static , 1 )
					for( int c=0 ; c<chunks ; c++ ) _ChunkFunction( omp_get_thread_num() , c );
				else if( schedule==DYNAMIC )
#pragma omp parallel for num_threads( threads ) schedule( dynamic , 1 )
					for( int c=0 ; c<chunks ; c++ ) _ChunkFunction( omp_get_thread_num() , c );
			}
#endif // _OPENMP
			else if( _ParallelType==ASYNC )
			{
				static std::vector< std::future< void > > futures;
				futures.resize( threads-1 );
				for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1] = std::async( std::launch::async , _ThreadFunction , t );
				_ThreadFunction( 0 );
				for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1].get();
			}
#ifdef SANITIZED_PR
#else // !SANITIZED_PR
			else if( _ParallelType==THREAD_POOL )
			{
				unsigned int targetTasks = 0;
				if( !SetAtomic( &_RemainingTasks , threads-1 , targetTasks ) )
				{
					WARN( "nested for loop, reverting to serial" );
					for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
				}
				else
				{
					_WaitingForWorkOrClose.notify_all();
					{
						std::unique_lock< std::mutex > lock( _Mutex );
						_DoneWithWork.wait( lock , [&]( void ){ return _RemainingTasks==0; } );
					}
				}
			}
#endif // SANITIZED_PR
		}

		static unsigned int NumThreads( void ){ return (unsigned int)_Threads.size()+1; }

		static void Init( ParallelType parallelType , unsigned int numThreads=std::thread::hardware_concurrency() )
		{
			_ParallelType = parallelType;
			if( _Threads.size() && !_Close )
			{
				_Close = true;
				_WaitingForWorkOrClose.notify_all();
				for( unsigned int t=0 ; t<_Threads.size() ; t++ ) _Threads[t].join();
			}
			_Close = true;
			numThreads--;
			_Threads.resize( numThreads );
#ifdef SANITIZED_PR
#else // !SANITIZED_PR
			if( _ParallelType==THREAD_POOL )
			{
				_RemainingTasks = 0;
				_Close = false;
				for( unsigned int t=0 ; t<numThreads ; t++ ) _Threads[t] = std::thread( _ThreadInitFunction , t );
			}
#endif // SANITIZED_PR
		}
		static void Terminate( void )
		{
			if( _Threads.size() && !_Close )
			{
				_Close = true;
				_WaitingForWorkOrClose.notify_all();
				for( unsigned int t=0 ; t<_Threads.size() ; t++ ) _Threads[t].join();
				_Threads.resize( 0 );
			}
		}
	private:
		ThreadPool( const ThreadPool & ){}
		ThreadPool &operator = ( const ThreadPool & ){ return *this; }

		template< typename Function >
		static void _ParallelSections( std::future< void > *futures , const Function &function ){ *futures = std::async( std::launch::async , function ); }
		template< typename Function , typename ... Functions >
		static void _ParallelSections( std::future< void > *futures , const Function &function , const Functions& ... functions )
		{
			*futures = std::async( std::launch::async , function );
			_ParallelSections( futures+1 , functions ... );
		}
		static void _ThreadInitFunction( unsigned int thread )
		{
			// Wait for the first job to come in
			std::unique_lock< std::mutex > lock( _Mutex );
			_WaitingForWorkOrClose.wait( lock );
			while( !_Close )
			{
				lock.unlock();
				// do the job
				_ThreadFunction( thread );

				// Notify and wait for the next job
				lock.lock();
				_RemainingTasks--;
				if( !_RemainingTasks ) _DoneWithWork.notify_all();
				_WaitingForWorkOrClose.wait( lock );
			}
		}

#ifdef SANITIZED_PR
		static std::atomic< bool > _Close;
		static std::atomic< unsigned int > _RemainingTasks;
#else // !SANITIZED_PR
		static bool _Close;
		static volatile unsigned int _RemainingTasks;
#endif // SANITIZED_PR
		static std::mutex _Mutex;
		static std::condition_variable _WaitingForWorkOrClose , _DoneWithWork;
		static std::vector< std::thread > _Threads;
		static std::function< void ( unsigned int ) > _ThreadFunction;
		static ParallelType _ParallelType;
	};

	size_t ThreadPool::DefaultChunkSize = 128;
	ThreadPool::ScheduleType ThreadPool::DefaultSchedule = ThreadPool::DYNAMIC;
#ifdef SANITIZED_PR
	std::atomic< bool > ThreadPool::_Close;
	std::atomic< unsigned int > ThreadPool::_RemainingTasks;
#else // !SANITIZED_PR
	bool ThreadPool::_Close;
	volatile unsigned int ThreadPool::_RemainingTasks;
#endif // SANITIZED_PR
	std::mutex ThreadPool::_Mutex;
	std::condition_variable ThreadPool::_WaitingForWorkOrClose;
	std::condition_variable ThreadPool::_DoneWithWork;
	std::vector< std::thread > ThreadPool::_Threads;
	std::function< void ( unsigned int ) > ThreadPool::_ThreadFunction;
	ThreadPool::ParallelType ThreadPool::_ParallelType;

	const std::vector< std::string >ThreadPool::ParallelNames =
	{
#ifdef _OPENMP
		"open mp" ,
#endif // _OPENMP
#ifdef SANITIZED_PR
#else // !SANITIZED_PR
		"thread pool" ,
#endif // SANITIZED_PR
		"async" ,
		"none"
	};
	const std::vector< std::string >ThreadPool::ScheduleNames = { "static" , "dynamic" };

	template< typename Value >
	Value ReadAtomic32( const volatile Value * value )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		long _value = InterlockedExchangeAdd( (long*)value , 0 );
		return *(Value*)(&_value);
#else // !_WIN32 && !_WIN64
		uint32_t _value =  __atomic_load_n( (uint32_t *)value , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
		return *(Value*)(&_value);
	}

	template< typename Value >
	Value ReadAtomic64( const volatile Value * value )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 _value = InterlockedExchangeAdd64( (__int64*)value , 0 );
#else // !_WIN32 && !_WIN64
		uint64_t _value = __atomic_load_n( (uint64_t *)value , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
		return *(Value*)(&_value);
	}

	template< typename Value >
	bool SetAtomic32( volatile Value *value , Value newValue , Value oldValue )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		long *_oldValue = (long *)&oldValue;
		long *_newValue = (long *)&newValue;
		return InterlockedCompareExchange( (long*)value , *_newValue , *_oldValue )==*_oldValue;
#else // !_WIN32 && !_WIN64
		uint32_t *_newValue = (uint32_t *)&newValue;
		return __atomic_compare_exchange_n( (uint32_t *)value , (uint32_t *)&oldValue , *_newValue , false , __ATOMIC_SEQ_CST , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
	}
	template< typename Value >
	bool SetAtomic64( volatile Value *value , Value newValue , Value oldValue )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 *_oldValue = (__int64 *)&oldValue;
		__int64 *_newValue = (__int64 *)&newValue;
		return InterlockedCompareExchange64( (__int64*)value , *_newValue , *_oldValue )==*_oldValue;
#else // !_WIN32 && !_WIN64
		uint64_t *_newValue = (uint64_t *)&newValue;
		return __atomic_compare_exchange_n( (uint64_t *)value , (uint64_t *)&oldValue , *_newValue , false , __ATOMIC_SEQ_CST , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
	}

	template< typename Number >
	void AddAtomic32( volatile Number &a , Number b )
	{
#ifdef SANITIZED_PR
		Number current = ReadAtomic32( &a );
#else // !SANITIZED_PR
		Number current = a;
#endif // SANITIZED_PR
		Number sum = current+b;
#if defined( _WIN32 ) || defined( _WIN64 )
		long *_current = (long *)&current;
		long *_sum = (long *)&sum;
#ifdef SANITIZED_PR
		while( InterlockedCompareExchange( (long*)&a , *_sum , *_current )!=*_current )
		{
			current = ReadAtomic32( &a );
			sum = current + b;
	}
#else // !SANITIZED_PR
		while( InterlockedCompareExchange( (long*)&a , *_sum , *_current )!=*_current ) current = a , sum = a+b;
#endif // SANITIZED_PR
#else // !_WIN32 && !_WIN64
		uint32_t *_current = (uint32_t *)&current;
		uint32_t *_sum = (uint32_t *)&sum;
#ifdef SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint32_t *)&a , *_current , *_sum )!=*_current )
		{
			current = ReadAtomic32( &a );
			sum = current+b;
		}
#else // !SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint32_t *)&a , *_current , *_sum )!=*_current ) current = a , sum = a+b;
#endif // SANITIZED_PR
#endif // _WIN32 || _WIN64
	}

	template< typename Number >
	void AddAtomic64( volatile Number &a , Number b )
	{
#if 1
#ifdef SANITIZED_PR
		Number current = ReadAtomic64( &a );
		Number sum = current+b;
		while( !SetAtomic64( &a , sum , current ) )
		{
			current = ReadAtomic64( &a );
			sum = current+b;
		}
#else // !SANITIZED_PR
		Number current = a;
		Number sum = current+b;
		while( !SetAtomic64( &a , sum , current ) ) current = a , sum = a+b;
#endif // SANITIZED_PR
#else
#ifdef SANITIZED_PR
		Number current = ReadAtomic64( &a );
#else // !SANITIZED_PR
		Number current = a;
#endif // SANITIZED_PR
		Number sum = current+b;
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 *_current = (__int64 *)&current;
		__int64 *_sum = (__int64 *)&sum;
#ifdef SANITIZED_PR
		while( InterlockedCompareExchange64( (__int64*)&a , *_sum , *_current )!=*_current )
		{
			current = ReadAtomic64( &a );
			sum = current+b;
	}
#else // !SANITIZED_PR
		while( InterlockedCompareExchange64( (__int64*)&a , *_sum , *_current )!=*_current ) current = a , sum = a+b;
#endif // SANITIZED_PR
#else // !_WIN32 && !_WIN64
		uint64_t *_current = (uint64_t *)&current;
		uint64_t *_sum = (uint64_t *)&sum;
#ifdef SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint64_t *)&a , *_current , *_sum )!=*_current )
		{
			current = ReadAtomic64( &a);
			sum = current+b;
		}
#else // !SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint64_t *)&a , *_current , *_sum )!=*_current ) current = a , sum = a+b;
#endif // SANITIZED_PR
#endif // _WIN32 || _WIN64
#endif
	}

	template< typename Value >
	bool SetAtomic( volatile Value *value , Value newValue , Value oldValue )
	{
		switch( sizeof(Value) )
		{
		case 4: return SetAtomic32( value , newValue , oldValue );
		case 8: return SetAtomic64( value , newValue , oldValue );
		default:
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex setAtomicMutex;
			std::lock_guard< std::mutex > lock( setAtomicMutex );
			if( *value==oldValue ){ *value = newValue ; return true; }
			else return false;
		}
	}

	template< typename Data >
	void AddAtomic( Data& a , Data b )
	{
		switch( sizeof(Data) )
		{
		case 4: return AddAtomic32( a , b );
		case 8: return AddAtomic64( a , b );
		default:
			WARN_ONCE( "should not use this function: " , sizeof(Data) );
			static std::mutex addAtomicMutex;
			std::lock_guard< std::mutex > lock( addAtomicMutex );
			a += b;
		}
	}

	template< typename Value >
	Value ReadAtomic( const volatile Value * value )
	{
		if constexpr( sizeof(Value)==4 ) return ReadAtomic32( value );
		else if constexpr( sizeof(Value)==8 ) return ReadAtomic64( value );
		else
		{
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex readAtomicMutex;
			std::lock_guard< std::mutex > lock( readAtomicMutex );
			return *value;
		}
	}

	//////////////////
	// Memory Stuff //
	//////////////////
	size_t getPeakRSS( void );
	size_t getCurrentRSS( void );

	struct Profiler
	{
		Profiler( unsigned int ms=0 )
		{
			_t = Time();
			_currentPeak = 0;
			_terminate = false;
			if( ms )
			{
				_thread = std::thread( &Profiler::_updatePeakMemoryFunction , std::ref( *this ) , ms );
				_spawnedSampler = true;
			}
			else _spawnedSampler = false;
		}

		~Profiler( void )
		{
			if( _spawnedSampler )
			{
				_terminate = true;
				_thread.join();
			}
		}

		void reset( void )
		{
			_t = Time();
			if( _spawnedSampler )
			{
				std::lock_guard< std::mutex > lock( _mutex );
				_currentPeak = 0;
			}
			else _currentPeak = 0;
		}

		void update( void )
		{
			size_t currentPeak = getCurrentRSS();
			if( _spawnedSampler )
			{
				std::lock_guard< std::mutex > lock( _mutex );
				if( currentPeak>_currentPeak ) _currentPeak = currentPeak;
			}
			else if( currentPeak>_currentPeak ) _currentPeak = currentPeak;
		}

		std::string operator()( bool showTime=true ) const
		{
			std::stringstream ss;
			double dt = Time()-_t;
			double  localPeakMB = ( (double)_currentPeak )/(1<<20);
			double globalPeakMB = ( (double)getPeakRSS() )/(1<<20);
			{
				StreamFloatPrecision sfp( ss , 1 );
				if( showTime ) ss << dt << " (s), ";
				ss << localPeakMB << " (MB) / " << globalPeakMB << " (MB)";
			}
			return ss.str();
		}

		friend std::ostream &operator << ( std::ostream &os , const Profiler &profiler ){ return os << profiler(); }

	protected:
		std::thread _thread;
		std::mutex _mutex;
		double _t;
		std::atomic< bool > _spawnedSampler;
		std::atomic< size_t > _currentPeak;
		std::atomic< bool > _terminate;

		void _updatePeakMemoryFunction( unsigned int ms )
		{
			while( true )
			{
				std::this_thread::sleep_for( std::chrono::milliseconds( ms ) );
				update();
				if( _terminate ) return;
			}
		};
	};

	struct MemoryInfo
	{
		static size_t Usage( void ){ return getCurrentRSS(); }
		static int PeakMemoryUsageMB( void ){ return (int)( getPeakRSS()>>20 ); }
	};

#if defined( _WIN32 ) || defined( _WIN64 )
	inline void SetPeakMemoryMB( size_t sz )
	{
		sz <<= 20;
		SIZE_T peakMemory = sz;
		HANDLE h = CreateJobObject( NULL , NULL );
		AssignProcessToJobObject( h , GetCurrentProcess() );

		JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
		jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_MEMORY;
		jeli.JobMemoryLimit = peakMemory;
		if( !SetInformationJobObject( h , JobObjectExtendedLimitInformation , &jeli , sizeof( jeli ) ) ) WARN( "Failed to set memory limit" );
	}
#else // !_WIN32 && !_WIN64
	inline void SetPeakMemoryMB( size_t sz )
	{
		sz <<= 20;
		struct rlimit rl;
		getrlimit( RLIMIT_AS , &rl );
		rl.rlim_cur = sz;
		setrlimit( RLIMIT_AS , &rl );
	}
#endif // _WIN32 || _WIN64

	/*
	* Author:  David Robert Nadeau
	* Site:    http://NadeauSoftware.com/
	* License: Creative Commons Attribution 3.0 Unported License
	*          http://creativecommons.org/licenses/by/3.0/deed.en_US
	*/

	/**
	* Returns the peak (maximum so far) resident set size (physical
	* memory use) measured in bytes, or zero if the value cannot be
	* determined on this OS.
	*/
	inline size_t getPeakRSS( )
	{
#if defined(_WIN32)
		/* Windows -------------------------------------------------- */
		PROCESS_MEMORY_COUNTERS info;
		GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
		return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
		/* AIX and Solaris ------------------------------------------ */
		struct psinfo psinfo;
		int fd = -1;
		if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
			return (size_t)0L;      /* Can't open? */
		if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
		{
			close( fd );
			return (size_t)0L;      /* Can't read? */
		}
		close( fd );
		return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
		/* BSD, Linux, and OSX -------------------------------------- */
		struct rusage rusage;
		getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
		return (size_t)rusage.ru_maxrss;
#else
		return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
		/* Unknown OS ----------------------------------------------- */
		return (size_t)0L;          /* Unsupported. */
#endif
	}





	/**
	* Returns the current resident set size (physical memory use) measured
	* in bytes, or zero if the value cannot be determined on this OS.
	*/
	inline size_t getCurrentRSS( )
	{
#if defined(_WIN32) || defined( _WIN64 )
		/* Windows -------------------------------------------------- */
		PROCESS_MEMORY_COUNTERS info;
		GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
		return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
		/* OSX ------------------------------------------------------ */
		struct mach_task_basic_info info;
		mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
		if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
			(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
			return (size_t)0L;      /* Can't access? */
		return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
		/* Linux ---------------------------------------------------- */
		long rss = 0L;
		FILE* fp = NULL;
		if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
			return (size_t)0L;      /* Can't open? */
		if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
		{
			fclose( fp );
			return (size_t)0L;      /* Can't read? */
		}
		fclose( fp );
		return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
		/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
		return (size_t)0L;          /* Unsupported. */
#endif
	}
}

#endif // MY_MISCELLANY_INCLUDED
