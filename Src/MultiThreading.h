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
#ifndef MULTI_THREADING_INCLUDED
#define MULTI_THREADING_INCLUDED

#include <thread>
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
#include <mutex>
#include <vector>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "Array.h"
#include "MyAtomic.h"

namespace PoissonRecon
{
#ifdef USE_HH_THREADS
	namespace hh
	{
		inline int get_max_threads( void ) { return std::max< int >( (int)std::thread::hardware_concurrency() , 1 ); }

		class ThreadPoolIndexedTask
		{
		public:
			using Task = std::function< void (int) >;
			ThreadPoolIndexedTask( void )
			{
				const int num_threads = get_max_threads()-1;
				_threads.reserve( num_threads );
				for( int i=0 ; i<num_threads ; i++ ) _threads.emplace_back( &ThreadPoolIndexedTask::worker_main , this );
			}
			~ThreadPoolIndexedTask( void )
			{
				if( _threads.size() )
				{
					std::unique_lock< std::mutex > lock( _mutex );
					if( !_running ) ERROR_OUT( "not running" );
					if( _num_remaining_tasks ) ERROR_OUT( "tasks remaining" );
					if( _task_index!=_num_tasks ) ERROR_OUT( "task index and num tasks don't match" );
					_running = false;
					_task_index = 0;
					_num_tasks = 1;
					_condition_variable_worker.notify_all();
				}
				for( int i=0 ; i<_threads.size() ; i++ ) _threads[i].join();
			}
			int num_threads( void ) const { return (int)_threads.size()+1; }
			bool already_active( void ) const { return _num_remaining_tasks!=0; }  // detect nested execution
			void execute( int num_tasks , const Task& task_function )
			{
				if( already_active() )
				{
					WARN( "Nested execution of ThreadPoolIndexedTask is run serially" );
					for( int i=0 ; i<num_tasks ; i++ ) task_function( i );
				}
				else
				{
					std::unique_lock< std::mutex > lock( _mutex );
					_task_function = task_function;
					_num_tasks = num_tasks;
					_num_remaining_tasks = num_tasks;
					_task_index = 0;
					_condition_variable_worker.notify_all();
					_condition_variable_master.wait( lock , [this]{ return !_num_remaining_tasks; } );
				}
			}
			static ThreadPoolIndexedTask& default_threadpool( void )
			{
				static std::unique_ptr< ThreadPoolIndexedTask > thread_pool;
				// This is safe because thread_pool is nullptr only in the main thread before any other thread is launched.
				if( !thread_pool ) thread_pool = std::unique_ptr< ThreadPoolIndexedTask >( new ThreadPoolIndexedTask() );
				return *thread_pool;
			}

		private:
			std::mutex _mutex;
			bool _running = true;
			std::vector< std::thread > _threads;
			Task _task_function;
			int _num_tasks = 0;
			int _num_remaining_tasks = 0;
			int _task_index = 0;
			std::condition_variable _condition_variable_worker;
			std::condition_variable _condition_variable_master;

			void worker_main( void )
			{
				std::unique_lock< std::mutex > lock( _mutex );
				// Consider: https://stackoverflow.com/questions/233127/how-can-i-propagate-exceptions-between-threads
				// However, rethrowing the exception in the main thread loses the stack state, so not useful for debugging.
				for (;;)
				{
					_condition_variable_worker.wait( lock , [this]{ return _task_index < _num_tasks; } );
					if( !_running ) break;
					while( _task_index<_num_tasks )
					{
						int i = _task_index++;
						lock.unlock();
						_task_function(i);
						lock.lock();
						if( _num_remaining_tasks<=0 ) ERROR_OUT( "num remaining tasks not greater than zero" );
						if (!--_num_remaining_tasks) _condition_variable_master.notify_all();
					}
				}
			}
		};
	}
#endif // USE_HH_THREADS

	struct ThreadPool
	{
		enum ParallelType
		{
#ifdef _OPENMP
			OPEN_MP ,
#endif // _OPENMP
			THREAD_POOL ,
#ifdef USE_HH_THREADS
			THREAD_POOL_HH ,
#endif // USE_HH_THREADS
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

		template< typename Function , typename ... Functions >
		static void ParallelSections( const Function &function , const Functions & ... functions )
		{
			std::vector< std::future< void > > futures;
			if constexpr( sizeof ... (Functions) )
			{
				futures.reserve( sizeof...(Functions) );
				_ParallelSections( futures , functions... );
			}
			function();
			for( unsigned int i=0 ; i<futures.size() ; i++ ) futures[i].get();
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

			std::function< void (unsigned int , size_t ) > _ChunkFunction = [ &iterationFunction , begin , end , chunkSize ]( unsigned int thread , size_t chunk )
				{
					const size_t _begin = begin + chunkSize*chunk;
					const size_t _end = std::min< size_t >( end , _begin+chunkSize );
					for( size_t i=_begin ; i<_end ; i++ ) iterationFunction( thread , i );
				};
			std::function< void (unsigned int ) > _StaticThreadFunction = [ &_ChunkFunction , chunks , threads ]( unsigned int thread )
				{
					for( size_t chunk=thread ; chunk<chunks ; chunk+=threads ) _ChunkFunction( thread , chunk );
				};
			std::function< void (unsigned int ) > _DynamicThreadFunction = [ &_ChunkFunction , chunks , &index ]( unsigned int thread )
				{
					size_t chunk;
					while( ( chunk=index.fetch_add(1) )<chunks ) _ChunkFunction( thread , chunk );
				};

			{
				std::lock_guard< std::mutex > lock( _Mutex );

				if     ( schedule==STATIC  ) _ThreadFunction = _StaticThreadFunction;
				else if( schedule==DYNAMIC ) _ThreadFunction = _DynamicThreadFunction;
			}

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
#ifdef USE_HH_THREADS
			else if( _ParallelType==THREAD_POOL_HH )
			{
				const int num_threads = (int)std::min< size_t >( threads , range );
				hh::ThreadPoolIndexedTask* const thread_pool = &hh::ThreadPoolIndexedTask::default_threadpool();
				if( !thread_pool || thread_pool->already_active() )
				{
					// Traverse the range elements sequentially.
					for( size_t index=begin ; index<end ; ++index ) iterationFunction( 0 , index );
				}
				else thread_pool->execute( num_threads , _ThreadFunction );
			}
#endif // USE_HH_THREADS
			else if( _ParallelType==ASYNC )
			{
				static std::vector< std::future< void > > futures;
				futures.resize( threads-1 );
				for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1] = std::async( std::launch::async , _ThreadFunction , t );
				_ThreadFunction( 0 );
				for( unsigned int t=1 ; t<threads ; t++ ) futures[t-1].get();
			}
			else if( _ParallelType==THREAD_POOL )
			{
				unsigned int targetTasks = 0;
#ifdef SANITIZED_PR
				if( !_RemainingTasks.compare_exchange_weak( targetTasks , threads-1 ) )
#else // !SANITIZED_PR
				if( !SetAtomic( _RemainingTasks , threads-1 , targetTasks ) )
#endif // SANITIZED_PR
				{
					WARN( "nested for loop, reverting to serial" );
					for( size_t i=begin ; i<end ; i++ ) iterationFunction( 0 , i );
				}
				else
				{
					_WaitingForWorkOrClose.notify_all();

					// The main thread runs the last batch
					_ThreadFunction( (unsigned int)_Threads.size() );

					{
						std::unique_lock< std::mutex > lock( _Mutex );
						_DoneWithWork.wait( lock , [&]( void ){ return _RemainingTasks==0; } );
					}
				}
			}
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
			if( _ParallelType==THREAD_POOL )
			{
				_RemainingTasks = 0;
				_Close = false;
				for( unsigned int t=0 ; t<numThreads ; t++ ) _Threads[t] = std::thread( _ThreadInitFunction , t );
			}
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

		template< typename Function , typename ... Functions >
		static void _ParallelSections( std::vector< std::future< void > > &futures , const Function &function , const Functions & ... functions )
		{
			futures.push_back( std::async( std::launch::async , function ) );
			if constexpr( sizeof...(Functions) ) _ParallelSections( futures , functions... );
		}

		static void _ThreadInitFunction( unsigned int thread )
		{
			// Wait for the first job to come in
			std::unique_lock< std::mutex > lock( _Mutex );
			_WaitingForWorkOrClose.wait( lock );

			while( !_Close )
			{
				// do the job
				_ThreadFunction( thread );

				// Notify and wait for the next job
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

	const std::vector< std::string > ThreadPool::ParallelNames =
	{
#ifdef _OPENMP
		"open mp" ,
#endif // _OPENMP
		"thread pool" ,
#ifdef USE_HH_THREADS
		"thread pool (hh)" ,
#endif // USE_HH_THREADS
		"async" ,
		"none"
	};
	const std::vector< std::string > ThreadPool::ScheduleNames = { "static" , "dynamic" };
}
#endif // MULTI_THREADING_INCLUDED
