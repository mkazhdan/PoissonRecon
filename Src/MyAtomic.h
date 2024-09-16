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
#ifndef MY_ATOMIC_INCLUDED
#define MY_ATOMIC_INCLUDED

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

namespace PoissonRecon
{
	template< typename Value > Value ReadAtomic32_( const volatile Value * value );
	template< typename Value > Value ReadAtomic64_( const volatile Value * value );

	template< typename Value > Value SetAtomic32_( volatile Value * value , Value newValue );
	template< typename Value > Value SetAtomic64_( volatile Value * value , Value newValue );

	template< typename Value > bool SetAtomic32_( volatile Value * value , Value newValue , Value oldValue );
	template< typename Value > bool SetAtomic64_( volatile Value * value , Value newValue , Value oldValue );

	template< typename Value > void AddAtomic32_( volatile Value * a , Value b );
	template< typename Value > void AddAtomic64_( volatile Value * a , Value b );

	template< typename Value >
	Value SetAtomic( volatile Value & value , Value newValue )
	{
		if      constexpr( sizeof(Value)==4 ) return SetAtomic32_( &value , newValue );
		else if constexpr( sizeof(Value)==8 ) return SetAtomic64_( &value , newValue );
		else
		{
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex setAtomicMutex;
			std::lock_guard< std::mutex > lock( setAtomicMutex );
			Value oldValue = *(Value*)&value;
			*(Value*)&value = newValue;
			return oldValue;
		}
	}

	template< typename Value >
	bool SetAtomic( volatile Value & value , Value newValue , Value oldValue )
	{
		if      constexpr( sizeof(Value)==4 ) return SetAtomic32_( &value , newValue , oldValue );
		else if constexpr( sizeof(Value)==8 ) return SetAtomic64_( &value , newValue , oldValue );
		else
		{
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex setAtomicMutex;
			std::lock_guard< std::mutex > lock( setAtomicMutex );
			if( value==oldValue ){ value = newValue ; return true; }
			else return false;
		}
	}

	template< typename Value >
	void AddAtomic( volatile Value & a , Value b )
	{
		if      constexpr( sizeof(Value)==4 ) return AddAtomic32_( &a , b );
		else if constexpr( sizeof(Value)==8 ) return AddAtomic64_( &a , b );
		else
		{
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex addAtomicMutex;
			std::lock_guard< std::mutex > lock( addAtomicMutex );
			*(Value*)&a += b;
		}
	}

	template< typename Value >
	Value ReadAtomic( const volatile Value & value )
	{
		if      constexpr( sizeof(Value)==4 ) return ReadAtomic32_( &value );
		else if constexpr( sizeof(Value)==8 ) return ReadAtomic64_( &value );
		else
		{
			WARN_ONCE( "should not use this function: " , sizeof(Value) );
			static std::mutex readAtomicMutex;
			std::lock_guard< std::mutex > lock( readAtomicMutex );
			return *(Value*)&value;
		}
	}

	///////////////////////////////////////////////
	///////////////////////////////////////////////
	///////////////////////////////////////////////

	template< typename Value >
	Value ReadAtomic32_( const volatile Value * value )
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
	Value ReadAtomic64_( const volatile Value * value )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 _value = InterlockedExchangeAdd64( (__int64*)value , 0 );
#else // !_WIN32 && !_WIN64
		uint64_t _value = __atomic_load_n( (uint64_t *)value , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
		return *(Value*)(&_value);
	}

	template< typename Value >
	Value SetAtomic32_( volatile Value *value , Value newValue )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		long *_newValue = (long *)&newValue;
		long oldValue = InterlockedExchange( (long*)value , *_newValue );
#else // !_WIN32 && !_WIN64
		uint32_t *_newValue = (uint32_t *)&newValue;
		long oldValue = __atomic_exchange_n( (uint32_t *)value , *_newValue , __ATOMIC_SEQ_CST );
#endif // _WIN32 || _WIN64
		return *(Value*)&oldValue;
	}

	template< typename Value >
	Value SetAtomic64_( volatile Value * value , Value newValue )
	{
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 *_newValue = (__int64 *)&newValue;
		__int64 oldValue = InterlockedExchange64( (__int64*)value , *_newValue );
#else // !_WIN32 && !_WIN64
		uint64_t *_newValue = (uint64_t *)&newValue;
		uint64_t oldValue = __atomic_exchange_n( (uint64_t *)value , *_newValue , __ATOMIC_SEQ_CST );;
#endif // _WIN32 || _WIN64
		return *(Value*)&oldValue;
	}

	template< typename Value >
	bool SetAtomic32_( volatile Value *value , Value newValue , Value oldValue )
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
	bool SetAtomic64_( volatile Value * value , Value newValue , Value oldValue )
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

	template< typename Value >
	void AddAtomic32_( volatile Value *a , Value b )
	{
#ifdef SANITIZED_PR
		Value current = ReadAtomic32_( a );
#else // !SANITIZED_PR
		Value current = *a;
#endif // SANITIZED_PR
		Value sum = current+b;
#if defined( _WIN32 ) || defined( _WIN64 )
		long *_current = (long *)&current;
		long *_sum = (long *)&sum;
#ifdef SANITIZED_PR
		while( InterlockedCompareExchange( (long*)a , *_sum , *_current )!=*_current )
		{
			current = ReadAtomic32_( a );
			sum = current + b;
		}
#else // !SANITIZED_PR
		while( InterlockedCompareExchange( (long*)a , *_sum , *_current )!=*_current ) current = *(Value*)a , sum = *(Value*)a+b;
#endif // SANITIZED_PR
#else // !_WIN32 && !_WIN64
		uint32_t *_current = (uint32_t *)&current;
		uint32_t *_sum = (uint32_t *)&sum;
#ifdef SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint32_t *)a , *_current , *_sum )!=*_current )
		{
			current = ReadAtomic32_( a );
			sum = current+b;
		}
#else // !SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint32_t *)a , *_current , *_sum )!=*_current ) current = *(Value*)a , sum = *(Value*)a+b;
#endif // SANITIZED_PR
#endif // _WIN32 || _WIN64
	}

	template< typename Value >
	void AddAtomic64_( volatile Value * a , Value b )
	{
#if 1
#ifdef SANITIZED_PR
		Value current = ReadAtomic64_( a );
		Value sum = current+b;
		while( !SetAtomic64_( a , sum , current ) )
		{
			current = ReadAtomic64_( a );
			sum = current+b;
		}
#else // !SANITIZED_PR
		Value current = *a;
		Value sum = current+b;
		while( !SetAtomic64_( a , sum , current ) ) current = *(Value*)a , sum = *(Value*)a+b;
#endif // SANITIZED_PR
#else
#ifdef SANITIZED_PR
		Value current = ReadAtomic64_( a );
#else // !SANITIZED_PR
		Value current = a;
#endif // SANITIZED_PR
		Value sum = current+b;
#if defined( _WIN32 ) || defined( _WIN64 )
		__int64 *_current = (__int64 *)&current;
		__int64 *_sum = (__int64 *)&sum;
#ifdef SANITIZED_PR
		while( InterlockedCompareExchange64( (__int64*)a , *_sum , *_current )!=*_current )
		{
			current = ReadAtomic64_( a );
			sum = current+b;
		}
#else // !SANITIZED_PR
		while( InterlockedCompareExchange64( (__int64*)a , *_sum , *_current )!=*_current ) current = *(Value*)a , sum = a+b;
#endif // SANITIZED_PR
#else // !_WIN32 && !_WIN64
		uint64_t *_current = (uint64_t *)&current;
		uint64_t *_sum = (uint64_t *)&sum;
#ifdef SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint64_t *)a , *_current , *_sum )!=*_current )
		{
			current = ReadAtomic64_( a);
			sum = current+b;
		}
#else // !SANITIZED_PR
		while( __sync_val_compare_and_swap( (uint64_t *)a , *_current , *_sum )!=*_current ) current = *(Value*)a , sum = a+b;
#endif // SANITIZED_PR
#endif // _WIN32 || _WIN64
#endif
	}
}
#endif // MY_ATOMIC_INCLUDED
