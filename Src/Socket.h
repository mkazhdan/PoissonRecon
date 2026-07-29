/*
Copyright (c) 2008, Michael Kazhdan
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

#ifndef SOCKET_INCLUDED
#define SOCKET_INCLUDED

#ifdef _WIN32
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif // _WIN32_WINNT
#endif // _WIN32

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment( lib , "ws2_32.lib" )
#else // !_WIN32
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <errno.h>
#endif // _WIN32

#include <iostream>
#include <string>
#include <string.h>
#include <stdarg.h>
#include <thread>
#include "Array.h"
#include "MyMiscellany.h"
#include "Streams.h"

namespace PoissonRecon
{

	static const unsigned int SOCKET_CONNECT_WAIT = 500;		// Default time to wait on a socket (in ms)

	////////////////////////////////////////////////////////////////////////////
	// Platform abstraction over the native socket API.
	//
	// Only blocking, synchronous TCP is used, so the native BSD/Winsock calls
	// are sufficient; there is no dependency on an external networking library.
	////////////////////////////////////////////////////////////////////////////
#ifdef _WIN32
	typedef SOCKET NativeSocket;
	static const NativeSocket _INVALID_NATIVE_SOCKET_ = INVALID_SOCKET;
	inline int  _CloseNativeSocket( NativeSocket s ){ return closesocket( s ); }
	inline int  _LastSocketErrorCode( void ){ return WSAGetLastError(); }
	inline bool _SocketErrorIsInterrupt( void ){ return WSAGetLastError()==WSAEINTR; }
#else // !_WIN32
	typedef int NativeSocket;
	static const NativeSocket _INVALID_NATIVE_SOCKET_ = -1;
	inline int  _CloseNativeSocket( NativeSocket s ){ return close( s ); }
	inline int  _LastSocketErrorCode( void ){ return errno; }
	inline bool _SocketErrorIsInterrupt( void ){ return errno==EINTR; }
#endif // _WIN32

	// Winsock requires per-process initialization. The function-local static
	// makes this happen exactly once, and is thread-safe under C++11.
	inline void _InitSockets( void )
	{
#ifdef _WIN32
		struct Initializer
		{
			Initializer( void )
			{
				WSADATA wsaData;
				if( WSAStartup( MAKEWORD( 2 , 2 ) , &wsaData )!=0 ) MK_THROW( "WSAStartup failed" );
			}
			~Initializer( void ){ WSACleanup(); }
		};
		static Initializer initializer;
#endif // _WIN32
	}

	// Sockets are handed around as pointers (and compared against NULL), so the
	// native descriptor is wrapped in a heap-allocated holder that closes on
	// destruction -- matching the ownership semantics the callers already use.
	struct _SocketHolder
	{
		NativeSocket fd;
		_SocketHolder( NativeSocket fd=_INVALID_NATIVE_SOCKET_ ) : fd(fd){}
		~_SocketHolder( void ){ if( fd!=_INVALID_NATIVE_SOCKET_ ) _CloseNativeSocket( fd ); }
		_SocketHolder( const _SocketHolder & ) = delete;
		_SocketHolder &operator = ( const _SocketHolder & ) = delete;
	};

	typedef _SocketHolder *Socket;
	typedef _SocketHolder *AcceptorSocket;
	const Socket _INVALID_SOCKET_ = (Socket)NULL;
	const AcceptorSocket _INVALID_ACCEPTOR_SOCKET_ = (AcceptorSocket)NULL;

	// IPv4 address value-type. (The original implementation likewise only ever
	// surfaced v4 addresses.)
	struct EndpointAddress
	{
		struct in_addr addr;

		EndpointAddress( void ){ memset( &addr , 0 , sizeof(addr) ); }
		explicit EndpointAddress( struct in_addr a ) : addr(a){}

		std::string to_string( void ) const
		{
			char buffer[ INET_ADDRSTRLEN ];
			if( !inet_ntop( AF_INET , (void *)&addr , buffer , INET_ADDRSTRLEN ) ) return std::string();
			return std::string( buffer );
		}
	};

	inline const char *LastSocketError( void )
	{
		static thread_local char buffer[512];
		int err = _LastSocketErrorCode();
#ifdef _WIN32
		if( !FormatMessageA( FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS , NULL , err , 0 , buffer , sizeof(buffer) , NULL ) )
			snprintf( buffer , sizeof(buffer) , "socket error %d" , err );
		else
		{
			// Trim the trailing CR/LF FormatMessage appends.
			size_t len = strlen( buffer );
			while( len && ( buffer[len-1]=='\n' || buffer[len-1]=='\r' ) ) buffer[--len] = 0;
		}
#else // !_WIN32
		snprintf( buffer , sizeof(buffer) , "%s" , strerror( err ) );
#endif // _WIN32
		return buffer;
	}

	// recv/send may transfer fewer bytes than requested, so both loop until the
	// full buffer has been moved -- the behavior the callers rely on.
	template< class C > int socket_receive( Socket &s , C *destination , size_t len )
	{
		if( s==_INVALID_SOCKET_ ) MK_THROW( "Failed to read from socket: invalid socket" );
		char *ptr = (char *)destination;
		size_t transferred = 0;
		while( transferred<len )
		{
			size_t remaining = len - transferred;
			int chunk = (int)( remaining>(size_t)(1<<30) ? (1<<30) : remaining );
			int ret = (int)recv( s->fd , ptr+transferred , chunk , 0 );
			if( ret>0 ) transferred += (size_t)ret;
			else if( ret==0 ) MK_THROW( "Failed to read from socket: connection closed" );
			else if( !_SocketErrorIsInterrupt() ) MK_THROW( "Failed to read from socket: " , LastSocketError() );
		}
		return (int)transferred;
	}

	template< class C > int socket_send( Socket& s , const C* source , size_t len )
	{
		if( s==_INVALID_SOCKET_ ) MK_THROW( "Failed to write to socket: invalid socket" );
		const char *ptr = (const char *)source;
		size_t transferred = 0;
		while( transferred<len )
		{
			size_t remaining = len - transferred;
			int chunk = (int)( remaining>(size_t)(1<<30) ? (1<<30) : remaining );
			int ret = (int)send( s->fd , ptr+transferred , chunk , 0 );
			if( ret>0 ) transferred += (size_t)ret;
			else if( !_SocketErrorIsInterrupt() ) MK_THROW( "Failed to write to socket: " , LastSocketError() );
		}
		return (int)transferred;
	}

	inline bool AddressesEqual( const EndpointAddress& a1 ,  const EndpointAddress& a2 ){ return a1.to_string()==a2.to_string(); }

	inline void PrintHostAddresses( FILE* fp )
	{
		_InitSockets();
		char hostName[256];
		if( gethostname( hostName , sizeof(hostName) )!=0 ) MK_THROW( "gethostname failed: " , LastSocketError() );

		struct addrinfo hints , *results = NULL;
		memset( &hints , 0 , sizeof(hints) );
		hints.ai_family = AF_INET;
		hints.ai_socktype = SOCK_STREAM;
		if( getaddrinfo( hostName , NULL , &hints , &results )!=0 ) MK_THROW( "getaddrinfo failed for host: " , hostName );

		int count = 0;
		for( struct addrinfo *it=results ; it ; it=it->ai_next )
			if( it->ai_family==AF_INET )
			{
				EndpointAddress address( ( (struct sockaddr_in *)it->ai_addr )->sin_addr );
				fprintf( fp , "%d]  %s\n" , count++ , address.to_string().c_str() );
			}
		freeaddrinfo( results );
	}

#ifdef ARRAY_DEBUG
	template< class C >
	int socket_receive( Socket& s , Array< C > destination , size_t len )
	{
		if( len>destination.maximum()*sizeof( C ) )
			MK_THROW( "Size of socket_receive exceeds destination maximum: " , len , " > " , destination.maximum()*sizeof( C ) );
		return socket_receive( s , (char*)&destination[0] , len );
	}
	template< class C >
	int socket_send( Socket s , ConstArray< C > source , size_t len )
	{
		if( len>source.maximum()*sizeof( C ) )
			MK_THROW( "Size of socket_send exceeds source maximum: " , len , " > " , source.maximum()*sizeof( C ) );
		return socket_send( s , (char*)&source[0] , len );
	}
#endif // ARRAY_DEBUG

	class ConnectionData
	{
	public:
		EndpointAddress localAddr , peerAddr;
		int localPort , peerPort;
	};


	template< class C > bool ReceiveOnSocket( Socket &s ,      Pointer( C ) data , size_t dataSize );
	template< class C > bool SendOnSocket   ( Socket &s , ConstPointer( C ) data , size_t dataSize );
	template< class C > bool SendOnSocket   ( Socket &s ,      Pointer( C ) data , size_t dataSize );
	template< class C > void ReceiveOnSocket( Socket &s ,      Pointer( C ) data , size_t dataSize , const char *errorMessage , ... );
	template< class C > void SendOnSocket   ( Socket &s , ConstPointer( C ) data , size_t dataSize , const char *errorMessage , ... );
	template< class C > void SendOnSocket   ( Socket &s ,      Pointer( C ) data , size_t dataSize , const char *errorMessage , ... );

	AcceptorSocket GetListenSocket( int& port );
	Socket AcceptSocket( AcceptorSocket listen );
	Socket GetConnectSocket( const char* address , int port , int ms=5 , bool progress=false );
	Socket GetConnectSocket( EndpointAddress , int port , int ms=5 , bool progress=false );
	void CloseSocket( Socket& s );
	void CloseAcceptorSocket( AcceptorSocket& s );
	EndpointAddress GetLocalSocketEndpointAddress( Socket& s );
	int             GetLocalSocketPort           ( Socket& s );
	EndpointAddress GetLocalSocketEndpointAddress( Socket& s );
	int             GetPeerSocketPort            ( Socket& s );
	bool GetHostAddress( char* address , const char* prefix = NULL );
	bool GetHostEndpointAddress( EndpointAddress* address , const char* prefix=NULL );
	void PrintHostAddress( void );

	struct SocketStream : public BinaryStream
	{
		SocketStream( Socket socket=_INVALID_SOCKET_ ) : _socket(socket){}
	protected:
		Socket _socket;
		bool  _read(      Pointer( unsigned char ) ptr , size_t sz ){ return socket_receive( _socket , ptr , sizeof(unsigned char)*sz )==sz; }
		bool _write( ConstPointer( unsigned char ) ptr , size_t sz ){ return socket_send   ( _socket , ptr , sizeof(unsigned char)*sz )==sz; }
	};


#include "Socket.inl"
}
#endif // SOCKET_INCLUDED
