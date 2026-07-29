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

template<class C>
bool ReceiveOnSocket( Socket& s , Pointer( C ) data , size_t dataSize )
{
#ifdef ARRAY_DEBUG
	if( dataSize>data.maximum()*sizeof( C ) ) MK_THROW( "Size of socket read exceeds source maximum: " , dataSize , " > " , data.maximum()*sizeof( C ) );
#endif // ARRAY_DEBUG
	unsigned long long rec=0;
	while( rec!=dataSize )
	{
		int tmp = socket_receive( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec );
		if( tmp<=0 )
		{
			if( !tmp ) MK_THROW( "Connection Closed" );
			else       MK_THROW( "socket_receive from client failed: " , LastSocketError() );
			return false;
		}
		rec+=tmp;
	}
	return true;
}

template<class C>
bool SendOnSocket( Socket& s , ConstPointer( C ) data , size_t dataSize )
{
#ifdef ARRAY_DEBUG
	if( dataSize>data.maximum()*sizeof( C ) ) MK_THROW( "Size of socket write exceeds source maximum: " , dataSize , " > " , data.maximum()*sizeof( C ) );
#endif // ARRAY_DEBUG
	if( socket_send( s , ( ConstPointer( char ) )data , dataSize )<0 )
	{
		MK_THROW( "socket_send to client failed (" , s , "): " , LastSocketError() );
		return false;
	}
	return true;
}

template<class C>
bool SendOnSocket( Socket& s , Pointer( C ) data , size_t dataSize ){ return SendOnSocket( ( ConstPointer( C ) )data , dataSize ); }

template<class C>
void ReceiveOnSocket( Socket& s , Pointer( C ) data , size_t dataSize , const char* errorMessage , ... )
{
#ifdef ARRAY_DEBUG
	if( dataSize>data.maximum()*sizeof( C ) ) MK_THROW( "Size of socket read exceeds source maximum: " , dataSize , " > " , data.maximum()*sizeof( C ) );
#endif // ARRAY_DEBUG
	unsigned long long rec=0;
	while( rec!=dataSize )
	{
		int tmp = socket_receive( s , ( ( Pointer( char ) ) data) + rec , dataSize-rec );
		if( tmp<=0 )
		{
			if( !tmp ) MK_THROW( "Connection Closed" );
			else       MK_THROW( "socket_receive from client failed: " , LastSocketError() );
			{
				fprintf( stderr , "\t" );
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit( EXIT_FAILURE );
		}
		rec+=tmp;
	}
}

template<class C>
void SendOnSocket( Socket& s , ConstPointer( C ) data , size_t dataSize , const char* errorMessage , ... )
{
#ifdef ARRAY_DEBUG
	if( dataSize>data.maximum()*sizeof( C ) ) MK_THROW( "Size of socket write exceeds source maximum: " , dataSize , " > " , data.maximum()*sizeof( C ) );
#endif // ARRAY_DEBUG
	if( socket_send( s , ( ConstPointer( char ) )data , dataSize )<0 )
		MK_THROW( "socket_send to client failed: " , LastSocketError() );
}

template<class C>
void SendOnSocket( Socket& s , Pointer( C ) data , size_t dataSize , const char* errorMessage , ... )
{
#ifdef ARRAY_DEBUG
	if( dataSize>data.maximum()*sizeof( C ) ) MK_THROW( "Size of socket write exceeds source maximum: " , dataSize , " > " , data.maximum()*sizeof( C ) );
#endif // ARRAY_DEBUG
	if( socket_send( s , ( ConstPointer( char ) )data , dataSize )<0 )
		MK_THROW( "socket_send to client failed: " , LastSocketError() );
}

inline bool GetHostEndpointAddress( EndpointAddress* address , const char* prefix )
{
	_InitSockets();
	char hostName[256];
	if( gethostname( hostName , sizeof(hostName) )!=0 ) MK_THROW( "gethostname failed: " , LastSocketError() );

	struct addrinfo hints , *results = NULL;
	memset( &hints , 0 , sizeof(hints) );
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	if( getaddrinfo( hostName , NULL , &hints , &results )!=0 ) return false;

	bool found = false;
	for( struct addrinfo *it=results ; it && !found ; it=it->ai_next )
		if( it->ai_family==AF_INET )
		{
			EndpointAddress _address( ( (struct sockaddr_in *)it->ai_addr )->sin_addr );
			std::string addressString = _address.to_string();
			if( !prefix || strstr( addressString.c_str() , prefix ) )
			{
				*address = _address;
				found = true;
			}
		}
	freeaddrinfo( results );
	return found;
}

inline bool GetHostAddress( char* address , const char* prefix )
{
	EndpointAddress _address;
	if( !GetHostEndpointAddress( &_address , prefix ) ) return false;
	strcpy( address , _address.to_string().c_str() );
	return true;
}

// Local/peer endpoint queries, via getsockname/getpeername.
inline int GetLocalSocketPort( Socket& s )
{
	struct sockaddr_in addr;
	socklen_t len = sizeof(addr);
	if( getsockname( s->fd , (struct sockaddr *)&addr , &len )!=0 ) MK_THROW( "getsockname failed: " , LastSocketError() );
	return (int)ntohs( addr.sin_port );
}

inline EndpointAddress GetLocalSocketEndpointAddress( Socket& s )
{
	struct sockaddr_in addr;
	socklen_t len = sizeof(addr);
	if( getsockname( s->fd , (struct sockaddr *)&addr , &len )!=0 ) MK_THROW( "getsockname failed: " , LastSocketError() );
	return EndpointAddress( addr.sin_addr );
}

inline int GetPeerSocketPort( Socket& s )
{
	struct sockaddr_in addr;
	socklen_t len = sizeof(addr);
	if( getpeername( s->fd , (struct sockaddr *)&addr , &len )!=0 ) MK_THROW( "getpeername failed: " , LastSocketError() );
	return (int)ntohs( addr.sin_port );
}

inline EndpointAddress GetPeerSocketEndpointAddress( Socket& s )
{
	struct sockaddr_in addr;
	socklen_t len = sizeof(addr);
	if( getpeername( s->fd , (struct sockaddr *)&addr , &len )!=0 ) MK_THROW( "getpeername failed: " , LastSocketError() );
	return EndpointAddress( addr.sin_addr );
}

inline Socket GetConnectSocket( const char* address , int port , int ms , bool progress )
{
	_InitSockets();
	char _port[128];
	snprintf( _port , sizeof(_port) , "%d" , port );

	struct addrinfo hints;
	memset( &hints , 0 , sizeof(hints) );
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;

	// Retry until the peer accepts, matching the original blocking behavior.
	// A failed connect() leaves the descriptor unusable, so every attempt gets
	// a freshly created socket.
	long long sleepCount = 0;
	while( true )
	{
		struct addrinfo *results = NULL;
		NativeSocket fd = _INVALID_NATIVE_SOCKET_;
		bool connected = false;

		if( getaddrinfo( address , _port , &hints , &results )==0 )
		{
			for( struct addrinfo *it=results ; it && !connected ; it=it->ai_next )
			{
				fd = socket( it->ai_family , it->ai_socktype , it->ai_protocol );
				if( fd==_INVALID_NATIVE_SOCKET_ ) continue;
				if( connect( fd , it->ai_addr , (int)it->ai_addrlen )==0 ) connected = true;
				else
				{
					_CloseNativeSocket( fd );
					fd = _INVALID_NATIVE_SOCKET_;
				}
			}
			freeaddrinfo( results );
		}

		if( connected )
		{
			if( progress ) printf( "\n" ) , fflush( stdout );
			return new _SocketHolder( fd );
		}

		sleepCount++;
		std::this_thread::sleep_for( std::chrono::milliseconds( 1 ) );
		if( progress && ms>0 && !(sleepCount%ms) ) printf( "." ) , fflush( stdout );
	}
}

inline Socket GetConnectSocket( EndpointAddress address , int port , int ms , bool progress )
{
	return GetConnectSocket( address.to_string().c_str() , port , ms , progress );
}

inline Socket AcceptSocket( AcceptorSocket listen )
{
	if( listen==_INVALID_ACCEPTOR_SOCKET_ ) MK_THROW( "Invalid acceptor socket" );
	NativeSocket fd;
	while( true )
	{
		fd = accept( listen->fd , NULL , NULL );
		if( fd!=_INVALID_NATIVE_SOCKET_ ) break;
		if( !_SocketErrorIsInterrupt() ) MK_THROW( "accept failed: " , LastSocketError() );
	}
	return new _SocketHolder( fd );
}

inline AcceptorSocket GetListenSocket( int &port )
{
	_InitSockets();
	NativeSocket fd = socket( AF_INET , SOCK_STREAM , IPPROTO_TCP );
	if( fd==_INVALID_NATIVE_SOCKET_ ) MK_THROW( "Failed to create listen socket: " , LastSocketError() );

	// Allow rebinding a port left in TIME_WAIT by a previous run (the behavior
	// the acceptor previously provided by default).
	int reuse = 1;
	setsockopt( fd , SOL_SOCKET , SO_REUSEADDR , (const char *)&reuse , sizeof(reuse) );

	struct sockaddr_in addr;
	memset( &addr , 0 , sizeof(addr) );
	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = INADDR_ANY;
	addr.sin_port = htons( (unsigned short)port );

	if( bind( fd , (struct sockaddr *)&addr , sizeof(addr) )!=0 )
	{
		std::string err = LastSocketError();
		_CloseNativeSocket( fd );
		MK_THROW( "Failed to bind listen socket on port " , port , ": " , err );
	}
	if( ::listen( fd , SOMAXCONN )!=0 )
	{
		std::string err = LastSocketError();
		_CloseNativeSocket( fd );
		MK_THROW( "Failed to listen on port " , port , ": " , err );
	}

	// Report the actual port, which matters when 0 was requested.
	socklen_t len = sizeof(addr);
	if( getsockname( fd , (struct sockaddr *)&addr , &len )!=0 )
	{
		std::string err = LastSocketError();
		_CloseNativeSocket( fd );
		MK_THROW( "getsockname failed on listen socket: " , err );
	}
	port = (int)ntohs( addr.sin_port );

	return new _SocketHolder( fd );
}

inline void CloseSocket( Socket& s )
{
	delete s;
	s = _INVALID_SOCKET_;
}

inline void CloseAcceptorSocket( AcceptorSocket& s )
{
	delete s;
	s = _INVALID_ACCEPTOR_SOCKET_;
}
