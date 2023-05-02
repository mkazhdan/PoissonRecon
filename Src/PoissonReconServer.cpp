/*
Copyright (c) 2023, Michael Kazhdan
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

#define BIG_DATA
#include "PreProcessor.h"

#include "Socket.h"
#include "PoissonReconClientServer.h"
#include "PointPartitionClientServer.h"
#include "MergePlyClientServer.h"

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PoissonRecon.h"


cmdLineParameter< std::string >
	AddressPrefix( "prefix" ) ,
	In( "in" ) ,
	TempDir( "tempDir" ) ,
	Out( "out" );

cmdLineParameter< int >
	ClientCount( "count" ) ,
	Port( "port" , 0 ) ,
	Verbose( "verbose" , 0 ) ,
#ifdef FAST_COMPILE
#else // !FAST_COMPILE
	Degree( "degree" , DEFAULT_FEM_DEGREE ) ,
	BType( "bType" , DEFAULT_FEM_BOUNDARY ) ,
#endif // FAST_COMPILE
	Iters( "iters" , 8 ) ,
	BaseVCycles( "vCycles" , 1 ) ,
	KernelDepth( "kernelDepth" ) ,
	BaseDepth( "baseDepth" , 5 ) ,
	SolveDepth( "solveDepth" ) ,
	PadSize( "pad" , 4 ) ,
	BufferSize( "buffer" , BUFFER_IO ) ,
	Depth( "depth" , 8 ) ,
	PartitionDepth( "pDepth" , 5 ) ,
	FilesPerDir( "filesPerDir" , -1 ) ,
	MaxMemoryGB( "maxMemory" , 0 ) ,
	PeakMemorySampleMS( "sampleMS" , 10 ) ,
#ifdef _OPENMP
	ParallelType( "parallel" , (int)ThreadPool::OPEN_MP ) ,
#else // !_OPENMP
	ParallelType( "parallel" , (int)ThreadPool::THREAD_POOL ) ,
#endif // _OPENMP
	ScheduleType( "schedule" , (int)ThreadPool::DefaultSchedule ) ,
	ThreadChunkSize( "chunkSize" , (int)ThreadPool::DefaultChunkSize ) ,
	Threads( "threads" , (int)std::thread::hardware_concurrency() ) ;

cmdLineReadable
	Performance( "performance" ) ,
	NoLoadBalance( "noLoadBalance" ) ,
	Density( "density" ) ,
	LinearFit( "linearFit" ) ,
	OutputVoxelGrid( "grid" ) ,
	OutputBoundarySlices( "boundary" ) ,
	NoFuse( "noFuse" ) ,
	ShowDiscontinuity( "showDiscontinuity" );

cmdLineParameter< float >
	Scale( "scale" , 1.1f ) ,
	Width( "width" , 0.f ) ,
	Confidence( "confidence" , 0.f ) ,
	ConfidenceBias( "confidenceBias" , 0.f ) ,
	SamplesPerNode( "samplesPerNode" , 1.5f ) ,
	DataX( "data" , 32.f ) ,
	PointWeight( "pointWeight" ) ,
	CGSolverAccuracy( "cgAccuracy" , 1e-3f );

cmdLineReadable* params[] =
{
	&Port , &ClientCount , &AddressPrefix , &Performance , &Verbose ,
	&In ,
	&Scale ,
	&PartitionDepth ,
	&FilesPerDir ,
	&TempDir ,
	&Out ,
#ifdef FAST_COMPILE
#else // !FAST_COMPILE
	&Degree , &BType ,
#endif // FAST_COMPILE
	&Iters , &BaseVCycles , &KernelDepth , &BaseDepth , &PadSize , &BufferSize , &Depth ,
	&SolveDepth ,
	&NoLoadBalance , &Density , &LinearFit ,
	&NoFuse ,
	&Width , &Confidence , &ConfidenceBias , &SamplesPerNode , &DataX , &PointWeight , &CGSolverAccuracy ,
	&MaxMemoryGB , &ParallelType , &ScheduleType , &ThreadChunkSize , &Threads ,
	&PeakMemorySampleMS ,
	&OutputVoxelGrid ,
	&OutputBoundarySlices ,
	&ShowDiscontinuity ,
	NULL
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input points>\n" , In.name );
	printf( "\t --%s <networked temporary directory>\n" , TempDir.name );
	printf( "\t --%s <output polygon mesh>\n" , Out.name );
	printf( "\t --%s <client count>\n" , ClientCount.name );

	printf( "\t[--%s <preferred address prefix>]\n" , AddressPrefix.name );
	printf( "\t[--%s <listen port>=%d]\n" , Port.name , Port.value );

#ifdef FAST_COMPILE
#else // !FAST_COMPILE
	printf( "\t[--%s <b-spline degree>=%d]\n" , Degree.name , Degree.value );
	printf( "\t[--%s <boundary type>=%d]\n" , BType.name , BType.value );
	for( int i=0 ; i<BOUNDARY_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , BoundaryNames[i] );
#endif // FAST_COMPILE
	printf( "\t[--%s <minimum number of samples per node>=%f]\n" , SamplesPerNode.name, SamplesPerNode.value );
	printf( "\t[--%s <base depth>=%d]\n" , BaseDepth.name , BaseDepth.value );
	printf( "\t[--%s <max reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s <partition depth>=%d]\n" , PartitionDepth.name , PartitionDepth.value );
	printf( "\t[--%s <kernel depth>]\n" , KernelDepth.name );
	printf( "\t[--%s <solver depth>]\n" , SolveDepth.name );
	printf( "\t[--%s <scale factor>=%f]\n" , Scale.name , Scale.value );
	printf( "\t[--%s <grid width>]\n" , Width.name );
	printf( "\t[--%s <iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s <base MG solver v-cycles>=%d]\n" , BaseVCycles.name , BaseVCycles.value );
	printf( "\t[--%s <cg solver accuracy>=%g]\n" , CGSolverAccuracy.name , CGSolverAccuracy.value );
	printf( "\t[--%s <interpolation weight>=%.3e * <b-spline degree>]\n" , PointWeight.name , DefaultPointWeightMultiplier );
	printf( "\t[--%s <normal confidence exponent>=%f]\n" , Confidence.name , Confidence.value );
	printf( "\t[--%s <normal confidence bias exponent>=%f]\n" , ConfidenceBias.name , ConfidenceBias.value );
	printf( "\t[--%s <pull factor>=%f]\n" , DataX.name , DataX.value );
	printf( "\t[--%s <pad size>=%d]\n" , PadSize.name , PadSize.value );
	printf( "\t[--%s <buffer size>=%d]\n" , BufferSize.name , BufferSize.value );
	printf( "\t[--%s <num threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s <parallel type>=%d]\n" , ParallelType.name , ParallelType.value );
	for( size_t i=0 ; i<ThreadPool::ParallelNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ParallelNames[i].c_str() );
	printf( "\t[--%s <schedue type>=%d]\n" , ScheduleType.name , ScheduleType.value );
	for( size_t i=0 ; i<ThreadPool::ScheduleNames.size() ; i++ ) printf( "\t\t%d] %s\n" , (int)i , ThreadPool::ScheduleNames[i].c_str() );
	printf( "\t[--%s <thread chunk size>=%d]\n" , ThreadChunkSize.name , ThreadChunkSize.value );
	printf( "\t[--%s <peak memory sampling rate (ms)>=%d]\n" , PeakMemorySampleMS.name , PeakMemorySampleMS.value );
	printf( "\t[--%s <maximum memory (in GB)>=%d]\n" , MaxMemoryGB.name , MaxMemoryGB.value );
	printf( "\t[--%s <slab files per directory>=%u]\n" , FilesPerDir.name , (unsigned int)FilesPerDir.value );
	printf( "\t[--%s <verbosity>=%d]\n" , Verbose.name , Verbose.value );
	printf( "\t[--%s]\n" , NoLoadBalance.name );
	printf( "\t[--%s]\n" , Density.name );
	printf( "\t[--%s]\n" , LinearFit.name );
	printf( "\t[--%s]\n" , OutputVoxelGrid.name );
	printf( "\t[--%s]\n" , OutputBoundarySlices.name );
	printf( "\t[--%s]\n" , ShowDiscontinuity.name );
	printf( "\t[--%s]\n" , NoFuse.name );

	printf( "\t[--%s]\n" , Performance.name );
}

template< typename Real , unsigned int Dim >
std::pair< PointPartition::PointSetInfo< Real , Dim > , PointPartition::Partition >
Partition
(
	std::vector< Socket > &clientSockets ,
	const PointPartitionClientServer::ClientPartitionInfo< Real > &clientPartitionInfo ,
	bool loadBalance ,
	bool performance
)
{
	Timer timer;

	std::pair< PointPartition::PointSetInfo< Real , Dim > , PointPartition::Partition > pointSetInfoAndPartition = PointPartitionClientServer::RunServer< Real , Dim >( clientSockets , clientPartitionInfo , loadBalance );
	unsigned int peakMem = 0;
	for( unsigned int c=0 ; c<clientSockets.size() ; c++ )
	{
		unsigned int _peakMem;
		SocketStream( clientSockets[c] ).read( _peakMem );
		peakMem = std::max< unsigned int >( peakMem , _peakMem );
	}

	if( performance )
	{
		StreamFloatPrecision sfp( std::cout , 2 );
		std::cout << "Partition performance:" << std::endl;
		std::cout << "\tPeak client memory: " << peakMem << " (MB)" << std::endl;
		std::cout << "\tPeak server memory: " << MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
		std::cout << "\tServer Time: " << timer.wallTime() << " (s)" << std::endl;
	}
	return pointSetInfoAndPartition;
}

template< typename Real , unsigned int Dim , BoundaryType BType , unsigned int Degree >
std::vector< unsigned int >
Reconstruct
(
	const PointPartition::PointSetInfo< Real , Dim > &pointSetInfo ,
	const PointPartition::Partition &pointPartition ,
	std::vector< Socket > &clientSockets ,
	const PoissonReconClientServer::ClientReconstructionInfo< Real , Dim > &clientReconInfo
)
{
	Timer timer;

	// Clean up files if they were not propertly cleaned up before.
	for( unsigned int i=0 ; i<clientSockets.size() ; i++ ) PoissonReconClientServer::ClientServerStream< false >::Reset( i , clientReconInfo );

	std::vector< unsigned int > sharedVertexCounts = PoissonReconClientServer::RunServer< Real , Dim , BType , Degree >( pointSetInfo , pointPartition , clientSockets , clientReconInfo , BaseVCycles.value , PeakMemorySampleMS.value<0 ? 0 : PeakMemorySampleMS.value , ShowDiscontinuity.set , OutputBoundarySlices.set );

	unsigned int peakMem = 0;
	for( unsigned int i=0 ; i<clientSockets.size() ; i++ )
	{
		unsigned int _peakMem;
		SocketStream( clientSockets[i] ).read( _peakMem );
		peakMem = std::max< unsigned int >( peakMem , _peakMem );
	}

	if( Performance.set )
	{
		StreamFloatPrecision sfp( std::cout , 2 );
		std::cout << "Reconstruction performance:" << std::endl;
		std::cout << "\tPeak client memory: " << peakMem << " (MB)" << std::endl;
		std::cout << "\tPeak server memory: " << MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
		std::cout << "\tServer Time: " << timer.wallTime() << " (s)" << std::endl;
	}

	return sharedVertexCounts;
}


#ifdef FAST_COMPILE
#else // !FAST_COMPILE
template< typename Real , unsigned int Dim , BoundaryType BType >
std::vector< unsigned int > Reconstruct( unsigned int degree , const PointPartition::PointSetInfo< Real , Dim > &pointSetInfo , const PointPartition::Partition &partition , std::vector< Socket > clientSockets , const PoissonReconClientServer::ClientReconstructionInfo< Real , Dim > &clientReconInfo )
{
	switch( degree )
	{
		case 1: return Reconstruct< Real , Dim , BType , 1 >( pointSetInfo , partition , clientSockets , clientReconInfo );
		case 2: return Reconstruct< Real , Dim , BType , 2 >( pointSetInfo , partition , clientSockets , clientReconInfo );
		default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
	}
	return std::vector< unsigned int >();
}

template< typename Real , unsigned int Dim >
std::vector< unsigned int > Reconstruct( BoundaryType bType , unsigned int degree , const PointPartition::PointSetInfo< Real , Dim > &pointSetInfo , const PointPartition::Partition &partition , std::vector< Socket > clientSockets , const PoissonReconClientServer::ClientReconstructionInfo< Real , Dim > &clientReconInfo )
{
	for( unsigned int i=0 ; i<clientSockets.size() ; i++ )
	{
		SocketStream( clientSockets[i] ).write( bType );
		SocketStream( clientSockets[i] ).write( degree );
	}
	switch( bType )
	{
		case BOUNDARY_FREE:      return Reconstruct< Real , Dim , BOUNDARY_FREE      >( degree , pointSetInfo , partition , clientSockets , clientReconInfo );
		case BOUNDARY_NEUMANN:   return Reconstruct< Real , Dim , BOUNDARY_NEUMANN   >( degree , pointSetInfo , partition , clientSockets , clientReconInfo );
		case BOUNDARY_DIRICHLET: return Reconstruct< Real , Dim , BOUNDARY_DIRICHLET >( degree , pointSetInfo , partition , clientSockets , clientReconInfo );
		default: ERROR_OUT( "Not a valid boundary type: " , bType );
	}
	return std::vector< unsigned int >();
}
#endif // FAST_COMPILE

template< typename Real , unsigned int Dim>
void Merge
(
	const std::vector< unsigned int > &sharedVertexCounts ,
	std::string header ,
	std::vector< Socket > &clientSockets ,
	const MergePlyClientServer::ClientMergePlyInfo &clientMergePlyInfo
)
{
	Timer timer;

	MergePlyClientServer::RunServer< Real , Dim >( TempDir.value , TempDir.value , header , Out.value , clientSockets , sharedVertexCounts , clientMergePlyInfo , PeakMemorySampleMS.value );

	unsigned int peakMem = 0;
	for( unsigned int i=0 ; i<clientSockets.size() ; i++ )
	{
		unsigned int _peakMem;
		SocketStream( clientSockets[i] ).read( _peakMem );
		peakMem = std::max< unsigned int >( peakMem , _peakMem );
	}

	if( Performance.set )
	{
		StreamFloatPrecision sfp( std::cout , 2 );
		std::cout << "Merge performance:" << std::endl;
		std::cout << "\tPeak client memory: " << peakMem << " (MB)" << std::endl;
		std::cout << "\tPeak server memory: " << MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
		std::cout << "\tServer Time: " << timer.wallTime() << " (s)" << std::endl;
	}
}

int main( int argc , char* argv[] )
{
#ifdef ARRAY_DEBUG
	WARN( "Array debugging enabled" );
#endif // ARRAY_DEBUG
#ifdef USE_DOUBLE
	typedef double Real;
#else // !USE_DOUBLE
	typedef float  Real;
#endif // USE_DOUBLE
	static const unsigned int Dim = DEFAULT_DIMENSION;

	cmdLineParse( argc-1 , &argv[1] , params );

	if( !In.set || !TempDir.set || !Out.set || !ClientCount.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}
	if( PadSize.value<0 ) ERROR_OUT( "Padding size cannot be negative" );

	if( Verbose.value>1 )
	{
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "** Running Poisson Reconstruction Server (Version " << VERSION << ") **" << std::endl;
		std::cout << "***********************************************************" << std::endl;
		std::cout << "***********************************************************" << std::endl;
	}

#ifdef FAST_COMPILE
	if( !PointWeight.set ) PointWeight.value = DefaultPointWeightMultiplier * DEFAULT_FEM_DEGREE;
#else // !FAST_COMPILE
	if( !PointWeight.set ) PointWeight.value = DefaultPointWeightMultiplier * Degree.value;
#endif // FAST_COMPILE

	if( Depth.set && Width.value>0 )
	{
		WARN( "Both --" , Depth.name  , " and --" , Width.name , " set, ignoring --" , Width.name );
		Width.value = 0;
	}

	if( MaxMemoryGB.value>0 ) SetPeakMemoryMB( MaxMemoryGB.value<<10 );
	ThreadPool::DefaultChunkSize = ThreadChunkSize.value;
	ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)ScheduleType.value;
	ThreadPool::Init( (ThreadPool::ParallelType)ParallelType.value , Threads.value );
	if( !Threads.set && Verbose.value>1 ) std::cout << "Running with " << Threads.value << " threads" << std::endl;
	std::string header;

	// Create the connections to the clients
	std::vector< Socket > clientSockets( ClientCount.value );
	{
		int port = Port.value;

		char address[512];
		GetHostAddress( address , AddressPrefix.value.c_str() );

		// Create a listening SOCKET for connecting to server
		AcceptorSocket listenSocket = GetListenSocket( port );
		if( listenSocket == _INVALID_ACCEPTOR_SOCKET_ ) ERROR_OUT( "Could not create listener socket" );
		std::cout << "Server Address: " << address << ":" << port << std::endl;
		{
			std::stringstream ss;
			ss << "PR_" << port;
			header = ss.str();
		}

		// Establish a connection to the clients
		for( unsigned int i=0 ; i<(unsigned int)ClientCount.value ; i++ )
		{
			clientSockets[i] = AcceptSocket( listenSocket );
			if( Verbose.value>1 ) std::cout << "Connected to process: " << (i+1) << " / " << ClientCount.value << std::endl;
		}
		CloseAcceptorSocket( listenSocket );
	}
	std::vector< unsigned int > sharedVertexCounts;


	std::pair< PointPartition::PointSetInfo< Real , Dim > , PointPartition::Partition > pointSetInfoAndPartition;

	// Get the partitioned points
	{
		PointPartition::CreatePointSlabDirs( PointPartition::FileDir( TempDir.value , header ) , 1<<PartitionDepth.value , FilesPerDir.value );
		PointPartitionClientServer::ClientPartitionInfo< Real > clientPartitionInfo;
		clientPartitionInfo.in = In.value;
		clientPartitionInfo.tempDir = TempDir.value;
		clientPartitionInfo.outDir = TempDir.value;
		clientPartitionInfo.outHeader = header;
		clientPartitionInfo.slabs = 1<<PartitionDepth.value;
		clientPartitionInfo.filesPerDir = FilesPerDir.value;
		clientPartitionInfo.bufferSize = BufferSize.value;
		clientPartitionInfo.scale = Scale.value;
		clientPartitionInfo.verbose = Verbose.value>1;
		pointSetInfoAndPartition = Partition< Real , Dim >( clientSockets , clientPartitionInfo , !NoLoadBalance.set , Performance.set );
	}

	// Reconstruct the slabs
	{
		if( Width.value>0 )
		{
			Real maxScale = 0;
			for( unsigned int i=0 ; i<Dim ; i++ ) maxScale = std::max< Real >( maxScale , (Real)1./pointSetInfoAndPartition.first.modelToUnitCube(i,i) );
			Depth.value = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/Width.value )/log(2.) ) );
		}
		PoissonReconClientServer::ClientReconstructionInfo< Real , Dim > clientReconInfo;
		clientReconInfo.inDir = TempDir.value;
		clientReconInfo.tempDir = TempDir.value;
		clientReconInfo.outDir = TempDir.value;
		clientReconInfo.header = header;
		clientReconInfo.bufferSize = BufferSize.value;
		clientReconInfo.iters = Iters.value;
		clientReconInfo.pointWeight = PointWeight.value;
		clientReconInfo.confidence = Confidence.value;
		clientReconInfo.confidenceBias = ConfidenceBias.value;
		clientReconInfo.kernelDepth = KernelDepth.value;
		clientReconInfo.solveDepth = SolveDepth.value;
		clientReconInfo.samplesPerNode = SamplesPerNode.value;
		clientReconInfo.dataX = DataX.value;
		clientReconInfo.density = Density.set;
		clientReconInfo.linearFit = LinearFit.set;
		PoissonReconClientServer::ClientReconstructionInfo< Real , Dim >::MergeType::TOPOLOGY_AND_FUNCTION;
		clientReconInfo.ouputVoxelGrid = OutputVoxelGrid.set;
		clientReconInfo.verbose = Verbose.value;
		clientReconInfo.filesPerDir = FilesPerDir.value;
		clientReconInfo.padSize = PadSize.value;

		clientReconInfo.reconstructionDepth = Depth.value;
		clientReconInfo.sharedDepth = 0;
		while( ((size_t)1<<clientReconInfo.sharedDepth) < pointSetInfoAndPartition.first.pointsPerSlab.size() ) clientReconInfo.sharedDepth++;
		if( ((size_t)1<<clientReconInfo.sharedDepth)!=pointSetInfoAndPartition.first.pointsPerSlab.size() ) ERROR_OUT( "Number of point slabs is not a power of two: " , pointSetInfoAndPartition.first.pointsPerSlab.size() );
		clientReconInfo.baseDepth = BaseDepth.value;

		if( clientReconInfo.pointWeight<0 ) ERROR_OUT( "Expected non-negative point-weight" );

		if( clientReconInfo.sharedDepth>clientReconInfo.reconstructionDepth ) ERROR_OUT( "Slab depth cannot exceed reconstruction depth: " , clientReconInfo.sharedDepth , " <= "  , clientReconInfo.reconstructionDepth );
		if( clientReconInfo.baseDepth>clientReconInfo.sharedDepth )
		{
			if( BaseDepth.set ) ERROR_OUT( "Base depth cannot exceed shared depth: " , clientReconInfo.baseDepth , " <="  , clientReconInfo.sharedDepth );
			else clientReconInfo.baseDepth = clientReconInfo.sharedDepth;
		}
		if( !KernelDepth.set ) KernelDepth.value = clientReconInfo.reconstructionDepth-2;
		clientReconInfo.kernelDepth = KernelDepth.value;

		if( clientReconInfo.kernelDepth>clientReconInfo.reconstructionDepth )
		{
			WARN( "Kernel depth should not exceed depth: " , clientReconInfo.kernelDepth , " <= " , clientReconInfo.reconstructionDepth );
			clientReconInfo.kernelDepth = clientReconInfo.reconstructionDepth;
		}

		clientReconInfo.solveDepth = SolveDepth.set ? SolveDepth.value : clientReconInfo.reconstructionDepth;
		if( clientReconInfo.solveDepth>clientReconInfo.reconstructionDepth )
		{
			WARN( "Solve depth cannot exceed reconstruction depth: " , clientReconInfo.solveDepth , " <= " , clientReconInfo.reconstructionDepth );
			clientReconInfo.solveDepth = clientReconInfo.reconstructionDepth;
		}
		if( clientReconInfo.solveDepth<clientReconInfo.baseDepth )
		{
			WARN( "Solve depth cannot be smaller than base depth: " , clientReconInfo.solveDepth , " >= " , clientReconInfo.baseDepth );
			clientReconInfo.solveDepth = clientReconInfo.baseDepth;
		}
#ifdef FAST_COMPILE
		sharedVertexCounts = Reconstruct< Real , Dim , DEFAULT_FEM_BOUNDARY , DEFAULT_FEM_DEGREE >( pointSetInfoAndPartition.first , pointSetInfoAndPartition.second , clientSockets , clientReconInfo );
#else // !FAST_COMPILE
		sharedVertexCounts = Reconstruct< Real , Dim >( (BoundaryType)BType.value , Degree.value , pointSetInfoAndPartition.first , pointSetInfoAndPartition.second , clientSockets , clientReconInfo );
#endif // FAST_COMPILE
	}

	if( Verbose.value>1 && !NoFuse.set ) 
		for( unsigned int i=0 ; i<sharedVertexCounts.size() ; i++ ) std::cout << "Vertices[" << (i+1) << "] " << sharedVertexCounts[i] << std::endl;

	PointPartition::RemovePointSlabDirs( PointPartition::FileDir( TempDir.value , header ) );

	{
		if( NoFuse.set ) for( unsigned int i=0 ; i<sharedVertexCounts.size() ; i++ ) sharedVertexCounts[i] = 0;
		MergePlyClientServer::ClientMergePlyInfo clientMergePlyInfo;
		clientMergePlyInfo.auxProperties = pointSetInfoAndPartition.first.auxiliaryProperties;
		clientMergePlyInfo.bufferSize = BufferSize.value;
		clientMergePlyInfo.verbose = Verbose.value!=0;
		Merge< Real , Dim >( sharedVertexCounts , header , clientSockets , clientMergePlyInfo );

		auto InFile = [&]( unsigned int idx )
		{
			std::stringstream ss;
			ss << header << "." << idx << ".ply";
			return PointPartition::FileDir( TempDir.value , ss.str() );
		};

		for( unsigned int i=0 ; i<clientSockets.size() ; i++ )
		{
			std::string fileName = InFile(i);
			std::remove( fileName.c_str() );
		}
	}

	for( unsigned int i=0 ; i<clientSockets.size() ; i++ ) CloseSocket( clientSockets[i] );

	ThreadPool::Terminate();

	return EXIT_SUCCESS;
}
