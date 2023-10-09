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

#include "PreProcessor.h"
#include "Reconstructors.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include "MyMiscellany.h"
#include "CmdLineParser.h"

cmdLineParameter< char* > Out( "out" );
cmdLineReadable SSDReconstruction( "ssd" ) , UseColor( "color" ) , Verbose( "verbose" );
cmdLineParameter< int >	Depth( "depth" , 8 ) , SampleNum( "samples" , 100000 );

cmdLineReadable* params[] = { &Out , &SSDReconstruction , &UseColor , &Verbose , &Depth , &SampleNum , nullptr };

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <number of samples>\n" , SampleNum.name );
	printf( "\t[--%s <ouput mesh>]\n" , Out.name );
	printf( "\t[--%s <reconstruction depth>=%d]\n" , Depth.name , Depth.value );
	printf( "\t[--%s]\n" , UseColor.name );
	printf( "\t[--%s]\n" , SSDReconstruction.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

// A simple structure for representing colors. 
// Assuming values are in the range [0,1].
template< typename Real >
struct RGBColor
{
	// The channels
	RGBColor( Real r=0 , Real g=0 , Real b=0 ) : r(r) , g(g) , b(b){}
	Real r,g,b;

	// Methods supporting affine re-combination
	RGBColor &operator += ( const RGBColor &c ){ r += c.r , g += c.g , b += c.b ; return *this; }
	RGBColor &operator *= ( Real s ){ r *= s , g *= s , b *= s ;  return *this; }
	RGBColor &operator /= ( Real s ){ return operator *= (1/s); }

	RGBColor operator + ( const RGBColor &c ) const { return RGBColor( r+c.r , g+c.g , b+c.b ); }
	RGBColor operator * ( Real s ) const { return RGBColor( r*s , g*s , b*s ); }
	RGBColor operator / ( Real s ) const { return operator * (1/s); }
};

// A stream for generating random samples on the sphere
template< typename Real , unsigned int Dim >
struct SphereSampleStream : public Reconstructor::InputSampleStream< Real , Dim >
{
	// from https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::random_device randomDevice;
	std::default_random_engine generator;
	std::uniform_real_distribution< Real > distribution;

	// Constructs a stream that contains the specified number of samples
	SphereSampleStream( unsigned int sz ) : _size(sz) , _current(0) , generator(0) , distribution((Real)-1.0,(Real)1.0) {}

	// Overrides the pure abstract method from InputSampleStream< Real , Dim >
	void reset( void ){ generator.seed(0) ; _current = 0; }

	// Overrides the pure abstract method from InputSampleStream< Real , Dim >
	bool base_read( Point< Real , Dim > &p , Point< Real , Dim > &n )
	{
		if( _current<_size )
		{
			p = n = RandomSpherePoint( generator , distribution );
			_current++;
			return true;
		}
		else return false;
	}

	static Point< Real , Dim > RandomSpherePoint( std::default_random_engine &generator , std::uniform_real_distribution< Real > &distribution )
	{
		while( true )
		{
			Point< Real , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = distribution( generator );
			if( Point< Real , Dim >::SquareNorm( p )<1 ) return p / (Real)sqrt( Point< Real , Dim >::SquareNorm(p) );
		}
	}
protected:
	unsigned int _size , _current;
};

// A stream for generating random samples with color on the sphere
template< typename Real , unsigned int Dim >
struct SphereSampleWithColorStream : public Reconstructor::InputSampleWithDataStream< Real , Dim , RGBColor< Real > >
{
	// from https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::random_device randomDevice;
	std::default_random_engine generator;
	std::uniform_real_distribution< Real > distribution;

	// Constructs a stream that contains the specified number of samples
	SphereSampleWithColorStream( unsigned int sz ) :
		_size(sz) , _current(0) , generator(0) , distribution((Real)-1.0,(Real)1.0) ,
		Reconstructor::InputSampleWithDataStream< Real , Dim , RGBColor< Real > >( RGBColor< Real >() ) {}

	// Overrides the pure abstract method from InputSampleWithDataStream< Real , Dim , RGBColor< Real > >
	void reset( void ){ generator.seed(0) ; _current = 0; }

	// Overrides the pure abstract method from InputSampleWithDataStream< Real , Dim , RGBColor< Real > >
	bool base_read( Point< Real , Dim > &p , Point< Real , Dim > &n , RGBColor< Real > &c )
	{
		if( _current<_size )
		{
			p = n = RandomSpherePoint( generator , distribution );
			_current++;
			c.r = c.g = c.b = 0;
			if     ( p[0]<-1.f/3 ) c.r = 1.f;
			else if( p[0]< 1.f/3 ) c.g = 1.f;
			else                   c.b = 1.f;
			return true;
		}
		else return false;
	}

	static Point< Real , Dim > RandomSpherePoint( std::default_random_engine &generator , std::uniform_real_distribution< Real > &distribution )
	{
		while( true )
		{
			Point< Real , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = distribution( generator );
			if( Point< Real , Dim >::SquareNorm( p )<1 ) return p / (Real)sqrt( Point< Real , Dim >::SquareNorm(p) );
		}
	}
protected:
	unsigned int _size , _current;
};

// A stream into which we can write polygons of the form std::vector< node_index_type >
template< typename Index >
struct PolygonStream : public Reconstructor::OutputPolygonStream
{
	// Construct a stream that adds polygons to the vector of polygons
	PolygonStream( std::vector< std::vector< Index > > &polygonStream ) : _polygons( polygonStream ) {}

	// Override the pure abstract method from OutputPolygonStream
	void base_write( const std::vector< node_index_type > &polygon )
	{
		std::vector< Index > poly( polygon.size() );
		for( unsigned int i=0 ; i<polygon.size() ; i++ ) poly[i] = (Index)polygon[i];
		_polygons.push_back( poly );
	}
protected:
	std::vector< std::vector< Index > > &_polygons;
};

// A stream into which we can write the output vertices of the extracted mesh
template< typename Real , unsigned int Dim >
struct VertexStream : public Reconstructor::OutputVertexStream< Real , Dim >
{
	// Construct a stream that adds vertices into the coordinates
	VertexStream( std::vector< Real > &vCoordinates ) : _vCoordinates( vCoordinates ) {}

	// Override the pure abstract method from Reconstructor::OutputVertexStream< Real , Dim >
	void base_write( Point< Real , Dim > p , Point< Real , Dim > , Real ){ for( unsigned int d=0 ; d<Dim ; d++ ) _vCoordinates.push_back( p[d] ); }
protected:
	std::vector< Real > &_vCoordinates;
};

// A stream into which we can write the output vertices and colors of the extracted mesh
template< typename Real , unsigned int Dim >
struct VertexWithColorStream : public Reconstructor::OutputVertexWithDataStream< Real , Dim , RGBColor< Real > >
{
	// Construct a stream that adds vertices into the coordinates
	VertexWithColorStream( std::vector< Real > &vCoordinates , std::vector< Real > &rgbCoordinates ) :
		_vCoordinates( vCoordinates ) , _rgbCoordinates( rgbCoordinates ) {}

	// Override the pure abstract methodfrom Reconstructor::OutputVertexWithColorStream< Real , Dim >
	void base_write( Point< Real , Dim > p , Point< Real , Dim > , Real , RGBColor< Real > c )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) _vCoordinates.push_back( p[d] );
		_rgbCoordinates.push_back( c.r );
		_rgbCoordinates.push_back( c.g );
		_rgbCoordinates.push_back( c.b );
	}
protected:
	std::vector< Real > &_vCoordinates;
	std::vector< Real > &_rgbCoordinates;
};

template< typename Real >
void WritePly( std::string fileName , size_t vNum , const Real *vCoordinates , const Real *rgbCoordinates , const std::vector< std::vector< int > > &polygons )
{
	std::fstream file( fileName , std::ios::out );
	file << "ply" << std::endl;
	file << "format ascii 1.0" << std::endl;
	file << "element vertex " << vNum << std::endl;
	file << "property float x" << std::endl << "property float y" << std::endl << "property float z" << std::endl;
	if( rgbCoordinates ) file << "property uchar red" << std::endl << "property uchar green" << std::endl << "property uchar blue" << std::endl;
	file << "element face " << polygons.size() << std::endl;
	file << "property list uchar int vertex_indices" << std::endl;
	file << "end_header" << std::endl;

	auto ColorChannel = []( Real v ){ return std::max<int>( 0 , std::min<int>( 255 , (int)floor(255*v+0.5) ) ); };

	for( size_t i=0 ; i<vNum ; i++ )
	{
		file << vCoordinates[3*i+0] << " " << vCoordinates[3*i+1] << " " << vCoordinates[3*i+2];
		if( rgbCoordinates ) file << " " << ColorChannel( rgbCoordinates[3*i+0] ) << " " << ColorChannel( rgbCoordinates[3*i+1] ) << " " << ColorChannel( rgbCoordinates[3*i+2] );
		file << std::endl;
	}
	for( const auto &polygon : polygons )
	{
		file << polygon.size();
		for( auto vIdx : polygon ) file << " " << vIdx;
		file << std::endl;
	}
}

template
<
	typename Real ,			// Arithmetic type (float or double)
	unsigned int Dim ,		// Dimensionality of the reconstruction (=3)
	typename ReconType ,	// Reconstructor type (Reconstructor::Poisson or Reconstructor::SSD)
	bool UseColor			// Should color be reconstructed as well?
>
void Execute( void )
{
	// Finite-elements signature
	static const unsigned int FEMSig = FEMDegreeAndBType< ReconType::DefaultFEMDegree , ReconType::DefaultFEMBoundary >::Signature;

	// Parameters for performing the reconstruction
	typename ReconType::template SolutionParameters< Real > solverParams;

	solverParams.verbose = Verbose.set;
	solverParams.depth = (unsigned int)Depth.value;

	// Parameters for exracting the level-set surface
	Reconstructor::LevelSetExtractionParameters extractionParams;
	extractionParams.linearFit = SSDReconstruction.set;		// Since the SSD solution approximates a TSDF, linear fitting works well
	extractionParams.verbose = Verbose.set;

	if constexpr( UseColor )
	{
		// The type of the reconstructor
		using Implicit = typename ReconType::template Implicit< Real , Dim , FEMSig , RGBColor< Real > >;

		// A stream generating random points on the sphere with color
		SphereSampleWithColorStream< Real , Dim > sampleStream( SampleNum.value );

		// Construct the implicit representation
		Implicit implicit( sampleStream , solverParams );

		// Scale the color information to give extrapolation preference to data at finer depths
		implicit.weightAuxDataByDepth( (Real)32. );


		// vectors for storing the polygons (specifically, triangles), the coordinates of the vertices, and the colors at the vertices
		std::vector< std::vector< int > > polygons;
		std::vector< Real > vCoordinates , rgbCoordinates;

		// Streams backed by these vectors
		VertexWithColorStream< Real , Dim > vStream( vCoordinates , rgbCoordinates );
		PolygonStream< int > pStream( polygons );

		// Extract the iso-surface
		implicit.extractLevelSet( vStream , pStream , extractionParams );

		if( Out.set ) WritePly( Out.value , vStream.size() , vCoordinates.data() , rgbCoordinates.data() , polygons );
	}
	else
	{
		// The type of the reconstructor
		using Implicit = typename ReconType::template Implicit< Real , Dim , FEMSig >;

		// A stream generating random points on the sphere
		SphereSampleStream< Real , Dim > sampleStream( SampleNum.value );

		// Construct the implicit representation
		Implicit implicit( sampleStream , solverParams );

		// vectors for storing the polygons (specifically, triangles) and the coordinates of the vertices
		std::vector< std::vector< int > > polygons;
		std::vector< Real > vCoordinates;

		// Streams backed by these vectors
		PolygonStream< int > pStream( polygons );
		VertexStream< Real , Dim > vStream( vCoordinates );

		// Extract the iso-surface
		implicit.extractLevelSet( vStream , pStream , extractionParams );

		if( Out.set ) WritePly( Out.value , vStream.size() , vCoordinates.data() , (Real*)nullptr , polygons );
	}
}

int main( int argc , char* argv[] )
{
	Timer timer;
	cmdLineParse( argc-1 , &argv[1] , params );
#ifdef _OPENMP
	ThreadPool::Init( ThreadPool::OPEN_MP , std::thread::hardware_concurrency() );
#else // !_OPENMP
	ThreadPool::Init( ThreadPool::THREAD_POOL , std::thread::hardware_concurrency() );
#endif // _OPENMP

	if( !SampleNum.set )
	{
		ShowUsage( argv[0] );
		return 0;
	}

	if( Verbose.set )
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "************************************************" << std::endl;
		std::cout << "** Running SSD Reconstruction (Version " << ADAPTIVE_SOLVERS_VERSION << ") **" << std::endl;
		std::cout << "************************************************" << std::endl;
		std::cout << "************************************************" << std::endl;
	}
	
	// Solve using single float precision, in dimension 3, w/ finite-elements of degree 2 for SSD and degree 1 for Poisson, and using Neumann boundaries
	if( SSDReconstruction.set )
		if( UseColor.set ) Execute< float , 3 , Reconstructor::   SSD  , true  >();
		else               Execute< float , 3 , Reconstructor::   SSD  , false >();
	else
		if( UseColor.set ) Execute< float , 3 , Reconstructor::Poisson , true  >();
		else               Execute< float , 3 , Reconstructor::Poisson , false >();

	if( Verbose.set )
	{
		printf( "Time (Wall/CPU): %.2f / %.2f\n" , timer.wallTime() , timer.cpuTime() );
		printf( "Peak Memory (MB): %d\n" , MemoryInfo::PeakMemoryUsageMB() );
	}

	ThreadPool::Terminate();
	return EXIT_SUCCESS;
}
