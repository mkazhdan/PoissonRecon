/*
Copyright (c) 2016, Michael Kazhdan
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

#undef FAST_COMPILE				// If enabled, only a single version of the reconstruction code is compiled
#undef ARRAY_DEBUG				// If enabled, array access is tested for validity
#define MAX_MEMORY_GB 15		// If non-zero, the maximum memory to be used by the application

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "Array.h"

cmdLineParameterArray< char* , 2 >
	In( "in" );
cmdLineParameter< char* >
	Out( "out" );
cmdLineParameter< float >
	Scale( "scale" , 1.f );

cmdLineReadable* params[] =
{
	&In , &Out , &Scale ,
	NULL
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input voxel grid1 , input voxel grid 2>\n" , In.name );
	printf( "\t[--%s <output voxel grid>]\n" , Out.name );
	printf( "\t[--%s <output scale>=%f]\n" , Scale.name , Scale.value );
}


int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , &argv[1] , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	auto ReadVoxel = []( const char* fileName , int& res )
	{
		FILE* fp = fopen( fileName , "rb" );
		if( !fp ) fprintf( stderr , "[ERROR] Failed to read voxel: %s\n" , fileName ) , exit( 0 );
		if( fread( &res , sizeof(int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to read restolution from file: %s\n" , fileName ) , exit( 0 );
		Pointer( float ) v = AllocPointer< float >( res*res*res );
		if( !v ) fprintf( stderr , "[ERROR] Failed to allocate voxel grid: %d x %d x %d\n" , res , res , res ) , exit( 0 );
		if( fread( v , sizeof(float) , res*res*res , fp )!=res*res*res ) fprintf( stderr , "[ERROR] Failed to read voxel values from file: %s\n" , fileName ) , exit( 0 );
		fclose( fp );
		return v;
	};

	int res;
	Pointer( float ) v1 ; Pointer( float ) v2;
	{
		int res1 , res2;
		v1 = ReadVoxel( In.values[0] , res1 );
		v2 = ReadVoxel( In.values[1] , res2 );
		if( res1!=res2 ) fprintf( stderr , "[ERROR] Voxel resolutions don't match: %d x %d\n" , res1 , res2 ) , exit( 0 );
		res = res1;
	}

	double l1Error = 0 , l2Error = 0;
	double l1Norm1 = 0 , l1Norm2 = 0 , l2Norm1 = 0 , l2Norm2 = 0;
#pragma omp parallel for reduction ( + : l1Error , l2Error , l1Norm1 , l1Norm2 , l2Norm1 , l2Norm2 )
	for( int i=0 ; i<res*res*res ; i++ )
	{
		l1Error += fabs( v1[i] - v2[i] ) , l1Norm1 += fabs( v1[i] ) , l1Norm2 += fabs( v2[i] );
		l2Error += ( v1[i] - v2[i] ) * ( v1[i] - v2[i] ) , l2Norm1 += v1[i] * v1[i] , l2Norm2 += v2[i] * v2[i];
	}
	l1Error /= res*res*res , l1Norm1 /= res*res*res , l1Norm2 /= res*res*res;
	l2Error /= res*res*res , l2Norm1 /= res*res*res , l2Norm2 /= res*res*res;
	l2Error = sqrt( l2Error ) , l2Norm1 = sqrt( l2Norm1 ) , l2Norm2 = sqrt( l2Norm2 );
	printf(  "L1 / L2 differences: %g %g\n" , 2*l1Error / ( l1Norm1 + l1Norm2 ) , 2*l2Error / ( l2Norm1 + l2Norm2 ) );

	if( Out.set )
	{
		if( Scale.value<0 ) Scale.value = 1. / ( ( l1Norm1 + l1Norm2 ) / 2 );
#pragma omp parallel for
		for( int i=0 ; i<res*res*res ; i++ ) v1[i] = ( v1[i] - v2[i] ) * Scale.value;

		FILE* fp = fopen( Out.value , "wb" );
		if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , Out.value ) , exit( 0 );

		if( fwrite( &res , sizeof(int) , 1 , fp )!=1 ) fprintf( stderr , "[ERROR] Failed to write voxel resolution\n" ) , fclose( fp ) , exit( 0 );
		fwrite( v1 , sizeof(float) , res*res*res , fp );
		fclose( fp );
	}
	FreePointer( v1 );
	FreePointer( v2 );

	return EXIT_SUCCESS;
}
