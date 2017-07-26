/*
Copyright (c) 2013, Michael Kazhdan
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

#undef ARRAY_DEBUG

#include <iostream>
#include <string>

#include "CmdLineParser.h"
#include "Ply.h"
#include "MyTime.h"

#include "SurfaceTrimmer.h"

/**
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include <algorithm>
#include "Geometry.h"
#include "MAT.h"
**/

cmdLineString In( "in" ) , Out( "out" );
cmdLineInt Smooth( "smooth" , 5 );
cmdLineFloat Trim( "trim" ) , IslandAreaRatio( "aRatio" , 0.001f );
cmdLineReadable PolygonMesh( "polygonMesh" );

typedef std::vector<int> Polygon;

cmdLineReadable* params[] =
{
	&In , &Out , &Trim , &PolygonMesh , &Smooth , &IslandAreaRatio
};

void ShowUsage( char* ex )
{
	printf( "Usage: %s\n" , ex );
	printf( "\t --%s <input polygon mesh>\n" , In.name );
	printf( "\t[--%s <trimming value>]\n" , Trim.name );
	printf( "\t[--%s <ouput polygon mesh>]\n" , Out.name );
	printf( "\t[--%s <smoothing iterations>=%d]\n" ,
        Smooth.name , Smooth.value );
	printf( "\t[--%s <relative area of islands>=%f]\n" ,
        IslandAreaRatio.name , IslandAreaRatio.value );
	printf( "\t[--%s]\n" , PolygonMesh.name );
}

template <typename Vertex>
void readData(std::vector<Vertex>& vertices, std::vector<Polygon>& polygons,
    std::vector<std::string> comments)
{
    int fileType = 0;
    char **cComments;
    int numComments = 0;
    PlyReadPolygons(In.value, vertices, polygons, Vertex::ReadProperties,
          Vertex::ReadComponents, fileType, &cComments, &numComments);
    for (auto i = 0; i < numComments; ++i)
    {
        std::string c = cComments[i];
        comments.push_back(c);
        free(cComments[i]);
    }
}

template <typename Vertex>
void writeData(std::vector<Vertex>& vertices, std::vector<Polygon>& polygons,
    std::vector<std::string> comments)
{
    if (!Out.set)
        return;

    int fileType = 0;

    char **cComments = new char *[comments.size()];
    for (size_t i = 0; i < comments.size(); ++i)
        cComments[i] = strdup(comments[i].data());

    PlyWritePolygons(Out.value, vertices, polygons, Vertex::WriteProperties,
          Vertex::WriteComponents, fileType, cComments, comments.size());

    for (size_t i = 0; i < comments.size(); ++i)
        free(cComments[i]);
    delete [] cComments;
}


void addComments(std::vector<std::string> comments)
{
    std::string leader("\t--");

    comments.push_back("Running Surface Trimmer (V5)");
    if (In.set)
        comments.push_back(leader + In.name + " " + In.value);
    if (Out.set)
        comments.push_back(leader + Out.name + " " + Out.value);
    if (Trim.set)
        comments.push_back(leader + Trim.name + " " +
            std::to_string(Trim.value));
    if (Smooth.set)
        comments.push_back(leader + IslandAreaRatio.name + " " +
            std::to_string(IslandAreaRatio.value));
    if (PolygonMesh.set)
        comments.push_back(leader + PolygonMesh.name);
}


template <typename Vertex>
void printVertexValues(std::vector<Vertex>& vertices)
{
    if (vertices.empty())
        return;

    float min = vertices[0].value;
    float max = vertices[0].value;

    for (Vertex& v : vertices)
    {
        min = std::min<float>(min, v.value);
        max = std::max<float>(min, v.value);
    }
    std::cout << "Value Range: [" << min << "," << max << "]\n";
}


template <typename Vertex>
int execute()
{
    float min;
    float max;

    using ST = SurfaceTrimmer<Vertex, float>;

    std::vector<Vertex> vertices;
    std::vector<typename ST::Polygon> polygons;
    std::vector<std::string> comments;

    readData(vertices, polygons, comments);
    addComments(comments);

    ST trimmer(vertices, polygons, Trim.value,
        IslandAreaRatio.value, PolygonMesh.set);

    while (Smooth.value--)
        trimmer.SmoothValues();
    printVertexValues(vertices);

    double t = Time();

    int result = trimmer.Execute();

    comments.push_back("#Trimmed In " +  std::to_string((Time() - t)) +
        " (s)");
    writeData(vertices, polygons, comments);
    return result;
}


int main( int argc , char* argv[] )
{
	int paramNum = sizeof(params)/sizeof(cmdLineReadable*);
	cmdLineParse( argc-1 , &argv[1] , paramNum , params , 0 );

	if( !In.set || !Trim.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	bool readFlags[ PlyColorAndValueVertex< float >::ReadComponents ];
	if (!PlyReadHeader(In.value, PlyColorAndValueVertex<float>::ReadProperties,
         PlyColorAndValueVertex<float>::ReadComponents, readFlags))
    {
        fprintf(stderr , "[ERROR] Failed to read ply header: %s\n" , In.value);
        exit( 0 );
    }

	bool hasValue = readFlags[3];
	bool hasColor = ( readFlags[4] || readFlags[7] ) &&
        ( readFlags[5] || readFlags[8] ) && ( readFlags[6] || readFlags[9] );

	if( !hasValue )
    {
        fprintf( stderr , "[ERROR] Ply file does not contain values\n" );
        exit( 0 );
    }

	if ( hasColor)
        return execute<PlyColorAndValueVertex<float>>();
	else
        return execute<PlyValueVertex<float>>();
}

