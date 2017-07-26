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
prior writften permission.

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

#pragma once

#include "PointSource.h"
#include "../Ply.h"

class PLYPointSource : public PointSource
{
    struct Vertex
    {
        Point3D<double> point;
        Point3D<double> normal;
    };

    std::string m_filename;
	PlyFile *m_ply;
	int m_numElements;
    int m_numInstances;
	char **m_elements;
	int m_current;

	void _free( void );
public:
	PLYPointSource(const std::string& filename);
	~PLYPointSource();
	void reset();
	bool nextPoint(Point& p);
};

class ColorPLYPointSource : public ColorPointSource
{
    struct Vertex
    {
        Point3D<double> point;
        Point3D<double> normal;
        Point3D<double> color;
    };

	std::string m_filename;
	PlyFile *m_ply;
	int m_numElements;
    int m_numInstances;
	char **m_elements;
	int m_current;
	void _free( void );
public:
	ColorPLYPointSource(const std::string& namefile);
	~ColorPLYPointSource();
	void reset();
	bool nextPoint(Point& p, Point3D<double>& color);
};

#include "PlyPointSource.inl"
