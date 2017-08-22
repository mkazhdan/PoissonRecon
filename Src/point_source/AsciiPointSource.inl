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

#include <iostream>

///////////////////////////////
// ASCIIPointSource          //
///////////////////////////////
ASCIIPointSource::ASCIIPointSource(const std::string& filename)
{
    m_fp = fopen(filename.data(), "r");
    if (!m_fp)
    {
        std::cerr << "Failed to open file for reading: " <<
            filename << std::endl;
        exit(0);
    }
}

ASCIIPointSource::~ASCIIPointSource()
{
    fclose(m_fp);
}

void ASCIIPointSource::reset()
{
    fseek(m_fp, SEEK_SET, 0);
}

bool ASCIIPointSource::nextPoint(Point& p)
{
    if (fscanf(m_fp , " %lf %lf %lf %lf %lf %lf ",
        &p.p[0], &p.p[1], &p.p[2], &p.n[0], &p.n[1], &p.n[2]) != 6)
        return false;
    return true;
}

///////////////////////////////
// ColorASCIIPointSource     //
///////////////////////////////
ColorASCIIPointSource::ColorASCIIPointSource(const std::string& filename)
{
    m_fp = fopen(filename.data(), "r");
    if (!m_fp)
    {
        std::cerr << "Failed to open file for reading: " <<
            filename << std::endl;
        exit(0);
    }
}

ColorASCIIPointSource::~ColorASCIIPointSource()
{
    fclose(m_fp);
}

void ColorASCIIPointSource::reset()
{
    fseek(m_fp, SEEK_SET, 0);
}

bool ColorASCIIPointSource::nextPoint(Point& p, Point3D<double>& color)
{
    if (fscanf(m_fp , " %lf %lf %lf %lf %lf %lf ",
        &p.p[0] , &p.p[1] , &p.p[2] , &p.n[0], &p.n[1], &p.n[2]) != 6)
         return false;
    char c[3];
    if (fscanf(m_fp, " %c %c %c ", &c[0], &c[1], &c[2]) != 3)
        return false;
    color.coords[0] = c[0];
    color.coords[1] = c[1];
    color.coords[2] = c[2];
    return true;
}


