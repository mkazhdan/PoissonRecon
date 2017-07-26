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

//////////////////////////
// BinaryPointSource    //
//////////////////////////
BinaryPointSource::BinaryPointSource(const std::string& filename) : m_current(0)
{
    m_fp = fopen(filename.data(), "rb");
    if (!m_fp)
    {
        std::cerr << "Failed to open file for reading: " <<
            filename << std::endl;
        exit(0);
    }
}

BinaryPointSource::~BinaryPointSource()
{
    fclose(m_fp);
}

void BinaryPointSource::reset()
{
    fseek(m_fp, SEEK_SET, 0);
    m_current = 0;
    m_points.clear();
}

bool BinaryPointSource::nextPoint(Point& p)
{
    if (m_current < m_points.size())
    {
        p = m_points[m_current++];
        return true;
    }
    else
    {
        m_current = 0;
        m_points.resize(1024);
        size_t count = fread(m_points.data(), sizeof(Point),
            m_points.size(), m_fp);
        m_points.resize(count);

        if (count == 0)
            return false;
        else
            return nextPoint(p);
    }
}

///////////////////////////////
// ColorBinaryPointSource    //
///////////////////////////////
ColorBinaryPointSource::ColorBinaryPointSource(const std::string& filename) :
    m_current(0)
{
    m_fp = fopen(filename.data(), "rb");
    if (!m_fp)
    {
        std::cerr << "Failed to open file for reading: " <<
            filename << std::endl;
        exit(0);
    }
}

ColorBinaryPointSource::~ColorBinaryPointSource()
{
    fclose(m_fp);
}

void ColorBinaryPointSource::reset()
{
    fseek(m_fp, SEEK_SET, 0);
    m_current = 0;
    m_points.clear();
}

bool ColorBinaryPointSource::nextPoint(Point& p, Point3D<double>& color)
{
    if (m_current < m_points.size())
    {
        p = m_points[m_current].first;
        color = m_points[m_current].second;
        m_current++;
        return true;
    }
    else
    {
        m_current = 0;
        m_points.resize(1024);
        size_t count = fread(m_points.data(), sizeof(Point),
            m_points.size(), m_fp);
        m_points.resize(count);

        if (count == 0)
            return false;
        else
            return nextPoint(p, color);
    }
}

