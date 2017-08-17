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

class MemoryPointSource : public PointSource
{
private:
    const Point *m_points;
    size_t m_pointCount;
    size_t m_current;
public:
    MemoryPointSource(size_t pointCount, const Point *points) :
        m_points(points), m_pointCount(pointCount), m_current(0)
    {}

    virtual void reset()
        { m_current = 0; }

    virtual bool nextPoint(Point& p)
    {
        if (m_current >= m_pointCount)
            return false;
        p = m_points[m_current++];
        return true;
    }
};

class ColorMemoryPointSource : public ColorPointSource
{
private:
    const PointPair *m_points;
    size_t m_pointCount;
    size_t m_current;
public:
    ColorMemoryPointSource(size_t pointCount, const PointPair *points) :
        m_points(points), m_pointCount(pointCount), m_current(0)
    {}

    virtual bool nextPoint(Point& p, Point3D<double>& color)
    {
        if (m_current >= m_pointCount)
            return false;
        const PointPair& pp(m_points[m_current++]);
        p = pp.first;
        color = pp.second;
        return true;
    }

	virtual void reset() 
	{
		m_current = 0;
	}
};


#include "PointSource.inl"
