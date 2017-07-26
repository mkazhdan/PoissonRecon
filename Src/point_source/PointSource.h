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

#include "../Geometry.h"

class PointSource
{
protected:
    using Point = OrientedPoint3D<double>; // Contains Normals

public:
    virtual ~PointSource()
    {}

    virtual void reset() = 0;
    virtual bool nextPoint(Point& point) = 0;
    void boundingBox(Point3D<double>& min, Point3D<double>& max)
	{
        for (size_t i = 0; i < 3; ++i)
        {
            min[i] = std::numeric_limits<double>::max();
            max[i] = std::numeric_limits<double>::lowest();
        }

        Point p;
		while (nextPoint(p))
		{
			for( int i=0 ; i<3 ; i++ )
			{
                min[i] = std::min(min[i], p.p[i]);
                max[i] = std::max(max[i], p.p[i]);
			}
		}
		reset();
	}
};
typedef std::unique_ptr<PointSource> PointSourcePtr;

class ColorPointSource : public PointSource
{
protected:
    using PointPair = std::pair<Point, Point3D<double>>;
public:
    virtual bool nextPoint(Point& p, Point3D<double>& color) = 0;
    virtual bool nextPoint(Point& p)
    {
        Point3D<double> color;
        return nextPoint(p, color);
    }
};
