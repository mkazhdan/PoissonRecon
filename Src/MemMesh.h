/*
Copyright (c) 2017, Andrew Bell
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

#include "Mesh.h"

namespace Kazhdan
{

class BaseMemMesh : public Mesh
{
protected:
    std::vector<Polygon> m_polys;
    int m_polyIter;
    int m_pointIter;

public:
    BaseMemMesh() : m_polyIter(0), m_pointIter(0)
    {}

    virtual int polygonCount() const
        { return (int)m_polys.size(); }
    virtual void newPolygon(std::vector<int>& poly)
        { 
#pragma omp critical(MESH_POLY)
			m_polys.push_back(poly); 
		}

    virtual void resetIterator()
    {
        m_polyIter = 0;
        m_pointIter = 0;
    }

    virtual bool nextPolygon(Polygon& p)
    {
        p = m_polys[m_polyIter++];
        return m_polyIter < m_polys.size();
    }
};

class MemMesh : public BaseMemMesh
{
    struct Point
    {
        double m_x;
        double m_y;
        double m_z;
    };

    std::vector<Point> m_points;

public:
    virtual int pointCount() const
        { return (int)m_points.size(); }
    virtual int newPoint(const std::array<double, 3>& position)
    {
#pragma omp critical(MESH_POINT)
        m_points.push_back(Point({position[0], position[1], position[2]}));
        return (int)m_points.size() - 1;
    }

    virtual bool nextPoint(Kazhdan::Point& p)
    {
        p.m_position[0] = m_points[m_pointIter].m_x;
        p.m_position[1] = m_points[m_pointIter].m_y;
        p.m_position[2] = m_points[m_pointIter].m_z;
        return ++m_pointIter < m_points.size();
    }
};

class ColorMemMesh : public MemMesh
{
    struct Point
    {
        double m_x;
        double m_y;
        double m_z;
        uint8_t m_red;
        uint8_t m_green;
        uint8_t m_blue;
    };

    std::vector<Point> m_points;

public:
    virtual int pointCount() const
        { return (int)m_points.size(); }
    virtual int newPoint(const std::array<double, 3>& position,
        const std::array<uint8_t, 3>& color)
    {
#pragma omp critical( MESH_POINT )
        m_points.push_back(
            Point({position[0], position[1], position[2],
                   color[0], color[1], color[2]}));
        return (int)m_points.size() - 1;
    }

    virtual bool nextPoint(Kazhdan::Point& p)
    {
        p.m_position[0] = m_points[m_pointIter].m_x;
        p.m_position[1] = m_points[m_pointIter].m_y;
        p.m_position[2] = m_points[m_pointIter].m_z;
        p.m_color[0] = m_points[m_pointIter].m_red;
        p.m_color[1] = m_points[m_pointIter].m_green;
        p.m_color[2] = m_points[m_pointIter].m_blue;
        return ++m_pointIter < m_points.size();
    }
};

class DensityMemMesh : public MemMesh
{
    struct Point
    {
        double m_x;
        double m_y;
        double m_z;
        double m_density;
    };

    std::vector<Point> m_points;

public:
    virtual int pointCount() const
        { return (int)m_points.size(); }
    virtual int newPoint(const std::array<double, 3>& position,
        double density)
    {
#pragma omp critical( MESH_POINT )
        m_points.push_back(
            Point({position[0], position[1], position[2], density}));
        return (int)m_points.size() - 1;
    }

    virtual bool nextPoint(Kazhdan::Point& p)
    {
        p.m_position[0] = m_points[m_pointIter].m_x;
        p.m_position[1] = m_points[m_pointIter].m_y;
        p.m_position[2] = m_points[m_pointIter].m_z;
        p.m_density = m_points[m_pointIter].m_density;
        return ++m_pointIter < m_points.size();
    }
};

class CompleteMemMesh : public MemMesh
{
    struct Point
    {
        double m_x;
        double m_y;
        double m_z;
        double m_density;
        uint8_t m_red;
        uint8_t m_green;
        uint8_t m_blue;
    };

    std::vector<Point> m_points;

public:
    virtual int pointCount() const
        { return (int)m_points.size(); }
    virtual int newPoint(const std::array<double, 3>& position,
        const std::array<uint8_t, 3>& color, double density)
    {
#pragma omp critical( MESH_POINT )
        m_points.push_back(
            Point({position[0], position[1], position[2], density,
                   color[0], color[1], color[2]}));
        return (int)m_points.size();
    }

    virtual bool nextPoint(Kazhdan::Point& p)
    {
        p.m_position[0] = m_points[m_pointIter].m_x;
        p.m_position[1] = m_points[m_pointIter].m_y;
        p.m_position[2] = m_points[m_pointIter].m_z;
        p.m_color[0] = m_points[m_pointIter].m_red;
        p.m_color[1] = m_points[m_pointIter].m_green;
        p.m_color[2] = m_points[m_pointIter].m_blue;
        p.m_density = m_points[m_pointIter].m_density;
        return ++m_pointIter < m_points.size();
    }
};

} // namespace
