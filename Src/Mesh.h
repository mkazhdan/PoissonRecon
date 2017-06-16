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

#include <array>
#include <vector>

namespace Kazhdan
{

struct Point
{
    std::array<double, 3> m_position;
    std::array<uint8_t, 3> m_color;
    double m_density;
};

using Polygon = std::vector<int>;

class Mesh
{
public:
    virtual int pointCount() const = 0;
    virtual int polygonCount() const = 0;
    virtual void newPolygon(std::vector<int>& poly) = 0;
    virtual int newPoint(const std::array<double, 3>& position) = 0;
    virtual int newPoint(const std::array<double, 3>& position,
        const std::array<uint8_t, 3>& color)
    {
        throw std::runtime_error("Mesh doesn't support color data.");
    }

    virtual int newPoint(const std::array<double, 3>& position,
        double density)
    {
        throw std::runtime_error("Mesh doesn't support density data.");
    }

    virtual int newPoint(const std::array<double, 3>& position,
        const std::array<uint8_t, 3>& color, double density)
    {
        throw std::runtime_error("Mesh doesn't support color "
            "and density data.");
    }

    virtual bool hasColor() const
        { return false; }
    virtual bool hasDensity() const
        { return false; }

    virtual void resetIterator() = 0;
    virtual bool nextPolygon(Polygon& poly) = 0;
    virtual bool nextPoint(Point& point) = 0;
};

} // namespace
