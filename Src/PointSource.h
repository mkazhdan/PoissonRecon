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

#include "Ply.h"
#include "Geometry.h"

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

class TransformedPointSource : public PointSource
{
    XForm4x4<double> m_xForm;
    XForm3x3<double> m_normalXForm;
    PointSource& m_stream;

public:
    TransformedPointSource(const XForm4x4<double>& xForm, PointSource& stream) :
        m_xForm(xForm), m_stream(stream)
    {
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                m_normalXForm(i, j) = m_xForm(i, j);
        m_normalXForm = m_normalXForm.transpose().inverse();
    }

    virtual void reset()
        { m_stream.reset(); }

	virtual bool nextPoint(Point& p)
	{
		bool ret = m_stream.nextPoint(p);
		p.p = m_xForm * p.p;
        p.n = m_normalXForm * p.n;
		return ret;
	}
};

class ColorTransformedPointSource : public ColorPointSource
{
private:
    XForm4x4<double> m_xForm;
    XForm3x3<double> m_normalXForm;
    ColorPointSource& m_stream;

public:
    ColorTransformedPointSource(const XForm4x4<double>& xForm,
        ColorPointSource& stream) : m_xForm(xForm), m_stream(stream)
    {
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                m_normalXForm(i, j) = m_xForm(i, j);
        m_normalXForm = m_normalXForm.transpose().inverse();
    }

    virtual void reset()
        { m_stream.reset(); }

	virtual bool nextPoint(Point& p, Point3D<double>& color)
	{
		bool ret = m_stream.nextPoint(p, color);
		p.p = m_xForm * p.p;
        p.n = m_normalXForm * p.n;
		return ret;
	}
};

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
};

class ASCIIPointSource : public PointSource
{
private:
	FILE *m_fp;
public:
	ASCIIPointSource(const std::string& filename);
	~ASCIIPointSource();

	void reset();
	bool nextPoint(Point& p);
};

class ColorASCIIPointSource : public ColorPointSource
{
private:
	FILE* m_fp;
public:
	ColorASCIIPointSource(const std::string& filename);
	~ColorASCIIPointSource();

	void reset();
	bool nextPoint(Point& point, Point3D<double>& color);
};

class BinaryPointSource : public PointSource
{
private:
    FILE *m_fp;
    std::vector<Point> m_points;
    size_t m_current;
public:
    BinaryPointSource(const std::string& filename);
    ~BinaryPointSource();

    void reset();
    bool nextPoint(Point& point);
};

class ColorBinaryPointSource : public ColorPointSource
{
private:
    FILE *m_fp;
    std::vector<PointPair> m_points;
    size_t m_current;
public:
    ColorBinaryPointSource(const std::string& filename);
    ~ColorBinaryPointSource();

    void reset();
    bool nextPoint(Point& point, Point3D<double>& color);
};

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

#include "PointSource.inl"
