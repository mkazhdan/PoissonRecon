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

class TransformedPointSource : public PointSource
{
    XForm4x4<double> m_xForm;
    XForm3x3<double> m_normalXForm;
    PointSource& m_stream;

public:
    TransformedPointSource(const XForm4x4<double>& xForm, PointSource& stream) :
        m_xForm(xForm), m_stream(stream)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
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
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
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
