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

//////////////////////
// PLYPointSource   //
//////////////////////
PLYPointSource::PLYPointSource(const std::string& filename) :
    m_filename(filename), m_ply(nullptr)
{
    reset();
}

PLYPointSource::~PLYPointSource()
{
    _free();
}

void PLYPointSource::reset()
{
    int fileType;
    float version;

    // The _strdup looks leaky, but I don't think I care.
    static PlyProperty requiredProps[] =
    {
        { _strdup("x") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("y") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("z") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("nx") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("ny") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("nz") , PLY_DOUBLE, PLY_DOUBLE }
    };
    const int numRequiredProps = sizeof(requiredProps) / sizeof(PlyProperty);

    Vertex v;
    requiredProps[0].offset = (int)((char *)&v.point.coords[0] - (char *)&v);
    requiredProps[1].offset = (int)((char *)&v.point.coords[1] - (char *)&v);
    requiredProps[2].offset = (int)((char *)&v.point.coords[2] - (char *)&v);
    requiredProps[3].offset = (int)((char *)&v.normal.coords[0] - (char *)&v);
    requiredProps[4].offset = (int)((char *)&v.normal.coords[1] - (char *)&v);
    requiredProps[5].offset = (int)((char *)&v.normal.coords[2] - (char *)&v);

    if (m_ply)
        _free();
    m_ply = ply_open_for_reading(const_cast<char *>(m_filename.data()),
        &m_numElements, &m_elements, &fileType, &version);
    if (!m_ply)
    {
        std::cerr << "[ERROR] Failed to open ply file for reading: " <<
            m_filename << std::endl;
        exit(0);
    }
    m_current = 0;

    PlyProperty** plist;
    bool foundVertices = false;
    for(int i = 0; i < m_numElements; ++i)
    {
        int numProperties;
        char* elementName = m_elements[i];
        plist = ply_get_element_description(m_ply , elementName,
            &m_numInstances, &numProperties);
        if (!plist)
        {
            std::cerr << "[ERROR] Failed to get element description: " <<
                elementName << std::endl;
            exit(0);
        }

        if (equal_strings("vertex", elementName))
        {
            foundVertices = true;
            for (int prop = 0; prop < numRequiredProps; ++prop)
                if (!ply_get_property(m_ply, elementName,
                    &(requiredProps[prop])))
                {
                    std::cerr << "[ERROR] Failed to find property in "
                        "ply file: " << requiredProps[prop].name << std::endl;
                    exit(0);
                }
        }
        for(size_t j = 0 ; j < numProperties; ++j)
        {
            free(plist[j]->name);
            free(plist[j]);
        }
        free(plist);
        if (foundVertices)
            break;
    }
    if (!foundVertices)
    {
        std::cerr << "[ERROR] Could not find vertices in ply file" << std::endl;
        exit(0);
    }
}

void PLYPointSource::_free()
{
    if (m_ply)
        ply_close(m_ply);
    m_ply = NULL;
    if (m_elements)
    {
        for (int i =0; i < m_numElements; ++i)
            free(m_elements[i]);
        free(m_elements);
    }
}

bool PLYPointSource::nextPoint(Point& p)
{
    if (m_current >= m_numInstances)
        return false;

    Vertex op;
    ply_get_element(m_ply, (void *)&op);
    p.p = op.point;
    p.n = op.normal;
    m_current++;
    return true;
}

///////////////////////////
// ColorPLYPointSource   //
///////////////////////////

ColorPLYPointSource::ColorPLYPointSource(const std::string& filename) :
    m_filename(filename), m_ply(nullptr)
{
    reset();
}

void ColorPLYPointSource::reset()
{
    int fileType;
    float version;
    PlyProperty** plist;

    // The _strdup looks leaky, but I don't think I care.
    static PlyProperty requiredProps[] =
    {
        { _strdup("x") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("y") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("z") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("nx") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("ny") , PLY_DOUBLE, PLY_DOUBLE },
        { _strdup("nz") , PLY_DOUBLE, PLY_DOUBLE }
    };
    const int numRequiredProps = sizeof(requiredProps) / sizeof(PlyProperty);

    Vertex v;
    requiredProps[0].offset = (int)((char *)&v.point.coords[0] - (char *)&v);
    requiredProps[1].offset = (int)((char *)&v.point.coords[1] - (char *)&v);
    requiredProps[2].offset = (int)((char *)&v.point.coords[2] - (char *)&v);
    requiredProps[3].offset = (int)((char *)&v.normal.coords[0] - (char *)&v);
    requiredProps[4].offset = (int)((char *)&v.normal.coords[1] - (char *)&v);
    requiredProps[5].offset = (int)((char *)&v.normal.coords[2] - (char *)&v);

    PlyProperty reqDataProps[] =
    {
        { "r"     ,  PLY_UCHAR, PLY_DOUBLE },
        { "g"     , PLY_UCHAR, PLY_DOUBLE },
        { "b"     , PLY_UCHAR, PLY_DOUBLE },
        { "red"   , PLY_UCHAR, PLY_DOUBLE },
        { "green" ,PLY_UCHAR, PLY_DOUBLE },
        { "blue"  , PLY_UCHAR, PLY_DOUBLE }
    };
    const int numReqDataProps = sizeof(reqDataProps) / sizeof(PlyProperty);

    reqDataProps[0].offset = (int)((char *)&v.color.coords[0] - (char *)&v);
    reqDataProps[1].offset = (int)((char *)&v.color.coords[1] - (char *)&v);
    reqDataProps[2].offset = (int)((char *)&v.color.coords[2] - (char *)&v);
    reqDataProps[3].offset = (int)((char *)&v.color.coords[0] - (char *)&v);
    reqDataProps[4].offset = (int)((char *)&v.color.coords[1] - (char *)&v);
    reqDataProps[5].offset = (int)((char *)&v.color.coords[2] - (char *)&v);

    if (m_ply)
        _free();
    m_ply = ply_open_for_reading(const_cast<char *>(m_filename.data()),
        &m_numElements, &m_elements, &fileType, &version);
    if (!m_ply)
    {
        std::cerr << "[ERROR] Failed to open ply file for reading: " <<
            m_filename << std::endl;
        exit(0);
    }

    m_current = 0;
    bool foundVertices = false;
    for (size_t i = 0; i < m_numElements; ++i)
    {
        int numProperties;
        char* elementName = m_elements[i];
        plist = ply_get_element_description(m_ply, elementName,
            &m_numInstances, &numProperties);
        if (!plist)
        {
            std::cerr << "[ERROR] Failed to get element description: " <<
                elementName << std::endl;
            exit(0);
        }

        if (equal_strings("vertex", elementName))
        {
            foundVertices = true;

            for (int prop = 0; prop < numRequiredProps; ++prop)
                if (!ply_get_property(m_ply, elementName,
                    &(requiredProps[prop])))
                {
                    std::cerr << "[ERROR] Failed to find property in "
                        "ply file: " << requiredProps[prop].name << std::endl;
                    exit( 0 );
                }

            int propIndicator = 0;
            for (size_t prop = 0; prop < numReqDataProps; ++prop)
            {
                // The properties r g b/red green blue are represented by
                // the low 3 bits.  When the propIndicator is 7, we know
                // we found all three.
                int bit = (int)pow(2, prop % 3);
                PlyProperty& property(reqDataProps[prop]);
                if (ply_get_property(m_ply, elementName, &property))
                    propIndicator |= bit;
            }
            if (propIndicator != 7)
            {
                std::cerr << "[ERROR] Missing at least one of red/gree/blue "
                    "in vertex property list." << std::endl;
                exit(0);
            }
        }
        for (int j=0; j < numProperties; ++j)
        {
            free(plist[j]->name);
            free(plist[j]);
        }
        free(plist);
        if (foundVertices)
            break;
    }
    if (!foundVertices)
    {
        std::cerr << "[ERROR] Could not find vertices in ply file" << std::endl;
        exit(0);
    }
}

void ColorPLYPointSource::_free()
{
    if (m_ply)
    {
        ply_close(m_ply);
        m_ply = nullptr;
    }
    if (m_elements)
    {
        for(int i = 0; i< m_numElements; ++i)
            free(m_elements[i]);
        free(m_elements);
        m_elements = nullptr;
    }
}

ColorPLYPointSource::~ColorPLYPointSource()
{
    _free();
}

bool ColorPLYPointSource::nextPoint(Point& p, Point3D<double>& color)
{
    if (m_current >= m_numInstances)
        return false;

    Vertex op;
    ply_get_element(m_ply, (void *)&op);
    p.p = op.point;
    p.n = op.normal;
    color = op.color;
    m_current++;
    return true;
}
