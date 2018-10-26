/*************************************************************************
zdvis: Lagrangian Visualization for Vector, Tensor, and Multifield Data.

Author: Zi'ang Ding

Copyright (c) 2016-2018, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/
#ifndef _ZD_LIB_GLYPH_HPP_
#define _ZD_LIB_GLYPH_HPP_

#include "utils/Base/ZD_Point.hpp"
#include "utils/define.hpp"
#include "utils/UnitSphere/ZD_UnitSphereLow.hpp"
#include "utils/UnitSphere/ZD_UnitSphereMedian.hpp"
#include "utils/UnitSphere/ZD_UnitSphereHigh.hpp"
#include "utils/ColorMap/ZD_ColorMap.hpp"

#include <teem/tijk.h>

namespace ZD {
    template <typename T>
    class CGlyph {
    public:
        CPoint<T, 3>   *m_pVertices;
        CPoint<T, 3>   *m_pNormals;
        CPoint<T, 3>   *m_pColors;
        CPoint<int, 3> *m_pFaces;        // triangles
        unsigned int m_vCount;            // vertex count
        unsigned int m_fCount;            // triangle count

    public:
        CGlyph();
        ~CGlyph();

    private:
        void CreateLowUnitSphere();
        void CreateMedianUnitSphere();
        void CreateHighUnitSphere();

        void CreateSHGlyph(T *sh, unsigned int order);
        void CreateSHSphere(T *sh, unsigned int order);
        void CreateDTIGlyph(T *dti);

    public:
        int CreateSHGlyphLow(T *sh, unsigned int order);    // create low resolution SH glyph
        int CreateSHGlyphMedian(T *sh, unsigned int order);    // create median resolution SH glyph
        int CreateSHGlyphHigh(T *sh, unsigned int order);    // create high resolution SH glyph
        int CreateSHSphereLow(T *sh, unsigned int order);
        int CreateSHSphereMedian(T *sh, unsigned int order);
        int CreateSHSphereHigh(T *sh, unsigned int order);
        int CreateDTIGlyphLow(T *dti);                        // create low resolution DTI glyph
        int CreateDTIGlyphMedian(T *dti);                        // create low resolution DTI glyph
        int CreateDTIGlyphHigh(T *dti);                        // create high resolution DTI glyph

        void Normalization();
    };
}

template<typename T>
ZD::CGlyph<T>::CGlyph()
{
    m_pVertices = nullptr;
    m_pNormals = nullptr;
    m_pColors = nullptr;
    m_pFaces = nullptr;
    m_vCount = 0;
    m_fCount = 0;
}


template<typename T>
ZD::CGlyph<T>::~CGlyph()
{
    SafeDelete(m_pVertices);
    SafeDelete(m_pNormals);
    SafeDelete(m_pColors);
    SafeDelete(m_pFaces);
    m_vCount = 0;
    m_fCount = 0;
}


template<typename T>
void ZD::CGlyph<T>::CreateLowUnitSphere()
{
    m_vCount = unit_sphere_low_vertex_count;
    m_fCount = unit_sphere_low_triangle_count;

    m_pVertices = new ZD::CPoint<T, 3>[m_vCount];
    m_pNormals = new ZD::CPoint<T, 3>[m_vCount];
    m_pColors = new ZD::CPoint<T, 3>[m_vCount];
    m_pFaces = new ZD::CPoint<int, 3>[m_fCount];

    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        m_pVertices[vid][0] = unit_sphere_low_vertices[vid*3+0];
        m_pVertices[vid][1] = unit_sphere_low_vertices[vid*3+1];
        m_pVertices[vid][2] = unit_sphere_low_vertices[vid*3+2];

        m_pNormals[vid] = m_pVertices[vid];

        m_pColors[vid][0] = std::abs(m_pVertices[vid][0]);
        m_pColors[vid][1] = std::abs(m_pVertices[vid][1]);
        m_pColors[vid][2] = std::abs(m_pVertices[vid][2]);
    }
    memcpy(m_pFaces, unit_sphere_low_triangles, sizeof(int)*m_fCount * 3);
}

template<typename T>
void ZD::CGlyph<T>::CreateMedianUnitSphere()
{
    m_vCount = unit_sphere_median_vertex_count;
    m_fCount = unit_sphere_median_triangle_count;

    m_pVertices = new ZD::CPoint<T, 3>[m_vCount];
    m_pNormals = new ZD::CPoint<T, 3>[m_vCount];
    m_pColors = new ZD::CPoint<T, 3>[m_vCount];
    m_pFaces = new ZD::CPoint<int, 3>[m_fCount];

    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        m_pVertices[vid][0] = unit_sphere_median_vertices[(vid * 3) + 0];
        m_pVertices[vid][1] = unit_sphere_median_vertices[(vid * 3) + 1];
        m_pVertices[vid][2] = unit_sphere_median_vertices[(vid * 3) + 2];

        m_pNormals[vid] = m_pVertices[vid];

        m_pColors[vid][0] = std::abs(m_pVertices[vid][0]);
        m_pColors[vid][1] = std::abs(m_pVertices[vid][1]);
        m_pColors[vid][2] = std::abs(m_pVertices[vid][2]);
    }
    memcpy(m_pFaces, unit_sphere_median_triangles, sizeof(int)*m_fCount * 3);
}

template<typename T>
void ZD::CGlyph<T>::CreateHighUnitSphere()
{
    m_vCount = unit_sphere_high_vertex_count;
    m_fCount = unit_sphere_high_triangle_count;

    m_pVertices = new ZD::CPoint<T, 3>[m_vCount];
    m_pNormals = new ZD::CPoint<T, 3>[m_vCount];
    m_pColors = new ZD::CPoint<T, 3>[m_vCount];
    m_pFaces = new ZD::CPoint<int, 3>[m_fCount];

    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        m_pVertices[vid][0] = unit_sphere_high_vertices[(vid*3)+0];
        m_pVertices[vid][1] = unit_sphere_high_vertices[(vid*3)+1];
        m_pVertices[vid][2] = unit_sphere_high_vertices[(vid*3)+2];

        m_pNormals[vid] = m_pVertices[vid];

        m_pColors[vid][0] = std::abs(m_pVertices[vid][0]);
        m_pColors[vid][1] = std::abs(m_pVertices[vid][1]);
        m_pColors[vid][2] = std::abs(m_pVertices[vid][2]);
    }
    memcpy(m_pFaces, unit_sphere_high_triangles, sizeof(int)*m_fCount*3);
}

namespace ZD {
template<>
void ZD::CGlyph<double>::CreateSHGlyph(double * sh, unsigned int order)
{
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        double phi = std::atan2(m_pVertices[vid][1], m_pVertices[vid][0]);
        double tmp = std::sqrt(m_pVertices[vid][0] * m_pVertices[vid][0] + m_pVertices[vid][1] * m_pVertices[vid][1]);
        double theta = std::atan2(tmp, m_pVertices[vid][2]);
        double length = tijk_eval_esh_d(sh, order, theta, phi);
        length = length < 0.0 ? 0.0 : length;
        m_pVertices[vid][0] *= length;
        m_pVertices[vid][1] *= length;
        m_pVertices[vid][2] *= length;
    }
}

template<>
void ZD::CGlyph<float>::CreateSHGlyph(float * sh, unsigned int order)
{
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        float phi = std::atan2(m_pVertices[vid][1], m_pVertices[vid][0]);
        float tmp = std::sqrt(m_pVertices[vid][0] * m_pVertices[vid][0] + m_pVertices[vid][1] * m_pVertices[vid][1]);
        float theta = std::atan2(tmp, m_pVertices[vid][2]);
        float length = tijk_eval_esh_f(sh, order, theta, phi);
        length = length < 0.0 ? 0.0 : length;
        m_pVertices[vid][0] *= length;
        m_pVertices[vid][1] *= length;
        m_pVertices[vid][2] *= length;
    }
}

template<>
void ZD::CGlyph<double>::CreateSHSphere(double *sh, unsigned int order)
{
    double maxLength = 0.0;
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        double phi = std::atan2(m_pVertices[vid][1], m_pVertices[vid][0]);
        double tmp = std::sqrt(m_pVertices[vid][0] * m_pVertices[vid][0] + m_pVertices[vid][1] * m_pVertices[vid][1]);
        double theta = std::atan2(tmp, m_pVertices[vid][2]);
        double length = tijk_eval_esh_d(sh, order, theta, phi);
        length = length < 0.0 ? 0.0 : length;
        maxLength = length > maxLength ? length : maxLength;
        m_pColors[vid][0] = length;
    }

    CColorMapRainbow<double> colorMap(0.0, maxLength);
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        m_pColors[vid] = colorMap.GetColor(m_pColors[vid][0]);
    }
}

template<>
void ZD::CGlyph<float>::CreateSHSphere(float *sh, unsigned int order)
{
    float maxLength = 0.0;
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        float phi = std::atan2(m_pVertices[vid][1], m_pVertices[vid][0]);
        float tmp = std::sqrt(m_pVertices[vid][0] * m_pVertices[vid][0] + m_pVertices[vid][1] * m_pVertices[vid][1]);
        float theta = std::atan2(tmp, m_pVertices[vid][2]);
        float length = tijk_eval_esh_f(sh, order, theta, phi);
        length = length < 0.0 ? 0.0 : length;
        maxLength = length > maxLength ? length : maxLength;
        m_pColors[vid][0] = length;
    }

    CColorMapRainbow<float> colorMap(0.0, maxLength);
    for (unsigned int vid = 0; vid < m_vCount; ++vid) {
        m_pColors[vid] = colorMap.GetColor(m_pColors[vid][0]);
    }
}
}

template<typename T>
void ZD::CGlyph<T>::CreateDTIGlyph(T * dti)
{

}

template<typename T>
int ZD::CGlyph<T>::CreateSHGlyphLow(T * sh, unsigned int order)
{
    CreateLowUnitSphere();

    CreateSHGlyph(sh, order);

    return 0;
}

template<typename T>
int ZD::CGlyph<T>::CreateSHGlyphMedian(T * sh, unsigned int order)
{
    CreateMedianUnitSphere();

    CreateSHGlyph(sh, order);

    return 0;
}

template<typename T>
int ZD::CGlyph<T>::CreateSHGlyphHigh(T * sh, unsigned int order)
{
    CreateHighUnitSphere();

    CreateSHGlyph(sh, order);

    return 0;
}

template<typename T>
int ZD::CGlyph<T>::CreateSHSphereLow(T * sh, unsigned int order)
{
    CreateLowUnitSphere();

    CreateSHSphere(sh, order);

    return 0;
}

template<typename T>
int ZD::CGlyph<T>::CreateSHSphereMedian(T * sh, unsigned int order)
{
    CreateMedianUnitSphere();

    CreateSHSphere(sh, order);

    return 0;
}


template<typename T>
int ZD::CGlyph<T>::CreateSHSphereHigh(T * sh, unsigned int order)
{
    CreateHighUnitSphere();

    CreateSHSphere(sh, order);

    return 0;
}


template<typename T>
int ZD::CGlyph<T>::CreateDTIGlyphLow(T * dti)
{
    CreateLowUnitSphere();
    CreateDTIGlyph(dti);
    return 0;
}

template<typename T>
int ZD::CGlyph<T>::CreateDTIGlyphMedian(T * dti)
{
    CreateMedianUnitSphere();
    CreateDTIGlyph(dti);
    return 0;
}


template<typename T>
int ZD::CGlyph<T>::CreateDTIGlyphHigh(T * dti)
{
    CreateHighUnitSphere();
    CreateDTIGlyph(dti);
    return 0;
}


template<typename T>
void ZD::CGlyph<T>::Normalization()
{
    if (m_pVertices != nullptr) {
        T maxLen = 0.0;
        for (unsigned int i = 0; i < m_vCount; ++i) {
            T len = m_pVertices[i].Length();
            if (len > maxLen)
                maxLen = len;
        }

        for (unsigned int i = 0; i < m_vCount; ++i) {
            m_pVertices[i] = m_pVertices[i] / maxLen;
        }
    }
}


#endif
