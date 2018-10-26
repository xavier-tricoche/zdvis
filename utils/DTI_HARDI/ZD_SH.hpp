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
#ifndef _ZD_LIB_SH_HPP_
#define _ZD_LIB_SH_HPP_

#include "utils/Base/ZD_Field.hpp"
#include "utils/define.hpp"
#include "utils/UnitSphere/ZD_UnitSphereHigh.hpp"

#include "ZD_HOT.hpp"

#include <teem/tijk.h>

namespace ZD {
    template <typename T> 
    class CSH {
    private:
        CField3<T, 15> *m_pSH4Volume;        // 4th order
        CField3<T, 28> *m_pSH6Volume;        // 6th order
        unsigned int m_order;

        static T m_weights[45];

    public:
        CSH(const int w, const int h, const int d, const unsigned int order);
        CSH(const char *shPathname);
        ~CSH();

    public:
        const unsigned int GetOrder() const;
        const int GetSize(const int k) const;

    public:
        void GetSH(const CPoint<int, 3> p, T *sh) const;
        void GetSH(const int x, const int y, const int z, T *sh) const;
        void GetSH(const CPoint<T, 3> p, T *sh) const;
        void GetSH(const T x, const T y, const T z, T *sh) const;

    public:
        void SetSH(const CPoint<int, 3> p, T *sh);
        void SetSH(const int x, const int y, const int z, T *sh);

    public:
        static const T Evaluate(T *sh, const int order, const T theta, const T phi);
        static const T Evaluate(T *sh, const int order, const CPoint<T, 3>& v);
        static const T MaxValue(T *sh, const int order, CPoint<T, 3> *dir = nullptr);
    };
} // namespace ZD

template<typename T>
ZD::CSH<T>::CSH(const int w, const int h, const int d, const unsigned int order)
{
    if (order == 4) {
        m_pSH4Volume = new CField3<T, 15>();
        m_pSH4Volume->CreateField(w, h, d, nullptr);
        m_order = 4;
    } else if (order == 6) {
        m_pSH6Volume = new CField3<T, 28>();
        m_pSH6Volume->CreateField(w, h, d, nullptr);
        m_order = 6;
    }
    else {
        ;
    }
}

template<typename T>
ZD::CSH<T>::CSH(const char * shPathname)
{
    m_pSH4Volume = nullptr;
    m_pSH6Volume = nullptr;
    m_order = 0;

    Nrrd *field = nrrdNew();
    if (nrrdLoad(field, shPathname, NULL)) {
        printf("cannot load file %s:\n%s\n", shPathname, biffGetDone(NRRD));
        exit(-1);
    }
    unsigned int n = field->axis[0].size;
    nrrdNuke(field);

    if (n == 15) {
        m_pSH4Volume = new CField3<T, 15>();
        m_pSH4Volume->OpenNrrdFile(shPathname);
        m_order = 4;
    }
    else if (n == 28) {
        m_pSH6Volume = new CField3<T, 28>();
        m_pSH6Volume->OpenNrrdFile(shPathname);
        m_order = 6;
    }
    else {
        printf("unknown order. (n=%d).\n", n);
        exit(-1);
    }
}

template<typename T>
ZD::CSH<T>::~CSH()
{
    SafeDelete(m_pSH4Volume);
    SafeDelete(m_pSH6Volume);
    m_order = 0;
}

template<typename T>
const unsigned int ZD::CSH<T>::GetOrder() const
{
    return m_order;
}

template<typename T>
const int ZD::CSH<T>::GetSize(const int k) const
{
    if (m_order == 4)
        return m_pSH4Volume->GetSize(k);
    else if (m_order == 6)
        return m_pSH6Volume->GetSize(k);
    else
        return -1;
}


template<typename T>
void ZD::CSH<T>::GetSH(const CPoint<int, 3> p, T * sh) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pSH4Volume->GetValue(p);
        memcpy(sh, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pSH6Volume->GetValue(p);
        memcpy(sh, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CSH<T>::GetSH(const int x, const int y, const int z, T * sh) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pSH4Volume->GetValue(x, y, z);
        memcpy(sh, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pSH6Volume->GetValue(x, y, z);
        memcpy(sh, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CSH<T>::GetSH(const CPoint<T, 3> p, T * sh) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pSH4Volume->GetValue(p);
        memcpy(sh, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pSH6Volume->GetValue(p);
        memcpy(sh, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CSH<T>::GetSH(const T x, const T y, const T z, T * sh) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pSH4Volume->GetValue(x, y, z);
        memcpy(sh, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pSH6Volume->GetValue(x, y, z);
        memcpy(sh, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}


template<typename T>
void ZD::CSH<T>::SetSH(const CPoint<int, 3> p, T * sh)
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        memcpy(tmp.m_data, sh, sizeof(T) * 15);
        m_pSH4Volume->SetValue(p, tmp);
    } else if (m_order == 6) {
        CPoint<T, 28> tmp;
        memcpy(tmp.m_data, sh, sizeof(T) * 28);
        m_pSH6Volume->SetValue(p, tmp);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CSH<T>::SetSH(const int x, const int y, const int z, T * sh)
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        memcpy(tmp.m_data, sh, sizeof(T) * 15);
        m_pSH4Volume->SetValue(x, y, z, tmp);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        memcpy(tmp.m_data, sh, sizeof(T) * 28);
        m_pSH6Volume->SetValue(x, y, z, tmp);
    }
    else {
        ;
    }
}

template<typename T>
const T ZD::CSH<T>::Evaluate(T * sh, const int order, const T theta, const T phi)
{
    T basis[45];
    T stheta = std::sin(theta);
    T ctheta = std::cos(theta);
    T cstheta = ctheta * stheta;
    T stheta2 = stheta * stheta, stheta4 = stheta2 * stheta2, stheta6 = stheta4 * stheta2, stheta8 = stheta4 * stheta4;
    T ctheta2 = ctheta * ctheta, ctheta4 = ctheta2 * ctheta2, ctheta6 = ctheta4 * ctheta2, ctheta8 = ctheta4 * ctheta4;
    T sphi = std::sin(phi);
    T cphi = std::cos(phi);
    T s2phi = 2.0 * cphi * sphi;
    T c2phi = cphi * cphi - sphi * sphi;
    T s3phi = s2phi * cphi + c2phi * sphi;
    T c3phi = c2phi * cphi - s2phi * sphi;
    T s4phi = 2.0 * c2phi * s2phi;
    T c4phi = c2phi * c2phi - s2phi * s2phi;
    T s5phi = s4phi * cphi + c4phi * sphi;
    T c5phi = c4phi * cphi - s4phi * sphi;
    T s6phi = s4phi * c2phi + c4phi * s2phi;
    T c6phi = c4phi * c2phi - s4phi * s2phi;
    T s7phi = s4phi * c3phi + c4phi * s3phi;
    T c7phi = c4phi * c3phi - s4phi * s3phi;
    T s8phi = 2.0 * c4phi * s4phi;
    T c8phi = c4phi * c4phi - s4phi * s4phi;


    basis[0] = m_weights[0];
    basis[5] = basis[1] = m_weights[1] * stheta2;
    basis[5] *= s2phi; basis[1] *= c2phi;
    basis[4] = basis[2] = m_weights[2] * cstheta;
    basis[4] *= sphi; basis[2] *= cphi;
    basis[3] = m_weights[3] * (3.0*ctheta2 - 1.0);
    basis[14] = basis[6] = m_weights[6] * stheta4;
    basis[14] *= s4phi; basis[6] *= c4phi;
    basis[13] = basis[7] = m_weights[7] * cstheta*stheta2;
    basis[13] *= s3phi; basis[7] *= c3phi;
    basis[12] = basis[8] = m_weights[8] * (7 * ctheta2 - 1)*stheta2;
    basis[12] *= s2phi; basis[8] *= c2phi;
    basis[11] = basis[9] = m_weights[9] * (7 * ctheta2*ctheta - 3 * ctheta)*stheta;
    basis[11] *= sphi; basis[9] *= cphi;
    basis[10] = m_weights[10] * (35.0*ctheta4 - 30.0*ctheta2 + 3.0);
    basis[27] = basis[15] = m_weights[15] * stheta6;
    basis[27] *= s6phi; basis[15] *= c6phi;
    basis[26] = basis[16] = m_weights[16] * stheta4*cstheta;
    basis[26] *= s5phi; basis[16] *= c5phi;
    basis[25] = basis[17] = m_weights[17] * stheta4*(11 * ctheta2 - 1.0);
    basis[25] *= s4phi; basis[17] *= c4phi;
    basis[24] = basis[18] = m_weights[18] * stheta2*stheta*(11 * ctheta2*ctheta - 3 * ctheta);
    basis[24] *= s3phi; basis[18] *= c3phi;
    basis[23] = basis[19] = m_weights[19] * stheta2*(33 * ctheta4 - 18 * ctheta2 + 1.0);
    basis[23] *= s2phi; basis[19] *= c2phi;
    basis[22] = basis[20] = m_weights[20] * stheta*(33 * ctheta4*ctheta - 30.0*ctheta2*ctheta + 5 * ctheta);
    basis[22] *= sphi; basis[20] *= cphi;
    basis[21] = m_weights[21] * (231 * ctheta6 - 315 * ctheta4 + 105 * ctheta2 - 5.0);
    basis[44] = basis[28] = m_weights[28] * stheta8;
    basis[44] *= s8phi; basis[28] *= c8phi;
    basis[43] = basis[29] = m_weights[29] * stheta6*cstheta;
    basis[43] *= s7phi; basis[29] *= c7phi;
    basis[42] = basis[30] = m_weights[30] * stheta6*(15 * ctheta2 - 1);
    basis[42] *= s6phi; basis[30] *= c6phi;
    basis[41] = basis[31] = m_weights[31] * stheta4*cstheta*(5 * ctheta2 - 1);
    basis[41] *= s5phi; basis[31] *= c5phi;
    basis[40] = basis[32] = m_weights[32] * stheta4*(65 * ctheta4 - 26 * ctheta2 + 1);
    basis[40] *= s4phi; basis[32] *= c4phi;
    basis[39] = basis[33] = m_weights[33] * stheta2*cstheta*(39 * ctheta4 - 26 * ctheta2 + 3);
    basis[39] *= s3phi; basis[33] *= c3phi;
    basis[38] = basis[34] = m_weights[34] * stheta2*(143 * ctheta6 - 143 * ctheta4 + 33 * ctheta2 - 1);
    basis[38] *= s2phi; basis[34] *= c2phi;
    basis[37] = basis[35] = m_weights[35] * cstheta*(715 * ctheta6 - 1001 * ctheta4 + 385 * ctheta2 - 35);
    basis[37] *= sphi; basis[35] *= cphi;
    basis[36] = m_weights[36] * (6435 * ctheta8 - 12012 * ctheta6 + 6930 * ctheta4 - 1260 * ctheta2 + 35);

    T res = 0.0;
    for (unsigned int i = 0; i < 45; i++)
        res += basis[i] * sh[i];

    return res;
}

template<typename T>
const T ZD::CSH<T>::Evaluate(T * sh, const int order, const CPoint<T, 3>& v)
{
    const T phi = std::atan2(v[1], v[0]);
    const T temp = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    const T theta = std::atan2(temp, v[2]);
    return Evaluate(sh, order, theta, phi);
}


template<typename T>
const T ZD::CSH<T>::MaxValue(T * sh, const int order, CPoint<T, 3>* dir)
{
    T maxValue = 0.0;
    CPoint<T, 3> maxDir;
    for (unsigned int i = 0; i < unit_sphere_high_vertex_count; ++i) {
        CPoint<T, 3> v(unit_sphere_high_vertices[i * 3 + 0], unit_sphere_high_vertices[i * 3 + 1], unit_sphere_high_vertices[i * 3 + 2]);
        T value = Evaluate(sh, order, v);
        if (value > maxValue) {
            maxValue = value;
            maxDir = v;
        }
    }
    if (dir)
        *dir = maxDir;
    return maxValue;
}

template<>
float ZD::CSH<float>::m_weights[45] = {
     static_cast<float>(1.0) / sqrt(static_cast<float>(4.0*AIR_PI)),

     std::sqrt(static_cast<float>(15.0 / (16.0*AIR_PI))),
    -std::sqrt(static_cast<float>(15.0 / (4.0*AIR_PI))),
     std::sqrt(static_cast<float>(5.0 / (16.0*AIR_PI))),
    -std::sqrt(static_cast<float>(15.0 / (4.0*AIR_PI))),
     std::sqrt(static_cast<float>(15.0 / (16.0*AIR_PI))),

     std::sqrt(static_cast<float>(315.0 / (256.0*AIR_PI))),
    -std::sqrt(static_cast<float>(315.0 / (32.0*AIR_PI))),
     std::sqrt(static_cast<float>(45.0 / (64.0*AIR_PI))),
    -std::sqrt(static_cast<float>(45.0 / (32.0*AIR_PI))),
     std::sqrt(static_cast<float>(9.0 / (256.0*AIR_PI))),
    -std::sqrt(static_cast<float>(45.0 / (32.0*AIR_PI))),
     std::sqrt(static_cast<float>(45.0 / (64.0*AIR_PI))),
    -std::sqrt(static_cast<float>(315.0 / (32.0*AIR_PI))),
     std::sqrt(static_cast<float>(315.0 / (256.0*AIR_PI))),

    static_cast<float>( 1.0 / 64.0*std::sqrt(6006.0 / AIR_PI)),
    static_cast<float>(-3.0 / 32.0*std::sqrt(2002.0 / AIR_PI)),
    static_cast<float>( 3.0 / 32.0*std::sqrt(91.0 / AIR_PI)),
    static_cast<float>(-1.0 / 32.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>( 1.0 / 64.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(-1.0 / 16.0*std::sqrt(273.0 / AIR_PI)),
    static_cast<float>( 1.0 / 32.0*std::sqrt(13.0 / AIR_PI)),
    static_cast<float>(-1.0 / 16.0*std::sqrt(273.0 / AIR_PI)),
    static_cast<float>( 1.0 / 64.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(-1.0 / 32.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>( 3.0 / 32.0*std::sqrt(91.0 / AIR_PI)),
    static_cast<float>(-3.0 / 32.0*std::sqrt(2002.0 / AIR_PI)),
    static_cast<float>( 1.0 / 64.0*std::sqrt(6006.0 / AIR_PI)),

    static_cast<float>( 3.0 / 256.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>( 1.0 / 128.0*sqrt(2 * 7293.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17017.0 / AIR_PI)),
    static_cast<float>( 3.0 / 128.0*sqrt(1309.0 / AIR_PI)),
    static_cast<float>(-1.0 / 64.0*sqrt(19635.0 / AIR_PI)),
    static_cast<float>( 3.0 / 128.0*sqrt(2 * 595.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>( 1.0 / 256.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>( 3.0 / 128.0*sqrt(2 * 595.0 / AIR_PI)),
    static_cast<float>(-1.0 / 64.0*sqrt(19635.0 / AIR_PI)),
    static_cast<float>( 3.0 / 128.0*sqrt(1309.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17017.0 / AIR_PI)),
    static_cast<float>( 1.0 / 128.0*sqrt(2 * 7293.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>( 3.0 / 256.0*sqrt(12155.0 / AIR_PI))
};

template<>
double ZD::CSH<double>::m_weights[45] = {
    static_cast<float>(1.0 / sqrt(4.0*AIR_PI)),

    static_cast<float>(std::sqrt(15.0 / (16.0*AIR_PI))),
    static_cast<float>(-std::sqrt(15.0 / (4.0*AIR_PI))),
    static_cast<float>(std::sqrt(5.0 / (16.0*AIR_PI))),
    static_cast<float>(-std::sqrt(15.0 / (4.0*AIR_PI))),
    static_cast<float>(std::sqrt(15.0 / (16.0*AIR_PI))),

    static_cast<float>(std::sqrt(315.0 / (256.0*AIR_PI))),
    static_cast<float>(-std::sqrt(315.0 / (32.0*AIR_PI))),
    static_cast<float>(std::sqrt(45.0 / (64.0*AIR_PI))),
    static_cast<float>(-std::sqrt(45.0 / (32.0*AIR_PI))),
    static_cast<float>(std::sqrt(9.0 / (256.0*AIR_PI))),
    static_cast<float>(-std::sqrt(45.0 / (32.0*AIR_PI))),
    static_cast<float>(std::sqrt(45.0 / (64.0*AIR_PI))),
    static_cast<float>(-std::sqrt(315.0 / (32.0*AIR_PI))),
    static_cast<float>(std::sqrt(315.0 / (256.0*AIR_PI))),

    static_cast<float>(1.0 / 64.0*std::sqrt(6006.0 / AIR_PI)),
    static_cast<float>(-3.0 / 32.0*std::sqrt(2002.0 / AIR_PI)),
    static_cast<float>(3.0 / 32.0*std::sqrt(91.0 / AIR_PI)),
    static_cast<float>(-1.0 / 32.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(1.0 / 64.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(-1.0 / 16.0*std::sqrt(273.0 / AIR_PI)),
    static_cast<float>(1.0 / 32.0*std::sqrt(13.0 / AIR_PI)),
    static_cast<float>(-1.0 / 16.0*std::sqrt(273.0 / AIR_PI)),
    static_cast<float>(1.0 / 64.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(-1.0 / 32.0*std::sqrt(2730.0 / AIR_PI)),
    static_cast<float>(3.0 / 32.0*std::sqrt(91.0 / AIR_PI)),
    static_cast<float>(-3.0 / 32.0*std::sqrt(2002.0 / AIR_PI)),
    static_cast<float>(1.0 / 64.0*std::sqrt(6006.0 / AIR_PI)),

    static_cast<float>(3.0 / 256.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>(1.0 / 128.0*sqrt(2 * 7293.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17017.0 / AIR_PI)),
    static_cast<float>(3.0 / 128.0*sqrt(1309.0 / AIR_PI)),
    static_cast<float>(-1.0 / 64.0*sqrt(19635.0 / AIR_PI)),
    static_cast<float>(3.0 / 128.0*sqrt(2 * 595.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>(1.0 / 256.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17.0 / AIR_PI)),
    static_cast<float>(3.0 / 128.0*sqrt(2 * 595.0 / AIR_PI)),
    static_cast<float>(-1.0 / 64.0*sqrt(19635.0 / AIR_PI)),
    static_cast<float>(3.0 / 128.0*sqrt(1309.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(17017.0 / AIR_PI)),
    static_cast<float>(1.0 / 128.0*sqrt(2 * 7293.0 / AIR_PI)),
    static_cast<float>(-3.0 / 64.0*sqrt(12155.0 / AIR_PI)),
    static_cast<float>(3.0 / 256.0*sqrt(12155.0 / AIR_PI))
};


#endif
