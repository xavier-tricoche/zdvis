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
#ifndef _ZD_LIB_POINT_HPP_
#define _ZD_LIB_POINT_HPP_

#include <stdlib.h>
#include <memory.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <typeinfo>

#include <stdarg.h> 

namespace ZD {
    /* point */
    template <typename T, unsigned int N>
    class CPoint {
    public:
        T m_data[N];

    public:
        CPoint() {
            memset(this->m_data, 0, sizeof(T)*N);
        }

        CPoint(const CPoint<T, N> &vec) {
            memcpy(this->m_data, vec.m_data, sizeof(T)*N);
        };

        CPoint(const T *vec) {
            memcpy(this->m_data, vec, sizeof(T)*N);
        };
        
        CPoint(const T& value) {
            std::fill(m_data, m_data+N, value);
        }
        
        CPoint(const T& value1, const T& value2) {
            static_assert(N >= 2, "Too many arguments in CPoint constructor");
            m_data[0] = value1;
            m_data[1] = value2;
            std::fill(m_data+2, m_data+N, static_cast<T>(0));
        }
        
        CPoint(const T& value1, const T& value2, const T& value3) {
            static_assert(N >= 3, "Too many arguments in CPoint constructor");
            m_data[0] = value1;
            m_data[1] = value2;
            m_data[2] = value3;
            std::fill(m_data+3, m_data+N, static_cast<T>(0));
        }
          
        CPoint(const T& value1, const T& value2, const T& value3, const T& value4) {
            static_assert(N >= 4, "Too many arguments in CPoint constructor");
            m_data[0] = value1;
            m_data[1] = value2;
            m_data[2] = value3;
            m_data[3] = value4;
            std::fill(m_data+4, m_data+N, static_cast<T>(0));
        }

        unsigned int Size() const { return N; }

        void SetZero() {
            memset(this->m_data, 0, sizeof(T)*N);
        }
        void SetNan() {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] = std::numeric_limits<T>::quiet_NaN();
            }
        }

        inline void Normalize() {
            T len = this->Length();
            if (len > 0.0) {
                for (unsigned int i = 0; i < N; ++i)
                    this->m_data[i] /= len;
            }
        }

        inline T Length() const {
            T len = 0.0;
            for (unsigned int i = 0; i < N; ++i)
                len += this->m_data[i] * this->m_data[i];
            len = std::sqrt(len);
            return len;
        }

        inline T Length2() const {
            T len = 0.0;
            for (unsigned int i = 0; i < N; ++i)
                len += this->m_data[i] * this->m_data[i];
            return len;
        }

        //
        T & operator[] (unsigned int n) {
            return this->m_data[n];
        }

        const T & operator[] (unsigned int n) const {
            return this->m_data[n];
        }

        //
        CPoint<T, N> & operator= (const CPoint<T, N> & rhs) {
            memcpy(this->m_data, rhs.m_data, sizeof(T)*N);
            return *this;
        }

        CPoint<T, N> & operator+= (const CPoint<T, N> & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] += rhs.m_data[i];
            }
            return *this;
        }
        
        CPoint<T, N> & operator-= (const CPoint<T, N> & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] -= rhs.m_data[i];
            }
            return *this;
        }
        
        CPoint<T, N> & operator*= (const CPoint<T, N> & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] *= rhs.m_data[i];
            }
            return *this;
        }

        CPoint<T, N> & operator*= (const T & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] *= rhs;
            }
            return *this;
        }
        
        CPoint<T, N> & operator/= (const CPoint<T, N> & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] /= rhs.m_data[i];
            }
            return *this;
        }

        CPoint<T, N> & operator/= (const T & rhs) {
            for (unsigned int i = 0; i < N; ++i) {
                this->m_data[i] /= rhs;
            }
            return *this;
        }
    };

    /* operators */
    template <typename T, unsigned int N> inline
    CPoint<T, N> operator+(const CPoint<T, N> & v1, const CPoint<T, N> & v2) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = v1.m_data[i] + v2.m_data[i];
        return v;
    }

    template <typename T, unsigned int N> inline
        CPoint<T, N> operator-(const CPoint<T, N> & v1) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = -v1.m_data[i];
        return v;
    }

    template <typename T, unsigned int N> inline
    CPoint<T, N> operator-(const CPoint<T, N> & v1, const CPoint<T, N> & v2) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = v1.m_data[i] - v2.m_data[i];
        return v;
    }

    template <typename S, typename T, unsigned int N> inline
    CPoint<T, N> operator*(const S & s1, const CPoint<T, N> & v2) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = s1 * v2.m_data[i];
        return v;
    }

    template <typename S, typename T, unsigned int N> inline
    CPoint<T, N> operator*(const CPoint<T, N> & v1, const S & s2) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = v1.m_data[i] * s2;
        return v;
    }

    template <typename S, typename T, unsigned int N> inline
    CPoint<T, N> operator/(const CPoint<T, N> & v1, const S & s2) {
        CPoint<T, N> v;
        for (unsigned int i = 0; i < N; ++i)
            v.m_data[i] = v1.m_data[i] / s2;
        return v;
    }

    /* functions */
    template <typename T, unsigned int N> inline
    T InnerProduct(const CPoint<T, N>& v1, const CPoint<T, N>& v2)
    {
        T res = 0.0; 
        for (unsigned int i = 0; i < N; ++i) {
            res += v1.m_data[i] * v2.m_data[i];
        }
        return res;
    }

    template <typename T, unsigned int N> inline
    T InnerProduct(const CPoint<T, N> *v1, const CPoint<T, N> *v2)
    {
        T res = 0.0;
        for (unsigned int i = 0; i < N; ++i) {
            res += v1->m_data[i] * v2->m_data[i];
        }
        return res;
    }

    template <typename T, unsigned int N> inline
    T Distance(const CPoint<T, N>& v1, const CPoint<T, N>& v2)
    {
        CPoint<T, N> v = v1 - v2;
        return v.Length();
    }

    template <typename T, unsigned int N> inline
    T Distance(const CPoint<T, N> *v1, const CPoint<T, N> *v2)
    {
        CPoint<T, N> v = *v1 - *v2;
        return v.Length();
    }

    template <typename T> inline
    CPoint<T, 3> CrossProduct(const CPoint<T, 3>& v1, const CPoint<T, 3>& v2)
    {
        CPoint<T, 3> c;
        c.m_data[0] = v1.m_data[1] * v2.m_data[2] - v2.m_data[1] * v1.m_data[2];
        c.m_data[1] = v1.m_data[2] * v2.m_data[0] - v2.m_data[2] * v1.m_data[0];
        c.m_data[2] = v1.m_data[0] * v2.m_data[1] - v2.m_data[0] * v1.m_data[1];
        return c;
    }

    template <typename T> inline
    CPoint<T, 3> CrossProduct(const CPoint<T, 3> *v1, const CPoint<T, 3> *v2)
    {
        CPoint<T, 3> c;
        c.m_data[0] = v1->m_data[1] * v2->m_data[2] - v2->m_data[1] * v1->m_data[2];
        c.m_data[1] = v1->m_data[2] * v2->m_data[0] - v2->m_data[2] * v1->m_data[0];
        c.m_data[2] = v1->m_data[0] * v2->m_data[1] - v2->m_data[0] * v1->m_data[1];
        return c;
    }

#if _MSC_VER != 1600
    template <typename T, unsigned int N>
    using CVector = CPoint<T, N>;
#endif
}


#endif
