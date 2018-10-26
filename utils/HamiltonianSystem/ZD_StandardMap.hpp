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
#ifndef _ZD_LIB_STANDARD_MAP_HPP_
#define _ZD_LIB_STANDARD_MAP_HPP_

#include "ZD_HamiltonianSystem.hpp"

namespace ZD {
    template <typename T>
    class CStandardMap : public CHamiltonianSystem<T, 2> {
    public:
        typedef CPoint<T, 2> vec2;

    private:
        const T m_k;

    public:
        CStandardMap(const T k);
        ~CStandardMap();

    public:
        /* integral for a fixed time */
        virtual T NextTime(vec2& p, const T startTime, const T time, const T direction) const;
        virtual T NextTime(vec2& p, const T startTime, const T time, 
            const T direction, std::vector<vec2>& buffer) const;
        
        /* integral for a fixed length */
        virtual T NextLength(vec2& p, const T startTime, const T length, const T direction) const;
        virtual T NextLength(vec2& p, const T startTime, const T length, const T direction,
            std::vector<vec2> &buffer) const;

        /* integral for one map iteration */
        virtual T NextMap(vec2& p, const T startTime, const T direction,
            const CPoincarePlane<T, 2> *plane) const;
        virtual T NextMap(vec2& p, const T startTime, const T direction,
            const CPoincarePlane<T, 2> *plane, std::vector<vec2> &buffer) const;

        virtual CIntegrator<T, 2> *CreateIntegrator() const;

        virtual const vec2 Velocity(const vec2& p, const T& time) const;

    public:
        inline T my_mod(T a, const T b) const;
        inline vec2 my_mod(vec2 p, const T b) const;
    };
}

template<typename T>
ZD::CStandardMap<T>::CStandardMap(const T k) : m_k(k)
{
}

template<typename T>
ZD::CStandardMap<T>::~CStandardMap()
{
}

template<typename T>
T ZD::CStandardMap<T>::NextTime(vec2 & p, const T startTime, const T time, const T direction) const
{
    return T();
}

template<typename T>
T ZD::CStandardMap<T>::NextTime(vec2 & p, const T startTime, const T time, const T direction, std::vector<vec2>& buffer) const
{
    return T();
}

template<typename T>
T ZD::CStandardMap<T>::NextLength(vec2 & p, const T startTime, const T length, const T direction) const
{
    return T();
}

template<typename T>
T ZD::CStandardMap<T>::NextLength(vec2 & p, const T startTime, const T length, const T direction, std::vector<vec2>& buffer) const
{
    return T();
}

template<typename T>
T ZD::CStandardMap<T>::NextMap(vec2 & p, const T startTime, const T direction, const CPoincarePlane<T, 2>* plane) const
{
    if (direction > 0.0) {
        T xn = p[0];
        T pn = p[1];

        p[1] = pn + m_k * std::sin(xn);            // p_{n+1} = p_{n} + k * sin(x_{n})
        p[0] = xn + p[1];                        // x_{n+1} = x_{x} + p_{n+1}

        return 1.0;
    }
    else {
        T xn = p[0];
        T pn = p[1];

        p[0] = xn - pn;                            // x_{n-1} = x_{n} - p_{n}
        p[1] = pn - m_k * std::sin(p[0]);        // p_{n-1} = p_{n} - k * sin(x_{n-1})

        return -1.0;
    }
}

template<typename T>
T ZD::CStandardMap<T>::NextMap(vec2 & p, const T startTime, const T direction, const CPoincarePlane<T, 2>* plane, std::vector<vec2>& buffer) const
{
    if (direction > 0.0) {
        T xn = p[0];
        T pn = p[1];

        p[1] = pn + m_k * std::sin(xn);            // p_{n+1} = p_{n} + k * sin(x_{n})
        p[0] = xn + p[1];                        // x_{n+1} = x_{x} + p_{n+1}
        
        buffer.push_back(p);

        return 1.0;
    }
    else {
        T xn = p[0];
        T pn = p[1];

        p[0] = xn - pn;                            // x_{n-1} = x_{n} - p_{n}
        p[1] = pn - m_k * std::sin(p[0]);        // p_{n-1} = p_{n} - k * sin(x_{n-1})

        buffer.push_back(p);
        return -1.0;
    }
}

template<typename T>
ZD::CIntegrator<T, 2>* ZD::CStandardMap<T>::CreateIntegrator() const
{
    return nullptr;
}


template<typename T>
const typename ZD::CStandardMap<T>::vec2 ZD::CStandardMap<T>::Velocity(const vec2 & p, const T & time) const
{
    return vec2(0.0, 0.0);
}

template<typename T>
inline T ZD::CStandardMap<T>::my_mod(T a, const T b) const
{
    while (a < 0.0)
        a += b;
    return std::fmod(a, b);
}

template<typename T>
inline typename ZD::CStandardMap<T>::vec2 ZD::CStandardMap<T>::my_mod(vec2 p, const T b) const
{
    vec2 pp;
    pp[0] = this->my_mod(p[0], b);
    pp[1] = this->my_mod(p[1], b);
    return pp;
}


#endif
