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
#ifndef _ZD_LIB_INTEGRATOR_RK4_HPP_
#define _ZD_LIB_INTEGRATOR_RK4_HPP_

#include "ZD_Integrator.hpp"

namespace ZD {
    template <typename T, unsigned int N>
    class CIntegratorRK4 : public CIntegrator<T, N> {
    private:
        T m_stepSize;

    public:
        CIntegratorRK4(const T stepSize);
        ~CIntegratorRK4();

    public:
        virtual T Next(const CPoint<T, N> &p, const T &t,
            const T direction, const CDynamicSystem<T, N> *sys);

        virtual CPoint<T, N> v(const T& t) const;
    };
}

template<typename T, unsigned int N>
ZD::CIntegratorRK4<T, N>::CIntegratorRK4(const T stepSize)
{
    this->m_stepSize = stepSize;
}

template<typename T, unsigned int N>
ZD::CIntegratorRK4<T, N>::~CIntegratorRK4()
{

}

template<typename T, unsigned int N>
T ZD::CIntegratorRK4<T, N>::Next(const CPoint<T, N> &p, const T &t,
    const T direction, const CDynamicSystem<T, N> *sys)
{
    CPoint<T, N> k1, k2, k3, k4;
    CPoint<T, N> p1, p2, p3, np;

    k1 = sys->Velocity(p, t);
    p1 = p + k1 * (m_stepSize * 0.5 * direction);
    k2 = sys->Velocity(p1, t + m_stepSize * 0.5 * direction);
    p2 = p + k2 * (m_stepSize * 0.5 * direction);
    k3 = sys->Velocity(p2, t + m_stepSize * 0.5 * direction);
    p3 = p + k3 * m_stepSize * direction;
    k4 = sys->Velocity(p3, t + m_stepSize * direction);
    
    np = p + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (m_stepSize * direction / 6.0);

    this->m_step.m_t0 = t;
    this->m_step.m_t1 = t + m_stepSize * direction;
    this->m_step.m_v0 = p;
    this->m_step.m_v1 = np;

    //printf("RK4: %f %f %f %f %f %f %f %f\n", k1[0], k1[1], k2[0], k2[1], k3[0], k3[1], k4[0], k4[1]);

    return m_stepSize;
}


template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CIntegratorRK4<T, N>::v(const T& t) const
{
    const T a = (t - this->m_step.m_t0) / (this->m_step.m_t1 - this->m_step.m_t0);

    return this->m_step.m_v0 * (1.0 - a) + this->m_step.m_v1 * a;
}


#endif
