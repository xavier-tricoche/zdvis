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
#ifndef _ZD_LIB_INTEGRATOR_EULER_HPP_
#define _ZD_LIB_INTEGRATOR_EULER_HPP_

#include "ZD_Integrator.hpp"

namespace ZD {
    template <typename T, unsigned int N>
    class CIntegratorEuler : public CIntegrator<T, N> {
    private:
        T m_stepSize;

    public:
        CIntegratorEuler(const T stepSize);
        ~CIntegratorEuler();

    public:
        virtual T Next(const CPoint<T, N> &p, const T &t,
            const T direction, const CDynamicSystem<T, N> *sys);

        virtual CPoint<T, N> v(const T& t) const;
    };
}

template<typename T, unsigned int N>
ZD::CIntegratorEuler<T, N>::CIntegratorEuler(const T stepSize)
{
    this->m_stepSize = stepSize;
}

template<typename T, unsigned int N>
ZD::CIntegratorEuler<T, N>::~CIntegratorEuler()
{
}

template<typename T, unsigned int N>
T ZD::CIntegratorEuler<T, N>::Next(const CPoint<T, N>& p, const T & t, const T direction, const CDynamicSystem<T, N>* sys)
{
    CPoint<T, N> v = sys->Velocity(p, t);
    CPoint<T, N> np = p + v * m_stepSize * direction;

    this->m_step.m_t0 = t;
    this->m_step.m_t1 = t + m_stepSize * direction;
    this->m_step.m_v0 = p;
    this->m_step.m_v1 = np;

    return m_stepSize;
}

template<typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CIntegratorEuler<T, N>::v(const T & t) const
{
    const T a = (t - this->m_step.m_t0) / (this->m_step.m_t1 - this->m_step.m_t0);

    return this->m_step.m_v0 * (1.0 - a) + this->m_step.m_v1 * a;
}

#endif
