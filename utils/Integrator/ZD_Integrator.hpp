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
#ifndef _ZD_LIB_INTEGRATOR_HPP_
#define _ZD_LIB_INTEGRATOR_HPP_

#include "../Base/ZD_Point.hpp"
#include "../Base/ZD_DynamicSystem.hpp"

#include <exception>

namespace ZD {
    enum INTEGRATION_METHOD {
        IM_NONE,
        IM_EULER,
        IM_RK2,
        IM_RK4,
        IM_RK45
    };

    INTEGRATION_METHOD SetIntegrationMethod(const char *str);

    template <typename T, unsigned int N>
    class CIntegrator {
    public:
        class CStep {
        public:
            CPoint<T, N> m_p[5];
            CPoint<T, N> m_v0, m_v1;
            T m_t0, m_t1;

        public:
            CStep() {
            }
            ~CStep() {
            }
        };

    protected:
        CStep m_step;

    public:
        // dynamic system
        virtual ~CIntegrator() {}
        virtual T Next(const CPoint<T, N> &p, const T &t, 
            const T direction, const CDynamicSystem<T, N> *sys) = 0;

        virtual CPoint<T, N> v(const T& t) const = 0;

        const T& t0() const { return m_step.m_t0; }
        const T& t1() const { return m_step.m_t1; }
        const CPoint<T, N>& v0() const { return m_step.m_v0; }
        const CPoint<T, N>& v1() const { return m_step.m_v1; }
    };
}

#endif
