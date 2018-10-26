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
#ifndef _ZD_LIB_FLOW_STEADY_DOUBLE_GYRE_HPP_
#define _ZD_LIB_FLOW_STEADY_DOUBLE_GYRE_HPP_

#include "ZD_Flow.hpp"
#include "utils/define.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"

#include <math.h>

namespace ZD {
    template<typename T>
    class CFlowSteadyDoubleGyre : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        const value_type m_A;

    public:
        CFlowSteadyDoubleGyre();
        CFlowSteadyDoubleGyre(const value_type A);

    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
    };
} // namespace ZD

template<typename T>
inline ZD::CFlowSteadyDoubleGyre<T>::CFlowSteadyDoubleGyre() : m_A(0.1)
{
}

template<typename T>
ZD::CFlowSteadyDoubleGyre<T>::CFlowSteadyDoubleGyre(const value_type A) : m_A(A)
{
}


template<typename T>
inline const bool ZD::CFlowSteadyDoubleGyre<T>::CheckPosition(const point_type& p) const
{
    return true;
}

template<typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowSteadyDoubleGyre<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowSteadyDoubleGyre<T>::Velocity(const point_type& p, const value_type & time) const
{
    const value_type x = p[0];
    const value_type y = p[1];

    value_type u = -m_A * ZD_PI * std::sin(x * ZD_PI) * std::cos(y * ZD_PI);
    value_type v =  m_A * ZD_PI * std::cos(x * ZD_PI) * std::sin(y * ZD_PI);

    return point_type(u, v);
}

template<typename T>
inline const T ZD::CFlowSteadyDoubleGyre<T>::Scalar(const point_type& p, const value_type & time) const
{
    return 0.0;
}

template<typename T>
inline void ZD::CFlowSteadyDoubleGyre<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = 0.0;
    min[1] = 0.0;
    max[0] = 2.0;
    max[1] = 2.0;

    return;
}

#endif
