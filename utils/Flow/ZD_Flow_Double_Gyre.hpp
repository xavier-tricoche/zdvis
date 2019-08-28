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
#ifndef _ZD_LIB_FLOW_DOUBLE_GYRE_HPP_
#define _ZD_LIB_FLOW_DOUBLE_GYRE_HPP_

#include "ZD_Flow.hpp"
#include "utils/define.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"

#include <math.h>

namespace ZD {
    template<typename T>
    class CFlowDoubleGyre : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        const value_type m_A;
        const value_type m_e;
        const value_type m_w;

    public:
        CFlowDoubleGyre();
        CFlowDoubleGyre(const value_type A, const value_type e, const value_type w);

    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& _time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& _time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
    };
} // namespace ZD

template<typename T>
inline ZD::CFlowDoubleGyre<T>::CFlowDoubleGyre() : m_A(0.1), m_e(0.25), m_w(0.6283185307179586476925286766559)
{
}

template<typename T>
ZD::CFlowDoubleGyre<T>::CFlowDoubleGyre(const value_type A, const value_type e, const value_type w) : m_A(A), m_e(e), m_w(w)
{
}

template<typename T>
inline const bool ZD::CFlowDoubleGyre<T>::CheckPosition(const point_type& p) const
{
    return true;
}

template<typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowDoubleGyre<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowDoubleGyre<T>::Velocity(const point_type& p, const value_type & _time) const
{
    const value_type& x = p[0];
    const value_type& y = p[1];
    
    value_type sin1 = std::sin(m_w*_time);
    value_type cos1 = std::cos(ZD_PI * y);

    value_type fxt = m_e * sin1 * x * x + (1.0 - 2.0 * m_e * sin1) * x;

    value_type uu = -1.0 * ZD_PI * m_A * std::sin(ZD_PI * fxt) * cos1;
    value_type vv = ZD_PI * m_A * sin1 * std::cos(ZD_PI * fxt);
    vv = vv * (2.0 * x * m_e * sin1 + 1.0 - 2.0 * m_e * sin1);

    return point_type(uu, vv);
}

template<typename T>
inline const T ZD::CFlowDoubleGyre<T>::Scalar(const point_type& p, const value_type & _time) const
{
    const value_type x = p[0];
    const value_type y = p[1];

    const value_type fxt = m_e * std::sin(m_w * _time) * x * x + (1.0 - 2.0 * m_e * std::sin(m_w * _time)) * x;
    const value_type gxt = 2.0 * m_e * std::sin(m_w * _time) * x + (1.0 - 2.0 * m_e * std::sin(m_w * _time));
    const value_type C = std::cos(ZD_PI * y);
    const value_type S = std::sin(ZD_PI * y);
    const value_type CP = std::cos(ZD_PI * fxt);
    const value_type SP = std::sin(ZD_PI * fxt);

    const value_type dudx = -1.0 * ZD_PI * ZD_PI * m_A * gxt * C * CP;
    const value_type dudy = ZD_PI * ZD_PI * m_A * S * SP;
    const value_type dvdx = ZD_PI * m_A * 2.0 * m_e * std::sin(m_w * _time) * CP * S - ZD_PI * ZD_PI * m_A * gxt * S * SP;
    const value_type dvdy = ZD_PI * ZD_PI * m_A * gxt * CP * S;

    const value_type s2 = (dudx - dvdy) * (dudx - dvdy) + (dudy + dvdx) * (dudy + dvdx);
    const value_type w2 = (dvdx - dudy) * (dvdx - dudy);
    return (s2 - w2);
}

template<typename T>
inline void ZD::CFlowDoubleGyre<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = 0.0;
    min[1] = 0.0;
    max[0] = 2.0;
    max[1] = 1.0;

    return;
}

#endif
