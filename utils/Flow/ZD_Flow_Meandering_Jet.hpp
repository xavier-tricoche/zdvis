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
#ifndef _ZD_LIB_FLOW_MEANDERING_JET_HPP_
#define _ZD_LIB_FLOW_MEANDERING_JET_HPP_

#include "ZD_Flow.hpp"
#include "utils/define.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"

#include <math.h>

namespace ZD {
    /* bounding box: [0.0 -- 10.0] x [-4.0 -- 4.0] */
    template <typename T>
    class CFlowMeanderingJet : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        const value_type m_B0;
        const value_type m_L;
        const value_type m_k;
        const value_type m_c;
        const value_type m_w;
        const value_type m_e;

    public:
        CFlowMeanderingJet();
        ~CFlowMeanderingJet();
        //CFlowMeanderingJet();

    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
    };
}


template <typename T>
ZD::CFlowMeanderingJet<T>::CFlowMeanderingJet() : m_B0(1.2), m_L(10.0), m_k(ZD_PI / 5.0), m_c(0.1), m_w(0.1), m_e(0.3)
{

}

template <typename T>
ZD::CFlowMeanderingJet<T>::~CFlowMeanderingJet()
{

}


template<typename T>
inline const bool ZD::CFlowMeanderingJet<T>::CheckPosition(const point_type& p) const
{
    return true;
}

template<typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowMeanderingJet<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}


template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowMeanderingJet<T>::Velocity(const point_type& p, const value_type & time) const
{
    const value_type x = p[0];
    const value_type y = p[1];

    const value_type t = time;
    const value_type A = std::sin(m_k * x);
    const value_type B = (m_B0 + m_e * cos(m_w * t));
    const value_type C = std::cos(m_k * x);
    const value_type temp_1 = m_k * m_k * A * A * B * B + 1.0;
    const value_type temp_2 = std::pow(std::tanh((y - C * B) / std::sqrt(temp_1)), 2) - 1.0;

    const value_type vv = (
        m_k * A * B / std::sqrt(temp_1) -
        (m_k * m_k * m_k * C * A * B * B * (y - C * B) / std::pow((double)temp_1, (double)1.5)
            )) * temp_2;

    const value_type uu = -m_c - temp_2 / std::sqrt(temp_1);

    return point_type(uu, vv);

    //const T x = p[0];
    //const T y = p[1];

    //T fxt = m_e * std::sin(m_w * time) * x * x + (1.0 - 2.0 * m_e * std::sin(m_w * time)) * x;

    //T uu = -1.0 * ZD_PI * m_A * std::sin(ZD_PI * fxt) * std::cos(ZD_PI * y);
    //T vv = ZD_PI * m_A * std::sin(ZD_PI * y) * std::cos(ZD_PI * fxt);
    //vv = vv * (2.0 * x * m_e * std::sin(m_w * time) + 1.0 - 2.0 * m_e * std::sin(m_w * time));

    //return CPoint<T, 2>(uu, vv);
}

template<typename T>
inline const T ZD::CFlowMeanderingJet<T>::Scalar(const point_type& p, const value_type & time) const
{
    return 0.0;
}

template<typename T>
inline void ZD::CFlowMeanderingJet<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = 0.0;
    min[1] = -4.0;
    max[0] = 10.0;
    max[1] = 4.0;

    return;
}

#endif
