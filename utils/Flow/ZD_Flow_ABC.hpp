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
#ifndef _ZD_LIB_FLOW_ABC_HPP_
#define _ZD_LIB_FLOW_ABC_HPP_

#include "ZD_Flow.hpp"
#include <utils/define.hpp>
#include <utils/Integrator/ZD_Integrator_RK4.hpp>

namespace ZD {
    template<typename T>
    class CFlowABC : public CFlow<T, 3> {
    private:
        const T m_A;
        const T m_B;
        const T m_C;
        const T m_q;
        const T m_w;

    public:
        CFlowABC();
        CFlowABC(const T A, const T B, const T C);
        ~CFlowABC();

    protected:
        virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        virtual inline const CPoint<T, 3> Velocity(const CPoint<T, 3>& p, const T& t) const;
        virtual inline const T Scalar(const CPoint<T, 3>& p, const T& time) const;
        virtual inline const bool CheckPosition(const CPoint<T, 3>& p) const;
        virtual inline void GetBBox(CPoint<T, 3>& min, CPoint<T, 3>& max) const;
    };
} // namespace ZD

template<typename T>
ZD::CFlowABC<T>::CFlowABC() 
    : m_A(1.7320508075688772935274463415059), 
      m_B(1.4142135623730950488016887242097), 
      m_C(1.0), m_q(0.1), m_w(2.0*ZD_PI)
{
}

template<typename T>
ZD::CFlowABC<T>::CFlowABC(const T A, const T B, const T C) 
    : m_A(A), m_B(B), m_C(C), m_q(0.1), m_w(2.0*ZD_PI)
{
}

template<typename T>
ZD::CFlowABC<T>::~CFlowABC()
{
}

template<typename T>
inline const bool ZD::CFlowABC<T>::CheckPosition(const CPoint<T, 3>& p) const
{
    return true;
}

template <typename T>
inline ZD::CIntegrator<T, 3> * ZD::CFlowABC<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 3>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 3> ZD::CFlowABC<T>::Velocity(const CPoint<T, 3>& p, const T & t) const
{
    T x = m_A * sin(p[2]) + m_C * cos(p[1]);
    T y = m_B * sin(p[0]) + m_A * cos(p[2]);
    T z = m_C * sin(p[1]) + m_B * cos(p[0]);

    return CPoint<T, 3>(x, y, z);
}

template<typename T>
inline const T ZD::CFlowABC<T>::Scalar(const CPoint<T, 3>& p, const T & time) const
{
    return 0.0;
}

template <typename T>
inline void ZD::CFlowABC<T>::GetBBox(CPoint<T, 3>& min, CPoint<T, 3>& max) const
{
    min[0] = 0.0;
    min[1] = 0.0;
    min[2] = 0.0;

    max[0] = 2.0 * ZD_PI;
    max[1] = 2.0 * ZD_PI;
    max[2] = 2.0 * ZD_PI;
}

#endif
