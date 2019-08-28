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
#ifndef _ZD_LIB_FLOW_TEST_HPP_
#define _ZD_LIB_FLOW_TEST_HPP_

#include "ZD_Flow.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"

namespace ZD {
    template <typename T>
    class CFlowTest2D : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& t) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
    };

    template <typename T>
    class CFlowTest3D : public CFlow<T, 3> {
    public:
        typedef T value_type;
        typedef CPoint<T, 3> point_type;
    protected:
        virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& t) const;
        virtual inline const value_type Scalar(const point_type& p, const T& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
    };
}

template <typename T>
inline const bool ZD::CFlowTest2D<T>::CheckPosition(const point_type& p) const
{
    return true;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowTest2D<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template <typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowTest2D<T>::Velocity(const point_type& p, const T& time) const
{
    point_type v;
    v[0] = 1.0;
    v[1] = 2.0;
    return v;
}

template <typename T>
inline const T ZD::CFlowTest2D<T>::Scalar(const point_type& p, const T& time) const
{
    return 1.0;
}

template <typename T>
inline void ZD::CFlowTest2D<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = 0.0;
    min[1] = 0.0;

    max[0] = 1.0;
    max[1] = 1.0;

    return;
}

template <typename T>
inline const bool ZD::CFlowTest3D<T>::CheckPosition(const point_type& p) const
{
    return true;
}

template <typename T>
inline ZD::CIntegrator<T, 3> * ZD::CFlowTest3D<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 3>(this->m_stepSize);
}

template <typename T>
inline const ZD::CPoint<T, 3> ZD::CFlowTest3D<T>::Velocity(const point_type& p, const T& time) const
{
    point_type v;
    v[0] = 1.0 * time;
    v[1] = 2.0 * time;
    v[2] = 3.0 * time;
    return v;
}

template <typename T>
inline const T ZD::CFlowTest3D<T>::Scalar(const point_type& p, const T& time) const
{
    return 1.0;
}

template <typename T>
inline void ZD::CFlowTest3D<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = 0.0;
    min[1] = 0.0;
    min[2] = 0.0;

    max[0] = 1.0;
    max[1] = 1.0;
    max[2] = 1.0;

    return;
}

#endif
