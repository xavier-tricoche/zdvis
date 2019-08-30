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
#ifndef _ZD_LIB_FLOW_HPP_
#define _ZD_LIB_FLOW_HPP_

#include <utils/Base/ZD_Point.hpp>
#include <utils/Base/ZD_Line.hpp>
#include <utils/Base/ZD_DynamicSystem.hpp>
#include <utils/Integrator/ZD_Integrator.hpp>
#include <utils/Integrator/ZD_Integrator_Euler.hpp>
#include <utils/Integrator/ZD_Integrator_RK4.hpp>
#include <utils/Integrator/ZD_Integrator_RK45.hpp>
#include <utils/define.hpp>
#include <exception>

#include <vector>

namespace ZD {
    template <typename T, unsigned int N>
    class CFlow : public CDynamicSystem<T, N> {
    public:
        typedef T value_type;
        typedef CPoint<T, N> point_type;
    protected:
        value_type m_stepSize;

    public:
        CFlow();
        virtual ~CFlow() {}

    public:
        const value_type NextTime(point_type& p, 
                                  const value_type startTime, 
                                  const value_type time, 
                                  const value_type direction) const;
        const value_type NextTime(point_type& p, 
                                  const value_type startTime, 
                                  const value_type time,
                                  const value_type direction, 
                                  std::vector<point_type>& buffer) const;
        const value_type NextTime(point_type& p, 
                                  const value_type startTime, 
                                  const value_type time, 
                                  const value_type direction,
                                  CIntegrator<T, N>* integrator) const;
        const value_type NextTime(point_type& p, 
                                  const value_type startTime, 
                                  const value_type time,
                                  const value_type direction, 
                                  std::vector<point_type>& buffer, 
                                  CIntegrator<T, N>* integrator) const;

        virtual const bool CheckPosition(const point_type& p) const = 0;

    protected:
        virtual CIntegrator<T, N> *CreateIntegrator() const = 0;
    public:
        virtual const value_type Scalar(const point_type& p, 
                                        const value_type& time) const = 0;
        virtual void GetBBox(point_type& min, point_type& max) const  = 0;

    public:
        void SetStepSize(const value_type stepSize);
        const value_type GetStepSize() const;


    };
} // namespace ZD

template<typename T, unsigned int N>
ZD::CFlow<T, N>::CFlow()
{
}

template<typename T, unsigned int N>
inline const T ZD::CFlow<T, N>::
NextTime(point_type& p, const value_type startTime, const value_type time,
         const value_type direction, CIntegrator<T, N>* integrator) const
{
    value_type integratedTime = 0.0;
    while (integratedTime < time) {
        if (CheckPosition(p) == false) {
            break;
        }
        try {
            value_type temp = integrator->Next(p, startTime + integratedTime*direction, direction, this);
            integratedTime += temp;
        }
        catch(std::runtime_error& e) {
            // something went wrong during ODE step, most probably a failed interpolation
            std::cout << "exception caught in NextTime: " << e.what() << '\n';
            break;
        }
        catch(...) {
            std::cout << "unknown exception caught\n";
            break;
        }
        if (integratedTime > time) {
            p = integrator->v(time*direction+startTime);
            integratedTime = time;
            break;
        }
        else {
            p = integrator->v1();
        }
    }

    SafeDelete(integrator);
    return integratedTime;
}

template<typename T, unsigned int N>
inline const T ZD::CFlow<T, N>::NextTime(point_type& p, 
                                         const value_type startTime,
                                         const value_type time, 
                                         const value_type direction) const
{
    CIntegrator<T, N>* intg = CreateIntegrator();
    return NextTime(p, startTime, time, direction, intg);
}

template<typename T, unsigned int N>
inline const T ZD::CFlow<T, N>::NextTime(point_type& p, 
                                         const value_type startTime,
                                         const value_type time,
                                         const value_type direction, 
                                         std::vector<point_type>& buffer,
                                         CIntegrator<T, N>* integrator) const
{
    value_type integratedTime = 0.0;
    while (integratedTime < time) {
        if (CheckPosition(p) == false)
            break;
        value_type temp = integrator->Next(p, startTime + integratedTime*direction, direction, this);
        integratedTime += temp;
        if (integratedTime > time) {
            p = integrator->v(time*direction+startTime);
            integratedTime = time;
            buffer.push_back(p);
            break;
        }
        else {
            p = integrator->v1();
            buffer.push_back(p);
        }
    }

    SafeDelete(integrator);
    return integratedTime;
}

template<typename T, unsigned int N>
inline const T ZD::CFlow<T, N>::NextTime(point_type& p, 
                                         const value_type startTime, 
                                         const value_type time, 
                                         const value_type direction, 
                                         std::vector<point_type>& buffer) const
{
    CIntegrator<T, N> *intg = CreateIntegrator();
    return NextTime(p, startTime, time, direction, buffer, intg);
}

template<typename T, unsigned int N>
void ZD::CFlow<T, N>::SetStepSize(const value_type stepSize)
{
    this->m_stepSize = stepSize;
}

template<typename T, unsigned int N>
const T ZD::CFlow<T, N>::GetStepSize() const
{
    return this->m_stepSize;
}

#endif
