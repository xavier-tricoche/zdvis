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
#ifndef _ZD_LIB_HAMILTONIAN_SYSTEM_HPP_
#define _ZD_LIB_HAMILTONIAN_SYSTEM_HPP_

#include "utils/Base/ZD_Point.hpp"
#include "utils/Base/ZD_DynamicSystem.hpp"
#include "utils/Integrator/ZD_Integrator.hpp"
#include "utils/Tool/ZD_StringTool.hpp"

#include <vector>

namespace ZD {
    enum HAMILTONIAN_SYSTEM_FLOWMAP_METHOD {
        HSFM_NONE,
        HSFM_FIXED_TIME,
        HSFM_FIXED_LENGTH,
        HSFM_MAP_ITER
    };

    HAMILTONIAN_SYSTEM_FLOWMAP_METHOD SetHamiltonianSystemFlowmapMethod(const char *str);

    /* class poincare plane */
    template <typename T, unsigned int N> 
    class CPoincarePlane{
    public:
        /* trajectory and Poincare plane intersection */
        virtual bool Intersect(CIntegrator<T, N> *integrator, CPoint<T, N> &v, T &t) const = 0;
        
        ///* project a point on Poincare plane to phase space */
        //virtual void Project(const CPoint<T, 2>& in, CPoint<T, N>& out) const = 0;
        //
        ///* unproject a point in phase space to Poincare plane */
        //virtual void Unproject(const CPoint<T, N>& in, CPoint<T, 2>& out) const = 0;
        
        T secant_method(const T& v0, const T& t0, const T& v1, const T& t1) const {
            return (v1*t0 - v0*t1) / (v1 - v0);
        }
        
    };

    /* class Hamiltonian system */
    template <typename T, unsigned int N>
    class CHamiltonianSystem : public CDynamicSystem<T, N> {
    public:
        CHamiltonianSystem();
        ~CHamiltonianSystem();

        /* integrate for next step */
        //virtual T Next(CPoint<T, N>& p, const T startTime, 
        //    const T direction, CIntegrator<T, N> *integrator) const = 0;

        /* integrate for a fixed time */
        virtual T NextTime(CPoint<T, N>& p, const T startTime, const T time, 
            const T direction) const = 0;
        virtual T NextTime(CPoint<T, N>& p, const T startTime, const T time, 
            const T direction, std::vector<CPoint<T, N>>& buffer) const = 0;

        /* integrate for a fixed trajectory length */
        virtual T NextLength(CPoint<T, N>& p, const T startTime, const T length, 
            const T direction) const = 0;
        virtual T NextLength(CPoint<T, N>& p, const T startTime, const T length, 
            const T direction, std::vector<CPoint<T, N>>& buffer) const = 0;

        /* integrate for one map iteration */
        virtual T NextMap(CPoint<T, N>& p, const T startTime, const T direction, 
            const CPoincarePlane<T, N> *plane) const = 0;
        virtual T NextMap(CPoint<T, N>& p, const T startTime, const T direction,
            const CPoincarePlane<T, N> *plane, std::vector<CPoint<T, N>>& buffer) const = 0;

        virtual CIntegrator<T, N> *CreateIntegrator() const = 0;

    public:
        //virtual const CPoint<T, N> GetVelocity(const CPoint<T, N>& p, const T& time) const = 0;
    };
}

ZD::HAMILTONIAN_SYSTEM_FLOWMAP_METHOD ZD::SetHamiltonianSystemFlowmapMethod(const char * str)
{
    std::string fmStr = CStringTool::ToLower(str);
    if (fmStr == "fixed_time") {
        return HAMILTONIAN_SYSTEM_FLOWMAP_METHOD::HSFM_FIXED_TIME;
    } else if (fmStr == "fixed_length") {
        return HAMILTONIAN_SYSTEM_FLOWMAP_METHOD::HSFM_FIXED_LENGTH;
    } else if (fmStr == "map_iter") {
        return HAMILTONIAN_SYSTEM_FLOWMAP_METHOD::HSFM_MAP_ITER;
    }
    else {
        return HAMILTONIAN_SYSTEM_FLOWMAP_METHOD::HSFM_FIXED_TIME;
    }
}

template <typename T, unsigned int N>
ZD::CHamiltonianSystem<T, N>::CHamiltonianSystem()
{
}

template <typename T, unsigned int N>
ZD::CHamiltonianSystem<T, N>::~CHamiltonianSystem()
{
}

#endif
