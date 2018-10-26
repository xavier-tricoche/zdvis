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
#ifndef _ZD_LIB_FLOW_CYLINDER_2D_HPP_
#define _ZD_LIB_FLOW_CYLINDER_2D_HPP_

#include "ZD_Flow.hpp"
#include <utils/Base/ZD_Field.hpp>
#include <utils/configure.hpp>
#include <utils/define.hpp>
#include <utils/Integrator/ZD_Integrator_RK4.hpp>

namespace ZD {
    template <typename T>
    class CFlowCylinder2D : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        CField2<T, 2> *m_pVelocities;
        value_type *m_pTimes;

        value_type m_boundingBox[4];
        value_type m_minTime;
        value_type m_maxTime;

        int m_gridSize[2];

    public:
        CFlowCylinder2D(const value_type minTime, const value_type maxTime);
        ~CFlowCylinder2D();

    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        bool CreateTimes();
        bool ReadVelocities(const value_type minTime, const value_type maxTime);
        inline void TimeToIndex(const value_type time, int *indices) const;

        inline const point_type P2C(const point_type& p) const;
    };
}

template <typename T>
ZD::CFlowCylinder2D<T>::CFlowCylinder2D(const value_type minTime, const value_type maxTime)
{
    m_boundingBox[0] = -0.5;        // min x
    m_boundingBox[1] = 7.5;            // max x
    m_boundingBox[2] = -0.5;        // min y
    m_boundingBox[3] = 0.5;            // max y

    m_minTime = 15.0;
    m_maxTime = 23.0;

    m_pTimes = new T[CYLINDER2D_TOTAL_COUNT];
    CreateTimes();

    m_pVelocities = new CField2<T, 2>[CYLINDER2D_TOTAL_COUNT];
    ReadVelocities(minTime, maxTime);
}

template <typename T>
ZD::CFlowCylinder2D<T>::~CFlowCylinder2D()
{
    SafeDeleteArray(m_pTimes);
    SafeDeleteArray(m_pVelocities);
}


template<typename T>
inline const bool ZD::CFlowCylinder2D<T>::CheckPosition(const point_type& p) const
{
    if (p[0] < m_boundingBox[0] || p[0] > m_boundingBox[1] ||
        p[1] < m_boundingBox[2] || p[1] > m_boundingBox[3])
        return false;
    else
        return true;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowCylinder2D<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowCylinder2D<T>::Velocity(const point_type& p, const value_type & time) const
{
    assert(m_pVelocities != nullptr);

    int indices[2];
    TimeToIndex(time, indices);
    value_type f = (time - m_pTimes[indices[0]]) / (m_pTimes[indices[1]] - m_pTimes[indices[0]]);
    point_type cp = P2C(p);
    point_type v0 = m_pVelocities[indices[0]].GetValue(cp);
    point_type v1 = m_pVelocities[indices[1]].GetValue(cp);

    point_type v = v0 * (1.0 - f) + v1 * f;
    return v;
}

template<typename T>
inline const T ZD::CFlowCylinder2D<T>::Scalar(const point_type& p, const value_type & time) const
{
    return 0.0;
}


template <typename T>
bool ZD::CFlowCylinder2D<T>::CreateTimes()
{
    const value_type deltaTime = 0.008;
    for (int i = 0; i < CYLINDER2D_TOTAL_COUNT; ++i) {
        m_pTimes[i] = m_minTime + (value_type)i * deltaTime;
    }
    return true;
}

template <typename T>
bool ZD::CFlowCylinder2D<T>::ReadVelocities(const value_type minTime, const value_type maxTime)
{
    int indices[2];
    TimeToIndex(minTime, indices);
    int minID = indices[0];
    TimeToIndex(maxTime, indices);
    int maxID = indices[1];

    for (int id = minID; id <= maxID; ++id) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%svec_t=%04d.nrrd", CYLINDER2D_VECTOR_PATH, id);
        if (m_pVelocities[id].OpenNrrdFile(pathname) != 0)
            return false;
    }

    memcpy(m_gridSize, m_pVelocities[minID].GetSize(), sizeof(int)*2);
    return true;
}

template<typename T>
inline void ZD::CFlowCylinder2D<T>::TimeToIndex(const value_type time, int * indices) const
{
    indices[0] = 0;
    indices[1] = 0;

    for (int i = 0; i < CYLINDER2D_TOTAL_COUNT - 1; ++i) {
        if (time >= m_pTimes[i] && time < m_pTimes[i + 1]) {
            indices[0] = i;
            indices[1] = i + 1;
            return;
        }
    }
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowCylinder2D<T>::P2C(const point_type& p) const
{
    point_type cp;
    cp[0] = (p[0] - m_boundingBox[0]) / (m_boundingBox[1] - m_boundingBox[0]) * T(m_gridSize[0] - 1);
    cp[1] = (p[1] - m_boundingBox[2]) / (m_boundingBox[3] - m_boundingBox[2]) * T(m_gridSize[1] - 1);
    return cp;
}

template<typename T>
inline void ZD::CFlowCylinder2D<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];

    return;
}

#endif
