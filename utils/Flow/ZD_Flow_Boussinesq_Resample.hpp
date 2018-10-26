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
#ifndef _ZD_LIB_FLOW_BOUSSINESQ_RESAMPLE_HPP_
#define _ZD_LIB_FLOW_BOUSSINESQ_RESAMPLE_HPP_

#include <vector>
#include <string>
#include <assert.h>

#include "ZD_Flow.hpp"
#include <utils/Base/ZD_Field.hpp>
#include <utils/Integrator/ZD_Integrator_RK4.hpp>
#include <utils/configure.hpp>
#include <utils/define.hpp>

namespace ZD {
    template <typename T>
    class CFlowBoussinesqResample : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
        typedef CIntegrator<T, 2> integ_type;
    private:
        CField2<T, 2> *m_pVelocities;

        std::vector<std::string> m_timesStr;
        value_type *m_pTimes;

        value_type m_boundingBox[4];
        value_type m_minTime;
        value_type m_maxTime;

        int m_gridSize[2];

    public:
        CFlowBoussinesqResample(const value_type minTime, const value_type maxTime);
        ~CFlowBoussinesqResample();

    protected:
        virtual inline integ_type *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;
        virtual inline const bool CheckPosition(const point_type& p) const;

    private:
        bool ReadTimeFile();
        bool ReadVelocities(const value_type minTime, const value_type maxTime);
        inline void TimeToIndex(const value_type time, int *indices) const;
        inline const point_type P2C(const point_type& p) const;
    };
} // namespace ZD

template <typename T>
ZD::CFlowBoussinesqResample<T>::CFlowBoussinesqResample(const value_type minTime, const value_type maxTime)
{
    this->m_pVelocities = nullptr;
    this->m_pTimes = nullptr;

    ReadTimeFile();

    m_minTime = this->m_pTimes[0];
    m_maxTime = this->m_pTimes[BOUSSINESQ_TOTAL_COUNT - 1];

    m_boundingBox[0] = -0.5;    // minx
    m_boundingBox[1] = 0.5;    // maxx
    m_boundingBox[2] = -0.5;    // miny
    m_boundingBox[3] = 2.5;    // maxy

    m_gridSize[0] = m_gridSize[1] = 0;

    ReadVelocities(minTime, maxTime);
}

template<typename T>
ZD::CFlowBoussinesqResample<T>::~CFlowBoussinesqResample()
{
    SafeDeleteArray(m_pVelocities);
    SafeDeleteArray(m_pTimes);
}

template <typename T>
inline const bool ZD::CFlowBoussinesqResample<T>::CheckPosition(const point_type& p) const
{
    if (p[0] > m_boundingBox[0] && p[0] < m_boundingBox[1] &&
        p[1] > m_boundingBox[2] && p[1] < m_boundingBox[3])
        return true;
    else
        return false;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowBoussinesqResample<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}


template <typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowBoussinesqResample<T>::Velocity(const point_type& p, const T& time) const
{
    int indices[2];
    TimeToIndex(time, indices);
    T f = (time - m_pTimes[indices[0]]) / (m_pTimes[indices[1]] - m_pTimes[indices[0]]);
    point_type cp = P2C(p);
    point_type v0 = this->m_pVelocities[indices[0]].GetValue(cp);
    point_type v1 = this->m_pVelocities[indices[1]].GetValue(cp);

    CPoint<T, 2> v = v0 * (1.0 - f) + v1 * f;
    return v;
}

template<typename T>
inline const T ZD::CFlowBoussinesqResample<T>::Scalar(const point_type& p, const value_type& time) const
{
    return 0.0;
}

template <typename T>
bool ZD::CFlowBoussinesqResample<T>::ReadTimeFile()
{
    assert(this->m_pTimes == nullptr);

    FILE *fp = fopen(BOUSSINESQ_TIME_PATHNAME, "r");
    if (fp == nullptr) {
        return false;
    }
    else {
        this->m_pTimes = new T[BOUSSINESQ_TOTAL_COUNT];
        for (int i = 0; i < BOUSSINESQ_TOTAL_COUNT; ++i) {
            char temp[256];
            fscanf(fp, "%s\n", temp);
            this->m_timesStr.push_back(temp);
            double t;
            sscanf(temp, "%lf", &t);
            this->m_pTimes[i] = t;
        }
        fclose(fp);
        return true;
    }
}

template<typename T>
bool ZD::CFlowBoussinesqResample<T>::ReadVelocities(const value_type minTime, const value_type maxTime)
{
    m_pVelocities = new CField2<T, 2>[BOUSSINESQ_TOTAL_COUNT];
    int indices[2];
    TimeToIndex(minTime-ZD_EPSILON, indices);
    int minID = indices[0];
    TimeToIndex(maxTime+ZD_EPSILON, indices);
    int maxID = indices[1];

    for (int id = minID; id <= maxID; ++id) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%svec_%04d.nrrd", BOUSSINESQ_VECTOR_PATH, id);
        if (m_pVelocities[id].OpenNrrdFile(pathname) != 0)
            return false;
    }

    memcpy(m_gridSize, m_pVelocities[minID].GetSize(), sizeof(int) * 2);
    return true;
}

template <typename T>
inline void ZD::CFlowBoussinesqResample<T>::TimeToIndex(const value_type time, int *indices) const
{
    indices[0] = 0;
    indices[1] = BOUSSINESQ_TOTAL_COUNT - 1;
    while ((indices[0] + 1) < indices[1]) {
        int mid = (indices[0] + indices[1]) / 2.0;
        if (this->m_pTimes[mid] > time) {
            indices[1] = mid;
        }
        else {
            indices[0] = mid;
        }
    }
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowBoussinesqResample<T>::P2C(const point_type& p) const
{
    point_type cp;
    cp[0] = (p[0] - m_boundingBox[0]) / (m_boundingBox[1] - m_boundingBox[0]) * T(m_gridSize[0] - 1);
    cp[1] = (p[1] - m_boundingBox[2]) / (m_boundingBox[3] - m_boundingBox[2]) * T(m_gridSize[1] - 1);
    return cp;
}

template<typename T>
inline void ZD::CFlowBoussinesqResample<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];

    return;
}

#endif
