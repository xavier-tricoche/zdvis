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
#ifndef _ZD_LIB_FLOW_CYLINDER_STEADY_HPP_
#define _ZD_LIB_FLOW_CYLINDER_STEADY_HPP_

#include "ZD_Flow.hpp"
#include "ZD_Flow_IO_Gerris.hpp"
#include <utils/Integrator/ZD_Integrator_RK4.hpp>
#include <utils/configure.hpp>
#include <utils/define.hpp>

#include <assert.h>
#include <vector>
#include <string>

namespace ZD {
    template <typename T>
    class CFlowCylinderSteady : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        CFlowIOGerris<T> *m_pDatasets;

        std::vector<std::string> m_timesStr;
        value_type *m_pTimes;

        value_type m_boundingBox[4];
        value_type m_minTime;
        value_type m_maxTime;
        value_type m_time;

    public:
        CFlowCylinderSteady(const value_type time, const char *scalarName = nullptr);
        ~CFlowCylinderSteady();

    protected:
        virtual inline CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        bool ReadTimeFile();
        bool ReadGerrisOutput(const value_type time, const char *scalarName);

        inline void TimeToIndex(const T time, int *indices) const;
    };
}

template <typename T>
ZD::CFlowCylinderSteady<T>::CFlowCylinderSteady(const value_type time, const char *scalarName)
{
    this->m_pDatasets = nullptr;
    this->m_pTimes = nullptr;

    ReadTimeFile();

    ReadGerrisOutput(time, scalarName);

    m_boundingBox[0] = -0.5;    // minx
    m_boundingBox[1] = 5.5;        // maxx
    m_boundingBox[2] = -0.5;    // miny
    m_boundingBox[3] = 0.5;        // maxy

    m_time = time;
}

template <typename T>
ZD::CFlowCylinderSteady<T>::~CFlowCylinderSteady()
{
    SafeDeleteArray(this->m_pDatasets);
    SafeDeleteArray(this->m_pTimes);
}

template <typename T>
inline const bool ZD::CFlowCylinderSteady<T>::CheckPosition(const point_type& p) const
{
    if (p[0] > m_boundingBox[0] && p[0] < m_boundingBox[1] &&
        p[1] > m_boundingBox[2] && p[1] < m_boundingBox[3])
        return true;
    else
        return false;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowCylinderSteady<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template <typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowCylinderSteady<T>::Velocity(const point_type& p, const value_type& time) const
{
    value_type f = (m_time - m_minTime) / (m_maxTime - m_minTime);
    point_type v0 = this->m_pDatasets[0].GetVelocity(p);
    point_type v1 = this->m_pDatasets[1].GetVelocity(p);
    point_type v = v0 * (1.0 - f) + v1 * f;
    return v;
}

template <typename T>
inline const T ZD::CFlowCylinderSteady<T>::Scalar(const point_type& p, const value_type& time) const
{
    value_type f = (m_time - m_minTime) / (m_maxTime - m_minTime);
    value_type s0 = this->m_pDatasets[0].GetScalar(p);
    value_type s1 = this->m_pDatasets[1].GetScalar(p);
    value_type s = s0 * (1.0 - f) + s1 * f;
    return s;
}

template <typename T>
bool ZD::CFlowCylinderSteady<T>::ReadTimeFile()
{
    assert(this->m_pTimes == nullptr);

    FILE *fp = fopen(CYLINDER_TIME_PATHNAME, "r");
    if (fp == nullptr) {
        return false;
    }
    else {
        this->m_pTimes = new T[CYLINDER_TOTAL_COUNT];
        for (int i = 0; i < CYLINDER_TOTAL_COUNT; ++i) {
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

template <typename T>
bool ZD::CFlowCylinderSteady<T>::ReadGerrisOutput(const value_type time, const char *scalarName)
{
    assert(this->m_pDatasets == nullptr);

    this->m_pDatasets = new CFlowIOGerris<T>[2];

    int tempID[2];
    TimeToIndex(time, tempID);
    m_minTime = m_pTimes[tempID[0]];
    m_maxTime = m_pTimes[tempID[1]];

    for (int i = 0; i < 2; ++i) {
        std::string timeStr = this->m_timesStr[tempID[i]];
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%ssim-%s.vtk", CYLINDER_DATA_PATH, timeStr.c_str());
        this->m_pDatasets[i].ReadVTK(pathname);

        int gridSize = 32;
        this->m_pDatasets[i].BuildGrid(gridSize);

        if (scalarName != nullptr) {
            std::string tmp = scalarName;
            m_pDatasets[i].SetScalarByName(tmp);
        }
    }

    return true;
}

template <typename T>
inline void ZD::CFlowCylinderSteady<T>::TimeToIndex(const value_type time, int *indices) const
{
    indices[0] = 0;
    indices[1] = CYLINDER_TOTAL_COUNT - 1;
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
inline void ZD::CFlowCylinderSteady<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];

    return;
}

#endif
