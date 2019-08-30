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
#ifndef _ZD_LIB_FLOW_GAUSSIAN_VORTICES_HPP_
#define _ZD_LIB_FLOW_GAUSSIAN_VORTICES_HPP_

#include "ZD_Flow.hpp"
#include "ZD_Flow_IO_Gerris.hpp"

#include "utils/Integrator/ZD_Integrator_RK4.hpp"
#include "utils/configure.hpp"
#include "utils/define.hpp"

#include <assert.h>
#include <vector>
#include <string>

namespace ZD {
    template <typename T>
    class CFlowGaussianVortices : public CFlow<T, 2> {
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

    public:
        CFlowGaussianVortices(const value_type minTime, 
                              const value_type maxTime, 
                              const char *scalarName = nullptr);
        ~CFlowGaussianVortices();

    protected:
        virtual CIntegrator<T, 2> *CreateIntegrator() const;

    public:
        virtual const point_type Velocity(const point_type& p,
                                          const value_type& time) const;
        virtual const value_type Scalar(const point_type& p, 
                                        const value_type& time) const;
        virtual const bool CheckPosition(const point_type& p) const;
        virtual void GetBBox(point_type& min, point_type& max) const;

    private:
        bool ReadTimeFile();
        bool ReadGerrisOutput(const value_type minTime, 
                              const value_type maxTime, 
                              const char *scalarName);

        inline void TimeToIndex(const T time, int *indices) const;
    };
}

template <typename T>
ZD::CFlowGaussianVortices<T>::
CFlowGaussianVortices(const value_type minTime, const value_type maxTime, 
                      const char *scalarName)
{
    this->m_pDatasets = nullptr;
    this->m_pTimes = nullptr;

    ReadTimeFile();

    ReadGerrisOutput(minTime, maxTime, scalarName);

    m_minTime = this->m_pTimes[0];
    m_maxTime = this->m_pTimes[GAUSSIAN_VORTICES_TOTAL_COUNT - 1];
    
    m_boundingBox[0] = -0.5;    // minx
    m_boundingBox[1] =  0.5;    // maxx
    m_boundingBox[2] = -0.5;    // miny
    m_boundingBox[3] =  0.5;    // maxy
}

template <typename T>
ZD::CFlowGaussianVortices<T>::~CFlowGaussianVortices()
{
    SafeDeleteArray(this->m_pDatasets);
    SafeDeleteArray(this->m_pTimes);
}

template <typename T>
inline const bool 
ZD::CFlowGaussianVortices<T>::CheckPosition(const point_type& p) const
{
    if (p[0] > m_boundingBox[0] && p[0] < m_boundingBox[1] &&
        p[1] > m_boundingBox[2] && p[1] < m_boundingBox[3])
        return true;
    else
        return false;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::
CFlowGaussianVortices<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template <typename T>
inline const ZD::CPoint<T, 2> 
ZD::CFlowGaussianVortices<T>::Velocity(const point_type& p, 
                                       const value_type& time) const
{
    int indices[2];
    TimeToIndex(time, indices);
    value_type f = (time - m_pTimes[indices[0]]) / 
                   (m_pTimes[indices[1]] - m_pTimes[indices[0]]);

    point_type v0 = this->m_pDatasets[indices[0]].GetVelocity(p);
    point_type v1 = this->m_pDatasets[indices[1]].GetVelocity(p);

    point_type v = v0 * (1.0 - f) + v1 * f;
    return v;
}

template <typename T>
inline const T 
ZD::CFlowGaussianVortices<T>::
Scalar(const point_type& p, const value_type& time) const
{
    int indices[2];
    TimeToIndex(time, indices);
    value_type f = (time - m_pTimes[indices[0]]) / (m_pTimes[indices[1]] - m_pTimes[indices[0]]);

    value_type s0 = this->m_pDatasets[indices[0]].GetScalar(p);
    value_type s1 = this->m_pDatasets[indices[1]].GetScalar(p);

    value_type s = s0 * (1.0 - f) + s1 * f;
    return s;
}

template <typename T>
bool ZD::CFlowGaussianVortices<T>::ReadTimeFile()
{
    assert(this->m_pTimes == nullptr);

    FILE *fp = fopen(GAUSSIAN_VORTICES_TIME_PATHNAME, "r");
    if (fp == nullptr) {
        return false;
    }
    else {
        this->m_pTimes = new T[GAUSSIAN_VORTICES_TOTAL_COUNT];
        for (int i = 0; i < GAUSSIAN_VORTICES_TOTAL_COUNT; ++i) {
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
bool ZD::CFlowGaussianVortices<T>::
ReadGerrisOutput(const value_type minTime, const value_type maxTime,
                 const char *scalarName)
{
    assert(this->m_pDatasets == nullptr);

    this->m_pDatasets = new CFlowIOGerris<T>[GAUSSIAN_VORTICES_TOTAL_COUNT];

    int tempID[2];
    TimeToIndex(minTime, tempID);
    int minID = tempID[0];
    TimeToIndex(maxTime, tempID);
    int maxID = tempID[1];

    for (int i = minID; i <= maxID; ++i) {
        std::string timeStr = this->m_timesStr[i];
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%ssim-%s.vtk", GAUSSIAN_VORTICES_DATA_PATH, timeStr.c_str());
        this->m_pDatasets[i].ReadVTK(pathname);

        int gridSize = 32;
        this->m_pDatasets[i].BuildGrid(gridSize);

        if (scalarName != nullptr) {
            std::string tmp = scalarName;
            m_pDatasets[i].SetScalarByName(tmp);
        }

        //this->m_pDatasets[i].AddNoise();
    }

    return true;
}

template <typename T>
inline void ZD::CFlowGaussianVortices<T>::
TimeToIndex(const value_type time, int *indices) const
{
    indices[0] = 0;
    indices[1] = GAUSSIAN_VORTICES_TOTAL_COUNT - 1;
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
inline void ZD::CFlowGaussianVortices<T>::
GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];

    return;
}

#endif
