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
#ifndef _ZD_LIB_FLOW_CONVECTION_HPP_
#define _ZD_LIB_FLOW_CONVECTION_HPP_

#include "ZD_Flow.hpp"
#include <utils/Base/ZD_Field.hpp>
#include <utils/configure.hpp>
#include <utils/define.hpp>
#include <utils/Integrator/ZD_Integrator_RK4.hpp>

#include <iostream>
#include <assert.h>

namespace ZD {
    template <typename T>
    class CFlowConvection : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2>      point_type;
        typedef CIntegrator<T, 2> integ_type;
    private:
        CField2<T, 2> *m_pVelocities;
        CField2<T, 1> *m_pScalars;
        
        value_type *m_pTimes;
        
        value_type m_boundingBox[4];
        value_type m_minTime;
        value_type m_maxTime;

        int m_gridSize[2];

    public:
        CFlowConvection(const value_type minTime, const value_type maxTime, 
            const bool readVelocities, const char *scalarName = nullptr);

        ~CFlowConvection();

    protected:
        virtual inline integ_type *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        bool CreateTimes();
        bool ReadVelocities(const value_type minTime, const value_type maxTime);
        bool ReadScalars(const value_type minTime, const value_type maxTime, const char *scalarName);

        inline void TimeToIndex(const value_type time, int *indices) const;

        inline const point_type P2C(const point_type& p) const;
    };
} // namespace ZD

template<typename T>
ZD::CFlowConvection<T>::CFlowConvection(const value_type minTime, const value_type maxTime, 
    const bool readVelocitiesFlag, const char * scalarName)
{
    m_boundingBox[0] = 0.0;            // min x
    m_boundingBox[1] = 0.5;            // max x
    m_boundingBox[2] = 0.0;            // min y
    m_boundingBox[3] = 0.25;        // max y

    m_minTime = 0.0;
    m_maxTime = 5.0;

    m_pVelocities = nullptr;
    m_pScalars = nullptr;
    m_pTimes = nullptr;

    m_gridSize[0] = m_gridSize[1] = 0;

    if (CreateTimes() == false) {
        std::cout << "Error in creating time array." << std::endl;
    }

    if (readVelocitiesFlag && ReadVelocities(minTime, maxTime) == false) {
        std::cout << "Error in reading velocity files." << std::endl;
    }

    if (scalarName != nullptr && ReadScalars(minTime, maxTime, scalarName) == false) {
        std::cout << "Error in reading scalar "<< scalarName << "." << std::endl;
    }
}

template<typename T>
ZD::CFlowConvection<T>::~CFlowConvection()
{
    SafeDeleteArray(m_pVelocities);
    SafeDeleteArray(m_pScalars);
    SafeDeleteArray(m_pTimes);
}

template<typename T>
inline const bool ZD::CFlowConvection<T>::CheckPosition(const point_type& p) const
{
    if (p[0] < m_boundingBox[0] || p[0] > m_boundingBox[1] ||
        p[1] < m_boundingBox[2] || p[1] > m_boundingBox[3])
        return false;
    else
        return true;
}

template <typename T>
inline ZD::CIntegrator<T, 2> * ZD::CFlowConvection<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowConvection<T>::Velocity(const point_type& p, const value_type & time) const
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
inline const T ZD::CFlowConvection<T>::Scalar(const point_type& p, const value_type & time) const
{
    assert(m_pScalars != nullptr);
    
    int indices[2];
    TimeToIndex(time, indices);
    value_type f = (time - m_pTimes[indices[0]]) / (m_pTimes[indices[1]] - m_pTimes[indices[0]]);
    point_type cp = P2C(p);
    CPoint<T, 1> s0 = m_pScalars[indices[0]].GetValue(cp);
    CPoint<T, 1> s1 = m_pScalars[indices[1]].GetValue(cp);
    
    CPoint<T, 1> s = s0 * (1.0 - f) + s1 * f;
    return s[0];
}

template<typename T>
bool ZD::CFlowConvection<T>::CreateTimes()
{
    m_pTimes = new value_type[CONVECTION_TOTAL_COUNT];
    for (int i = 0; i < CONVECTION_TOTAL_COUNT; ++i) {
        m_pTimes[i] = m_minTime + (m_maxTime - m_minTime) * T(i) / T(CONVECTION_TOTAL_COUNT - 1);
    }
    return true;
}

template<typename T>
bool ZD::CFlowConvection<T>::ReadVelocities(const value_type minTime, const value_type maxTime)
{
    m_pVelocities = new CField2<T, 2>[CONVECTION_TOTAL_COUNT];
    int indices[2];
    TimeToIndex(minTime, indices);
    int minID = indices[0];
    TimeToIndex(maxTime, indices);
    int maxID = indices[1];

    for (int id = minID; id <= maxID; ++id) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%svec_%04d.nhdr", CONVECTION_VECTOR_PATH, id);
        if (m_pVelocities[id].OpenNrrdFile(pathname) != 0)
            return false;
    }

    memcpy(m_gridSize, m_pVelocities[minID].GetSize(), sizeof(int) * 2);
    return true;
}

template<typename T>
bool ZD::CFlowConvection<T>::ReadScalars(const value_type minTime, const value_type maxTime, const char * scalarName)
{
    std::string path;
    std::string name;
    if (strcmp(scalarName, "temperature") == 0) {
        path = CONVECTION_TEMPERATURE_PATH;
        name = "t_";
    }
    else {
        return false;
    }

    m_pScalars = new CField2<T, 1>[CONVECTION_TOTAL_COUNT];
    
    int indices[2];
    TimeToIndex(minTime, indices);
    int minID = indices[0];
    TimeToIndex(maxTime, indices);
    int maxID = indices[1];

    for (int id = minID; id <= maxID; ++id) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%s%s%04d.nhdr", path.c_str(), name.c_str(), id);
        if (m_pScalars[id].OpenNrrdFile(pathname) != 0)
            return false;
    }

    memcpy(m_gridSize, m_pScalars[minID].GetSize(), sizeof(int) * 2);
    return true;
}

template<typename T>
inline void ZD::CFlowConvection<T>::TimeToIndex(const value_type time, int * indices) const
{
    indices[0] = 0;
    indices[1] = 0;
    
    for (int i = 0; i < CONVECTION_TOTAL_COUNT-1; ++i) {
        if (time >= m_pTimes[i] && time < m_pTimes[i + 1]) {
            indices[0] = i;
            indices[1] = i + 1;
            return;
        }
    }
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowConvection<T>::P2C(const point_type& p) const
{
    point_type cp;
    cp[0] = (p[0] - m_boundingBox[0]) / (m_boundingBox[1] - m_boundingBox[0]) * T(m_gridSize[0] - 1);
    cp[1] = (p[1] - m_boundingBox[2]) / (m_boundingBox[3] - m_boundingBox[2]) * T(m_gridSize[1] - 1);
    return cp;
}

template<typename T>
inline void ZD::CFlowConvection<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];
    
    return;
}

#endif
