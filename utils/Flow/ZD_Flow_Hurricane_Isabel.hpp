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
#ifndef _ZD_LIB_FLOW_HURRICANE_ISBEL_HPP_
#define _ZD_LIB_FLOW_HURRICANE_ISBEL_HPP_

#include "ZD_Flow.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "utils/define.hpp"
#include "utils/configure.hpp"
#include "utils/Tool/ZD_StringTool.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"
#include "utils/Integrator/ZD_Integrator_RK45.hpp"

#include <assert.h>

namespace ZD {
    template <typename T>
    class CFlowHurricaneIsabel : public CFlow<T, 3> {
    public:
        typedef T value_type;
        typedef CPoint<T, 3> point_type;
    private:
        CField3<T, 3> *m_pVectorFields;
        CField3<T, 1> *m_pScalarFields;
        
        value_type *m_pTimes;
        
        CField2<T, 1> *m_pHeightMap;

        value_type m_boundingBox[6];
        value_type m_minTime, m_maxTime;

        point_type m_dispersion;

    public:
        CFlowHurricaneIsabel(const value_type minTime, const value_type maxTime,
            const bool readVelocities, const char *scalarName = nullptr);
        ~CFlowHurricaneIsabel();

    protected:
        virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        inline void TimeToIndex(const value_type& time, int *indices) const;
        bool ReadVelocity(const value_type minTime, const value_type maxTime);
        bool ReadScalar(const value_type minTime, const value_type maxTime, const char *name);
    };
}

template <typename T>
ZD::CFlowHurricaneIsabel<T>::CFlowHurricaneIsabel(const value_type minTime, const value_type maxTime,
    const bool readVelocities, const char *scalarName)
{
    m_boundingBox[0] = 41.7;        // Latitude
    m_boundingBox[1] = 23.7;
    m_boundingBox[2] = 83.0;        // Longitude
    m_boundingBox[3] = 62.0;
    m_boundingBox[4] = 0.035;        // Vertical, km
    m_boundingBox[5] = 19.835;

    m_pTimes = new T[HURRICANE_ISABEL_TOTAL_COUNT];
    for (int i = 0; i < HURRICANE_ISABEL_TOTAL_COUNT; ++i) {
        m_pTimes[i] = i;        // hour
    }

    m_pHeightMap = new CField2<T, 1>;
    m_pHeightMap->OpenNrrdFile(HURRICANE_ISABEL_HEIGHT_MAP_PATHNAME);

    if (readVelocities) {
        m_pVectorFields = new CField3<T, 3>[HURRICANE_ISABEL_TOTAL_COUNT];
        ReadVelocity(minTime, maxTime);
    }
    else {
        m_pVectorFields = nullptr;
    }

    if (scalarName != nullptr) {
        m_pScalarFields = new CField3<T, 1>[HURRICANE_ISABEL_TOTAL_COUNT];
        ReadScalar(minTime, maxTime, scalarName);
    }
    else {
        m_pScalarFields = nullptr;
    }

    m_dispersion[0] = (m_boundingBox[1] - m_boundingBox[0]) / (500.0 - 1.0);
    m_dispersion[1] = (m_boundingBox[3] - m_boundingBox[2]) / (500.0 - 1.0);
    m_dispersion[2] = (m_boundingBox[5] - m_boundingBox[4]) / (100.0 - 1.0);
}

template <typename T>
ZD::CFlowHurricaneIsabel<T>::~CFlowHurricaneIsabel()
{
    SafeDeleteArray(m_pVectorFields);
    SafeDeleteArray(m_pScalarFields);
    SafeDeleteArray(m_pTimes);
    SafeDelete(m_pHeightMap);
}

template <typename T>
inline const bool ZD::CFlowHurricaneIsabel<T>::CheckPosition(const point_type& p) const
{
    if (p[0] > this->m_boundingBox[0] || p[0] < this->m_boundingBox[1] ||
        p[1] > this->m_boundingBox[2] || p[1] < this->m_boundingBox[3] ||
        p[2] < this->m_boundingBox[4] || p[2] > this->m_boundingBox[5]) {
        return false;
    }
    else {
        /* check height map */
        value_type fx = (p[0] - this->m_boundingBox[0]) / m_dispersion[0];
        value_type fy = (p[1] - this->m_boundingBox[2]) / m_dispersion[1];
        value_type height = m_pHeightMap->GetValue(fx, fy)[0] * 0.001; /* m to km */
        if (p[2] > height)
            return true;
        else
            return false;
    }
}

template <typename T>
inline ZD::CIntegrator<T, 3> * ZD::CFlowHurricaneIsabel<T>::CreateIntegrator() const
{
    //return new CIntegratorRK45<T, 3>(1e-7, this->m_stepSize);
    return new CIntegratorRK4<T, 3>(this->m_stepSize);
}

template <typename T>
inline const ZD::CPoint<T, 3> ZD::CFlowHurricaneIsabel<T>::Velocity(const point_type& p, const T& time) const
{
    int id[2];
    TimeToIndex(time, id);

    value_type f = (time - m_pTimes[id[0]]) / (m_pTimes[id[1]] - m_pTimes[id[0]]);

    value_type fx = (p[0] - this->m_boundingBox[0]) / m_dispersion[0];
    value_type fy = (p[1] - this->m_boundingBox[2]) / m_dispersion[1];
    value_type fz = (p[2] - this->m_boundingBox[4]) / m_dispersion[2];

    point_type v0 = m_pVectorFields[id[0]].GetValue(fx, fy, fz);
    point_type v1 = m_pVectorFields[id[1]].GetValue(fx, fy, fz);
    point_type dir = v0 * (1.0 - f) + v1 * f;
    dir = dir * 3.6;

    value_type r = ZD_EARTH_R + p[2];

    value_type lonr = r * std::cos(p[0] / 180.0 * ZD_PI);
    value_type latr = r;
    value_type dlon = lonr / 360.0;
    value_type dlat = latr / 360.0;

    dir[0] = -dir[0] / dlat;
    dir[1] = -dir[1] / dlon;

    return dir;
}

template <typename T>
inline const T ZD::CFlowHurricaneIsabel<T>::Scalar(const point_type& p, const value_type& time) const
{
    assert(m_pScalarFields != nullptr);

    int id[2];
    TimeToIndex(time, id);

    value_type f = (time - m_pTimes[id[0]]) / (m_pTimes[id[1]] - m_pTimes[id[0]]);

    value_type fx = (p[0] - this->m_boundingBox[0]) / m_dispersion[0];
    value_type fy = (p[1] - this->m_boundingBox[2]) / m_dispersion[1];
    value_type fz = (p[2] - this->m_boundingBox[4]) / m_dispersion[2];

    value_type s0 = m_pScalarFields[id[0]].GetValue(fx, fy, fz)[0];
    value_type s1 = m_pScalarFields[id[1]].GetValue(fx, fy, fz)[0];

    value_type sp = s0 * (1.0 - f) + s1 * f;

    return sp;
}

template <typename T>
bool ZD::CFlowHurricaneIsabel<T>::ReadVelocity(const value_type minTime, const value_type maxTime)
{
    int temp[2];
    TimeToIndex(minTime, temp);
    int minID = temp[0];
    TimeToIndex(maxTime, temp);
    int maxID = temp[1];

    for (int i = minID; i <= maxID; ++i) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%svelocity_%02d.nrrd", HURRICANE_ISABEL_VECTOR_PATH, i + 1);
        m_pVectorFields[i].OpenNrrdFile(pathname);
    }

    return true;
}

template <typename T>
bool ZD::CFlowHurricaneIsabel<T>::ReadScalar(const value_type minTime, const value_type maxTime, const char *name)
{
    std::string tmp = CStringTool::ToUpper(name);
    std::string path;
    std::string scalarName;
    if (strcmp(tmp.c_str(), "CLOUD") == 0) {
        path = HURRICANE_ISABEL_CLOUD_PATH;
        scalarName = "CLOUD";
    }
    else if (strcmp(tmp.c_str(), "P") == 0) {
        path = HURRICANE_ISABEL_PRESSURE_PATH;
        scalarName = "P";
    }
    else if (strcmp(tmp.c_str(), "TC") == 0) {
        path = HURRICANE_ISABEL_TEMPERATURE_PATH;
        scalarName = "TC";
    }
    else if (strcmp(tmp.c_str(), "PRECIP") == 0) {
        path = HURRICANE_ISABEL_PRECIP_PATH;
        scalarName = "PRECIP";
    }
    else if (strcmp(tmp.c_str(), "SNOW") == 0) {
        path = HURRICANE_ISABEL_SNOW_PATH;
        scalarName = "QSNOW";
    }
    else {
        return false;
    }

    int temp[2];
    TimeToIndex(minTime, temp);
    int minID = temp[0];
    TimeToIndex(maxTime, temp);
    int maxID = temp[1];

    for (int i = minID; i <= maxID; ++i) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%s%s_%02d.nrrd", path.c_str(), scalarName.c_str(), i + 1);
        m_pScalarFields[i].OpenNrrdFile(pathname);
    }

    return true;
}

template <typename T>
inline void ZD::CFlowHurricaneIsabel<T>::TimeToIndex(const value_type& time, int *indices) const
{
    indices[0] = 0;
    indices[1] = 0;
    for (int i = 0; i < HURRICANE_ISABEL_TOTAL_COUNT - 1; ++i) {
        if (time >= m_pTimes[i] && time < m_pTimes[i + 1]) {
            indices[0] = i;
            indices[1] = i + 1;
            return;
        }
    }
}

template <typename T>
inline void ZD::CFlowHurricaneIsabel<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    min[2] = m_boundingBox[4];

    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];
    max[2] = m_boundingBox[5];
}

#endif
