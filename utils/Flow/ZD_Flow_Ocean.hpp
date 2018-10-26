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
#ifndef _ZD_LIB_FLOW_OCEAN_HPP_
#define _ZD_LIB_FLOW_OCEAN_HPP_

#include "ZD_Flow.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "utils/configure.hpp"
#include "utils/define.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"

#include <vector>
#include <string>
#include <fstream>

namespace ZD{
    template <typename T>
    class CFlowOcean : public CFlow<T, 2> {
    public:
        typedef T value_type;
        typedef CPoint<T, 2> point_type;
    private:
        CField2<T, 2> *m_pVelocities;
        T *m_pTimes;

        value_type m_boundingBox[4];
        value_type m_minTime;
        value_type m_maxTime;

        value_type m_gridSize[2];

    public:
        CFlowOcean(const value_type minTime, const value_type maxTime);
        ~CFlowOcean();

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
} // namespace ZD

template<typename T>
ZD::CFlowOcean<T>::CFlowOcean(const value_type minTime, const value_type maxTime)
{
    m_boundingBox[0] = 262.04998779296875;            // min x (longitude)
    m_boundingBox[1] = 279.74999999998760;            // max x (longitude)
    m_boundingBox[2] = 18.049999237060547;            // min y (latitude )
    m_boundingBox[3] = 30.789993286135107;            // max y (latitude )

    m_minTime = 0.0;                // from 2013/04/01 to 2013/08/31
    m_maxTime = 3670.0;                // about 3670 hours

    m_pTimes = new T[OCEAN_TOTAL_COUNT];
    CreateTimes();

    m_pVelocities = new CField2<T, 2>[OCEAN_TOTAL_COUNT];
    ReadVelocities(minTime, maxTime);
}

template <typename T>
ZD::CFlowOcean<T>::~CFlowOcean()
{
    SafeDeleteArray(m_pTimes);
    SafeDeleteArray(m_pVelocities);
}


template<typename T>
inline ZD::CIntegrator<T, 2>* ZD::CFlowOcean<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 2>(this->m_stepSize);
}

template<typename T>
inline const bool ZD::CFlowOcean<T>::CheckPosition(const point_type& p) const
{
    if (p[0] < m_boundingBox[0] || p[0] > m_boundingBox[1] ||
        p[1] < m_boundingBox[2] || p[1] > m_boundingBox[3])
        return false;
    else
        return true;
}

template<typename T>
inline void ZD::CFlowOcean<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];

    return;
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowOcean<T>::Velocity(const point_type& p, const value_type & time) const
{
    assert(m_pVelocities != nullptr);

    int indices[2];
    TimeToIndex(time, indices);
    value_type f = (time - m_pTimes[indices[0]]) / (m_pTimes[indices[1]] - m_pTimes[indices[0]]);
    point_type cp = P2C(p);
    point_type v0 = m_pVelocities[indices[0]].GetValue(cp);
    point_type v1 = m_pVelocities[indices[1]].GetValue(cp);

    point_type v = v0 * (1.0 - f) + v1 * f;

    // convert from degree / second to degree / hour
    v = v * 3600.0;
    //v[0] = -v[0];
    return v;
}

template<typename T>
inline const T ZD::CFlowOcean<T>::Scalar(const point_type& p, const value_type & time) const
{
    return 0.0;
}

template<typename T>
bool ZD::CFlowOcean<T>::CreateTimes()
{
    const value_type deltaTime = 3.0;    // 3 hours
    for (int i = 0; i < OCEAN_TOTAL_COUNT; ++i) {
        m_pTimes[i] = m_minTime + (value_type)i * deltaTime;
    }
    return true;
}

template<typename T>
bool ZD::CFlowOcean<T>::ReadVelocities(const value_type minTime, const value_type maxTime)
{
    // read all file names
    std::vector<std::string> pathnames(OCEAN_TOTAL_COUNT);
    std::ifstream ifs;
    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sall_files.txt", OCEAN_VECTOR_PATH);
    ifs.open(pathname);
    if (ifs.is_open()) {
        for (int i = 0; i < OCEAN_TOTAL_COUNT; ++i) {
            std::string temp;
            int tmp;
            ifs >> pathnames[i] >> temp >> tmp;
        }
        ifs.close();
    } else {
        std::cout << "Cannot open all_files.txt!" << std::endl;
        return false;
    }


    int indices[2];
    TimeToIndex(minTime, indices);
    int minID = indices[0];
    TimeToIndex(maxTime, indices);
    int maxID = indices[1];

    for (int id = minID; id <= maxID; ++id) {
        sprintf(pathname, "%s%s", OCEAN_VECTOR_PATH, pathnames[id].c_str());
        if (m_pVelocities[id].OpenNrrdFile(pathname) != 0)
            return false;
    }

    memcpy(m_gridSize, m_pVelocities[minID].GetSize(), sizeof(int) * 2);
    return true;
}

template<typename T>
inline void ZD::CFlowOcean<T>::TimeToIndex(const value_type time, int * indices) const
{
    indices[0] = 0;
    indices[1] = 0;

    for (int i = 0; i < OCEAN_TOTAL_COUNT - 1; ++i) {
        if (time >= m_pTimes[i] && time < m_pTimes[i + 1]) {
            indices[0] = i;
            indices[1] = i + 1;
            return;
        }
    }
}

template<typename T>
inline const ZD::CPoint<T, 2> ZD::CFlowOcean<T>::P2C(const point_type& p) const
{
    const value_type dx = 0.0098388061184095932;
    const value_type dy = 0.0089781494355705149;
    point_type cp;
    cp[0] = (p[0] - m_boundingBox[0]) / dx;
    cp[1] = (p[1] - m_boundingBox[2]) / dy;
    return cp;
}

#endif
