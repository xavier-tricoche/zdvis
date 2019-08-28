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
#ifndef _ZD_LIB_FLOW_IONIZATION_FRONT_HPP_
#define _ZD_LIB_FLOW_IONIZATION_FRONT_HPP_

#include "ZD_Flow.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "utils/define.hpp"
#include "utils/configure.hpp"
#include "utils/Tool/ZD_StringTool.hpp"
#include "utils/Integrator/ZD_Integrator_RK4.hpp"
#include "utils/Integrator/ZD_Integrator_RK45.hpp"

namespace ZD {
    template <typename T>
    class CFlowIonizationFront : public CFlow<T, 3> {
    public:
        typedef T value_type;
        typedef CPoint<T, 3> point_type;
    private:
        CField3<T, 3> *m_pVectorFields;
        CField3<T, 1> *m_pScalarFields;

        value_type *m_pTimes;
        
        value_type m_boundingBox[6];
        value_type m_minTime, m_maxTime;

        point_type m_dispersion;

    public:
        CFlowIonizationFront(const value_type minTime, const value_type maxTime,
            const bool readVelocities, const char *scalarName = nullptr);
        ~CFlowIonizationFront();

    protected:
        virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        inline void TimeToIndex(const value_type& time, int *indices) const;
        bool ReadTime();
        bool ReadVelocity(const value_type minTime, const value_type maxTime);
        bool ReadScalar(const value_type minTime, const value_type maxTime, const char *name);
    };
}

template<typename T>
inline ZD::CFlowIonizationFront<T>::CFlowIonizationFront(const value_type minTime, const value_type maxTime, 
    const bool readVelocities, const char * scalarName)
{
    // convert all unit to cm and cm/s
    const double unitScale = 3.08567758149137 * 1e18;
    m_boundingBox[0] = 0.0;
    //m_boundingBox[1] = 0.599 * unitScale;
    m_boundingBox[1] = 0.599;
    m_boundingBox[2] = 0.0;
    //m_boundingBox[3] = 0.247 * unitScale;
    m_boundingBox[3] = 0.247;
    m_boundingBox[4] = 0.0;
    //m_boundingBox[5] = 0.247 * unitScale;
    m_boundingBox[5] = 0.247;

    m_pTimes = new T[IONIZATION_FRONT_TOTAL_COUNT];
    ReadTime();

    if (readVelocities) {
        m_pVectorFields = new CField3<T, 3>[IONIZATION_FRONT_TOTAL_COUNT];
        ReadVelocity(minTime, maxTime);
    }
    else {
        m_pVectorFields = nullptr;
    }

    if (scalarName != nullptr) {
        m_pScalarFields = new CField3<T, 1>[IONIZATION_FRONT_TOTAL_COUNT];
        ReadScalar(minTime, maxTime, scalarName);
    }
    else {
        m_pScalarFields = nullptr;
    }

    // grid size 600 x 248 x 248
    m_dispersion[0] = (m_boundingBox[1] - m_boundingBox[0]) / (600.0 - 1.0);
    m_dispersion[1] = (m_boundingBox[3] - m_boundingBox[2]) / (248.0 - 1.0);
    m_dispersion[2] = (m_boundingBox[5] - m_boundingBox[4]) / (248.0 - 1.0);
}

template <typename T>
ZD::CFlowIonizationFront<T>::~CFlowIonizationFront()
{
    SafeDeleteArray(m_pVectorFields);
    SafeDeleteArray(m_pScalarFields);
    SafeDeleteArray(m_pTimes);
}

template<typename T>
inline ZD::CIntegrator<T, 3>* ZD::CFlowIonizationFront<T>::CreateIntegrator() const
{
    return new CIntegratorRK4<T, 3>(this->m_stepSize);
}

template<typename T>
inline const ZD::CPoint<T, 3> ZD::CFlowIonizationFront<T>::Velocity(const point_type& p, const value_type & time) const
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

    return dir;
}

template<typename T>
inline const T ZD::CFlowIonizationFront<T>::Scalar(const point_type& p, const value_type & time) const
{
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

template<typename T>
inline const bool ZD::CFlowIonizationFront<T>::CheckPosition(const point_type& p) const
{
    if (p[0] < this->m_boundingBox[0] || p[0] > this->m_boundingBox[1] ||
        p[1] < this->m_boundingBox[2] || p[1] > this->m_boundingBox[3] ||
        p[2] < this->m_boundingBox[4] || p[2] > this->m_boundingBox[5]) {
        return false;
    }
    else {
        return true;
    }
}


template<typename T>
inline void ZD::CFlowIonizationFront<T>::TimeToIndex(const value_type & time, int * indices) const
{
    indices[0] = 0;
    indices[1] = 0;
    for (int i = 0; i < IONIZATION_FRONT_TOTAL_COUNT - 1; ++i) {
        if (time >= m_pTimes[i] && time < m_pTimes[i + 1]) {
            indices[0] = i;
            indices[1] = i + 1;
            return;
        }
    }
}

template<typename T>
inline bool ZD::CFlowIonizationFront<T>::ReadTime()
{
    FILE *fp = fopen(IONIZATION_FRONT_TIME_PATHNAME, "r");
    if (fp != nullptr) {
        char temp[256];
        fscanf(fp, "%s %s\n", temp, temp);
        for (int i = 0; i < IONIZATION_FRONT_TOTAL_COUNT; ++i) {
            double time;
            fscanf(fp, "%s %lf\n", temp, &time);
            m_pTimes[i] = time;
        }
        return true;
    }
    else {
        fprintf(stderr, "Error! Cannot open time file.\n");
        return false;
    }
}

template<typename T>
bool ZD::CFlowIonizationFront<T>::ReadVelocity(const value_type minTime, const value_type maxTime)
{
    int temp[2];
    TimeToIndex(minTime, temp);
    int minID = temp[0];
    TimeToIndex(maxTime, temp);
    int maxID = temp[1];

    for (int i = minID; i <= maxID; ++i) {
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%svelocity.%04d.nrrd", IONIZATION_FRONT_VECTOR_PATH, i);
        m_pVectorFields[i].OpenNrrdFile(pathname);
    }

    return true;
}

template<typename T>
bool ZD::CFlowIonizationFront<T>::ReadScalar(const value_type minTime, const value_type maxTime, const char * name)
{
    std::string tmp = CStringTool::ToUpper(name);
    std::string path;
    std::string scalarName;
    if (strcmp(tmp.c_str(), "DENSITY") == 0) {
        path = IONIZATION_FRONT_TOTAL_PARTICLE_DENSITY_PATH;
        scalarName = "total_particle_density";
    }
    else if (strcmp(tmp.c_str(), "TEMPERATURE") == 0) {
        path = IONIZATION_FRONT_GAS_TEMPERATURE_PATH;
        scalarName = "gas_temperature";
    }
    else if (strcmp(tmp.c_str(), "H_MASS") == 0) {
        path = IONIZATION_FRONT_H_MASS_PATH;
        scalarName = "h_mass";
    }
    else if (strcmp(tmp.c_str(), "H+_MASS") == 0) {
        path = IONIZATION_FRONT_HP_MASS_PATH;
        scalarName = "h+_mass";
    }
    else if (strcmp(tmp.c_str(), "HE_MASS") == 0) {
        path = IONIZATION_FRONT_HE_MASS_PATH;
        scalarName = "he_mass";
    }
    else if (strcmp(tmp.c_str(), "HE+_MASS") == 0) {
        path = IONIZATION_FRONT_HEP_MASS_PATH;
        scalarName = "he+_mass";
    }
    else if (strcmp(tmp.c_str(), "HE++_MASS") == 0) {
        path = IONIZATION_FRONT_HEPP_MASS_PATH;
        scalarName = "he++_mass";
    }
    else if (strcmp(tmp.c_str(), "H-_MASS") == 0) {
        path = IONIZATION_FRONT_HM_MASS_PATH;
        scalarName = "h-_mass";
    }
    else if (strcmp(tmp.c_str(), "H2_MASS") == 0) {
        path = IONIZATION_FRONT_H_2_MASS_PATH;
        scalarName = "h_2_mass";
    }
    else if (strcmp(tmp.c_str(), "H2+_MASS") == 0) {
        path = IONIZATION_FRONT_H_2P_MASS_PATH;
        scalarName = "h_2+_mass";
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
        sprintf(pathname, "%s%s.%04d.nrrd", path.c_str(), scalarName.c_str(), i);
        m_pScalarFields[i].OpenNrrdFile(pathname);
    }

    return true;
}

template <typename T>
inline void ZD::CFlowIonizationFront<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    min[2] = m_boundingBox[4];

    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];
    max[2] = m_boundingBox[5];
}

#endif
