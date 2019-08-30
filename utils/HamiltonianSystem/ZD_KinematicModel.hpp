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
#ifndef _ZD_LIB_KINEMATIC_MODEL_HPP
#define _ZD_LIB_KINEMATIC_MODEL_HPP

#include <vector>

#include "utils/Base/ZD_Point.hpp"
#include "utils/Integrator/ZD_Integrator.hpp"
#include "utils/define.hpp"

#include <limits>

using namespace ZD;

namespace ZD {
    /* assume the poincare section in Kinematic model is the x-z plane */
    template <typename T>
    class CKinematicModel {
    public:
        typedef CPoint<T, 3> vec3;
    private:
        T m_a;
        T m_alpha;
        T m_eps;
        T m_maxStepSize;

    public:
        CKinematicModel(const T a, const T alpha, const T maxStepSize, const T eps=1e-12);
        ~CKinematicModel();

    public:
        T Next(vec3& p, const T time, const T maxStepSize, CIntegrator<T, 3> *integrator) const;

        /* integral for a fixed time */
        T NextTime(vec3& p, const T startTime, const T time, const T direction,
            CIntegrator<T, 3> *integrator) const;
        T NextTime(vec3& p, const T startTime, const T time, const T direction,
            CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const;

        /* integral for a fixed length */
        T NextLength(vec3& p, const T startTime, const T length, const T direction,
            CIntegrator<T, 3> *integrator) const;
        T NextLength(vec3& p, const T startTime, const T length, const T direction,
            CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const;

        /* integral for one map iteration */
        T NextMap(vec3& p, const T startTime, const T direction,
            CIntegrator<T, 3> *integrator) const;
        T NextMap(vec3& p, const T startTime, const T direction,
            CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const;

    public:
        /* input and output use a cylinderical coordinate system */
        /*  input: p[0] -- radial distance                       */
        /*         p[1] -- angular coordinate                    */
        /*         p[2] -- height                                */
        /* output: p[0] -- radial velocity component             */
        /*         p[1] -- azimuthal velocity component          */
        /*         p[2] -- vertical velocity component           */
        inline const vec3 GetVelocity(const vec3& p, const T time) const;

    public:
        /* project Cartesian coordinate system into cylinderical coordniate system */
        inline void ProjectCoordinate(const vec3& src, vec3& dst) const;
        inline vec3 ProjectCoordinate(const vec3& src) const;
        /* unproject cylinderical coordniate system inot Cartesian coordinate system */
        inline void UnprojectCoordinate(const vec3& src, vec3& dst) const;
        inline vec3 UnprojectCoordinate(const vec3& src) const;

    public:
        /* input uses a cylinderical coordinate system */
        inline bool CheckPoint(const vec3& p) const;
    };
}

template <typename T>
ZD::CKinematicModel<T>::CKinematicModel(const T a, const T alpha, const T maxStepSize, const T eps)
{
    this->m_a = a;
    this->m_alpha = alpha;
    this->m_eps = eps;
    this->m_maxStepSize = std::abs(maxStepSize);
}

template <typename T>
ZD::CKinematicModel<T>::~CKinematicModel()
{

}

template <typename T>
T ZD::CKinematicModel<T>::Next(vec3& p, const T time, const T maxStepSize, CIntegrator<T, 3> *integrator) const
{
    if (CheckPoint(p) == false)
        return -1.0;


    T dt = integrator->Next(p, time, maxStepSize);
    p = integrator->v1();
    return dt;
}

template <typename T>
T ZD::CKinematicModel<T>::NextTime(vec3& p, const T startTime, const T time, const T direction,
    CIntegrator<T, 3> *integrator) const
{
    T integrationTime = 0.0;
    T remainTime = time;
    while (remainTime > 0.0) {
        if (CheckPoint(p) == false)
            return -1.0;

        T t = Next(p, startTime + integrationTime, direction*this->m_maxStepSize, integrator);
        integrationTime += t;
        if (t > remainTime) {
            p = integrator->v(startTime + time);
            return startTime + time;
        }
        else {
            remainTime -= t;
        }
    }
}

template <typename T>
T ZD::CKinematicModel<T>::NextTime(vec3& p, const T startTime, const T time, const T direction,
                               CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const
{
    T integrationTime = 0.0;
    T remainTime = time;
    while (remainTime > 0.0) {
        if (CheckPoint(p) == false)
            return -1.0;

        T t = Next(p, startTime + integrationTime, direction*this->m_maxStepSize, integrator);
        integrationTime += t;
        if (t > remainTime) {
            p = integrator->v(startTime + time);
            buffer.push_back(this->UnprojectCoordinate(p));
            return startTime + time;
        }
        else {
            remainTime -= t;
            buffer.push_back(this->UnprojectCoordinate(p));
        }
    }
}

template <typename T>
T ZD::CKinematicModel<T>::NextLength(vec3& p, const T startTime, const T length, const T direction,
    CIntegrator<T, 3> *integrator) const
{
    T integrationTime = 0.0;
    T remainLength = length;
    while (remainLength > 0.0) {
        if (CheckPoint(p) == false)
            return -1.0;

        vec3 pp = p;
        T t = Next(p, startTime + integrationTime, direction*this->m_maxStepSize, integrator);
        integrationTime += t * direction;
        T length = (this->UnprojectCoordinate(p) - this->UnprojectCoordinate(pp)).Length();
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec3 mp = integrator->v(tt);
            length = (this->UnprojectCoordinate(mp) - this->UnprojectCoordinate(pp)).Length();
            while (std::abs(length - remainLength) > this->m_eps) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = (this->UnprojectCoordinate(mp) - this->UnprojectCoordinate(pp)).Length();
            }
            p = mp;
            return tt;
        }
        else {
            remainLength -= length;
        }
    }
}

template <typename T>
T ZD::CKinematicModel<T>::NextLength(vec3& p, const T startTime, const T length, const T direction,
                                 CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const
{
    T integrationTime = 0.0;
    T remainLength = length;
    while (remainLength > 0.0) {
        if (CheckPoint(p) == false)
            return -1.0;

        vec3 pp = p;
        T t = Next(p, startTime + integrationTime, direction*this->m_maxStepSize, integrator);
        integrationTime += t;
        T length = (this->UnprojectCoordinate(p) - this->UnprojectCoordinate(pp)).Length();
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec3 mp = integrator->v(tt);
            length = (this->UnprojectCoordinate(mp) - this->UnprojectCoordinate(pp)).Length();
            while (std::abs(length - remainLength) > this->m_eps) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = (this->UnprojectCoordinate(mp) - this->UnprojectCoordinate(pp)).Length();
            }
            p = mp;
            buffer.push_back(this->UnprojectCoordinate(p));
            return tt;
        }
        else {
            remainLength -= length;
            buffer.push_back(this->UnprojectCoordinate(p));
        }
    }
}

template <typename T>
T ZD::CKinematicModel<T>::NextMap(vec3& p, const T startTime, const T direction,
    CIntegrator<T, 3> *integrator) const
{
    return 0.0;
}

template <typename T>
T ZD::CKinematicModel<T>::NextMap(vec3& p, const T startTime, const T direction,
                              CIntegrator<T, 3> *integrator, std::vector<vec3> &buffer) const
{
    return 0.0;
}


template <typename T>
inline const typename ZD::CKinematicModel<T>::vec3 ZD::CKinematicModel<T>::GetVelocity(const vec3& p, const T time) const
{
    const T& r = p[0];
    const T& z = p[2];

    vec3 res;
    res[0] = (2.0 * z - 1.0) / 3.0 * r * (this->m_a - r);
    res[1] = this->m_alpha * r;
    res[2] = z * (1.0 - z) * (2.0 * this->m_a / 3.0 - r);

    return res;
}

template<typename T>
inline void ZD::CKinematicModel<T>::ProjectCoordinate(const vec3 & src, vec3 & dst) const
{
    dst[0] = std::sqrt(src[0] * src[0] + src[1] * src[1]);
    dst[1] = std::atan2(src[0], src[1]);
    if (dst[1] < 0.0)
        dst[1] += 2.0 * ZD_PI;
    dst[2] = src[2];
}

template <typename T>
inline typename ZD::CKinematicModel<T>::vec3 ZD::CKinematicModel<T>::ProjectCoordinate(const vec3& src) const
{
    vec3 dst;
    this->ProjectCoordinate(src, dst);
    return dst;
}

template <typename T>
void ZD::CKinematicModel<T>::UnprojectCoordinate(const vec3& src, vec3& dst) const
{
    dst[0] = std::cos(src[1]) * src[0];
    dst[1] = std::sin(src[1]) * src[0];
    dst[2] = src[2];
}

template <typename T>
inline typename ZD::CKinematicModel<T>::vec3 ZD::CKinematicModel<T>::UnprojectCoordinate(const vec3& src) const
{
    vec3 dst;
    this->UnprojectCoordinate(src, dst);
    return dst;
}

template <typename T>
inline bool ZD::CKinematicModel<T>::CheckPoint(const vec3& p) const
{
    if (p[0] < 0.0 || p[0] > this->m_a)
        return false;
    if (p[2] < 0.0 || p[2] > 1.0)
        return false;

    return true;
}


#endif
