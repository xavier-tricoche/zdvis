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
#ifndef _ZD_LIB_KINEMATIC_CHAOS_MODEL_HPP_
#define _ZD_LIB_KINEMATIC_CHAOS_MODEL_HPP_

#include "ZD_HamiltonianSystem.hpp"
#include "utils/define.hpp"
#include "utils/Integrator/ZD_Integrator_RK45.hpp"
#include <limits>

namespace ZD {
    /* assume the poincare section in Kinematic chaos model is the x-z plane */
    template <typename T>
    class CKinematicChaosModel : public CHamiltonianSystem<T, 3> {
    public:
        typedef CPoint<T, 3> vec3;

    public:
        /* x-z plane */
        class CPoincarePlaneXZ : public CPoincarePlane <T, 3> {
        public:
            virtual bool Intersect(CIntegrator<T, 3> *integrator, vec3 &v, T &t) const;
            //virtual void Project(const CPoint<T, 2>& in, vec3& out) const;
            //virtual void Unproject(const vec3& in, CPoint<T, 2>& out) const;
        };

    
    private:
        T m_a;
        T m_alpha, m_beta;
        T m_x0;
        T m_e1, m_e2;
        T m_n;
        T m_eps;

    public:
        CKinematicChaosModel(const T a, const T alpha, const T beta, const T x0, 
            const T e1, const T e2, const T eps = 1e-12);
        ~CKinematicChaosModel();

    public:
        //virtual T Next(vec3& p, const T time, const T direction, CIntegrator<T, 3> *integrator) const;

        /* integral for a fixed time */
        virtual T NextTime(vec3& p, const T startTime, const T time, const T direction) const;
        virtual T NextTime(vec3& p, const T startTime, const T time, const T direction,
            std::vector<vec3> &buffer) const;

        /* integral for a fixed length */
        virtual T NextLength(vec3& p, const T startTime, const T length, const T direction) const;
        virtual T NextLength(vec3& p, const T startTime, const T length, const T direction,
            std::vector<vec3> &buffer) const;

        /* integral for one map iteration */
        virtual T NextMap(vec3& p, const T startTime, const T direction,
            const CPoincarePlane<T, 3> *plane) const;
        virtual T NextMap(vec3& p, const T startTime, const T direction,
            const CPoincarePlane<T, 3> *plane, std::vector<vec3> &buffer) const;

    public:
        /* input and output use a Cartesian coordinate system */
        virtual inline const vec3 Velocity(const vec3& p, const T& t) const;

        virtual CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        /* input uses a Cartesian coordinate system */
        inline bool CheckPoint(const vec3& p) const;
    };
}

template<typename T>
bool ZD::CKinematicChaosModel<T>::CPoincarePlaneXZ::Intersect(CIntegrator<T, 3> *integrator, vec3 &v, T &t) const
{
    T t0 = integrator->t0();
    T t1 = integrator->t1();

    T y0 = integrator->v0()[1];
    T y1 = integrator->v1()[1];

    if (y0 * y1 >= 0.0)
        return false;        // no intersection

    T tt = (t0 + t1) / 2.0;
    T yf = integrator->v(tt)[1];

    const T eps = 1e-15;

    while (std::fabs(yf) > eps && std::fabs(t1-t0) > eps * std::fabs(t1)) {
        if (yf * y0 > 0.0) {
            y0 = yf;
            t0 = tt;
        }
        else {
            y1 = yf;
            t1 = tt;
        }
        tt = (t0 + t1) / 2.0;
        yf = integrator->v(tt)[1];
    }

    t = tt;
    v = integrator->v(tt);
    v[1] = 0.0;

    return true;
}

template<typename T>
ZD::CKinematicChaosModel<T>::CKinematicChaosModel(const T a, const T alpha, const T beta, 
    const T x0, const T e1, const T e2, const T eps)
{
    this->m_a     = a;
    this->m_alpha = alpha;
    this->m_beta  = beta;
    this->m_x0    = x0;
    this->m_e1    = e1;
    this->m_e2    = e2;

    this->m_n = 1;

    this->m_eps = eps;
}

template<typename T>
ZD::CKinematicChaosModel<T>::~CKinematicChaosModel()
{

}

template<typename T>
T ZD::CKinematicChaosModel<T>::NextTime(vec3 & p, const T startTime, const T time, const T direction) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainTime = time;

    while (remainTime > 0.0) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        integrationTime += t * direction;
        if (std::fabs(t) > std::fabs(remainTime)) {
            p = integrator->v(startTime + time * direction);
            SafeDelete(integrator);
            return startTime + time * direction;
        }
        else {
            p = integrator->v1();
            remainTime -= std::fabs(t);
        }
    }
}

template<typename T>
T ZD::CKinematicChaosModel<T>::NextTime(vec3 & p, const T startTime, const T time, const T direction, std::vector<vec3>& buffer) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainTime = time;

    while (remainTime > 0.0) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        if (std::fabs(t) > std::fabs(remainTime)) {
            p = integrator->v(startTime + time * direction);
            buffer.push_back(p);
            SafeDelete(integrator);
            return startTime + time * direction;
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
            remainTime -= std::fabs(t);
            buffer.push_back(p);
        }
    }
}

template<typename T>
T ZD::CKinematicChaosModel<T>::NextLength(vec3 & p, const T startTime, const T length, const T direction) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainLength = length;

    while (remainLength > 0.0) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        T length = (integrator->v1() - integrator->v0()).Length();
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec3 mp = integrator->v(tt);
            length = (mp - integrator->v0()).Length();
            while (std::abs(length - remainLength) > this->m_eps) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = (mp - integrator->v0()).Length();
            }
            p = mp;
            SafeDelete(integrator);
            return tt;
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
            remainLength -= length;
        }
    }
}


template<typename T>
T ZD::CKinematicChaosModel<T>::NextLength(vec3 & p, const T startTime, const T length, const T direction, std::vector<vec3>& buffer) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainLength = length;

    while (remainLength > 0.0) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        T length = (integrator->v1() - integrator->v0()).Length();
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec3 mp = integrator->v(tt);
            length = (mp - integrator->v0()).Length();
            while (std::abs(length - remainLength) > this->m_eps) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = (mp - integrator->v0()).Length();
            }
            p = mp;
            buffer.push_back(p);
            SafeDelete(integrator);
            return tt;
        }
        else {
            integrationTime += t * direction;
            remainLength -= length;
            p = integrator->v1();
            buffer.push_back(p);
        }
    }
}

template<typename T>
T ZD::CKinematicChaosModel<T>::NextMap(vec3 & p, const T startTime, const T direction, 
    const CPoincarePlane<T, 3> *plane) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    vec3 ip;
    T it;

    while (true) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        
        if (plane->Intersect(integrator, ip, it)) {
            p = ip;
            SafeDelete(integrator);
            return it;
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
        }
    }
}

template<typename T>
T ZD::CKinematicChaosModel<T>::NextMap(vec3 & p, const T startTime, const T direction, 
    const CPoincarePlane<T, 3> *plane, std::vector<vec3>& buffer) const
{
    CIntegrator<T, 3> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    vec3 ip;
    T it;

    while (true) {
        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        
        if (plane->Intersect(integrator, ip, it)) {
            p = ip;
            buffer.push_back(p);
            SafeDelete(integrator);
            return it;
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
            buffer.push_back(p);
        }
    }
}

template<typename T>
inline const typename ZD::CKinematicChaosModel<T>::vec3 ZD::CKinematicChaosModel<T>::Velocity(const vec3 & p, const T& t) const
{
    const T& x = p[0];
    const T& y = p[1];
    const T& z = p[2];

    const T r = std::sqrt(x * x + y * y);
    const T theta = std::atan2(x, y);
    
    const T A = (1.0 - 2.0 * z) * (this->m_a - r);
    const T B = 1.0 - this->m_beta * z;
    const T C = (this->m_a * this->m_a - 2.0 * r * r) * std::cos(this->m_n * theta);
    const T D = (this->m_a * this->m_a - r * r) * std::sin(this->m_n * theta);

    vec3 res;
    res[0] = -x / 3.0 * A - this->m_alpha * y + this->m_e1 * B * y * (x - this->m_x0)                                           - 0.5 * this->m_e2 * B * (2.0 * y * C - this->m_n * x * D);
    res[1] = -y / 3.0 * A + this->m_alpha * x + this->m_e1 * B * (0.5 * (this->m_a * this->m_a - r * r) - x * (x - this->m_x0)) + 0.5 * this->m_e2 * B * (2.0 * x * C + this->m_n * y * D);
    res[2] = z * (1.0 - z) * (2.0 * this->m_a / 3.0 - r);

    //res[0] = -x / 3.0 * A - this->m_alpha * y;
    //res[1] = -y / 3.0 * A + this->m_alpha * x;
    //res[2] = z * (1.0 - z) * (2.0 * this->m_a / 3.0 - r);

    return res;
}

template<typename T>
ZD::CIntegrator<T, 3>* ZD::CKinematicChaosModel<T>::CreateIntegrator() const
{
    return new CIntegratorRK45<T, 3>(this->m_eps, 1.0);
}

template <typename T>
inline bool ZD::CKinematicChaosModel<T>::CheckPoint(const vec3& src) const
{
    return true;
}


#endif
