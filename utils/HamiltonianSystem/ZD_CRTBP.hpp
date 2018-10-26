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
#ifndef _ZD_LIB_CRTBP_HPP_
#define _ZD_LIB_CRTBP_HPP_

#include "ZD_HamiltonianSystem.hpp"
#include "../Integrator/ZD_Integrator_RK45.hpp"
#include "../define.hpp"

#include <limits>

namespace ZD {
    const double Jupiter_Europa_mu = 2.528017705e-5;
    const double Earth_Moon_mu = 0.012150571430596;
    const double Sun_Saturn_mu = 2.85804e-4;
    const double Test_Mu = 0.5;

    const double Earth_Radius = 0.01657178735872750104646295521139;
    const double Moon_Radius = 0.00451653411469233989116785265802;

    template <typename T>
    class CCRTBP : public CHamiltonianSystem<T, 6> {
        enum CRTBP_SYSTEM {
            CRTBP_EARTH_MOON,
            CRTBP_JUPITER_EUROPA,
            CRTBP_SUN_SATURN,
            CRTBP_TEST
        };

    public:
        typedef CPoint<T, 6> vec6;

    public:
        /* x-xd plane */
        class CPoincarePlanePxVx : public CPoincarePlane<T, 6> {
        public:
            virtual bool Intersect(CIntegrator<T, 6> *integrator, vec6& v, T& t) const;
        };

        /* x-y plane */
        class CPoincarePlanePxPy : public CPoincarePlane<T, 6> {
        public:
            virtual bool Intersect(CIntegrator<T, 6> *integrator, vec6 &v, T &t) const;
        };

    private:
        T m_mu;
        T m_C;
        T m_eps;

    public:
        CCRTBP(const T mu, const T C, const T eps);
        CCRTBP(const CRTBP_SYSTEM sys, const T C, const T eps);
        CCRTBP(const char *name, const T C, const T eps);
        ~CCRTBP();

    public:
        //virtual T Next(vec6& p, const T startTime,
        //    const T direction, CIntegrator<T, 6> *integrator) const;

        /**/
        virtual T NextTime(vec6& p, const T startTime, const T time, const T direction) const;
        virtual T NextTime(vec6& p, const T startTime, const T time, const T direction,
            std::vector<vec6> &buffer) const;

        /**/
        virtual T NextLength(vec6& p, const T startTime, const T length, const T direction) const;
        virtual T NextLength(vec6& p, const T startTime, const T length, const T direction,
            std::vector<vec6> &buffer) const;

        /**/
        virtual T NextMap(vec6& p, const T startTime, const T direction,
            const CPoincarePlane<T, 6> *plane) const;
        virtual T NextMap(vec6& p, const T startTime, const T direction,
            const CPoincarePlane<T, 6> *plane, std::vector<vec6> &buffer) const;

        /* get velocity and acceleration,
        input: 6x1 vector, position followed by velocity,
        return: 6x1 vector, velocity followed by acceleration */
        //vec6 Interpolate(const vec6& in, const T& time) const;

        /* unproject a point in poincare map to phase space */
        bool Unproject(const CPoint<T, 2>& p, vec6& out) const;
        /* project a point in phase space to poincare map */
        void Project(const vec6& in, CPoint<T, 2>& p) const;

    public:
        /* get velocity and acceleration,
        input: 6x1 vector, position followed by velocity,
        return: 6x1 vector, velocity followed by acceleration */
        virtual inline const vec6 Velocity(const vec6& p, const T& t) const;

        virtual CIntegrator<T, 6> *CreateIntegrator() const;

    private:
        /* compute energy */
        T Energy(const vec6 &v) const;

        /* check distance to two massive body */
        bool CheckDistance(const vec6 &v) const;

        /* check periapse condition */
        bool CheckPeriapse(const vec6 &v) const;

        T JacobiConstant(const vec6 &v) const;

        bool yd(const T& x, const T& xd, T &yd) const;

        inline T Distance(const vec6& src, const vec6& dst) const;

    public:
        T Mu() const { return this->m_mu; }
        T C() const { return this->m_C; }
    };
}

template<typename T>
bool ZD::CCRTBP<T>::CPoincarePlanePxVx::Intersect(CIntegrator<T, 6>* integrator, typename CCRTBP<T>::vec6& v, T& t) const
{
    T t0 = integrator->t0();
    T t1 = integrator->t1();

    T y0 = integrator->v0()[1];
    T y1 = integrator->v1()[1];

    if (y0 * y1 >= 0.0)
        return false;        // no intersection

    T tt = (t0 + t1) / 2.0;
    T yf = integrator->v(tt)[1];
    T lastyf;
    
    const T eps = 1e-15;

    while ((std::fabs(yf) > eps && std::fabs(t1 - t0) > eps * std::fabs(t1))) {
        if (yf * y0 > 0.0) {
            y0 = yf;
            t0 = tt;
        }
        else {
            y1 = yf;
            t1 = tt;
        }
        tt = t0 + (t1 - t0) / 2.0;
        lastyf = yf;
        yf = integrator->v(tt)[1];
    }

    t = tt;
    v = integrator->v(tt);
    v[1] = 0.0;

    return true;
}

template<typename T>
bool ZD::CCRTBP<T>::CPoincarePlanePxPy::Intersect(CIntegrator<T, 6>* integrator, typename CCRTBP<T>::vec6& v, T& t) const
{
    return false;
}

template<typename T>
ZD::CCRTBP<T>::CCRTBP(const T mu, const T C, const T eps)
{
    this->m_mu = mu;
    this->m_C = C;
    this->m_eps = eps;
}

template<typename T>
ZD::CCRTBP<T>::CCRTBP(const CRTBP_SYSTEM sys, const T C, const T eps)
{
    if (CRTBP_SYSTEM::CRTBP_EARTH_MOON == sys)
        this->m_mu = Earth_Moon_mu;
    else if (CRTBP_SYSTEM::CRTBP_JUPITER_EUROPA == sys)
        this->m_mu = Jupiter_Europa_mu;
    else if (CRTBP_SYSTEM::CRTBP_SUN_SATURN == sys)
        this->m_mu = Sun_Saturn_mu;
    else
        this->m_mu = Test_Mu;

    this->m_C = C;
    this->m_eps = eps;
}

template<typename T>
ZD::CCRTBP<T>::CCRTBP(const char * name, const T C, const T eps)
{
    if (strcmp(name, "earth") == 0 || strcmp(name, "earth_moon") == 0 ||
        strcmp(name, "Earth") == 0 || strcmp(name, "Earth_Monn") == 0) {
        this->m_mu = Earth_Moon_mu;
    }
    else if (strcmp(name, "jupiter") == 0 || strcmp(name, "jupiter_europa") == 0 ||
        strcmp(name, "Jupiter") == 0 || strcmp(name, "Jupiter_Europa") == 0) {
        this->m_mu = Jupiter_Europa_mu;
    }
    else if (strcmp(name, "sun") == 0 || strcmp(name, "sun_saturn") == 0 ||
        strcmp(name, "Sun") == 0 || strcmp(name, "Sun_Saturn") == 0) {
        this->m_mu = Sun_Saturn_mu;
    }
    else if (strcmp(name, "test") == 0 || strcmp(name, "Test") == 0) {
        this->m_mu = Test_Mu;
    }
    else {
        std::cout << "Unknow system name: " << name << std::endl;
        exit(EXIT_FAILURE);
    }

    this->m_C = C;
    this->m_eps = eps;
}

template<typename T>
ZD::CCRTBP<T>::~CCRTBP()
{
}

template<typename T>
T ZD::CCRTBP<T>::NextTime(typename CCRTBP<T>::vec6& p, const T startTime, const T time, const T direction) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainTime = time;
    while (remainTime > 0.0) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        
        if (t > remainTime) {
            p = integrator->v(startTime + time);
            SafeDelete(integrator);
            return startTime + time;
        }
        else {
            integrationTime += t;
            remainTime -= t;
            p = integrator->v1();
        }
    }
}

template<typename T>
T ZD::CCRTBP<T>::NextTime(typename CCRTBP<T>::vec6& p, const T startTime, const T time, const T direction, std::vector<vec6>& buffer) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainTime = time;
    while (remainTime > 0.0) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        
        if (t > remainTime) {
            p = integrator->v(startTime + time);
            buffer.push_back(p);
            return startTime + time;
        }
        else {
            integrationTime += t;
            remainTime -= t;
            p = integrator->v1();
            buffer.push_back(p);
        }
    }
}

template<typename T>
T ZD::CCRTBP<T>::NextLength(typename CCRTBP<T>::vec6& p, const T startTime, const T length, const T direction) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainLength = length;

    const T eps = 1e-15;

    while (remainLength > 0.0) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        T length = this->Distance(integrator->v0(), integrator->v1());
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec6 mp = integrator->v(tt);
            length = this->Distance(mp, integrator->v0());
            while (std::abs(length - remainLength) > eps && std::abs(t1 - t0) > eps * std::abs(t1)) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = this->Distance(mp, integrator->v0());
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
T ZD::CCRTBP<T>::NextLength(typename CCRTBP<T>::vec6& p, const T startTime, const T length, const T direction, std::vector<vec6>& buffer) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    T remainLength = length;

    const T eps = 1e-15;

    while (remainLength > 0.0) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, startTime + integrationTime, direction, this);
        T length = this->Distance(integrator->v0(), integrator->v1());
        if (length <= 0.0)
            return std::numeric_limits<T>::quiet_NaN();
        if (length > remainLength) {
            T t0 = integrator->t0();
            T t1 = integrator->t1();

            T tt = (t0 + t1) / 2.0;
            vec6 mp = integrator->v(tt);
            length = this->Distance(mp, integrator->v0());
            while (std::abs(length - remainLength) > eps && std::abs(t1 - t0) > eps * std::abs(t1)) {
                if (length > remainLength) {
                    t1 = tt;
                }
                else {
                    t0 = tt;
                }
                tt = (t0 + t1) / 2.0;
                mp = integrator->v(tt);
                length = this->Distance(mp, integrator->v0());
            }
            p = mp;
            buffer.push_back(p);
            SafeDelete(integrator);
            return tt;
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
            remainLength -= length;
            buffer.push_back(p);
        }
    }
}

template<typename T>
T ZD::CCRTBP<T>::NextMap(typename CCRTBP<T>::vec6& p, const T startTime, const T direction, const CPoincarePlane<T, 6>* plane) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    vec6 ip;
    T it;

    while (true) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, 0.0, direction, this);

        if (integrator->v1()[4] >= 0.0 && plane->Intersect(integrator, ip, it)) {
            p = ip;
            SafeDelete(integrator);
            integrationTime += it * direction;
            return (startTime + integrationTime);
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
        }
    }
}

template<typename T>
T ZD::CCRTBP<T>::NextMap(typename CCRTBP<T>::vec6& p, const T startTime, const T direction, const CPoincarePlane<T, 6>* plane, std::vector<vec6>& buffer) const
{
    CIntegrator<T, 6> *integrator = this->CreateIntegrator();
    T integrationTime = 0.0;
    vec6 ip;
    T it;

    while (true) {
        if (CheckDistance(p) == false) {
            SafeDelete(integrator);
            return std::numeric_limits<T>::quiet_NaN();
        }

        T t = integrator->Next(p, 0.0, direction, this);

        if (integrator->v1()[4] >= 0.0 && plane->Intersect(integrator, ip, it)) {
            p = ip;
            buffer.push_back(p);
            SafeDelete(integrator);
            integrationTime += it * direction;
            return (startTime + integrationTime);
        }
        else {
            integrationTime += t * direction;
            p = integrator->v1();
            buffer.push_back(p);
        }
    }
}

template<typename T>
bool ZD::CCRTBP<T>::Unproject(const CPoint<T, 2>& p, typename CCRTBP<T>::vec6& out) const
{
    memset(out.m_data, 0, sizeof(vec6));
    out.m_data[0] = p[0];        /* x */
    out.m_data[3] = p[1];        /* dx */
    return yd(out.m_data[0], out.m_data[3], out.m_data[4]);
}

template<typename T>
void ZD::CCRTBP<T>::Project(const typename CCRTBP<T>::vec6& in, CPoint<T, 2>& p) const
{
    p[0] = in.m_data[0];        /* x */
    p[1] = in.m_data[3];        /* dx */
}

template<typename T>
inline const typename ZD::CCRTBP<T>::vec6 ZD::CCRTBP<T>::Velocity(const typename CCRTBP<T>::vec6& in, const T & t) const
{
    T x1 = in[0] + m_mu;
    T x2 = in[0] - 1.0 + m_mu;
    T r1 = std::sqrt(x1 * x1 + in[1] * in[1] + in[2] * in[2]);
    T r2 = std::sqrt(x2 * x2 + in[1] * in[1] + in[2] * in[2]);
    T r13 = r1 * r1 * r1;
    T r23 = r2 * r2 * r2;

    vec6 out;
    // velocity
    out[0] = in[3];
    out[1] = in[4];
    out[2] = in[5];

    // acceleration
    out[3] = 2.0 * in[4] + in[0] - (1.0 - m_mu) * x1 / r13 - m_mu * x2 / r23;
    out[4] = -2.0 * in[3] + in[1] - (1.0 - m_mu) * in[1] / r13 - m_mu * in[1] / r23;
    out[5] = -(1.0 - m_mu) * in[2] / r13 - m_mu * in[2] / r23;

    return out;
}

template<typename T>
ZD::CIntegrator<T, 6>* ZD::CCRTBP<T>::CreateIntegrator() const
{
    return new CIntegratorRK45<T, 6>(this->m_eps, 0.1);
}

template<typename T>
bool ZD::CCRTBP<T>::CheckDistance(const typename CCRTBP<T>::vec6& v) const
{
    T x1 = v[0] + m_mu;
    T x2 = v[0] - 1.0 + m_mu;
    T r1s = x1 * x1 + v[1] * v[1] + v[2] * v[2];
    T r2s = x2 * x2 + v[1] * v[1] + v[2] * v[2];

    const T epsilon = 1e-16;
    if (r1s < epsilon || r2s < epsilon)
        return false;
    else
        return true;
}

template<typename T>
bool ZD::CCRTBP<T>::yd(const T & x, const T & xd, T & yd) const
{
    T tmp = x * x + 2.0 * (1.0 - m_mu) / std::abs(x + m_mu) + 2.0 * m_mu / std::abs(x - 1.0 + m_mu) - xd * xd - m_C;
    if (tmp < 0.0)
        return false;
    else {
        yd = std::sqrt(tmp);
        return true;
    }
}

template<typename T>
inline T ZD::CCRTBP<T>::Distance(const vec6& src, const vec6& dst) const
{
    T dx = dst[0] - src[0];
    T dy = dst[1] - src[1];
    T dz = dst[2] - src[2];

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

#endif
