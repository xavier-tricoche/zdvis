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
#ifndef _ZD_LIB_INTEGRATOR_RK56_HPP_
#define _ZD_LIB_INTEGRATOR_RK56_HPP_

#include "ZD_Integrator.hpp"
#include <algorithm>

namespace ZD {
    template <typename T, unsigned int N>
    class CIntegratorRK45 : public CIntegrator<T, N> {
    private:
        T m_reltol;
        T m_abstol;
        T m_facold;

        CPoint<T, N> m_k1;

        int m_steps;
        int m_acceptedSteps;
        int m_rejectedSteps;
        int m_totalEvals;

        T m_totalError;

        T m_stepSize;
        T m_maxStepSize;

        /* constants */
        static const T m_safe;
        static const T m_epsilon;
        static const T m_facl;
        static const T m_facr;
        static const T m_beta;

        static const T A21;
        static const T A31, A32;
        static const T A41, A42, A43;
        static const T A51, A52, A53, A54;
        static const T A61, A62, A63, A64, A65;
        static const T A71, A73, A74, A75, A76;
        static const T C2, C3, C4, C5;
        static const T D1, D3, D4, D5, D6, D7;
        static const T E1, E3, E4, E5, E6, E7;

    public:
        CIntegratorRK45();
        CIntegratorRK45(const T eps, const T maxStepSize);
        ~CIntegratorRK45();

        T   TotalError() const { return m_totalError; }
        int TotalEvals() const { return m_totalEvals; }

        //void Init(const T eps);

        virtual T Next(const CPoint<T, N>& p, const T& t, const T direction,
            const CDynamicSystem<T, N> *sys);

        virtual CPoint<T, N> v(const T& t) const;


        //virtual void Reset();
    };
}


/* class CIntegration3DRK45 */
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::m_safe = 0.9;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::m_epsilon = std::numeric_limits<T>::epsilon();
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::m_facl = 0.2;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::m_facr = 10.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::m_beta = 0.04;

template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A21 = 0.2;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A31 = 3.0 / 40.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A32 = 9.0 / 40.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A41 = 44.0 / 45.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A42 = -56.0 / 15.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A43 = 32.0 / 9.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A51 = 19372.0 / 6561.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A52 = -25360.0 / 2187.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A53 = 64448.0 / 6561.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A54 = -212.0 / 729.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A61 = 9017.0 / 3168.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A62 = -355.0 / 33.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A63 = 46732.0 / 5247.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A64 = 49.0 / 176.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A65 = -5103.0 / 18656.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A71 = 35.0 / 384.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A73 = 500.0 / 1113.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A74 = 125.0 / 192.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A75 = -2187.0 / 6784.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::A76 = 11.0 / 84.0;

template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::C2 = 0.2;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::C3 = 0.3;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::C4 = 0.8;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::C5 = 8.0 / 9.0;

template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D1 = -12715105075.0 / 11282082432.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D3 = 87487479700.0 / 32700410799.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D4 = -10690763975.0 / 1880347072.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D5 = 701980252875.0 / 199316789632.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D6 = -1453857185.0 / 822651844.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::D7 = 69997945.0 / 29380423.0;

template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E1 = 71.0 / 57600.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E3 = -71.0 / 16695.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E4 = 71.0 / 1920.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E5 = -17253.0 / 339200.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E6 = 22.0 / 525.0;
template <typename T, unsigned int N> const T ZD::CIntegratorRK45<T, N>::E7 = -1.0 / 40.0;

template <typename T, unsigned int N>
ZD::CIntegratorRK45<T, N>::CIntegratorRK45()
{
    m_reltol = 1e-10;
    m_abstol = 1e-10;

    m_steps = 0;
    m_acceptedSteps = 0;
    m_rejectedSteps = 0;
    m_totalEvals = 0;
    m_totalError = 0;

    m_facold = 0;

    m_stepSize = 1e-12;
    m_maxStepSize = 0.01;
}

template <typename T, unsigned int N>
ZD::CIntegratorRK45<T, N>::CIntegratorRK45(const T eps, const T maxStepSize)
{
    m_reltol = eps;
    m_abstol = eps;

    m_steps = 0;
    m_acceptedSteps = 0;
    m_rejectedSteps = 0;
    m_totalEvals = 0;
    m_totalError = 0;

    m_facold = 0;

    m_stepSize = 1e-12;
    m_maxStepSize = maxStepSize;
}

template <typename T, unsigned int N>
ZD::CIntegratorRK45<T, N>::~CIntegratorRK45()
{

}

//template <typename T, unsigned int N>
//void CIntegratorRK45<T, N>::Init(const T eps)
//{
//    m_reltol = eps;
//    m_abstol = eps;
//}

//template <typename T, unsigned int N>
//void CIntegratorRK45<T, N>::Reset()
//{
//    m_steps = 0;
//    m_acceptedSteps = 0;
//    m_rejectedSteps = 0;
//    m_totalEvals = 0;
//    m_totalError = 0;
//
//    m_facold = 0;
//
//    m_stepSize = 1e-12;
//    m_maxStepSize = 1.0;
//}

template <typename T, unsigned int N>
T ZD::CIntegratorRK45<T, N>::Next(const CPoint<T, N>& p, const T& t,
    const T direction, const CDynamicSystem<T, N> *sys)
{
    CPoint<T, N> k2, k3, k4, k5, k6, k7;

    bool reject = false;

    m_stepSize = std::fabs(m_stepSize) * direction;

    // first step? initialize k1 and step size;
    if (m_steps == 0) {
        if (std::fabs(m_stepSize) < std::fabs(m_epsilon))
            m_stepSize = this->m_maxStepSize * direction;
        m_k1 = sys->Velocity(p, t);
        m_totalEvals = 1;
    }

    // integration step loop
    while (true) {
        CPoint<T, N> np;
        bool underflowFlag = false;

        // stepsize underflow?
        if (0.1 * std::fabs(m_stepSize) <= m_epsilon) {
            //std::cerr << "step size too small (h = " << params.stepSize << "), increase the tolerance." << std::endl;
            m_stepSize = 1e-12 * direction;
            underflowFlag = true;
        }

        if (std::fabs(m_stepSize) > std::fabs(this->m_maxStepSize)) {
            m_stepSize = std::fabs(this->m_maxStepSize) * direction;
        }

        m_steps++;

        np = p + m_stepSize * A21 * m_k1;
        k2 = sys->Velocity(np, t + C2 * m_stepSize);

        np = p + m_stepSize * (A31 * m_k1 + A32 * k2);
        k3 = sys->Velocity(np, t + C3 * m_stepSize);

        np = p + m_stepSize * (A41 * m_k1 + A42 * k2 + A43 * k3);
        k4 = sys->Velocity(np, t + C4 * m_stepSize);

        np = p + m_stepSize * (A51 * m_k1 + A52 * k2 + A53 * k3 + A54 * k4);
        k5 = sys->Velocity(np, t + C5 * m_stepSize);

        np = p + m_stepSize * (A61 * m_k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5);
        k6 = sys->Velocity(np, t + m_stepSize);

        np = p + m_stepSize * (A71 * m_k1 + A73 * k3 + A74 * k4 + A75 * k5 + A76 * k6);
        k7 = sys->Velocity(np, t + m_stepSize);

        m_totalEvals += 6;

        T err = 0.0, newStepSize, fac11;
        // error estimation
        CPoint<T, N> ee = m_stepSize * (E1 * m_k1 + E3 * k3 + E4 * k4 + E5 * k5 + E6 * k6 + E7 * k7);
        for (size_t i = 0; i < p.Size(); ++i) {
            T sk = m_abstol + m_reltol * std::max(std::fabs(p[i]), std::fabs(np[i]));
            T sqr = ee[i] / sk;
            err += sqr * sqr;
        }
        err = sqrt(err / (T)p.Size());

        // compute next (potential) stepsize
        fac11 = pow(err, 0.2 - m_beta * 0.75);

        // Lund-stabilization
        T fac = fac11 / pow(m_facold, m_beta);

        // we require facl <= h_new/h <= facr
        fac = std::max(T(1.0) / m_facr, std::min(T(1.0) / m_facl, fac / m_safe));
        newStepSize = m_stepSize / fac;

        if (err <= 1.0 || underflowFlag == true) {
            T successStepSize = m_stepSize;

            // step accepted
            m_facold = std::max(err, (T)1.0e-4);

            m_acceptedSteps++;
            for (size_t i = 0; i < p.Size(); ++i)
                m_totalError += std::fabs(ee[i]);

            if (reject)
                newStepSize = direction * std::min(std::fabs(newStepSize), std::fabs(m_stepSize));

            this->m_step.m_p[0] = p;
            this->m_step.m_p[1] = np - p;
            this->m_step.m_p[2] = m_stepSize * m_k1 - this->m_step.m_p[1];
            this->m_step.m_p[3] = -m_stepSize * k7 + this->m_step.m_p[1] - this->m_step.m_p[2];
            this->m_step.m_p[4] = m_stepSize * (D1 * m_k1 + D3 * k3 + D4 * k4 + D5 * k5 + D6 * k6 + D7 * k7);

            this->m_step.m_t0 = t;
            this->m_step.m_v0 = p;
            this->m_step.m_t1 = t + std::fabs(m_stepSize) * direction;
            this->m_step.m_v1 = np;

            //p = np;
            m_k1 = k7;
            m_stepSize = newStepSize;
            return std::fabs(successStepSize);
        }
        else {
            // step rejected
            newStepSize = m_stepSize / std::min(T(1.0) / m_facl, fac11 / m_safe);
            reject = true;

            if (m_acceptedSteps >= 1)
                m_rejectedSteps++;

            m_stepSize = newStepSize;
        }
    }
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CIntegratorRK45<T, N>::v(const T& t) const
{
    const T a = (t - this->m_step.m_t0) / (this->m_step.m_t1 - this->m_step.m_t0);
    const T b = 1.0 - a;

    return this->m_step.m_p[0] + a*(this->m_step.m_p[1] + b*(this->m_step.m_p[2] + a*(this->m_step.m_p[3] + b*this->m_step.m_p[4])));
}


#endif
