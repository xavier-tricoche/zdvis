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
#ifndef _ZD_LIB_INTEGRATOR_RKDP45_HPP_
#define _ZD_LIB_INTEGRATOR_RKDP45_HPP_

#include "ZD_Integrator.hpp"
#include <math/dopri5.hpp>
#include <algorithm>

namespace ZD {
    template <typename T, unsigned int N>
    class CIntegratorRKDP45 : public CIntegrator<T, N> {
    public:
        typedef T scalar_type;
        typedef nvis::fixed_vector<T, N+2>     state_type;
        typedef nvis::fixed_vector<T, N>       point_type;
        typedef nvis::fixed_vector<T, N>       vector_type;
        typedef nvis::dopri5<state_type>       integrator_type;
        typedef typename integrator_type::step step_type;
        typedef ZD::CDynamicSystem<T, N>       system_type;
        typedef ZD::CPoint<T, N>               c_point_type;
        
        enum exit_state {
           OK = 0,
           LEFT_DOMAIN,
           ERROR,
           OTHER, 
        };
        
    private:
        scalar_type m_abs_tol, m_rel_tol;
        scalar_type m_step_size, m_current_time;
        std::vector<step_type> m_all_steps;
        size_t m_total_evals;
        integrator_type m_intg;
        scalar_type m_dir;
        
        struct rhs_type {
            const system_type* m_sys;
            mutable size_t m_n_evals;
            scalar_type m_direction;
            rhs_type(const system_type* sys, scalar_type dir) 
                : m_sys(sys), m_n_evals(0), m_direction(dir)
            {}
            
            state_type operator()(const double& t, const state_type& x) const {
                vector_type r = m_direction* *(vector_type*)&m_sys->Velocity(*(c_point_type*)&x[0], t)[0];
                state_type v;
                std::copy(&r[0], &r[N], &v[0]);
                v[N] = nvis::norm(r);
                v[N+1] = 0;
                ++m_n_evals;
                return v; 
            }
            
            size_t n_evals() const { return m_n_evals; }
        };
    public:
        CIntegratorRKDP45() 
            : m_intg(), m_abs_tol(1.0e-7), m_rel_tol(1.0e-7),
              m_all_steps(), m_total_evals(0),
              m_step_size(1.0e-12) {
                  m_intg.abstol = m_abs_tol;
                  m_intg.reltol = m_rel_tol;
                  m_intg.h = 0;
              }
        CIntegratorRKDP45(const T eps, const T maxStepSize) 
            : m_intg(), m_abs_tol(eps), m_rel_tol(eps), 
              m_step_size(maxStepSize), m_all_steps(), m_total_evals(0) {
                  m_intg.abstol = m_abs_tol;
                  m_intg.reltol = m_rel_tol;
                  m_intg.h = 0;
              }
        ~CIntegratorRKDP45() {}

        T   TotalError() const { return 0; } // not calculated
        int TotalEvals() const { return m_total_evals; }

        //void Init(const T eps);

        virtual T Next(const c_point_type& p, const scalar_type& t, 
                       const scalar_type direction,
                       const system_type *sys);

        virtual c_point_type v(const T& t) const;


        //virtual void Reset();
    };
}

template <typename T, unsigned int N>
T ZD::CIntegratorRKDP45<T, N>::Next(const c_point_type& p, const scalar_type& t,
    const scalar_type direction, const system_type *sys)
{
    if (m_all_steps.empty()) {
        m_total_evals = 0;
        std::copy(&p[0], &p[N], &m_intg.y[0]);
        m_intg.y[N] = 0;
        m_intg.y[N+1] = 0;
        m_intg.t = t;
        m_all_steps.clear();
    }
    m_dir = direction;
    
    rhs_type rhs(sys, direction);
    step_type s;
    typename integrator_type::result res = m_intg.do_step(rhs, s);
    m_all_steps.push_back(s);
    m_total_evals += rhs.m_n_evals;
    return s.t1() - s.t0();
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CIntegratorRKDP45<T, N>::v(const T& t) const
{
    return *(ZD::CPoint<T, N>*)&m_all_steps.back().y(t)[0];
}


#endif
