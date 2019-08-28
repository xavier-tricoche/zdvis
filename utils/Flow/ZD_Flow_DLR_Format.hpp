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
#ifndef _ZD_LIB_FLOW_DLR_FORMAT_HPP_
#define _ZD_LIB_FLOW_DLR_FORMAT_HPP_

#include "ZD_Flow.hpp"

// nvis
#include "math/bounding_box.hpp"
#include "math/fixed_vector.hpp"

// celltree
#include "celltree/celltree.hpp"
#include "celltree/celltree_builder.hpp"
#include "celltree/dataset.hpp"
#include "celltree/interpolator.hpp"
#include "celltree/mesh.hpp"

#include "utils/Tool/ZD_TimeTool.hpp"
#include "utils/Integrator/ZD_Integrator_RK45.hpp"
#include "utils/configure.hpp"
#include "utils/define.hpp"

#include <assert.h>
#include <exception>


namespace ZD {
template< typename T>
class CFlowDLRFormat : public CFlow<T, 3> {
public:
    typedef T value_type;
    typedef CPoint<T, 3> point_type;
private:
    
    typedef std::map<value_type, int> map_type;
    typedef typename map_type::const_iterator map_iterator_type;
    
    std::shared_ptr<CellTree::dataset> m_dataset;
    std::shared_ptr<CellTree::mesh> m_mesh;
    CellTree::celltree m_celltree;
    
    std::vector< std::shared_ptr<celltree::variable> > m_vectors;
    std::vector< std::shared_ptr<celltree::variable> > m_scalars;
    std::vector< std::shared_ptr<CellTree::interpolator> > m_interpolators;
    std::map<value_type, int> m_times_to_id;
    value_type m_boundingBox[6];
    value_type m_minTime;
    value_type m_maxTime;

    int m_testID;

public:
    CFlowDLRFormat(const std::string& filename, 
                   const value_type minTime, const value_type maxTime, 
                   const std::string& scalarName="");
    ~CFlowDLRFormat();

protected:
    virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

public:
    virtual inline const point_type Velocity(const point_type& p, const value_type& _time) const;
    virtual inline const value_type Scalar(const point_type& p, const value_type& _time) const;
    virtual inline const bool CheckPosition(const point_type& p) const;
    virtual inline void GetBBox(point_type& min, point_type& max) const;

};

template <typename T>
CFlowDLRFormat<T>::CFlowDLRFormat(const std::string& filename, const value_type minTime, 
                                  const value_type maxTime, const std::string& scalarName)
{
    m_testID = -1;
        
    m_dataset.reset(dataset::create(filename));
    size_t nb_ts = m_dataset->get_num_timesteps();
    for (size_t i=0; i<nb_ts; ++i) {
        m_times_to_id[m_dataset->get_timestep(i)]= i;
    }
    m_vectors.resize(nb_ts);
    m_interpolators.resize(nb_ts);  
    
    m_mesh.reset(m_dataset->read_mesh());
    celltree:::celltree_builder ct_builder;
    ct_builder.m_buckets = 5;
    ct_builder.m_leafsize = 8;
    ct_builder.build(m_celltree, *m_mesh);

    CellTree::interpolator::coord_type CellTree::mesh_min[3];
    CellTree::interpolator::coord_type CellTree::mesh_max[3];
    CellTree::mesh_traits<CellTree::mesh>::extents(*m_mesh, CellTree::mesh_min, CellTree::mesh_max);     
    m_boundingBox[0] = CellTree::mesh_min[0];
    m_boundingBox[1] = CellTree::mesh_max[0];
    m_boundingBox[2] = CellTree::mesh_min[1];
    m_boundingBox[3] = CellTree::mesh_max[1];
    m_boundingBox[4] = CellTree::mesh_min[2];
    m_boundingBox[5] = CellTree::mesh_max[2];
    
    for (size_t i=0; i<nb_ts; ++i) {
        m_vectors[i].reset(m_dataset->read_vector_variable(i, "velocity"));
        m_interpolators[i].reset(new CellTree::interpolator(const_cast<const CellTree::mesh*>(m_mesh.get()), 
        const_cast<const celltree::variable*>(m_vectors[i].get()), m_celltree));
    }
    
    if (!scalarName.empty()) {
        m_scalars.resize(nb_ts);
        for (size_t i=0; i<nb_ts; ++i) {
            m_scalars[i].reset(m_dataset->read_scalar_variable(i, scalarName));
        }
    }
}

template <typename T>
CFlowDLRFormat<T>::~CFlowDLRFormat()
{
}

template <typename T>
inline const bool CFlowDLRFormat<T>::CheckPosition(const point_type& p) const
{
    if (p[0] > m_boundingBox[1] || p[0] < m_boundingBox[0] ||
        p[1] > m_boundingBox[3] || p[1] < m_boundingBox[2] ||
        p[2] > m_boundingBox[5] || p[2] < m_boundingBox[4]) return false;
    
    point_type v;
    return m_interpolators[0]->operator()(0, p.m_data, v.m_data);
}

template <typename T>
inline CIntegrator<T, 3> * CFlowDLRFormat<T>::CreateIntegrator() const
{
    return new CIntegratorRK45<T, 3>(1e-8, this->m_stepSize);
}
// cubic 
template <typename T>
inline const CPoint<T, 3> CFlowDLRFormat<T>::Velocity(const point_type& p, const value_type& _time) const
{
    map_iterator_type it_next = m_times_to_id.lower_bound(_time);
    
    if ( (it_next == m_times_to_id.begin() && it_next->first > _time) ||
         (it_next == m_times_to_id.end()) )
        throw std::runtime_error("invalid time coordinate");
    
    map_iterator_type it_prev = it_next;
    --it_prev;
    // requested time lies between time coordinate of it_prev and it_next
    double u = (_time-it_prev->first)/(it_next->first - it_prev->first);
    
    point_type vprev, vnext;
    bool inside = m_interpolators[it_prev->second]->operator()(0, p.m_data, vprev.m_data);
    if (!inside) throw std::runtime_error("invalid position");
    
    inside = m_interpolators[it_prev->second]->operator()(0, p.m_data, vnext.m_data);
    return (1.0-u)*vprev + u*vnext;
}

// cubic
template <typename T>
inline const T CFlowDLRFormat<T>::Scalar(const point_type& p, const value_type& _time) const
{
    map_iterator_type it_next = m_times_to_id.lower_bound(_time);
    
    if ( (it_next == m_times_to_id.begin() && it_next->first > _time) ||
         (it_next == m_times_to_id.end()) )
        throw std::runtime_error("invalid position");
    
    map_iterator_type it_prev = it_next;
    --it_prev;
    // requested time lies between time coordinate of it_prev and it_next
    double u = (_time-it_prev->first)/(it_next->first - it_prev->first);
    
    value_type fprev, fnext;
    bool inside = m_interpolators[it_prev->second]->operator()(0, p.m_data, &fprev);
    if (!inside) throw std::runtime_error("invalid position");
    
    inside = m_interpolators[it_prev->second]->operator()(0, p.m_data, &fnext);
    return (1.0-u)*fprev + u*fnext;
}

template <typename T>
inline void CFlowDLRFormat<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    min[2] = m_boundingBox[4];

    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];
    max[2] = m_boundingBox[5];
}

} // namespace ZD

#endif
