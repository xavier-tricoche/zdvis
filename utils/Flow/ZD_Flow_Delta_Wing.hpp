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
#ifndef _ZD_LIB_FLOW_DELTA_WING_HPP_
#define _ZD_LIB_FLOW_DELTA_WING_HPP_

#include "ZD_Flow.hpp"

#include "math/bounding_box.hpp"
#include "math/fixed_vector.hpp"

#include "celltree/celltree.hpp"
#include "celltree/celltree_builder.hpp"
#include "celltree/dataset.hpp"
#include "celltree/datasetDLR.hpp"
#include "celltree/interpolator.hpp"
#include "celltree/mesh.hpp"
#include "celltree/mesh_traits.hpp"

#include "utils/Tool/ZD_TimeTool.hpp"
#include "utils/Integrator/ZD_Integrator_RK45.hpp"
#include "utils/configure.hpp"
#include "utils/define.hpp"

#include <assert.h>


namespace ZD {
    template< typename T>
    class CFlowDeltaWing : public CFlow<T, 3> {
    public:
        typedef T value_type;
        typedef CPoint<T, 3> point_type;
        typedef CellTree::mesh_traits<CellTree::mesh> mtraits;
    private:
        struct DeltaWingDataset {
            CellTree::dataset *ds;
            CellTree::celltree ct;
            CellTree::interpolator *intp;
        };

    private:
        CellTree::dataset* m_dataset;
        CellTree::mesh* m_mesh;
        CellTree::celltree m_celltree;
        std::vector<CellTree::interpolator*> m_vector_fields, m_scalar_fields;

        std::vector<value_type> m_pTimes;

        value_type m_boundingBox[6];
        value_type m_minTime;
        value_type m_maxTime;

        int m_testID;

    public:
        CFlowDeltaWing(const value_type minTime, const value_type maxTime, const bool readVelocity, const char *scalarName = nullptr, const std::string& celltree_name = std::string());
        ~CFlowDeltaWing();

    protected:
        virtual inline CIntegrator<T, 3> *CreateIntegrator() const;

    public:
        virtual inline const point_type Velocity(const point_type& p, const value_type& time) const;
        virtual inline const value_type Scalar(const point_type& p, const value_type& time) const;
        virtual inline const bool CheckPosition(const point_type& p) const;
        virtual inline void GetBBox(point_type& min, point_type& max) const;

    private:
        bool InitTimes();
        bool ReadVelocity(const value_type minTime, const value_type maxTime);
        bool ReadScalar(const value_type minTime, const value_type maxTime, const char *name);

        inline void TimeToIndex(const value_type time, int *indices) const;
    };
}

template <typename T>
ZD::CFlowDeltaWing<T>::CFlowDeltaWing(const value_type minTime, const value_type maxTime, const bool readVelocity, const char *scalarName,const std::string& celltree_name)
{
    m_testID = -1;

    InitTimes();

    m_dataset = new CellTree::datasetDLR(DELTA_WING_DLR_FILE);

    m_minTime = m_pTimes[0];
    m_maxTime = m_pTimes.back();

    m_boundingBox[0] = -20.5;
    m_boundingBox[1] = 20.5;
    m_boundingBox[2] = -10.0;
    m_boundingBox[3] = 10.0;
    m_boundingBox[4] = -10.0;
    m_boundingBox[5] = 10.0;


    m_mesh = m_dataset->read_mesh();

    bool read = false;

    if (!celltree_name.empty()) {
        std::cout << "Checking if provided celltree filename (" << celltree_name << ") exists... ";
        read = m_celltree.read(celltree_name);
        if (!read) std::cout << "nope\n";
        else std::cout << "yes: celltree successfully imported\n";
    }

    if (!read) {
        CellTree::celltree_builder builder;
        builder.m_buckets = 5;
        builder.m_leafsize = 8;
        builder.build(m_celltree, *m_mesh);
        if (!celltree_name.empty()) {
            std::cout << "Attempting to export celltree to " << celltree_name << "... ";
            bool written = m_celltree.write(celltree_name);
            if (!written) std::cout << "failed\n";
            else std::cout << "done\n";
        }
    }

    if (readVelocity == true)
        ReadVelocity(minTime, maxTime);

    if (scalarName != nullptr)
        ReadScalar(minTime, maxTime, scalarName);
}

template <typename T>
ZD::CFlowDeltaWing<T>::~CFlowDeltaWing()
{
    for (size_t i=0; i<m_vector_fields.size(); ++i) {
        SafeDelete(m_vector_fields[i]);
    }
    for (size_t i=0; i<m_scalar_fields.size(); ++i) {
        SafeDelete(m_scalar_fields[i]);
    }
    SafeDelete(m_mesh);
    SafeDelete(m_dataset);
}

template <typename T>
inline const bool ZD::CFlowDeltaWing<T>::CheckPosition(const point_type& p) const
{
    if (p[0]>m_boundingBox[0] && p[0] < m_boundingBox[1] &&
        p[1]>m_boundingBox[2] && p[1] < m_boundingBox[3] &&
        p[2]>m_boundingBox[4] && p[2] < m_boundingBox[5]) {
        nvis::vec3 pp;
        pp[0] = p[0];
        pp[1] = p[1];
        pp[2] = p[2];
        double dir[3] = { 0.0, 0.0, 0.0 };
        mtraits::size_type cell;
        double coord[3];
        if (m_vector_fields.size()) {
	        // std::cout << "Checkposition called for (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
            bool ok = m_vector_fields[0]->find(pp.begin(), cell, coord);
	        // if (ok)
                // std::cout << "search successful: local coordinates are (" << coord[0] << "," << coord[1] << "," << coord[2] << ") in cell #" << cell << '\n';
	        // else std::cout << "check failed\n";
	        return ok;
	    }
        else if (m_scalar_fields.size())
            return m_scalar_fields[0]->find(pp.begin(), cell, coord);
        else
            return false;
    }
    else {
        return false;
    }
}

template <typename T>
inline ZD::CIntegrator<T, 3>*ZD::CFlowDeltaWing<T>::CreateIntegrator() const
{
    return new CIntegratorRK45<T, 3>(1e-8, this->m_stepSize);
}
// cubic
template <typename T>
inline const ZD::CPoint<T, 3> ZD::CFlowDeltaWing<T>::Velocity(const point_type& p, const value_type& t) const
{
    if (m_pTimes.empty() || t < m_pTimes[0] || t>m_pTimes.back())
        return point_type(0,0,0);

    int id[2];
    TimeToIndex(t, id);

    nvis::vec3 pp;
    pp[0] = p[0];
    pp[1] = p[1];
    pp[2] = p[2];

    double dir_0[3] = { 0.0, 0.0, 0.0 };
    double dir_1[3] = { 0.0, 0.0, 0.0 };
    double dir_2[3] = { 0.0, 0.0, 0.0 };
    double dir_3[3] = { 0.0, 0.0, 0.0 };

    double local[3];
    mtraits::size_type cell;
    // std::cout << "Velocity query at (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
    if (m_vector_fields[id[0]]->find(pp.begin(), cell, local)) {
	// std::cout << "point location successful\n";
        if (id[0]-1 >= 0)
            m_vector_fields[id[0]-1]->evaluate(cell, local, dir_0);
        m_vector_fields[id[0]]->evaluate(cell, local, dir_1);
        m_vector_fields[id[1]]->evaluate(cell, local, dir_2);
        if (id[1]+1 < m_vector_fields.size())
            m_vector_fields[id[1]+1]->evaluate(cell, local, dir_3);
    }
    else {
        // std::cout << "point location failed\n";
        return point_type(0,0,0);
    }

    // Centripetal Catmull-Rom spline interpolation
    const point_type v0 = point_type(dir_0[0], dir_0[1], dir_0[2]);
    const point_type v1 = point_type(dir_1[0], dir_1[1], dir_1[2]);
    const point_type v2 = point_type(dir_2[0], dir_2[1], dir_2[2]);
    const point_type v3 = point_type(dir_3[0], dir_3[1], dir_3[2]);

    const value_type t0 = m_pTimes[std::max(id[0]-1, (int)0)];
    const value_type t1 = m_pTimes[id[0]];
    const value_type t2 = m_pTimes[id[1]];
    const value_type t3 = m_pTimes[std::min(id[1]+1, (int)(m_pTimes.size()-1))];

    const point_type a1 = (t1>t0) ? v0*(t1-t) / (t1-t0) + v1*(t-t0) / (t1-t0) : point_type(0,0,0);
    const point_type a2 = v1*(t2-t) / (t2-t1) + v2*(t-t1) / (t2-t1);
    const point_type a3 = (t3>t2) ? v2*(t3-t) / (t3-t2) + v3*(t-t2) / (t3-t2) : point_type(0,0,0);

    const point_type b1 = a1*(t2-t) / (t2-t0) + a2*(t-t0) / (t2-t0);
    const point_type b2 = a2*(t3-t) / (t3-t1) + a3*(t-t1) / (t3-t1);

    const point_type dir = b1*(t2-t) / (t2-t1) + b2*(t-t1) / (t2-t1);

    return dir;
}

// cubic
template <typename T>
inline const T ZD::CFlowDeltaWing<T>::Scalar(const point_type& p, const value_type& t) const
{
    if (m_pTimes.empty() || t < m_pTimes[0] || t>m_pTimes.back())
        return T(0);

    int id[2];
    TimeToIndex(t, id);

    nvis::vec3 pp;
    pp[0] = p[0];
    pp[1] = p[1];
    pp[2] = p[2];

    double s_0[1] = { 0.0 };
    double s_1[1] = { 0.0 };
    double s_2[1] = { 0.0 };
    double s_3[1] = { 0.0 };

    double local[3];
    mtraits::size_type cell;
    if (m_scalar_fields[id[0]]->find(pp.begin(), cell, local)) {
        if (id[0]-1 >= 0)
            m_scalar_fields[id[0]-1]->evaluate(cell, local, s_0);
        m_scalar_fields[id[0]]->evaluate(cell, local, s_1);
        m_scalar_fields[id[1]]->evaluate(cell, local, s_2);
        if (id[1]+1 < m_scalar_fields.size())
            m_scalar_fields[id[1]+1]->evaluate(cell, local, s_3);
    }

    const value_type t0 = m_pTimes[std::max(id[0]-1, (int)0)];
    const value_type t1 = m_pTimes[id[0]];
    const value_type t2 = m_pTimes[id[1]];
    const value_type t3 = m_pTimes[std::min(id[1] + 1, (int)(m_pTimes.size()-1))];

    const value_type a1 = (t1>t0) ? s_0[0]*(t1-t) / (t1-t0) + s_1[0]*(t-t0) / (t1-t0) : T(0);
    const value_type a2 = s_1[0]*(t2-t) / (t2-t1) + s_2[0]*(t-t1) / (t2-t1);
    const value_type a3 = (t3>t2) ? s_2[0]*(t3-t) / (t3-t2) + s_3[0]*(t-t2) / (t3-t2) : T(0);

    const value_type b1 = a1*(t2-t) / (t2-t0) + a2*(t-t0) / (t2-t0);
    const value_type b2 = a2*(t3-t) / (t3-t1) + a3*(t-t1) / (t3-t1);

    const value_type s = b1*(t2-t) / (t2-t1) + b2*(t-t1) / (t2-t1);

    return s;
}

template <typename T>
bool ZD::CFlowDeltaWing<T>::InitTimes()
{
    assert(m_pTimes.empty());
    m_pTimes.resize(DELTA_WING_TOTAL_COUNT);

    m_pTimes[0] = 0.01;
    for (int i = 1; i < DELTA_WING_TOTAL_COUNT; ++i) {
        m_pTimes[i] = (i + 1)*0.01 + 0.001;
    }

    return true;
}

template <typename T>
bool ZD::CFlowDeltaWing<T>::ReadVelocity(const value_type minTime, const value_type maxTime)
{
    assert(m_vector_fields.empty());
    m_vector_fields.resize(DELTA_WING_TOTAL_COUNT);
    for (int i = 0; i < DELTA_WING_TOTAL_COUNT; ++i) {
        m_vector_fields[i]= nullptr;
    }

    int temp[2];
    TimeToIndex(minTime, temp);
    int minID = temp[0];
    TimeToIndex(maxTime, temp);
    int maxID = temp[1];

    for (int i = minID; i <= maxID; ++i) {
        CTimeTool<T> timer;
        timer.StartTimer();
        std::cout << "Read vector field on time = " << m_pTimes[i] << std::endl;

        const CellTree::variable* v = m_dataset->read_vector_variable(i, "velocity");
        m_vector_fields[i] = new CellTree::interpolator(m_mesh, v, m_celltree);
    }

    m_testID = minID;

    return true;
}

template <typename T>
bool ZD::CFlowDeltaWing<T>::ReadScalar(const value_type minTime, const value_type maxTime, const char *name)
{
    assert(m_scalar_fields.empty());
    m_scalar_fields.resize(DELTA_WING_TOTAL_COUNT);
    for (int i = 0; i < DELTA_WING_TOTAL_COUNT; ++i) {
        m_scalar_fields[i] = nullptr;
    }

    int temp[2];
    TimeToIndex(minTime, temp);
    int minID = temp[0];
    TimeToIndex(maxTime, temp);
    int maxID = temp[1];

    for (int i = minID-1; i <= maxID+1; ++i) {
        CTimeTool<T> timer;
        timer.StartTimer();
        std::cout << "Read Scalar field " << name << " on time = " << m_pTimes[i] << std::endl;
        const CellTree::variable* v = m_dataset->read_scalar_variable(i, name);
        m_scalar_fields[i] = new CellTree::interpolator(m_mesh, v, m_celltree);
    }

    m_testID = minID;

    return true;
}

template <typename T>
inline void ZD::CFlowDeltaWing<T>::TimeToIndex(const value_type t, int *indices) const
{
    indices[0] = 0;
    indices[1] = m_pTimes.size()-1;
    if (t <= this->m_pTimes[0]) {
        indices[1] = 1;
    }
    else if (t >= this->m_pTimes[indices[1]]) {
        indices[0] = indices[1]-1;
    }
    else {
        // binary search for interval containing requested time
        while ((indices[0] + 1) < indices[1]) {
            int mid = (indices[0] + indices[1]) / 2.0;
            if (this->m_pTimes[mid]>t) {
                indices[1] = mid;
            }
            else {
                indices[0] = mid;
            }
        }
    }
}

template <typename T>
inline void ZD::CFlowDeltaWing<T>::GetBBox(point_type& min, point_type& max) const
{
    min[0] = m_boundingBox[0];
    min[1] = m_boundingBox[2];
    min[2] = m_boundingBox[4];

    max[0] = m_boundingBox[1];
    max[1] = m_boundingBox[3];
    max[2] = m_boundingBox[5];
}

#endif
