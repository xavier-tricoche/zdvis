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
#ifndef _ZD_LIB_FLOW_IO_GERRIS_HPP_
#define _ZD_LIB_FLOW_IO_GERRIS_HPP_

#include "utils/Base/ZD_Point.hpp"

#include <string>
#include <vector>
#include <assert.h>
#include <stdio.h>

#include "utils/define.hpp"
#include "utils/Tool/ZD_FileTool.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSet.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"

namespace ZD {
    template <typename T>
    class CFlowIOGerris {
    private:
        CPoint<T, 2> *m_pVertices;
        int m_vertexCount;

        CPoint<int, 4> *m_pCells;
        int m_cellCount;

        CPoint<T, 1> **m_pScalars;
        int m_scalarCount;
        std::vector<std::string> m_scalarNames;

        int m_scalarID;

    public:
        CFlowIOGerris();
        ~CFlowIOGerris();

    public:
        /* read output file as vtk unstructured grid */
        int ReadVTK(const char *pathname);

        inline int GetVertexCount() const;
        inline CPoint<T, 2> GetVertex(const int id) const;
        inline int GetCellCount() const;
        inline CPoint<int, 4> GetCell(const int id) const;
        inline int GetScalarCount() const;
        inline std::string GetScalarName(const int id) const;
        inline T GetScalar(const int id, const int k) const;

        const CPoint<T, 2>* GetVertices() const;
        const CPoint<int, 4>* GetCells() const;
        CPoint<T, 1>** GetScalars() const;

    private:
        std::vector<std::vector<int>> m_grids;

        int m_gridSize;
        T m_minX, m_minY, m_maxX, m_maxY;
        T m_dx, m_dy;

    public:
        void BuildGrid(const int size);
        const int GetCellID(const CPoint<T, 2>& p) const;

        CPoint<T, 2> GetVelocity(const CPoint<T, 2>& p) const;

        void SetScalarByName(const std::string& name);
        T GetScalar(const CPoint<T, 2>& p) const;

        void AddNoise();
    };
}

template<typename T>
ZD::CFlowIOGerris<T>::CFlowIOGerris()
{
    m_pVertices = nullptr;
    m_vertexCount = 0;

    m_pCells = nullptr;
    m_cellCount = 0;

    m_pScalars = nullptr;
    m_scalarCount = 0;

    m_scalarID = -1;
}

template<typename T>
ZD::CFlowIOGerris<T>::~CFlowIOGerris()
{
    SafeDeleteArray(m_pVertices);
    SafeDeleteArray(m_pCells);

    for (int i = 0; i < m_scalarCount; ++i)
        SafeDeleteArray(m_pScalars[i]);
    SafeDeleteArray(m_pScalars);
}


template<typename T>
int ZD::CFlowIOGerris<T>::ReadVTK(const char * pathname)
{
    printf("Reading VTK file %s ... \n", pathname);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(pathname);
    reader->Update();
    // reader->ReadAllScalarsOn();
    vtkSmartPointer<vtkUnstructuredGrid> data = reader->GetOutput();
    m_scalarCount = data->GetPointData()->GetNumberOfArrays();

    m_vertexCount = data->GetNumberOfPoints();
    m_cellCount = data->GetNumberOfCells();
    double x[3];
    m_pVertices = new CPoint<T, 2>[m_vertexCount];
    for (size_t i=0; i<m_vertexCount; ++i) {
        data->GetPoint(i, x);
        m_pVertices[i][0] = x[0];
        m_pVertices[i][1] = x[1];
    }

    m_pCells = new CPoint<int, 4>[m_cellCount];
    vtkSmartPointer<vtkIdList> vertIds = vtkSmartPointer<vtkIdList>::New();
    for (size_t i=0; i<m_cellCount; ++i) {
        data->GetCellPoints(i, vertIds);
        m_pCells[i][0] = vertIds->GetId(0);
        m_pCells[i][1] = vertIds->GetId(1);
        m_pCells[i][2] = vertIds->GetId(2);
        m_pCells[i][3] = vertIds->GetId(3);
    }

    m_pScalars = new CPoint<T, 1>*[m_scalarCount];
    for (size_t i=0; i<m_scalarCount; ++i) {
        m_scalarNames.push_back(data->GetPointData()->GetAbstractArray(i)->GetName());
        m_pScalars[i] = new CPoint<T, 1>[m_vertexCount];
        vtkSmartPointer<vtkDataArray> values =  vtkDataArray::SafeDownCast(data->GetPointData()->GetAbstractArray(i));
        for (size_t n=0; n<m_vertexCount; ++n) {
            m_pScalars[i][n] = values->GetTuple(n)[0];
        }
    }

    printf("done.\n");
    return 0;
}

template<typename T>
inline int ZD::CFlowIOGerris<T>::GetVertexCount() const
{
    return m_vertexCount;
}

template<typename T>
inline ZD::CPoint<T, 2> ZD::CFlowIOGerris<T>::GetVertex(const int id) const
{
    assert(m_pVertices != nullptr);
    return m_pVertices[id];
}

template<typename T>
inline int ZD::CFlowIOGerris<T>::GetCellCount() const
{
    return m_cellCount;
}

template<typename T>
inline ZD::CPoint<int, 4> ZD::CFlowIOGerris<T>::GetCell(const int id) const
{
    assert(m_pCells != nullptr);
    return m_pCells[id];
}

template<typename T>
inline int ZD::CFlowIOGerris<T>::GetScalarCount() const
{
    return m_scalarCount;
}

template<typename T>
inline std::string ZD::CFlowIOGerris<T>::GetScalarName(const int id) const
{
    return m_scalarNames[id];
}

template<typename T>
inline T ZD::CFlowIOGerris<T>::GetScalar(const int id, const int k) const
{
    assert(m_pScalars != nullptr);
    return m_pScalars[k][id][0];
}

template<typename T>
const ZD::CPoint<T, 2>* ZD::CFlowIOGerris<T>::GetVertices() const
{
    return m_pVertices;
}

template<typename T>
const ZD::CPoint<int, 4>* ZD::CFlowIOGerris<T>::GetCells() const
{
    return m_pCells;
}

template<typename T>
ZD::CPoint<T, 1>** ZD::CFlowIOGerris<T>::GetScalars() const
{
    return m_pScalars;
}


template<typename T>
void ZD::CFlowIOGerris<T>::BuildGrid(const int size)
{
    // min & max
    m_minX = m_maxX = m_pVertices[0][0];
    m_minY = m_maxY = m_pVertices[0][1];
    for (int i = 0; i < m_vertexCount; ++i) {
        if (m_pVertices[i][0] < m_minX)
            m_minX = m_pVertices[i][0];
        if (m_pVertices[i][0] > m_maxX)
            m_maxX = m_pVertices[i][0];
        if (m_pVertices[i][1] < m_minY)
            m_minY = m_pVertices[i][1];
        if (m_pVertices[i][1] > m_maxY)
            m_maxY = m_pVertices[i][1];
    }

    // build grid
    m_gridSize = size;
    m_grids.resize((m_gridSize + 1) * (m_gridSize + 1));
    m_dx = (m_maxX - m_minX) / T(m_gridSize);
    m_dy = (m_maxY - m_minY) / T(m_gridSize);

    for (int i = 0; i < m_cellCount; ++i) {
        int v0 = m_pCells[i][0];
        int v1 = m_pCells[i][1];
        int v2 = m_pCells[i][2];
        int v3 = m_pCells[i][3];

        T min_x = m_pVertices[v0][0];
        T max_x = m_pVertices[v1][0];
        T min_y = m_pVertices[v0][1];
        T max_y = m_pVertices[v2][1];

        int min_x_id = int(std::floor((min_x - m_minX) / m_dx));
        int max_x_id =  int(std::ceil((max_x - m_minX) / m_dx));
        int min_y_id = int(std::floor((min_y - m_minY) / m_dy));
        int max_y_id =  int(std::ceil((max_y - m_minY) / m_dy));
        for (int y = min_y_id; y <= max_y_id; ++y) {
            for (int x = min_x_id; x <= max_x_id; ++x) {
                m_grids[y*(m_gridSize + 1) + x].push_back(i);
            }
        }
    }
}

template<typename T>
const int ZD::CFlowIOGerris<T>::GetCellID(const CPoint<T, 2>& p) const
{
    int x_id = int(std::floor((p[0] - m_minX) / m_dx));
    int y_id = int(std::floor((p[1] - m_minY) / m_dy));

    if (x_id < 0 || x_id > m_gridSize ||
        y_id < 0 || y_id > m_gridSize)
        return -1;

    int grid_id = y_id * (m_gridSize + 1) + x_id;
    int cell_id = -1;
    int v0, v1, v2, v3;
    T min_x, max_x, min_y, max_y;
    for (int i = 0; i < m_grids[grid_id].size(); ++i) {
        int tmp = m_grids[grid_id][i];
        v0 = m_pCells[tmp][0];
        v1 = m_pCells[tmp][1];
        v2 = m_pCells[tmp][2];
        v3 = m_pCells[tmp][3];

        min_x = m_pVertices[v0][0];
        max_x = m_pVertices[v1][0];
        min_y = m_pVertices[v0][1];
        max_y = m_pVertices[v2][1];

        if (p[0] >= min_x && p[0] <= max_x && p[1] >= min_y && p[1] <= max_y) {
            cell_id = tmp;
            break;
        }
    }

    return cell_id;
}

template<typename T>
ZD::CPoint<T, 2> ZD::CFlowIOGerris<T>::GetVelocity(const CPoint<T, 2>& p) const
{
    int cellID = GetCellID(p);
    if (cellID == -1) {
        return CPoint<T, 2>(0.0, 0.0);
    } else {
        int i0 = m_pCells[cellID][0];
        int i1 = m_pCells[cellID][1];
        int i2 = m_pCells[cellID][2];
        int i3 = m_pCells[cellID][3];

        T min_x = m_pVertices[i0][0];
        T max_x = m_pVertices[i1][0];
        T min_y = m_pVertices[i0][1];
        T max_y = m_pVertices[i2][1];

        T x_factor = (p[0] - min_x) / (max_x - min_x);
        T y_factor = (p[1] - min_y) / (max_y - min_y);

        T u0 = m_pScalars[2][i0][0];
        T u1 = m_pScalars[2][i1][0];
        T u2 = m_pScalars[2][i2][0];
        T u3 = m_pScalars[2][i3][0];
        T u01 = u0 * (1.0 - x_factor) + u1 * x_factor;
        T u23 = u2 * (1.0 - x_factor) + u3 * x_factor;
        T uu = u01 * (1.0 - y_factor) + u23 * y_factor;

        T v0 = m_pScalars[3][i0][0];
        T v1 = m_pScalars[3][i1][0];
        T v2 = m_pScalars[3][i2][0];
        T v3 = m_pScalars[3][i3][0];
        T v01 = v0 * (1.0 - x_factor) + v1 * x_factor;
        T v23 = v2 * (1.0 - x_factor) + v3 * x_factor;
        T vv = v01 * (1.0 - y_factor) + v23 * y_factor;

        return CPoint<T, 2>(uu, vv);
    }
}

template<typename T>
void ZD::CFlowIOGerris<T>::SetScalarByName(const std::string& name)
{
    for (int i = 0; i < m_scalarNames.size(); ++i) {
        if (name == m_scalarNames[i])
            m_scalarID = i;
    }
}

template<typename T>
T ZD::CFlowIOGerris<T>::GetScalar(const CPoint<T, 2>& p) const
{
    assert(m_scalarID != -1);

    int cellID = GetCellID(p);
    if (cellID == -1) {
        return 0.0;
    }
    else {
        int i0 = m_pCells[cellID][0];
        int i1 = m_pCells[cellID][1];
        int i2 = m_pCells[cellID][2];
        int i3 = m_pCells[cellID][3];

        T min_x = m_pVertices[i0][0];
        T max_x = m_pVertices[i1][0];
        T min_y = m_pVertices[i0][1];
        T max_y = m_pVertices[i2][1];

        T x_factor = (p[0] - min_x) / (max_x - min_x);
        T y_factor = (p[1] - min_y) / (max_y - min_y);

        T s0 = m_pScalars[m_scalarID][i0][0];
        T s1 = m_pScalars[m_scalarID][i1][0];
        T s2 = m_pScalars[m_scalarID][i2][0];
        T s3 = m_pScalars[m_scalarID][i3][0];
        T s01 = s0 * (1.0 - x_factor) + s1 * x_factor;
        T s23 = s2 * (1.0 - x_factor) + s3 * x_factor;
        T ss = s01 * (1.0 - y_factor) + s23 * y_factor;

        return ss;
    }
}

template <typename T>
void ZD::CFlowIOGerris<T>::AddNoise()
{
    for (int i = 0; i < m_scalarCount; ++i) {
        //T minS = m_pScalars[i][0][0];
        //T maxS = m_pScalars[i][0][0];
        //for (int k = 0; k < m_vertexCount; ++k) {
        //    if (m_pScalars[i][k][0] > maxS)
        //        maxS = m_pScalars[i][k][0];
        //    if (m_pScalars[i][k][0] < minS)
        //        minS = m_pScalars[i][k][0];
        //}

        CRandomTool<T> rt;
        for (int k = 0; k < m_vertexCount; ++k) {
            m_pScalars[i][k][0] *= (1.0 + (rt.GetRandomNumer() - 0.5) * 0.1);
        }
    }
}

#endif
