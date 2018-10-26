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
#ifndef _ZD_LIB_MESH_HPP_
#define _ZD_LIB_MESH_HPP_

#include "ZD_Point.hpp"

#include "../define.hpp"

#include <algorithm>

#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>

#include <cmath>
#include <cfloat>

#include <vector>

namespace ZD {
    namespace MeshMisc {
        template <typename T>
        void ComputeAABB(const std::vector<CPoint<T, 3>>& vertices, CPoint<T, 3> *aabb);
        template <typename T>
        int FindVertex(const std::vector<CPoint<T, 3>>& vertices, const CPoint<T, 3>& v);
        template <typename T>
        int FindVertex(const std::vector<CPoint<T, 3>>& vertices, std::vector<int>& grid, const CPoint<T, 3>& v);
        template <typename T>
        int ComputeGridID(const CPoint<T, 3> *aabb, const CPoint<T, 3>& dis, const int grid_size, const CPoint<T, 3>& v);
    }

    template <typename T>
    class CMesh {
    public:
        CPoint<T, 3>   *m_pVertices;
        CPoint<int, 3> *m_pFaces;
        int m_vNum;
        int m_fNum;

        CPoint<T, 3> *m_pColors;            // colors;

    private:
        CPoint<T, 3> m_AABB[2];        /* bounding box */

    public:
        CMesh();
        CMesh(const std::vector<CPoint<T, 3>>& triangles);
        CMesh(const std::vector<CPoint<int, 3>>& triangles, const std::vector<CPoint<T, 3>>& vertices);
        ~CMesh();

    public:
        int OpenMesh(const char *pathname, FILE *fp = stdout);
        int SaveMesh(const char *pathname, FILE *fp = stdout) const;
        int OpenMeshVTK(const char *pathname, FILE *fp = stdout);
        int SaveMeshVTK(const char *pathname, FILE *fp = stdout) const;

    public:
        void AABB(CPoint<T, 3> *aabb);
        void Normalize();
        void Scale(const CPoint<T, 3>& scale);
        void Translate(const CPoint<T, 3>& translate);
        void SetColors(const std::vector<CPoint<T, 3>>& colors);

    public:
        void AddMesh(const CMesh *mesh);
        void AddTriangles(const std::vector<CPoint<T, 3>>& triangles);

    private:
        void CreateMesh(const std::vector<CPoint<T, 3>>& triangles);
        void CreateMesh(const std::vector<CPoint<int, 3>>& triangles, const std::vector<CPoint<T, 3>>& vertices);
        void ComputeAABB();
    };
} // namespace ZD

template<typename T>
int ZD::MeshMisc::FindVertex(const std::vector<CPoint<T, 3>>& vertices, const CPoint<T, 3>& v)
{
    for (int i = 0; i < vertices.size(); ++i) {
        if ((vertices[i] - v).Length() < ZD_EPSILON)
            return i;
    }
    return -1;
}

template<typename T>
int ZD::MeshMisc::FindVertex(const std::vector<CPoint<T, 3>>& vertices, std::vector<int>& grid, const CPoint<T, 3>& v)
{
    for (int i = 0; i < grid.size(); ++i) {
        if ((vertices[grid[i]] - v).Length() < ZD_EPSILON)
            return grid[i];
    }
    return -1;
}


template<typename T>
void ZD::MeshMisc::ComputeAABB(const std::vector<CPoint<T, 3>>& vertices, CPoint<T, 3>* aabb)
{
    aabb[0] = vertices[0];
    aabb[1] = vertices[0];

    for (int i = 1; i < vertices.size(); ++i) {
        if (vertices[i][0] < aabb[0][0]) aabb[0][0] = vertices[i][0];
        if (vertices[i][1] < aabb[0][1]) aabb[0][1] = vertices[i][1];
        if (vertices[i][2] < aabb[0][2]) aabb[0][2] = vertices[i][2];

        if (vertices[i][0] > aabb[1][0]) aabb[1][0] = vertices[i][0];
        if (vertices[i][1] > aabb[1][1]) aabb[1][1] = vertices[i][1];
        if (vertices[i][2] > aabb[1][2]) aabb[1][2] = vertices[i][2];
    }
}

template<typename T>
int ZD::MeshMisc::ComputeGridID(const CPoint<T, 3> *aabb, const CPoint<T, 3>& dis, const int grid_size, const CPoint<T, 3>& v)
{
    int ix = (v[0] - aabb[0][0]) / dis[0];
    ix = std::min(ix, grid_size - 1);
    ix = std::max(ix, 0);

    int iy = (v[1] - aabb[0][1]) / dis[1];
    iy = std::min(iy, grid_size - 1);
    iy = std::max(iy, 0);

    int iz = (v[2] - aabb[0][2]) / dis[2];
    iz = std::min(iz, grid_size - 1);
    iz = std::max(iz, 0);

    return ((iz * grid_size + iy) * grid_size + ix);
}

template<typename T>
ZD::CMesh<T>::CMesh()
{
    m_pVertices = nullptr;
    m_pFaces = nullptr;
    m_vNum = 0;
    m_fNum = 0;

    m_AABB[0].SetNan();
    m_AABB[1].SetNan();

    m_pColors = nullptr;
}

template<typename T>
ZD::CMesh<T>::CMesh(const std::vector<CPoint<T, 3>>& triangles)
{
    m_pVertices = nullptr;
    m_pFaces = nullptr;
    m_vNum = 0;
    m_fNum = 0;

    m_AABB[0].SetNan();
    m_AABB[1].SetNan();

    CreateMesh(triangles);

    ComputeAABB();

    m_pColors = nullptr;
}

template<typename T>
ZD::CMesh<T>::CMesh(const std::vector<CPoint<int, 3>>& triangles, const std::vector<CPoint<T, 3>>& vertices)
{
    m_pVertices = nullptr;
    m_pFaces = nullptr;
    m_vNum = 0;
    m_fNum = 0;

    m_AABB[0].SetNan();
    m_AABB[1].SetNan();

    CreateMesh(triangles, vertices);

    ComputeAABB();

    m_pColors = nullptr;
}

template<typename T>
ZD::CMesh<T>::~CMesh()
{
    SafeDeleteArray(m_pVertices);
    SafeDeleteArray(m_pFaces);
}


template<typename T>
int ZD::CMesh<T>::OpenMesh(const char * pathname, FILE * fp)
{
    if (m_pFaces != nullptr || m_pVertices != nullptr) {
        fprintf(fp, "error, already has data.\n");
        return -3;
    }

    if (pathname == nullptr) {
        fprintf(fp, "error, no input file!\n");
        return -1;
    }
    FILE *in = fopen(pathname, "rb");
    if (nullptr == in) {
        fprintf(fp, "error, cannot open input file %s.\n", pathname);
        return -2;
    }

    if (fread(&m_vNum, sizeof(int), 1, in) != 1) {
        fprintf(fp, "error, cannot read vertex number.\n");
        fclose(in);
        return -5;
    }

    try {
        m_pVertices = new CPoint<T, 3>[m_vNum];
    }
    catch (...) {
        fprintf(fp, "Cannot alloc memory for vertices!\n");
        SafeDeleteArray(m_pVertices);
        fclose(in);
        return -4;
    }

    if (fread(m_pVertices, sizeof(CPoint<T, 3>), m_vNum, in) != m_vNum) {
        fprintf(fp, "error, cannot read vertices.\n");
        fclose(in);
        return -6;
    }

    if (fread(&m_fNum, sizeof(int), 1, in) != 1) {
        fprintf(fp, "error, cannot read face number.\n");
        fclose(in);
        return -7;
    }

    try {
        m_pFaces = new CPoint<int, 3>[m_fNum];
    }
    catch (...) {
        fprintf(fp, "Cannot alloc memory for faces!\n");
        SafeDeleteArray(m_pVertices);
        SafeDeleteArray(m_pFaces);
        fclose(in);
        return -4;
    }

    if (fread(m_pFaces, sizeof(CPoint<int, 3>), m_fNum, in) != m_fNum) {
        fprintf(fp, "error, cannot read faces.\n");
        fclose(in);
        return -8;
    }

    fclose(in);

    return 0;
}

template<typename T>
int ZD::CMesh<T>::SaveMesh(const char * pathname, FILE * fp) const
{
    if (m_pFaces == nullptr || m_pVertices == nullptr) {
        fprintf(fp, "error, no data.\n");
        return -3;
    }

    if (pathname == nullptr) {
        fprintf(fp, "error, no output file!\n");
        return -1;
    }

    FILE *out = fopen(pathname, "wb");
    if (nullptr == out) {
        fprintf(fp, "error, cannot open output file %s.\n", pathname);
        return -2;
    }

    fwrite(&m_vNum, sizeof(int), 1, out);
    fwrite(m_pVertices, sizeof(CPoint<T, 3>), m_vNum, out);
    fwrite(&m_fNum, sizeof(int), 1, out);
    fwrite(m_pFaces, sizeof(CPoint<int, 3>), m_fNum, out);

    fclose(out);

    return 0;
}

template <typename T>
int ZD::CMesh<T>::OpenMeshVTK(const char *pathname, FILE *fp)
{
    if (m_pFaces != NULL || m_pVertices != NULL) {
        fprintf(fp, "error, already has data.\n");
        return -3;
    }

    if (pathname == NULL) {
        fprintf(fp, "error, no input file!\n");
        return -1;
    }
    FILE *in = fopen(pathname, "r");
    if (NULL == in) {
        fprintf(fp, "error, cannot open input file %s.\n", pathname);
        return -2;
    }

    /* only read ascii format */
    char token[256];
    do {
        fscanf(in, "%s", token);
    } while (strcmp(token, "POINTS") != 0);
    fscanf(in, "%d", &m_vNum);
    fscanf(in, "%s", token);

    try {
        m_pVertices = new CPoint<T, 3>[m_vNum];
    }
    catch (...) {
        fprintf(fp, "Cannot alloc memory for vertices!\n");
        SafeDeleteArray(m_pVertices);
        fclose(in);
        return -3;
    }

    for (int i = 0; i < m_vNum; ++i) {
        double x, y, z;
        fscanf(in, "%lf %lf %lf\n", &x, &y, &z);
        m_pVertices[i] = CPoint<T, 3>(x, y, z);
    }

    do {
        fscanf(in, "%s", token);
    } while (strcmp(token, "POLYGONS") != 0);
    fscanf(in, "%d", &m_fNum);
    fscanf(in, "%s", token);

    try {
        m_pFaces = new CPoint<int, 3>[m_fNum];
    }
    catch (...) {
        fprintf(fp, "Cannot alloc memory for faces!\n");
        SafeDeleteArray(m_pVertices);
        SafeDeleteArray(m_pFaces);
        fclose(in);
        return -4;
    }

    for (int i = 0; i < m_fNum; ++i) {
        int n, x, y, z;
        fscanf(in, "%d %d %d %d\n", &n, &m_pFaces[i][0], &m_pFaces[i][1], &m_pFaces[i][2]);
    }

    fclose(in);

    return 0;
}

template <typename T>
int ZD::CMesh<T>::SaveMeshVTK(const char *pathname, FILE *fp) const
{
    if (m_pFaces == nullptr || m_pVertices == nullptr) {
        fprintf(fp, "error, no data.\n");
        return -3;
    }

    if (pathname == nullptr) {
        fprintf(fp, "error, no input file!\n");
        return -1;
    }

    vtkSmartPointer<vtkPoints> pPoints = vtkSmartPointer<vtkPoints>::New();
    if (sizeof(T) == sizeof(double))
        pPoints->SetDataTypeToDouble();
    else
        pPoints->SetDataTypeToFloat();
    pPoints->SetNumberOfPoints(m_vNum);
    for (int vID = 0; vID < m_vNum; ++vID) {
        pPoints->SetPoint(vID, m_pVertices[vID][0], m_pVertices[vID][1], m_pVertices[vID][2]);
    }

    vtkSmartPointer<vtkCellArray> pTriangles = vtkSmartPointer<vtkCellArray>::New();
    for (int fID = 0; fID < m_fNum; ++fID) {
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0, m_pFaces[fID][0]);
        triangle->GetPointIds()->SetId(1, m_pFaces[fID][1]);
        triangle->GetPointIds()->SetId(2, m_pFaces[fID][2]);
        pTriangles->InsertNextCell(triangle);
    }

    vtkSmartPointer<vtkUnsignedCharArray> pColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    if (m_pColors) {
        pColors->SetNumberOfComponents(3);
        pColors->SetName("Colors");
        for (int vID = 0; vID < m_vNum; ++vID) {
            unsigned char color[3];
            color[0] = (unsigned char)(m_pColors[vID][0] * 255.0);
            color[1] = (unsigned char)(m_pColors[vID][1] * 255.0);
            color[2] = (unsigned char)(m_pColors[vID][2] * 255.0);
            pColors->InsertTuple3(vID, color[0], color[1], color[2]);
        }
    }

    vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
    pPolyData->SetPoints(pPoints);
    pPolyData->SetPolys(pTriangles);
    if (m_pColors) {
        pPolyData->GetPointData()->SetScalars(pColors);
    }

    vtkSmartPointer<vtkPolyDataWriter> pWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    pWriter->SetInputData(pPolyData);
    //pWriter->SetFileTypeToBinary();
    pWriter->SetFileTypeToASCII();
    pWriter->SetFileName(pathname);
    pWriter->Write();

    return 0;
}

template<typename T>
void ZD::CMesh<T>::AABB(CPoint<T, 3>* aabb)
{
    //if (std::isnan(m_AABB[0][0]))
    if (std::isnan(m_AABB[0][0]))
        ComputeAABB();

    aabb[0] = m_AABB[0];
    aabb[1] = m_AABB[1];
}


template<typename T>
void ZD::CMesh<T>::CreateMesh(const std::vector<CPoint<T, 3>>& triangles)
{
    CPoint<T, 3> aabb[2];
    MeshMisc::ComputeAABB(triangles, aabb);

    const int grid_size = 8;
    CPoint<T, 3> dis;
    dis[0] = (aabb[1][0] - aabb[0][0]) / (T)grid_size;
    dis[1] = (aabb[1][1] - aabb[0][1]) / (T)grid_size;
    dis[2] = (aabb[1][2] - aabb[0][2]) / (T)grid_size;

    std::vector<std::vector<int>> grids;
    grids.resize(grid_size*grid_size*grid_size);

    std::vector<CPoint<T, 3>> vertices;
    std::vector<CPoint<int, 3>> faces;

    for (int i = 0; i < triangles.size() / 3; ++i) {
        CPoint<T, 3> v0 = triangles.at(i*3+0);
        CPoint<T, 3> v1 = triangles.at(i*3+1);
        CPoint<T, 3> v2 = triangles.at(i*3+2);

        int g0 = MeshMisc::ComputeGridID(aabb, dis, grid_size, v0);
        int g1 = MeshMisc::ComputeGridID(aabb, dis, grid_size, v1);
        int g2 = MeshMisc::ComputeGridID(aabb, dis, grid_size, v2);

        int i0 = MeshMisc::FindVertex(vertices, grids[g0], v0);
        int i1 = MeshMisc::FindVertex(vertices, grids[g1], v1);
        int i2 = MeshMisc::FindVertex(vertices, grids[g2], v2);

        if (i0 == -1) {
            i0 = vertices.size();
            vertices.push_back(v0);
            grids[g0].push_back(i0);
        }
        if (i1 == -1) {
            i1 = vertices.size();
            vertices.push_back(v1);
            grids[g1].push_back(i1);
        }
        if (i2 == -1) {
            i2 = vertices.size();
            vertices.push_back(v2);
            grids[g2].push_back(i2);
        }

        CPoint<int, 3> f(i0, i1, i2);
        faces.push_back(f);

        //if (i % 1000 == 0) {
        //    printf("%d / %d\n", i, triangles.size() / 3);
        //}
    }

    m_vNum = vertices.size();
    m_pVertices = new CPoint<T, 3>[m_vNum];
    for (int i = 0; i < m_vNum; ++i)
        m_pVertices[i] = vertices[i];

    m_fNum = faces.size();
    m_pFaces = new CPoint<int, 3>[m_fNum];
    for (int i = 0; i < m_fNum; ++i)
        m_pFaces[i] = faces[i];
}

template<typename T>
void ZD::CMesh<T>::CreateMesh(const std::vector<CPoint<int, 3>>& triangles, const std::vector<CPoint<T, 3>>& vertices)
{
    m_vNum = vertices.size();
    m_pVertices = new CPoint<T, 3>[m_vNum];
    for (int i = 0; i < m_vNum; ++i)
        m_pVertices[i] = vertices[i];

    m_fNum = triangles.size();
    m_pFaces = new CPoint<int, 3>[m_fNum];
    for (int i = 0; i < m_fNum; ++i)
        m_pFaces[i] = triangles[i];
}

template<typename T>
void ZD::CMesh<T>::ComputeAABB()
{
    m_AABB[0] = m_pVertices[0];
    m_AABB[1] = m_pVertices[0];

    for (int i = 1; i < m_vNum; ++i) {
        if (m_pVertices[i][0] < m_AABB[0][0])
            m_AABB[0][0] = m_pVertices[i][0];
        if (m_pVertices[i][1] < m_AABB[0][1])
            m_AABB[0][1] = m_pVertices[i][1];
        if (m_pVertices[i][2] < m_AABB[0][2])
            m_AABB[0][2] = m_pVertices[i][2];

        if (m_pVertices[i][0] > m_AABB[1][0])
            m_AABB[1][0] = m_pVertices[i][0];
        if (m_pVertices[i][1] > m_AABB[1][1])
            m_AABB[1][1] = m_pVertices[i][1];
        if (m_pVertices[i][2] > m_AABB[1][2])
            m_AABB[1][2] = m_pVertices[i][2];
    }
}

template<typename T>
void ZD::CMesh<T>::SetColors(const std::vector<CPoint<T, 3>>& colors)
{
    SafeDeleteArray(m_pColors);
    m_pColors = new CPoint<T, 3>[m_vNum];
    for (int i = 0; i < m_vNum; ++i)
        m_pColors[i] = colors[i];
}

#endif
