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
#ifndef _ZD_LIB_LINE_HPP_
#define _ZD_LIB_LINE_HPP_

#include <vector>

#include "ZD_Point.hpp"

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

#include <assert.h>
#include <Eigen/Eigenvalues>


namespace ZD {
    template <typename T, unsigned int N>
    class CLine {
    public:
        CPoint<T, N> *m_pVertices;
        int m_count;
        CPoint<T, N> m_seed;
        CPoint<unsigned char, 3> *m_pColors;

    public:
        CLine(void);
        ~CLine(void);

    public:
        void CreateLine(const CPoint<T, N> *pts, const int count, const CPoint<T, N>& seed);
        void CreateLine(const std::vector<CPoint<T, N>>& pts);
        void CreateLine(const CLine<T, N> *line);
        void CreateLine(const CLine<T, N>& line);

        void ComputeColor();
        void SetColor(const T *color);

        const CPoint<T, 9> Parameterization3() const;

        ///* fixed number vertices parameterization */
        //void ParameterizeFixed(CPoint<T, N> *param, const int num) const;

    public:
        friend void ReadLines(const char *pathname, CLine<T, N> **lines, int &size)
        {
            FILE *fp = fopen(pathname, "rb");
            if (fp == nullptr)
                return;

            fread(&size, sizeof(int), 1, fp);
            *lines = new CLine<T, N>[size];

            for (int i = 0; i < size; ++i) {
                int count;
                fread(&count, sizeof(int), 1, fp);
                (*lines)[i].m_count = count;
                fread(&((*lines)[i].m_seed), sizeof(CPoint<T, N>), 1, fp);

                (*lines)[i].m_pVertices = new CPoint<T, N>[count];
                fread((*lines)[i].m_pVertices, sizeof(CPoint<T, N>), count, fp);
            }

            fclose(fp);
        }
        friend void SaveLines(const char *pathname, const CLine<T, N> *lines, const int size)
        {
            FILE *fp = fopen(pathname, "wb");
            if (fp == nullptr)
                return;

            fwrite(&size, sizeof(int), 1, fp);
            for (int i = 0; i < size; ++i) {
                fwrite(&(lines[i].m_count), sizeof(int), 1, fp);
                fwrite(&(lines[i].m_seed), sizeof(CPoint<T, N>), 1, fp);
                fwrite(lines[i].m_pVertices, sizeof(CPoint<T, N>), lines[i].m_count, fp);
            }
            fclose(fp);
        }
        friend void ReadLinesVTK(const char *pathname, CLine<T, N> **lines, int &size)
        {
            vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
            reader->SetFileName(pathname);
            reader->Update();

            int cell_size = reader->GetOutput()->GetNumberOfCells();
            vtkCellArray *cells = reader->GetOutput()->GetLines();
            size = 0;
            while (cell_size > 0) {
                vtkIdType *pts;
                vtkIdType npts;
                cells->GetNextCell(npts, pts);
                size++;
                cell_size -= (npts + 1);
            }

            *lines = new CLine<T, N>[size];
            cells->InitTraversal();
            for (int i = 0; i < size; ++i) {
                vtkIdType *pts;
                vtkIdType npts;
                cells->GetNextCell(npts, pts);

                (*lines)[i].m_count = npts;
                (*lines)[i].m_pVertices = new CPoint<T, N>[npts];
                (*lines)[i].m_pColors = new CPoint<unsigned char, 3>[npts];

                for (int ipt = 0; ipt < npts; ++ipt) {
                    double *pt = reader->GetOutput()->GetPoint(pts[ipt]);
                    if (N == 2) {
                        (*lines)[i].m_pVertices[ipt][0] = pt[0];
                        (*lines)[i].m_pVertices[ipt][1] = pt[1];
                    }
                    else if (N == 3) {
                        (*lines)[i].m_pVertices[ipt][0] = pt[0];
                        (*lines)[i].m_pVertices[ipt][1] = pt[1];
                        (*lines)[i].m_pVertices[ipt][2] = pt[2];
                    }
                    else {
                        ;
                    }
                    
                    if (reader->GetOutput()->GetPointData()->GetScalars("color")) {
                        double cr = reader->GetOutput()->GetPointData()->GetScalars("color")->GetComponent(pts[ipt], 0);
                        double cg = reader->GetOutput()->GetPointData()->GetScalars("color")->GetComponent(pts[ipt], 1);
                        double cb = reader->GetOutput()->GetPointData()->GetScalars("color")->GetComponent(pts[ipt], 2);
                        (*lines)[i].m_pColors[ipt][0] = cr;
                        (*lines)[i].m_pColors[ipt][1] = cg;
                        (*lines)[i].m_pColors[ipt][2] = cb;
                    }
                    else {
                        (*lines)[i].m_pColors[ipt][0] = 255;
                        (*lines)[i].m_pColors[ipt][1] = 255;
                        (*lines)[i].m_pColors[ipt][2] = 255;
                    }
                }
                (*lines)[i].m_seed = (*lines)[i].m_pVertices[0];
            }
        }
        friend void SaveLinesVTK(const char *pathname, const CLine<T, N> *lines, const int size)
        {
            vtkIdType npts = 0;
            for (int i = 0; i < size; ++i)
                npts += lines[i].m_count;

            vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
            pts->SetDataTypeToDouble();
            pts->SetNumberOfPoints(npts);
            vtkIdType offset = 0;
            if (N == 2) {
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < lines[i].m_count; ++j) {
                        pts->SetPoint(j + offset, lines[i].m_pVertices[j][0], lines[i].m_pVertices[j][1], 0.0);
                    }
                    offset += lines[i].m_count;
                }
            }
            else if (N == 3) {
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < lines[i].m_count; ++j) {
                        pts->SetPoint(j + offset, lines[i].m_pVertices[j][0], lines[i].m_pVertices[j][1], lines[i].m_pVertices[j][2]);
                    }
                    offset += lines[i].m_count;
                }
            }
            else {
                ;
            }
            
            vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
            offset = 0;
            for (int i = 0; i < size; ++i) {
                cells->InsertNextCell(lines[i].m_count);
                for (int j = 0; j < lines[i].m_count; ++j) {
                    cells->InsertNextCell(j + offset);
                }
                offset += lines[i].m_count;
            }

            vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
            colors->SetNumberOfComponents(3);
            colors->SetName("color");
            offset = 0;
            for (int i = 0; i < size; ++i) {
                if (lines[i].m_pColors == nullptr) {
                    for (int j = 0; j < lines[i].m_count; ++j) {
                        unsigned char cc[3] = { 255, 255, 255 };
                        //colors->InsertNextTupleValue(cc);
                        colors->InsertTuple3(j+offset, cc[0], cc[1], cc[2]);
                    }
                }
                else {
                    for (int j = 0; j < lines[i].m_count; ++j) {
                        unsigned char cc[3];
                        cc[0] = lines[i].m_pColors[j][0];
                        cc[1] = lines[i].m_pColors[j][1];
                        cc[2] = lines[i].m_pColors[j][2];
                        //colors->InsertNextTupleValue(cc);
                        colors->InsertTuple3(j+offset, cc[0], cc[1], cc[2]);
                    }
                }
                offset += lines[i].m_count;
            }

            vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
            polydata->SetPoints(pts);
            polydata->SetLines(cells);
            polydata->GetPointData()->SetScalars(colors);

            vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
            writer->SetInputData(polydata);
            writer->SetFileName(pathname);
            //writer->SetFileTypeToBinary();
            writer->Write();
        }
    };
} // namespace ZD


template <typename T, unsigned int N>
ZD::CLine<T, N>::CLine()
{
    m_pVertices = nullptr;
    m_count = 0;
    m_seed.SetZero();
    m_pColors = nullptr;
}

template <typename T, unsigned int N>
ZD::CLine<T, N>::~CLine()
{
    SafeDeleteArray(m_pVertices);
    m_count = 0;
    m_seed.SetZero();
    SafeDeleteArray(m_pColors);
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::CreateLine(const CPoint<T, N> *pts, const int count, const CPoint<T, N>& seed)
{
    SafeDelete(this->m_pVertices);
    SafeDelete(this->m_pColors);
    
    this->m_count = count;
    this->m_seed = seed;
    this->m_pVertices = new CPoint<T, N>[this->m_count];
    memcpy(this->m_pVertices, pts, sizeof(CPoint<T, N>)*count);
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::CreateLine(const std::vector<CPoint<T, N>>& pts)
{
    SafeDelete(this->m_pVertices);
    SafeDelete(this->m_pColors);

    if (pts.size() > 0) {
        this->m_count = pts.size();
        this->m_seed = pts[0];
        this->m_pVertices = new CPoint<T, N>[this->m_count];
        for (size_t i = 0; i < pts.size(); ++i) {
            this->m_pVertices[i] = pts[i];
        }
    }
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::CreateLine(const CLine<T, N> *line)
{
    this->CreateLine(line->m_pVertices, line->m_count, line->m_seed);
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::CreateLine(const CLine<T, N>& line)
{
    this->CreateLine(line.m_pVertices, line.m_count, line.m_seed);
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::ComputeColor()
{
    SafeDelete(this->m_pColors);
    this->m_pColors = new CPoint<unsigned char, 3>[this->m_count];

    if (this->m_count == 1) {
        this->m_pColors[0][0] = 128;
        this->m_pColors[0][1] = 128;
        this->m_pColors[0][2] = 128;
    } else {
        CPoint<T, N> dir;
        for (int i = 1; i < this->m_count - 1; ++i) {
            dir = m_pVertices[i + 1] - m_pVertices[i - 1];
            dir.Normalize();
            if (N == 2) {
                this->m_pColors[i][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
                this->m_pColors[i][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
                this->m_pColors[i][2] = 0.0;
            }
            else if (N == 3) {
                this->m_pColors[i][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
                this->m_pColors[i][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
                this->m_pColors[i][2] = (unsigned char)(std::abs(dir[2]) * 255.0);
            }
            else {
                ;
            }
        }

        if (N == 2) {
            dir = m_pVertices[1] - m_pVertices[0];
            dir.Normalize();
            this->m_pColors[0][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            this->m_pColors[0][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            this->m_pColors[0][2] = 0.0;

            dir = m_pVertices[this->m_count - 2] - m_pVertices[this->m_count - 1];
            dir.Normalize();
            this->m_pColors[this->m_count - 1][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            this->m_pColors[this->m_count - 1][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            this->m_pColors[this->m_count - 1][2] = 0.0;
        }
        else if (N == 3) {
            dir = m_pVertices[1] - m_pVertices[0];
            dir.Normalize();
            this->m_pColors[0][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            this->m_pColors[0][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            this->m_pColors[0][2] = (unsigned char)(std::abs(dir[2]) * 255.0);

            dir = m_pVertices[this->m_count - 2] - m_pVertices[this->m_count - 1];
            dir.Normalize();
            this->m_pColors[this->m_count - 1][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            this->m_pColors[this->m_count - 1][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            this->m_pColors[this->m_count - 1][2] = (unsigned char)(std::abs(dir[2]) * 255.0);
        }
        else {
            ;
        }
    }
}

template <typename T, unsigned int N>
void ZD::CLine<T, N>::SetColor(const T *color)
{
    SafeDelete(this->m_pColors);
    this->m_pColors = new CPoint<unsigned char, 3>[this->m_count];

    for (int i = 0; i < this->m_count; ++i) {
        this->m_pColors[i][0] = (unsigned char)(color[0] * 255.0);
        this->m_pColors[i][1] = (unsigned char)(color[1] * 255.0);
        this->m_pColors[i][2] = (unsigned char)(color[2] * 255.0);
    }
}

template <typename T, unsigned int N>
const ZD::CPoint<T, 9> ZD::CLine<T, N>::Parameterization3() const
{
    assert(N == 3);

    CPoint<T, 9> params;
    params.SetZero();
    if (m_count == 0) {
        ;
    }
    else if (m_count == 1) {
        params[0] = m_pVertices[0][0];
        params[1] = m_pVertices[0][1];
        params[2] = m_pVertices[0][2];
    }
    else {
        CPoint<T, N> mean;
        mean.SetZero();
        for (int j = 0; j < m_count; ++j) {
            mean += m_pVertices[j];
        }
        mean /= (T)m_count;
        params[0] = mean[0]; params[1] = mean[1]; params[2] = mean[2];

        Eigen::Matrix3d cov;
        cov.setZero();
        for (int j = 0; j < m_count; ++j) {
            cov(0, 0) += (m_pVertices[j][0] - mean[0]) * (m_pVertices[j][0] - mean[0]);
            cov(0, 1) += (m_pVertices[j][0] - mean[0]) * (m_pVertices[j][1] - mean[1]);
            cov(0, 2) += (m_pVertices[j][0] - mean[0]) * (m_pVertices[j][2] - mean[2]);
            cov(1, 0) += (m_pVertices[j][1] - mean[1]) * (m_pVertices[j][0] - mean[0]);
            cov(1, 1) += (m_pVertices[j][1] - mean[1]) * (m_pVertices[j][1] - mean[1]);
            cov(1, 2) += (m_pVertices[j][1] - mean[1]) * (m_pVertices[j][2] - mean[2]);
            cov(2, 0) += (m_pVertices[j][2] - mean[2]) * (m_pVertices[j][0] - mean[0]);
            cov(2, 1) += (m_pVertices[j][2] - mean[2]) * (m_pVertices[j][1] - mean[1]);
            cov(2, 2) += (m_pVertices[j][2] - mean[2]) * (m_pVertices[j][2] - mean[2]);
        }
        cov /= (T)m_count;

        Eigen::EigenSolver<Eigen::Matrix3d> es(cov);
        Eigen::Matrix3cd D = es.eigenvalues().asDiagonal();
        Eigen::Matrix3cd V = es.eigenvectors();

        D(0, 0) = std::sqrt(D(0, 0));
        D(1, 1) = std::sqrt(D(1, 1));
        D(2, 2) = std::sqrt(D(2, 2));
        Eigen::Matrix3cd temp = V * D * V.inverse();

        params[3] = temp(0, 0).real(); params[4] = temp(0, 1).real(); params[5] = temp(0, 2).real();
        params[6] = temp(1, 1).real(); params[7] = temp(1, 2).real(); params[8] = temp(2, 2).real();
    }
    return params;
}

#endif
