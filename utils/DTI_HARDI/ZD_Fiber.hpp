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
#ifndef _ZD_LIB_FIBER_HPP_
#define _ZD_LIB_FIBER_HPP_

#include <vector>


#include "utils/define.hpp"
#include "utils/Base/ZD_Point.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "utils/Base/ZD_Line.hpp"
#include "utils/Tool/ZD_ParameterizationTool.hpp"
#include "utils/Tool/ZD_ColorSpaceTool.hpp"

#include <Eigen/Eigenvalues>

#include "ZD_DTI.hpp"

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

namespace ZD {
    template<typename T>
    class CFiber {
    public:
        CPoint<T, 3> *m_pForward, *m_pBackward;
        int m_fCount, m_bCount;
        CPoint<T, 3> m_seedPos;
        CPoint<T, 3> m_seedDir;

        CPoint<unsigned char, 3> *m_pFColor, *m_pBColor;

    public:
        CFiber();
        CFiber(const CPoint<T, 3> *pf, const int fc, const CPoint<T, 3> *pb, const int bc,
            const CPoint<T, 3>& pos, const CPoint<T, 3>& dir);
        CFiber(const std::vector<CPoint<T, 3>>& forward, const std::vector<CPoint<T, 3>>& backward,
            const CPoint<T, 3>& pos, const CPoint<T, 3>& dir);
        CFiber(const CFiber<T> *fiber);
        ~CFiber();

    public:
        void CreateFiber(const CPoint<T, 3> *pf, const int fc, const CPoint<T, 3> *pb, const int bc,
            const CPoint<T, 3>& pos, const CPoint<T, 3>& dir);
        void CreateFiber(const std::vector<CPoint<T, 3>>& forward, const std::vector<CPoint<T, 3>>& backward,
            const CPoint<T, 3>& pos, const CPoint<T, 3>& dir);
        void CreateFiber(const CFiber<T> *fiber);

        void ComputeColor();

        void ComputeColorByFA(const CDTI<T> *dti);

        /* parameterization used in fiber function method */
        void ParameterizeFF(CPoint<T, 3>& param, const int length);        // use 1st order moment
        void ParameterizeFF(CPoint<T, 9>& param, const int length);        // use 1st and 2nd order moment
        /* parameterization used in FSR method */
        void ParameterizeFSR(CPoint<T, 9>& param, const int length);
        /* parameterization used in resample method */
        void ParameterizeResample(CPoint<T, 3> *param, const int length, const int k);
        /* scalar quantity parameterization */
        CPoint<T, 2> ParameterizeScalar(CField3<T, 1> *scalar, const int length);
        //CPoint<T, 1> ParameterizeScalar(CField3<T, 1> *scalar, const int length);

        const CPoint<T, 3> GetOrientation() const;

    public:
        friend int OpenDTIFibers(const char *pathname, CFiber<T>*& fibers, unsigned int& count)
        {
            FILE *fp = fopen(pathname, "rb");
            if (nullptr == fp)
                return 0;

            fread(&count, sizeof(int), 1, fp);
            fibers = new CFiber<T>[count];
            for (int i = 0; i < count; ++i) {
                CPoint<T, 3> pos, dir;
                int fc, bc;
                CPoint<T, 3> *forward, *backward;
                // seed pos
                fread(&pos, sizeof(CPoint<T, 3>), 1, fp);
                // seed dir
                fread(&dir, sizeof(CPoint<T, 3>), 1, fp);
                // forward count
                fread(&fc, sizeof(int), 1, fp);
                // forward vertices
                forward = new CPoint<T, 3>[fc];
                fread(forward, sizeof(CPoint<T, 3>), fc, fp);
                // backward count
                fread(&bc, sizeof(int), 1, fp);
                // backward vertices
                backward = new CPoint<T, 3>[bc];
                fread(backward, sizeof(CPoint<T, 3>), bc, fp);
                fibers[i].CreateFiber(forward, fc, backward, bc, pos, dir);
                SafeDeleteArray(forward);
                SafeDeleteArray(backward);
            }
            fclose(fp);

            return count;
        }
        friend int OpenHARDIFibers(const char *pathname, CFiber<T>** fibers, unsigned int *fiberCount, unsigned int& count) 
        {
            FILE *fp = fopen(pathname, "rb");
            if (nullptr == fp)
                return 0;

            // get the size of seeds
            fread(&count, sizeof(unsigned int), 1, fp);
            
            // get the fiber count for each seed
            fiberCount = new unsigned int[count];
            fread(fiberCount, sizeof(unsigned int), count, fp);

            // read fibers
            fibers = new CFiber<T>*[count];
            for (unsigned int i = 0; i < count; ++i) {
                fibers[i] = new CFiber<T>[fiberCount[i]];
                for (unsigned int j = 0; j < fiberCount[i]; ++j) {
                    CPoint<T, 3> pos, dir;
                    int fc, bc;
                    CPoint<T, 3> *forward, *backward;
                    // seed pos
                    fread(&pos, sizeof(CPoint<T, 3>), 1, fp);
                    // seed dir
                    fread(&dir, sizeof(CPoint<T, 3>), 1, fp);
                    // forward count
                    fread(&fc, sizeof(int), 1, fp);
                    // forward vertices
                    forward = new CPoint<T, 3>[fc];
                    fread(forward, sizeof(CPoint<T, 3>), fc, fp);
                    // backward count
                    fread(&bc, sizeof(int), 1, fp);
                    // backward vertices
                    backward = new CPoint<T, 3>[bc];
                    fread(backward, sizeof(CPoint<T, 3>), bc, fp);
                    fibers[i][j].CreateFiber(forward, fc, backward, bc, pos, dir);
                    SafeDeleteArray(forward);
                    SafeDeleteArray(backward);
                }
            }
            fclose(fp);

            return count;
        }

        friend void SaveDTIFibers(const char *pathname, const CFiber<T> *fibers, const unsigned int count)
        {
            FILE *fp = fopen(pathname, "wb");
            if (nullptr == fp)
                return;

            fwrite(&count, sizeof(int), 1, fp);
            for (unsigned int i = 0; i < count; ++i) {
                // seed pos
                fwrite(&(fibers[i].m_seedPos), sizeof(CPoint<T, 3>), 1, fp);
                // seed dir
                fwrite(&(fibers[i].m_seedDir), sizeof(CPoint<T, 3>), 1, fp);
                // forward count
                fwrite(&(fibers[i].m_fCount), sizeof(int), 1, fp);
                // forward vertices
                fwrite(fibers[i].m_pForward, sizeof(CPoint<T, 3>), fibers[i].m_fCount, fp);
                // backward count
                fwrite(&(fibers[i].m_bCount), sizeof(int), 1, fp);
                // backward vertices
                fwrite(fibers[i].m_pBackward, sizeof(CPoint<T, 3>), fibers[i].m_bCount, fp);
            }
            fclose(fp);
        }

        friend void SaveHARDIFibers(const char *pathname, const CFiber<T> **fibers, const unsigned int *fiberCount, const unsigned int count)
        {
            FILE *fp = fopen(pathname, "wb");
            if (nullptr == fp)
                return;

            // save the size of seeds
            fwrite(&count, sizeof(unsigned int), 1, fp);

            // save the number of fibers of each seed
            fwrite(fiberCount, sizeof(unsigned int), count, fp);

            // save fibers
            for (unsigned int i = 0; i < count; ++i) {
                for (unsigned int j = 0; j < fiberCount[i]; ++j) {
                    // seed pos
                    fwrite(&(fibers[i][j].m_seedPos), sizeof(CPoint<T, 3>), 1, fp);
                    // seed dir
                    fwrite(&(fibers[i][j].m_seedDir), sizeof(CPoint<T, 3>), 1, fp);
                    // forward count
                    fwrite(&(fibers[i][j].m_fCount), sizeof(int), 1, fp);
                    // forward vertices
                    fwrite(fibers[i][j].m_pForward, sizeof(CPoint<T, 3>), fibers[i][j].m_fCount, fp);
                    // backward count
                    fwrite(&(fibers[i][j].m_bCount), sizeof(int), 1, fp);
                    // backward vertices
                    fwrite(fibers[i][j].m_pBackward, sizeof(CPoint<T, 3>), fibers[i][j].m_bCount, fp);
                }
            }
            fclose(fp);
        }

        friend void SaveFibersAsVTK(const char *pathname, const CFiber<T> *fibers, const unsigned int count)
        {
            vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
            pts->SetDataTypeToDouble();
            pts->Initialize();
            vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
            colors->SetNumberOfComponents(3);
            int offset = 0;
            for (unsigned int i = 0; i < count; ++i) {
                cells->InsertNextCell(fibers[i].m_fCount + fibers[i].m_bCount - 1);
                for (int k = fibers[i].m_bCount; k > 0; --k) {
                    pts->InsertNextPoint(fibers[i].m_pBackward[k - 1].m_data);
                    cells->InsertNextCell(offset++);
                    unsigned char cc[3] = { 0, 0, 0 };
                    if (fibers[i].m_pBColor == nullptr)
                        //colors->InsertNextTupleValue(cc);
                        colors->InsertTuple3(offset, cc[0], cc[1], cc[2]);
                    else
                        //colors->InsertNextTupleValue(fibers[i].m_pBColor[k - 1].m_data);
                        colors->InsertTuple3(offset, fibers[i].m_pBColor[k-1].m_data[0], fibers[i].m_pBColor[k-1].m_data[1], fibers[i].m_pBColor[k-1].m_data[2]);
                }
                for (int k = 1; k < fibers[i].m_fCount; ++k) {
                    pts->InsertNextPoint(fibers[i].m_pForward[k].m_data);
                    cells->InsertNextCell(offset++);
                    unsigned char cc[3] = { 0, 0, 0 };
                    if (fibers[i].m_pFColor == nullptr)
                        //colors->InsertNextTupleValue(cc);
                        colors->InsertTuple3(offset, cc[0], cc[1], cc[2]);
                    else
                        //colors->InsertNextTupleValue(fibers[i].m_pFColor[k - 1].m_data);
                        colors->InsertTuple3(offset, fibers[i].m_pFColor[k-1].m_data[0], fibers[i].m_pFColor[k-1].m_data[1], fibers[i].m_pFColor[k-1].m_data[2]);
                }                
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
}

template<typename T>
ZD::CFiber<T>::CFiber()
{
    m_pForward = nullptr;
    m_pBackward = nullptr;

    m_fCount = m_bCount = 0;
    m_pFColor = m_pBColor = nullptr;
}

template<typename T>
ZD::CFiber<T>::CFiber(const CPoint<T, 3> *pf, const int fc, const CPoint<T, 3> *pb,
    const int bc, const CPoint<T, 3>& pos, const CPoint<T, 3>& dir)
{
    this->CreateFiber(pf, fc, pb, bc, pos, dir);
    m_pFColor = m_pBColor = nullptr;
}

template<typename T>
ZD::CFiber<T>::CFiber(const std::vector<CPoint<T, 3>>& forward, const std::vector<CPoint<T, 3>>& backward,
    const CPoint<T, 3>& pos, const CPoint<T, 3>& dir)
{
    this->CreateFiber(forward, backward, pos, dir);
    m_pFColor = m_pBColor = nullptr;
}

template<typename T>
ZD::CFiber<T>::CFiber(const CFiber<T> *fiber)
{
    this->CreateFiber(fiber);
    m_pFColor = m_pBColor = nullptr;
}

template<typename T>
ZD::CFiber<T>::~CFiber()
{
    SafeDeleteArray(m_pForward);
    SafeDeleteArray(m_pBackward);
    
    m_fCount = m_bCount = 0;

    SafeDeleteArray(m_pFColor);
    SafeDeleteArray(m_pBColor);
}


template<typename T>
void ZD::CFiber<T>::CreateFiber(const CPoint<T, 3> *pf, const int fc, const CPoint<T, 3> *pb, 
    const int bc, const CPoint<T, 3>& pos, const CPoint<T, 3>& dir)
{
    m_fCount = fc;
    m_pForward = new CPoint<T, 3>[m_fCount];
    memcpy(m_pForward, pf, sizeof(CPoint<T, 3>)*m_fCount);

    m_bCount = bc;
    m_pBackward = new CPoint<T, 3>[m_bCount];
    memcpy(m_pBackward, pb, sizeof(CPoint<T, 3>)*m_bCount);

    m_seedPos = pos;
    m_seedDir = dir;
}

template<typename T>
void ZD::CFiber<T>::CreateFiber(const std::vector<CPoint<T, 3>>& forward, 
    const std::vector<CPoint<T, 3>>& backward, const CPoint<T, 3>& pos, const CPoint<T, 3>& dir)
{
    m_fCount = forward.size();
    m_pForward = new CPoint<T, 3>[m_fCount];
    int i = 0;
    for (int i = 0; i < m_fCount; ++i)
        m_pForward[i] = forward[i];

    m_bCount = backward.size();
    m_pBackward = new CPoint<T, 3>[m_bCount];
    for (int i = 0; i < m_bCount; ++i)
        m_pBackward[i] = backward[i];

    m_seedPos = pos;
    m_seedDir = dir;
}

template<typename T>
void ZD::CFiber<T>::CreateFiber(const CFiber<T> *fiber)
{
    this->m_fCount = fiber->m_fCount;
    this->m_pForward = new CPoint<T, 3>[this->m_fCount];
    memcpy(this->m_pForward, fiber->m_pForward, sizeof(CPoint<T, 3>)*this->m_fCount);

    this->m_bCount = fiber->m_bCount;
    this->m_pBackward = new CPoint<T, 3>[this->m_bCount];
    memcpy(this->m_pBackward, fiber->m_pBackward, sizeof(CPoint<T, 3>)*this->m_bCount);

    this->m_seedPos = fiber->m_seedPos;
    this->m_seedDir = fiber->m_seedDir;
}

template<typename T>
void ZD::CFiber<T>::ComputeColor()
{
    assert(m_pForward != nullptr && m_pBackward != nullptr);

    m_pFColor = new CPoint<unsigned char, 3>[m_fCount];
    if (m_fCount == 1) {
        m_pFColor[0][0] = 1.0;
        m_pFColor[0][1] = 1.0;
        m_pFColor[0][2] = 1.0;
    }
    else {
        CPoint<T, 3> dir;
        for (int i = 1; i < this->m_fCount - 1; ++i) {
            dir = m_pForward[i + 1] - m_pForward[i - 1];
            dir.Normalize();
            m_pFColor[i][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            m_pFColor[i][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            m_pFColor[i][2] = (unsigned char)(std::abs(dir[2]) * 255.0);
        }
        dir = m_pForward[1] - m_pForward[0];
        dir.Normalize();
        m_pFColor[0][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
        m_pFColor[0][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
        m_pFColor[0][2] = (unsigned char)(std::abs(dir[2]) * 255.0);

        dir = m_pForward[this->m_fCount - 2] - m_pForward[this->m_fCount - 1];
        dir.Normalize();
        m_pFColor[this->m_fCount - 1][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
        m_pFColor[this->m_fCount - 1][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
        m_pFColor[this->m_fCount - 1][2] = (unsigned char)(std::abs(dir[2]) * 255.0);

    }
    
    m_pBColor = new CPoint<unsigned char, 3>[m_bCount];
    if (m_bCount == 1) {
        m_pBColor[0][0] = 1.0;
        m_pBColor[0][1] = 1.0;
        m_pBColor[0][2] = 1.0;
    }
    else {
        CPoint<T, 3> dir;
        for (int i = 1; i < this->m_bCount - 1; ++i) {
            dir = m_pBackward[i + 1] - m_pBackward[i - 1];
            dir.Normalize();
            m_pBColor[i][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
            m_pBColor[i][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
            m_pBColor[i][2] = (unsigned char)(std::abs(dir[2]) * 255.0);
        }
        dir = m_pBackward[1] - m_pBackward[0];
        dir.Normalize();
        m_pBColor[0][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
        m_pBColor[0][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
        m_pBColor[0][2] = (unsigned char)(std::abs(dir[2]) * 255.0);

        dir = m_pBackward[this->m_bCount - 2] - m_pBackward[this->m_bCount - 1];
        dir.Normalize();
        m_pBColor[this->m_bCount - 1][0] = (unsigned char)(std::abs(dir[0]) * 255.0);
        m_pBColor[this->m_bCount - 1][1] = (unsigned char)(std::abs(dir[1]) * 255.0);
        m_pBColor[this->m_bCount - 1][2] = (unsigned char)(std::abs(dir[2]) * 255.0);
    }

    if (m_fCount == 1 && m_bCount > 1) {
        m_pFColor[0] = m_pBColor[0];
    }
    if (m_fCount > 1 && m_bCount == 1) {
        m_pBColor[0] = m_pFColor[0];
    }
}

template<typename T>
void ZD::CFiber<T>::ComputeColorByFA(const CDTI<T> *dti)
{
    assert(m_pForward != nullptr && m_pBackward != nullptr);

    m_pFColor = new CPoint<unsigned char, 3>[m_fCount];
    m_pBColor = new CPoint<unsigned char, 3>[m_bCount];

    for (int i = 0; i < m_fCount; ++i) {
        T fa = dti->GetFA(m_pForward[i]);
        fa = (fa - 0.15) / 0.517;
        CPoint<T, 3> hsv;
        hsv[0] = (1.0 - fa) * 240.0;
        hsv[1] = 1.0;
        hsv[2] = 1.0;
        CColorSpaceTool<T>::HSV2RGB(hsv, m_pFColor[i]);
    }

    for (int i = 0; i < m_bCount; ++i) {
        T fa = dti->GetFA(m_pBackward[i]);
        fa = (fa - 0.15) / 0.517;
        CPoint<T, 3> hsv;
        hsv[0] = (1.0 - fa) * 240.0;
        hsv[1] = 1.0;
        hsv[2] = 1.0;
        CColorSpaceTool<T>::HSV2RGB(hsv, m_pBColor[i]);
    }
}

template <typename T>
void ZD::CFiber<T>::ParameterizeFF(CPoint<T, 3>& param, const int length)
{
    if ((m_fCount + m_bCount) == 2) {
        param[0] = m_seedPos[0]; param[1] = m_seedPos[1]; param[2] = m_seedPos[2];
    }
    else {
        int fc, bc;
        if (length < 0) {
            fc = m_fCount - 1;
            bc = m_bCount - 1;
        }
        else {
            bc = length < (m_bCount - 1) ? length : (m_bCount - 1);
            fc = length < (m_fCount - 1) ? length : (m_fCount - 1);
        }

        CPoint<T, 3> mean = CPoint<T, 3>((T)0.0, (T)0.0, (T)0.0);
        for (int j = 0; j <= bc; ++j) {
            mean += m_pBackward[j];
        }
        for (int j = 1; j <= fc; ++j) {
            mean += m_pForward[j];
        }
        mean /= T(bc + fc + 1);
        param[0] = mean[0]; param[1] = mean[1]; param[2] = mean[2];
    }
}

template<typename T>
void ZD::CFiber<T>::ParameterizeFF(CPoint<T, 9>& param, const int length)
{
    if ((m_fCount + m_bCount) == 2) {
        param[0] = m_seedPos[0]; param[1] = m_seedPos[1]; param[2] = m_seedPos[2];
        param[3] = 0.0; param[4] = 0.0; param[5] = 0.0;    param[6] = 0.0; param[7] = 0.0;    param[8] = 0.0;
    }
    else {
        int fc, bc;
        if (length < 0) {
            fc = m_fCount - 1;
            bc = m_bCount - 1;
        } else {
            bc = length < (m_bCount - 1) ? length : (m_bCount - 1);
            fc = length < (m_fCount - 1) ? length : (m_fCount - 1);
        }
        
        CPoint<T, 3> mean = CPoint<T, 3>((T)0.0, (T)0.0, (T)0.0);
        for (int j = 0; j <= bc; ++j) {
            mean += m_pBackward[j];
        }
        for (int j = 1; j <= fc; ++j) {
            mean += m_pForward[j];
        }
        mean /= T(bc + fc + 1);
        param[0] = mean[0]; param[1] = mean[1]; param[2] = mean[2];

        Eigen::Matrix3d cov;
        cov.setZero();
        for (int j = 0; j <= bc; ++j) {
            cov(0, 0) += (m_pBackward[j][0] - mean[0]) * (m_pBackward[j][0] - mean[0]);
            cov(0, 1) += (m_pBackward[j][0] - mean[0]) * (m_pBackward[j][1] - mean[1]);
            cov(0, 2) += (m_pBackward[j][0] - mean[0]) * (m_pBackward[j][2] - mean[2]);
            cov(1, 0) += (m_pBackward[j][1] - mean[1]) * (m_pBackward[j][0] - mean[0]);
            cov(1, 1) += (m_pBackward[j][1] - mean[1]) * (m_pBackward[j][1] - mean[1]);
            cov(1, 2) += (m_pBackward[j][1] - mean[1]) * (m_pBackward[j][2] - mean[2]);
            cov(2, 0) += (m_pBackward[j][2] - mean[2]) * (m_pBackward[j][0] - mean[0]);
            cov(2, 1) += (m_pBackward[j][2] - mean[2]) * (m_pBackward[j][1] - mean[1]);
            cov(2, 2) += (m_pBackward[j][2] - mean[2]) * (m_pBackward[j][2] - mean[2]);
        }
        for (int j = 1; j <= fc; ++j) {
            cov(0, 0) += (m_pForward[j][0] - mean[0]) * (m_pForward[j][0] - mean[0]);
            cov(0, 1) += (m_pForward[j][0] - mean[0]) * (m_pForward[j][1] - mean[1]);
            cov(0, 2) += (m_pForward[j][0] - mean[0]) * (m_pForward[j][2] - mean[2]);
            cov(1, 0) += (m_pForward[j][1] - mean[1]) * (m_pForward[j][0] - mean[0]);
            cov(1, 1) += (m_pForward[j][1] - mean[1]) * (m_pForward[j][1] - mean[1]);
            cov(1, 2) += (m_pForward[j][1] - mean[1]) * (m_pForward[j][2] - mean[2]);
            cov(2, 0) += (m_pForward[j][2] - mean[2]) * (m_pForward[j][0] - mean[0]);
            cov(2, 1) += (m_pForward[j][2] - mean[2]) * (m_pForward[j][1] - mean[1]);
            cov(2, 2) += (m_pForward[j][2] - mean[2]) * (m_pForward[j][2] - mean[2]);
        }
        cov /= T(bc + fc + 1);

        Eigen::EigenSolver<Eigen::Matrix3d> es(cov);
        Eigen::Matrix3cd D = es.eigenvalues().asDiagonal();
        Eigen::Matrix3cd V = es.eigenvectors();

        D(0, 0) = std::sqrt(D(0, 0));
        D(1, 1) = std::sqrt(D(1, 1));
        D(2, 2) = std::sqrt(D(2, 2));
        Eigen::Matrix3cd temp = V * D * V.inverse();

        param[3] = temp(0, 0).real(); param[4] = temp(0, 1).real(); param[5] = temp(0, 2).real();
        param[6] = temp(1, 1).real(); param[7] = temp(1, 2).real(); param[8] = temp(2, 2).real();
    }
}

template<typename T>
void ZD::CFiber<T>::ParameterizeFSR(CPoint<T, 9>& param, const int length)
{
    param[0] = this->m_seedDir[0];
    param[1] = this->m_seedDir[1];
    param[2] = this->m_seedDir[2];

    int fc, bc;
    if (length < 0) {
        fc = m_fCount - 1;
        bc = m_bCount - 1;
    }
    else {
        bc = length < (m_bCount - 1) ? length : (m_bCount - 1);
        fc = length < (m_fCount - 1) ? length : (m_fCount - 1);
    }

    param[3] = this->m_pForward[fc][0];
    param[4] = this->m_pForward[fc][1];
    param[5] = this->m_pForward[fc][2];

    param[6] = this->m_pBackward[bc][0];
    param[7] = this->m_pBackward[bc][1];
    param[8] = this->m_pBackward[bc][2];
}

template <typename T>
void ZD::CFiber<T>::ParameterizeResample(CPoint<T, 3> *param, const int length, const int k)
{
    if ((m_fCount + m_bCount) == 2) {
        for (int i = 0; i < k; ++i) {
            param[i] = m_pForward[0];
        }
    }
    else {
        int fc, bc;
        if (length < 0) {
            fc = m_fCount - 1;
            bc = m_bCount - 1;
        }
        else {
            bc = length < (m_bCount - 1) ? length : (m_bCount - 1);
            fc = length < (m_fCount - 1) ? length : (m_fCount - 1);
        }

        T length = 0.0;
        for (int j = 1; j <= bc; ++j) {
            CPoint<T, 3> vec = m_pBackward[j] - m_pBackward[j - 1];
            length += vec.Length();
        }
        for (int j = 1; j <= fc; ++j) {
            CPoint<T, 3> vec = m_pForward[j] - m_pForward[j - 1];
            length += vec.Length();
        }

        T stepSize = length / T(k - 1);
        param[0] = m_pBackward[bc];
        int id = 0;
        T remainLength = stepSize;
        ZD::CPoint<T, 3> lastPoint = param[0];
        for (int j = bc; j > 0; --j) {
            CPoint<T, 3> vec = m_pBackward[j-1] - lastPoint;
            T thisLength = vec.Length();
            if (thisLength < remainLength)
            {
                remainLength -= thisLength;
                lastPoint = m_pBackward[j - 1];
            }
            else
            {
                T factor = remainLength / thisLength;
                param[++id] = factor * m_pBackward[j - 1] + (1.0 - factor) * lastPoint;
                lastPoint = param[id];
                remainLength = stepSize;
                ++j;
            }
        }

        for (int j = 0; j < fc; ++j) {
            CPoint<T, 3> vec = m_pForward[j + 1] - lastPoint;
            T thisLength = vec.Length();
            if (thisLength < remainLength)
            {
                remainLength -= thisLength;
                lastPoint = m_pForward[j + 1];
            }
            else
            {
                T factor = remainLength / thisLength;
                param[++id] = factor * m_pForward[j + 1] + (1.0 - factor) * lastPoint;
                lastPoint = param[id];
                remainLength = stepSize;
                --j;
            }
        }
        param[k - 1] = m_pForward[fc];
    }
}

template<typename T>
ZD::CPoint<T, 2> ZD::CFiber<T>::ParameterizeScalar(CField3<T, 1> *scalar, const int length)
{
    int fc, bc;
    if (length < 0) {
        fc = m_fCount - 1;
        bc = m_bCount - 1;
    }
    else {
        bc = length < (m_bCount - 1) ? length : (m_bCount - 1);
        fc = length < (m_fCount - 1) ? length : (m_fCount - 1);
    }

    T *temp = new T[bc + fc + 1];
    for (int j = 0; j <= bc; ++j) {
        temp[j] = scalar->GetValue(m_pBackward[j])[0];
    }
    for (int j = 1; j <= fc; ++j) {
        temp[j+bc] = scalar->GetValue(m_pForward[j])[0];
    }
    CPoint<T, 2> p = CParameterizationTool<T, 2>::ParameterizeSingleScalar(temp, bc+fc+1);
    SafeDeleteArray(temp);
    return p;
}

template<typename T>
const ZD::CPoint<T, 3> ZD::CFiber<T>::GetOrientation() const
{
    return m_seedDir;
}


#endif
