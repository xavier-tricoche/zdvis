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
#ifndef _ZD_LIB_DTI_HPP_
#define _ZD_LIB_DTI_HPP_

#include "utils/Base/ZD_Field.hpp"
#include "utils/define.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"

#include <Eigen/Eigenvalues>

#include <teem/ten.h>


namespace ZD {
    /* flag, xx, xy, xz, yy, yz, zz */
    template<typename T>
    class CDTI {
    private:
        CField3<T, 7> *m_pDTIVolume;
        CField3<T, 1> *m_pFAVolume;

        const T m_angleThreshold;
        const T m_faThreshold;

        CPoint<T, 3> *m_pMLSPoints;
        T *m_pMLSWeights;
        const int m_MLSPointCount;
        const T m_MLSSize;
        
    public:
        const CPoint<T, 7> * GetDTIVolume() const;
        const CPoint<T, 1> * GetFAVolume() const;

        const CPoint<T, 3> * GetMLSPoints() const;
        const T * GetMLSWeights() const;

    public:
        CDTI(const char *dtiPathname, const char *faPathname, 
            const T at, const T ft, const int mlsCount, const T mlsSize);
        CDTI(const char *dtiPathname, const char *faPathname,
            const T at, const T ft);
        ~CDTI();

    public:
        const int *GetSize() const;
        CPoint<T, 7> GetDTI(const int x, const int y, const int z) const;
        //CPoint<T, 7> GetDTI(const T x, const T y, const T z) const;
        CPoint<T, 7> GetDTI(const CPoint<T, 3>& pos) const;

        CPoint<T, 7> GetDTIMLS(const CPoint<T, 3>& pos, const CPoint<T, 7>& kernel) const;

        T GetFA(const CPoint<T, 3>& pos) const;
        //CPoint<T, 1> GetFA(const CPoint<T, 3>& pos) const;
        T GetFAMLS(const CPoint<T, 3>& pos) const;
        //CPoint<T, 1> GetFAMLS(const CPoint<T, 3>& pos) const;



    public:
        inline bool NextLength(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, 
            const T length, const T stepSize, const bool normalize = false) const;
        inline bool NextLengthMLS(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, 
            const T length, CPoint<T, 7>& kernel, const T stepSize, const bool normalize = false) const;
        
        inline bool Next(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, const T step_size, const bool normalize = false) const;
        inline bool NextMLS(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, CPoint<T, 7>& kernel, const T step_size, const bool normalize = false) const;

    public:
        static T DTI2FA(const T *dti);
        static T DTI2CL(const T *dti);
        static void DTI2Eigens(const T *dti, T *evals, CPoint<T, 3> *evecs);
        static CPoint<T, 3> DTI2Direction(const T *dti);

    private:
        inline bool GetDirection(const CPoint<T, 3>& pos, 
            const CPoint<T, 3>& last_dir, CPoint<T, 3>& dir) const;
        inline bool GetDirectionMLS(const CPoint<T, 3>& pos, 
            const CPoint<T, 3>& last_dir, CPoint<T, 3>& dir, const CPoint<T, 7>& kernel) const;

        void CreateMLSSamples();

        void ComputeFAVolume();
        void ComputeFAVolumeMLS();
    };
}

template<typename T>
ZD::CDTI<T>::CDTI(const char *dtiPathname, const char *faPathname,
    const T at, const T ft, const int count, const T size)
    : m_angleThreshold(at), m_faThreshold(ft), m_MLSPointCount(count), m_MLSSize(size)
{
    m_pDTIVolume = new CField3<T, 7>();
    m_pDTIVolume->OpenNrrdFile(dtiPathname);
    
    m_pMLSPoints = nullptr;
    m_pMLSWeights = nullptr;

    if (m_MLSPointCount > 0)
        CreateMLSSamples();

    m_pFAVolume = new CField3<T, 1>();
    if (faPathname == nullptr) {
        if (m_MLSPointCount > 0)
            ComputeFAVolumeMLS();
        else
            ComputeFAVolume();
    }
    else
        m_pFAVolume->OpenNrrdFile(faPathname);
}

template<typename T>
ZD::CDTI<T>::CDTI(const char *dtiPathname, const char *faPathname,
    const T at, const T ft)
    : m_angleThreshold(at), m_faThreshold(ft), m_MLSPointCount(0), m_MLSSize(0.0)
{
    m_pDTIVolume = new CField3<T, 7>();
    m_pDTIVolume->OpenNrrdFile(dtiPathname);

    m_pMLSPoints = nullptr;
    m_pMLSWeights = nullptr;

    m_pFAVolume = new CField3<T, 1>();
    if (faPathname == nullptr)
        ComputeFAVolume();
    else
        m_pFAVolume->OpenNrrdFile(faPathname);
}


template<typename T>
ZD::CDTI<T>::~CDTI()
{
    SafeDelete(this->m_pDTIVolume);
    SafeDelete(this->m_pFAVolume);

    SafeDeleteArray(m_pMLSPoints);
    SafeDeleteArray(m_pMLSWeights);
}

template <typename T>
const ZD::CPoint<T, 7> * ZD::CDTI<T>::GetDTIVolume() const
{
    return this->m_pDTIVolume->GetData();
}

template <typename T>
const ZD::CPoint<T, 1> * ZD::CDTI<T>::GetFAVolume() const
{
    return this->m_pFAVolume->GetData();
}

template <typename T>
const ZD::CPoint<T, 3> * ZD::CDTI<T>::GetMLSPoints() const
{
    return this->m_pMLSPoints;
}

template <typename T>
const T * ZD::CDTI<T>::GetMLSWeights() const
{
    return this->m_pMLSWeights;
}


template<typename T>
const int * ZD::CDTI<T>::GetSize() const
{
    return this->m_pDTIVolume->GetSize();
}

template<typename T>
ZD::CPoint<T, 7> ZD::CDTI<T>::GetDTI(const int x, const int y, const int z) const
{
    return this->m_pDTIVolume->GetValue(x, y, z);
}

template<typename T>
ZD::CPoint<T, 7> ZD::CDTI<T>::GetDTI(const CPoint<T, 3>& pos) const
{
    return this->m_pDTIVolume->GetValue(pos);
}

template<typename T>
ZD::CPoint<T, 7> ZD::CDTI<T>::GetDTIMLS(const CPoint<T, 3>& pos, const CPoint<T, 7>& kernel) const
{
    assert(m_pMLSPoints != nullptr && m_pMLSWeights != nullptr);

    T m[9];
    TEN_T2M(m, kernel.m_data);

    CPoint<T, 3> *points = new CPoint<T, 3>[m_MLSPointCount];
    T *length = new T[m_MLSPointCount];
    T max_length = 0.0;
    for (int i = 0; i < m_MLSPointCount; ++i) {
        points[i][0] = m_pMLSPoints[i][0] * m[0] + m_pMLSPoints[i][1] * m[1] + m_pMLSPoints[i][2] * m[2];
        points[i][1] = m_pMLSPoints[i][0] * m[3] + m_pMLSPoints[i][1] * m[4] + m_pMLSPoints[i][2] * m[5];
        points[i][2] = m_pMLSPoints[i][0] * m[6] + m_pMLSPoints[i][1] * m[7] + m_pMLSPoints[i][2] * m[8];
        length[i] = std::sqrt(InnerProduct(points[i], points[i]));
        if (length[i] > max_length)
            max_length = length[i];
    }

    if (max_length > ZD_EPSILON) {
        for (int i = 0; i < m_MLSPointCount; ++i) {
            points[i][0] = points[i][0] / max_length * m_MLSSize + pos[0];
            points[i][1] = points[i][1] / max_length * m_MLSSize + pos[1];
            points[i][2] = points[i][2] / max_length * m_MLSSize + pos[2];
        }
    }

    CPoint<T, 7> smoothed_dti = CPoint<T, 7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    T weight = 0.0;
    for (int i = 0; i < m_MLSPointCount; ++i) {
        if (GetFA(points[i]) < m_faThreshold) {
            continue;
        }
        CPoint<T, 7> temp = GetDTI(points[i]);
        smoothed_dti += m_pMLSWeights[i] * temp;
        weight += m_pMLSWeights[i];
    }

    SafeDeleteArray(points);
    SafeDeleteArray(length);

    if (weight < ZD_EPSILON) {
        smoothed_dti = GetDTI(pos);
    }
    else {
        smoothed_dti /= weight;
    }

    return smoothed_dti;
}

template<typename T>
T ZD::CDTI<T>::GetFA(const CPoint<T, 3>& pos) const
{
    return (this->m_pFAVolume->GetValue(pos))[0];
}

template<typename T>
T ZD::CDTI<T>::GetFAMLS(const CPoint<T, 3>& pos) const
{
    assert(m_pMLSPoints != nullptr && m_pMLSWeights != nullptr);
    
    CPoint<T, 7> kernel = this->GetDTI(pos);
    
    T m[9];
    TEN_T2M(m, kernel.m_data);

    CPoint<T, 3> *points = new CPoint<T, 3>[m_MLSPointCount];
    T *length = new T[m_MLSPointCount];
    T max_length = 0.0;
    for (int i = 0; i < m_MLSPointCount; ++i) {
        points[i][0] = m_pMLSPoints[i][0] * m[0] + m_pMLSPoints[i][1] * m[1] + m_pMLSPoints[i][2] * m[2];
        points[i][1] = m_pMLSPoints[i][0] * m[3] + m_pMLSPoints[i][1] * m[4] + m_pMLSPoints[i][2] * m[5];
        points[i][2] = m_pMLSPoints[i][0] * m[6] + m_pMLSPoints[i][1] * m[7] + m_pMLSPoints[i][2] * m[8];
        length[i] = std::sqrt(InnerProduct(points[i], points[i]));
        if (length[i] > max_length)
            max_length = length[i];
    }

    if (max_length > ZD_EPSILON) {
        for (int i = 0; i < m_MLSPointCount; ++i) {
            points[i][0] = points[i][0] / max_length * m_MLSSize + pos[0];
            points[i][1] = points[i][1] / max_length * m_MLSSize + pos[1];
            points[i][2] = points[i][2] / max_length * m_MLSSize + pos[2];
        }
    }

    T smoothed_fa = 0.0;
    T weight = 0.0;
    for (int i = 0; i < m_MLSPointCount; ++i) {
        T temp = this->GetFA(points[i]);
        smoothed_fa += m_pMLSWeights[i] * temp;
        weight += m_pMLSWeights[i];
    }

    SafeDeleteArray(points);
    SafeDeleteArray(length);
    
    smoothed_fa /= weight;

    return smoothed_fa;
}

template<typename T>
inline bool ZD::CDTI<T>::NextLength(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, 
    const T length, const T stepSize, const bool normalize) const
{
    T left_length = length;
    while (left_length > 0.0) {
        T ss = left_length < stepSize ? left_length : stepSize;
        if (Next(pos, last_dir, ss, normalize) == false)
            return false;
        else
            left_length -= ss;
    }
    return true;
}

template<typename T>
inline bool ZD::CDTI<T>::NextLengthMLS(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, 
    const T length, CPoint<T, 7>& kernel, const T stepSize, const bool normalize) const
{
    T left_length = length;
    while (left_length > 0.0) {
        T ss = left_length < stepSize ? left_length : stepSize;
        if (NextMLS(pos, last_dir, kernel, ss, normalize) == false)
            return false;
        else
            left_length -= ss;
    }
    return true;
}

template<typename T>
inline bool ZD::CDTI<T>::Next(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, const T step_size, const bool normalize) const
{
    CPoint<T, 3> k1, k2, temp;
    T ss;

    if (normalize) {
        T fa = this->GetFA(pos);
        ss = fa * step_size;
    }
    else {
        ss = step_size;
    }

    if (GetDirection(pos, last_dir, k1) == false)
        return false;

    temp = pos + k1 * (ss * 0.5);
    if (GetDirection(temp, last_dir, k2) == false)
        return false;

    last_dir = k2;                    // update last_dir
    pos = pos + ss * last_dir;        // update pos
    return true;
}

template<typename T>
inline bool ZD::CDTI<T>::NextMLS(CPoint<T, 3>& pos, CPoint<T, 3>& last_dir, CPoint<T, 7>& kernel, const T step_size, const bool normalize) const
{
    CPoint<T, 3> k1, k2, temp;
    T ss;

    if (normalize) {
        T fa = this->GetFAMLS(pos);
        ss = fa * step_size;
    }
    else {
        ss = step_size;
    }

    if (GetDirectionMLS(pos, last_dir, k1, kernel) == false)
        return false;

    temp = pos + k1 * (ss * 0.5);
    if (GetDirectionMLS(temp, last_dir, k2, kernel) == false)
        return false;

    last_dir = k2;                        // update last_dir
    pos = pos + ss * last_dir;            // update pos
    kernel = GetDTIMLS(pos, kernel);    // update kernel
    return true;
}


template<typename T>
void ZD::CDTI<T>::DTI2Eigens(const T *dti, T *evals, CPoint<T, 3> *evecs)
{
    Eigen::Matrix3d mat;
    mat(0, 0) = dti[1]; mat(0, 1) = dti[2]; mat(0, 2) = dti[3];
    mat(1, 0) = dti[2]; mat(1, 1) = dti[4]; mat(1, 2) = dti[5];
    mat(2, 0) = dti[3]; mat(2, 1) = dti[5]; mat(2, 2) = dti[6];

    Eigen::EigenSolver<Eigen::Matrix3d> es(mat);
    
    T temp[3];
    int maxID, midID, minID;
    
    temp[0] = es.eigenvalues()(0).real();
    temp[1] = es.eigenvalues()(1).real();
    temp[2] = es.eigenvalues()(2).real();
    if (temp[0] > temp[1] && temp[0] > temp[2] && temp[1] > temp[2]) {
        maxID = 0; midID = 1; minID = 2;
    } else if (temp[0] > temp[1] && temp[0] > temp[2] && temp[2] > temp[1]) {
        maxID = 0; midID = 2; minID = 1;
    } else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[0] > temp[2]) {
        maxID = 1; midID = 0; minID = 2;
    } else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[2] > temp[0]) {
        maxID = 1; midID = 2; minID = 0;
    } else if (temp[2] > temp[0] && temp[2] > temp[1] && temp[0] > temp[1]) {
        maxID = 2; midID = 0; minID = 1;
    } else {
        maxID = 2; midID = 1; minID = 0;
    }

    if (evals != nullptr) {
        evals[0] = es.eigenvalues()(maxID).real();
        evals[1] = es.eigenvalues()(midID).real();
        evals[2] = es.eigenvalues()(minID).real();
    }
    if (evecs != nullptr) {
        evecs[0] = CPoint<T, 3>(es.eigenvectors()(0, maxID).real(), es.eigenvectors()(1, maxID).real(), es.eigenvectors()(2, maxID).real());
        evecs[1] = CPoint<T, 3>(es.eigenvectors()(0, midID).real(), es.eigenvectors()(1, midID).real(), es.eigenvectors()(2, midID).real());
        evecs[2] = CPoint<T, 3>(es.eigenvectors()(0, minID).real(), es.eigenvectors()(1, minID).real(), es.eigenvectors()(2, minID).real());
    }
}

template<typename T>
T ZD::CDTI<T>::DTI2FA(const T *dti)
{
    T t[6];
    t[0] = dti[1]; t[1] = dti[2]; t[2] = dti[3]; t[3] = dti[4]; t[4] = dti[5]; t[5] = dti[6];

    T mean = (t[0] + t[1] + t[2] + t[3] + t[4] + t[5]) / 6.0;
    if (mean == 0.0)
        return 0.0;

    t[0] = t[0] / mean;
    t[1] = t[1] / mean;
    t[2] = t[2] / mean;
    t[3] = t[3] / mean;
    t[4] = t[4] / mean;
    t[5] = t[5] / mean;

    T cross = t[1] * t[1] + t[2] * t[2] + t[4] * t[4];
    T j2 = t[0] * t[3] + t[3] * t[5] + t[5] * t[0] - cross;
    T j4 = t[0] * t[0] + t[3] * t[3] + t[5] * t[5] + 2.0 * cross;
    T fa = sqrt((j4 - j2) / j4);
    return fa;
}

template<typename T>
T ZD::CDTI<T>::DTI2CL(const T *dti)
{
    T *evals;
    DTI2Eigens(dti, evals, nullptr);

    T cl = (evals[0] - evals[1]) / (evals[0] + evals[1] + evals[2]);
    return cl;
}

template<typename T>
ZD::CPoint<T, 3> ZD::CDTI<T>::DTI2Direction(const T *dti)
{
    CPoint<T, 3> evecs[3];
    DTI2Eigens(dti, nullptr, evecs);
    return evecs[0];
}


template<typename T>
inline bool ZD::CDTI<T>::GetDirection(const CPoint<T, 3>& pos, const CPoint<T, 3>& last_dir, CPoint<T, 3>& dir) const
{
    CPoint<T, 7> dti = GetDTI(pos);
    //if (DTI2FA(dti) < m_faThreshold)
    //    return false;
    if (GetFA(pos) < m_faThreshold)
        return false;
    dir = DTI2Direction(dti.m_data);
    T ipv = InnerProduct(dir, last_dir);
    if (std::abs(ipv) < m_angleThreshold)
        return false;
    if (ipv < 0.0)
        dir = -dir;
    return true;
}

template<typename T>
inline bool ZD::CDTI<T>::GetDirectionMLS(const CPoint<T, 3>& pos, 
    const CPoint<T, 3>& last_dir, CPoint<T, 3>& dir, const CPoint<T, 7>& kernel) const
{
    CPoint<T, 7> dti = GetDTIMLS(pos, kernel);
    //if (DTI2FA(dti) < m_faThreshold)
    //    return false;
    if (GetFAMLS(pos) < m_faThreshold)
        return false;
    dir = DTI2Direction(dti.m_data);
    T ipv = InnerProduct(dir, last_dir);
    if (std::abs(ipv) < m_angleThreshold)
        return false;
    if (ipv < 0.0)
        dir = -dir;
    return true;
}

template<typename T>
void ZD::CDTI<T>::CreateMLSSamples()
{
    m_pMLSPoints = new CPoint<T, 3>[m_MLSPointCount];
    m_pMLSWeights = new T[m_MLSPointCount];

    m_pMLSPoints[0] = CPoint<T, 3>(0.0, 0.0, 0.0);
    m_pMLSWeights[0] = 1.0;

    CRandomTool<T> rt;

    for (int i = 1; i < m_MLSPointCount; ++i) {
        bool flag = false;
        CPoint<T, 3> pt;
        while (!flag) {
            T a = rt.GetRandomNumer() * 2.0 * ZD_PI;
            T b = rt.GetRandomNumer() * ZD_PI;
            T r = rt.GetRandomNumer();

            pt[0] = r * sin(a) * sin(b);
            pt[1] = r * sin(a) * cos(b);
            pt[2] = r * cos(a);

            flag = true;
            for (int j = 0; j < i; ++j) {
                T dis = Distance(m_pMLSPoints[j], pt);
                if (dis < 0.3) {
                    flag = false;
                    break;
                }
            }
        }

        m_pMLSPoints[i] = pt;
        m_pMLSWeights[i] = InnerProduct(m_pMLSPoints[i], m_pMLSPoints[i]) * 4.0;
        m_pMLSWeights[i] = exp(m_pMLSWeights[i] / -2.5);
    }
}

template<typename T>
void ZD::CDTI<T>::ComputeFAVolume()
{
    const int *size = m_pDTIVolume->GetSize();
    m_pFAVolume->CreateField(size);

    for (int z = 0; z < size[2]; ++z) {
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int y = 0; y < size[1]; ++y) {
            for (int x = 0; x < size[0]; ++x) {
                CPoint<T, 3> pos = CPoint<T, 3>((T)(x), (T)(y), (T)(z));
                CPoint<T, 7> tensor = GetDTI(pos);
                m_pFAVolume->SetValue(x, y, z, CPoint<T, 1>(DTI2FA(tensor.m_data)));
            }
        }
    }
}

template<typename T>
void ZD::CDTI<T>::ComputeFAVolumeMLS()
{
    const int *size = m_pDTIVolume->GetSize();
    m_pFAVolume->CreateField(size);

    for (int z = 0; z < size[2]; ++z) {
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int y = 0; y < size[1]; ++y) {
            for (int x = 0; x < size[0]; ++x) {
                CPoint<T, 7> kernel = ZD::CPoint<T, 7>(1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
                CPoint<T, 3> pos = CPoint<T, 3>((T)(x), (T)(y), (T)(z));
                CPoint<T, 7> tensor = GetDTIMLS(pos, kernel);
                m_pFAVolume->SetValue(x, y, z, CPoint<T, 1>(DTI2FA(tensor.m_data)));
            }
        }
    }
}


#endif
