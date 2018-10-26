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
#ifndef _ZD_LIB_TENSOR_HPP_
#define _ZD_LIB_TENSOR_HPP_

#include "utils/Base/ZD_Field.hpp"
#include "utils/DTI_HARDI/ZD_Fiber.hpp"
#include "utils/define.hpp"

#include <vector>
#include <Eigen/Eigenvalues>

namespace ZD {
    template <typename T>
    using CTensorLine = CFiber<T>;

    template <typename T>
    class CTensor {
    public:
        CTensor();
        ~CTensor();

    public:
        inline CPoint<T, 3> Direction(const CPoint<T, 3>& p, const CPoint<T, 3>& v) const;
        virtual inline CPoint<T, 3> Direction(const CPoint<T, 3>& p) const = 0;
    
    public:
        inline void IntegrateTensorLine(const CPoint<T, 3>& seed, 
            const T length, const T stepSize, CTensorLine<T> *tensorline) const;

        inline bool NextLength(CPoint<T, 3>& p, CPoint<T, 3>& v, const T length, const T stepSize) const;

    protected:
        inline bool Next(CPoint<T, 3>& p, CPoint<T, 3>& v, const T stepSize) const;
    };

    template <typename T>
    class CTensor7 : public CTensor<T> {
    protected:
        CField3<T, 7> *m_pTensorData;
    public:
        CTensor7(const char *pathname);
        ~CTensor7();
        virtual inline CPoint<T, 3> Direction(const CPoint<T, 3>& p) const;
    };

    template <typename T>
    class CTensor9 : public CTensor<T> {
    protected:
        CField3<T, 9> *m_pTensorData;
    public:
        CTensor9(const char *pathname);
        ~CTensor9();
        virtual inline CPoint<T, 3> Direction(const CPoint<T, 3>& p) const;
    };
}

template<typename T>
ZD::CTensor<T>::CTensor()
{
}

template<typename T>
ZD::CTensor<T>::~CTensor()
{
}

template<typename T>
inline void ZD::CTensor<T>::IntegrateTensorLine(const CPoint<T, 3>& seed,
    const T totalLength, const T stepSize, CTensorLine<T> *tensorline) const
{
    CPoint<T, 3> seed_dir = Direction(seed);

    CPoint<T, 3> p, dir;

    const T deltaLength = 1.0;

    // forward
    p = seed;
    dir = seed_dir;
    std::vector<CPoint<T, 3>> forward;
    forward.push_back(p);
    for (T length = 0.0; length < totalLength * (1.0 - ZD_EPSILON); length += deltaLength) {
        if (NextLength(p, dir, deltaLength, stepSize) == true)
            forward.push_back(p);
        else
            break;
    }

    // backward
    p = seed;
    dir = -seed_dir;
    std::vector<CPoint<T, 3>> backward;
    backward.push_back(p);
    for (T length = 0.0; length < totalLength * (1.0 - ZD_EPSILON); length += deltaLength) {
        if (NextLength(p, dir, deltaLength, stepSize) == true)
            backward.push_back(p);
        else
            break;
    }

    tensorline->CreateFiber(forward, backward, seed, Direction(seed));
}

template<typename T>
inline bool ZD::CTensor<T>::NextLength(CPoint<T, 3>& p, CPoint<T, 3>& v, const T length, const T stepSize) const
{
    T left_length = length;
    while (left_length > 0.0) {
        T ss = left_length < stepSize ? left_length : stepSize;
        if (Next(p, v, ss) == false)
            return false;
        else
            left_length -= ss;
    }
    return true;
}

template<typename T>
inline bool ZD::CTensor<T>::Next(CPoint<T, 3>& p, CPoint<T, 3>& v, const T stepSize) const
{
    CPoint<T, 3> k1;
    CPoint<T, 3> np;

    k1 = this->Direction(p, v);
    if (std::isnan(k1[0]))
        return false;
    
    np = p + k1 * stepSize;
    
    //v = np - p;
    //v.Normalize();
    v = k1;

    p = np;

    return true;
}

template<typename T>
inline ZD::CPoint<T, 3> ZD::CTensor<T>::Direction(const CPoint<T, 3>& p, const CPoint<T, 3>& v) const
{
    CPoint<T, 3> dir = Direction(p);
    if (InnerProduct(dir, v) < 0.0)
        dir = -dir;
    return dir;
}


template<typename T>
ZD::CTensor7<T>::CTensor7(const char * pathname)
{
    m_pTensorData = new CField3<T, 7>();
    m_pTensorData->OpenNrrdFile(pathname);
}

template<typename T>
ZD::CTensor7<T>::~CTensor7()
{
    SafeDelete(m_pTensorData);
}

template<typename T>
inline ZD::CPoint<T, 3> ZD::CTensor7<T>::Direction(const CPoint<T, 3>& p) const
{
    CPoint<T, 7> tensor;
    tensor = this->m_pTensorData->GetValue(p);
    if (tensor.Length2() < 1e-10) {
        CPoint<T, 3> dir;
        dir.SetNan();
        return dir;
    }

    Eigen::Matrix3d mat;
    mat(0, 0) = tensor[1]; mat(0, 1) = tensor[2]; mat(0, 2) = tensor[3];
    mat(1, 0) = tensor[2]; mat(1, 1) = tensor[4]; mat(1, 2) = tensor[5];
    mat(2, 0) = tensor[3]; mat(2, 1) = tensor[5]; mat(2, 2) = tensor[6];

    Eigen::EigenSolver<Eigen::Matrix3d> es(mat);

    T temp[3];
    int maxID, midID, minID;

    temp[0] = es.eigenvalues()(0).real();
    temp[1] = es.eigenvalues()(1).real();
    temp[2] = es.eigenvalues()(2).real();
    if (temp[0] > temp[1] && temp[0] > temp[2] && temp[1] > temp[2]) {
        maxID = 0; midID = 1; minID = 2;
    }
    else if (temp[0] > temp[1] && temp[0] > temp[2] && temp[2] > temp[1]) {
        maxID = 0; midID = 2; minID = 1;
    }
    else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[0] > temp[2]) {
        maxID = 1; midID = 0; minID = 2;
    }
    else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[2] > temp[0]) {
        maxID = 1; midID = 2; minID = 0;
    }
    else if (temp[2] > temp[0] && temp[2] > temp[1] && temp[0] > temp[1]) {
        maxID = 2; midID = 0; minID = 1;
    }
    else {
        maxID = 2; midID = 1; minID = 0;
    }

    int id = maxID;
    return CPoint<T, 3>(es.eigenvectors()(0, id).real(), es.eigenvectors()(1, id).real(), es.eigenvectors()(2, id).real());
}

template<typename T>
ZD::CTensor9<T>::CTensor9(const char * pathname)
{
    m_pTensorData = new CField3<T, 9>();
    m_pTensorData->OpenNrrdFile(pathname);
}

template<typename T>
ZD::CTensor9<T>::~CTensor9()
{
    SafeDelete(m_pTensorData);
}

template<typename T>
inline ZD::CPoint<T, 3> ZD::CTensor9<T>::Direction(const CPoint<T, 3>& p) const
{
    CPoint<T, 9> tensor;
    tensor = this->m_pTensorData->GetValue(p);
    if (tensor.Length2() < 1e-10) {
        CPoint<T, 3> dir;
        dir.SetNan();
        return dir;
    }

    Eigen::Matrix3d mat;
    mat(0, 0) = tensor[0]; mat(0, 1) = tensor[1]; mat(0, 2) = tensor[2];
    mat(1, 0) = tensor[3]; mat(1, 1) = tensor[4]; mat(1, 2) = tensor[5];
    mat(2, 0) = tensor[6]; mat(2, 1) = tensor[7]; mat(2, 2) = tensor[8];

    Eigen::EigenSolver<Eigen::Matrix3d> es(mat);

    T temp[3];
    int maxID, midID, minID;

    temp[0] = es.eigenvalues()(0).real();
    temp[1] = es.eigenvalues()(1).real();
    temp[2] = es.eigenvalues()(2).real();
    if (temp[0] > temp[1] && temp[0] > temp[2] && temp[1] > temp[2]) {
        maxID = 0; midID = 1; minID = 2;
    }
    else if (temp[0] > temp[1] && temp[0] > temp[2] && temp[2] > temp[1]) {
        maxID = 0; midID = 2; minID = 1;
    }
    else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[0] > temp[2]) {
        maxID = 1; midID = 0; minID = 2;
    }
    else if (temp[1] > temp[0] && temp[1] > temp[2] && temp[2] > temp[0]) {
        maxID = 1; midID = 2; minID = 0;
    }
    else if (temp[2] > temp[0] && temp[2] > temp[1] && temp[0] > temp[1]) {
        maxID = 2; midID = 0; minID = 1;
    }
    else {
        maxID = 2; midID = 1; minID = 0;
    }

    int id = maxID;
    return CPoint<T, 3>(es.eigenvectors()(0, id).real(), es.eigenvectors()(1, id).real(), es.eigenvectors()(2, id).real());
}


#endif
