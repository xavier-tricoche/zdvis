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
#ifndef _ZD_LIB_DOUBLE_POINT_LOAD_HPP
#define _ZD_LIB_DOUBLE_POINT_LOAD_HPP

#include "ZD_Tensor.hpp"
#include "utils/define.hpp"
#include <Eigen/Eigenvalues>

namespace ZD {
    template<typename T>
    class CDoublePointLoad : public CTensor<T> {
    private:
        CVector<T, 3> m_min, m_max;
        CPoint<T, 3> m_centers[2];
        T m_load, m_nu, m_twoPi;

    public:
        CDoublePointLoad(const T w, const T h, const T d, const T l = 1.0);
        ~CDoublePointLoad();

    public:
        virtual inline CPoint<T, 3> Direction(const CPoint<T, 3>& p) const;

    private:
        CPoint<T, 9> GetTensor(const CPoint<T, 3>& p) const;
    };
}

template<typename T>
ZD::CDoublePointLoad<T>::CDoublePointLoad(const T w, const T h, const T d, const T l)
{
    m_min = CVector<T, 3>(-0.5*w, -0.5*h, -d);
    m_max = CVector<T, 3>( 0.5*w,  0.5*h, 0.0);

    m_centers[0] = CPoint<T, 3>(-0.25*w, 0.0, 0.0);
    m_centers[1] = CPoint<T, 3>( 0.25*w, 0.0, 0.0);

    m_twoPi = 2.0 * ZD_PI;
    m_load = l;
    m_nu = 0.4;
}

template<typename T>
ZD::CDoublePointLoad<T>::~CDoublePointLoad()
{

}

template<typename T>
inline ZD::CPoint<T, 3> ZD::CDoublePointLoad<T>::Direction(const CPoint<T, 3>& p) const
{
    CPoint<T, 9> tensor = GetTensor(p);
    if (std::isnan(tensor[0]) == true) {
        CPoint<T, 3> dir;
        dir.SetNan();
        return dir;
    }
    else {
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
}

template<typename T>
ZD::CPoint<T, 9> ZD::CDoublePointLoad<T>::GetTensor(const CPoint<T, 3>& pos) const
{
    T s[3];
    T P, rho, rho2, rho3, rho5;
    T rhoPlusz2, zPlus2rho, txy, txz, tyz;

    CPoint<T, 9> tensor;
    tensor.SetZero();

    for (int n = 0; n < 2; ++n) {
        P = -m_load;
        CPoint<T, 3> dis;
        dis = m_centers[n] - pos;
        dis[0] *= -1.0;
        rho = dis.Length();
        if (rho < 1.0e-10) {
            tensor.SetNan();
            return tensor;            // singular tensor
        }

        rho2 = rho * rho;
        rho3 = rho * rho2;
        rho5 = rho2 * rho3;

        T x = dis[0];
        T y = dis[1];
        T z = dis[2];

        T x2 = x * x;
        T y2 = y * y;
        T z2 = z * z;

        rhoPlusz2 = (rho + z) * (rho + z);
        zPlus2rho = 2.0 * rho + z;

        // normal stresses
        s[0] = P / (m_twoPi * rho2) * (3.0 * z * x2 / rho3 - m_nu * (z / rho - rho / (rho + z) +
            x2 * (zPlus2rho) / (rho * rhoPlusz2)));
        s[1] = P / (m_twoPi * rho2) * (3.0 * z * y2 / rho3 - m_nu * (z / rho - rho / (rho + z) +
            y2 * (zPlus2rho) / (rho * rhoPlusz2)));
        s[2] = 3.0 * P * z2 * z / (m_twoPi * rho5);

        //shear stresses - negative signs are coordinate transformations
        //that is, equations (in text) are in different coordinate system
        //than volume is in.
        txy = -(P / (m_twoPi * rho2) * (3.0 * x * y * z / rho3 -
            m_nu * x * y * (zPlus2rho) / (rho * rhoPlusz2)));
        txz = -(3.0 * P * x * z2 / (m_twoPi * rho5));
        tyz = 3.0 * P * y * z2 / (m_twoPi * rho5);

        tensor[0] += s[0];  // Component(0,0);
        tensor[4] += s[1];  // Component(1,1);
        tensor[8] += s[2];  // Component(2,2);
        tensor[1] += txy;   // Component(0,1);  real symmetric matrix
        tensor[2] += txz;   // Component(0,2);
        tensor[5] += tyz;   // Component(1,2);
    }

    tensor[3] = tensor[1];
    tensor[6] = tensor[2];
    tensor[7] = tensor[5];

    return tensor;
}


#endif
