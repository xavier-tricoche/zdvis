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
#ifndef _ZD_LIB_FTLE_TOOL_HPP_
#define _ZD_LIB_FTLE_TOOL_HPP_

#include "ZD_EigenTool.hpp"
#include "../Base/ZD_Point.hpp"

#include <Eigen/Eigenvalues> 

namespace ZD {
    CPoint<int, 2> neighbors_D2N4[4] = {
        ZD::CPoint<int, 2>(-1,  0),    ZD::CPoint<int, 2>( 1,  0),
        ZD::CPoint<int, 2>( 0, -1),    ZD::CPoint<int, 2>( 0,  1)
    };

    CPoint<int, 2> neighbors_D2N9[9] = {
        ZD::CPoint<int, 2>(-1, -1), ZD::CPoint<int, 2>( 0, -1), ZD::CPoint<int, 2>( 1, -1),
        ZD::CPoint<int, 2>(-1,  0), ZD::CPoint<int, 2>( 0,  0), ZD::CPoint<int, 2>( 1,  0),
        ZD::CPoint<int, 2>(-1,  1), ZD::CPoint<int, 2>( 0,  1), ZD::CPoint<int, 2>( 1,  1)
    };

    CPoint<int, 3> neighbors_D3N6[6] = {
        ZD::CPoint<int, 3>(-1,  0,  0), ZD::CPoint<int, 3>(1, 0, 0),
        ZD::CPoint<int, 3>( 0, -1,  0), ZD::CPoint<int, 3>(0, 1, 0),
        ZD::CPoint<int, 3>( 0,  0, -1), ZD::CPoint<int, 3>(0, 0, 1)
    };

    CPoint<int, 3> neighbors_D3N27[27] = {
        ZD::CPoint<int, 3>(-1, -1, -1), ZD::CPoint<int, 3>(0, -1, -1), ZD::CPoint<int, 3>(1, -1, -1),
        ZD::CPoint<int, 3>(-1,  0, -1), ZD::CPoint<int, 3>(0,  0, -1), ZD::CPoint<int, 3>(1,  0, -1),
        ZD::CPoint<int, 3>(-1,  1, -1), ZD::CPoint<int, 3>(0,  1, -1), ZD::CPoint<int, 3>(1,  1, -1),
        ZD::CPoint<int, 3>(-1, -1,  0), ZD::CPoint<int, 3>(0, -1,  0), ZD::CPoint<int, 3>(1, -1,  0),
        ZD::CPoint<int, 3>(-1,  0,  0), ZD::CPoint<int, 3>(0,  0,  0), ZD::CPoint<int, 3>(1,  0,  0),
        ZD::CPoint<int, 3>(-1,  1,  0), ZD::CPoint<int, 3>(0,  1,  0), ZD::CPoint<int, 3>(1,  1,  0),
        ZD::CPoint<int, 3>(-1, -1,  1), ZD::CPoint<int, 3>(0, -1,  1), ZD::CPoint<int, 3>(1, -1,  1),
        ZD::CPoint<int, 3>(-1,  0,  1), ZD::CPoint<int, 3>(0,  0,  1), ZD::CPoint<int, 3>(1,  0,  1),
        ZD::CPoint<int, 3>(-1,  1,  1), ZD::CPoint<int, 3>(0,  1,  1), ZD::CPoint<int, 3>(1,  1,  1)
    };

    double kernel_2_x[] = {
         -3.0 / 32.0, 0.0,  3 / 32.0,
        -10.0 / 32.0, 0.0, 10 / 32.0,
         -3.0 / 32.0, 0.0,  3 / 32.0
    };

    double kernel_2_y[] = {
        -3.0 / 32.0, -10.0 / 32.0, -3 / 32.0,
                0.0,          0.0,       0.0,
         3.0 / 32.0,  10.0 / 32.0,  3 / 32.0,
    };

    double kernel_3_x[] = {
        -1.0 / 42.0, 0.0 / 42.0, 1.0 / 42.0,
        -3.0 / 42.0, 0.0 / 42.0, 3.0 / 42.0,
        -1.0 / 42.0, 0.0 / 42.0, 1.0 / 42.0,
        -3.0 / 42.0, 0.0 / 42.0, 3.0 / 42.0,
        -6.0 / 42.0, 0.0 / 42.0, 6.0 / 42.0,
        -3.0 / 42.0, 0.0 / 42.0, 3.0 / 42.0,
        -1.0 / 42.0, 0.0 / 42.0, 1.0 / 42.0,
        -3.0 / 42.0, 0.0 / 42.0, 3.0 / 42.0,
        -1.0 / 42.0, 0.0 / 42.0, 1.0 / 42.0,
    };

    double kernel_3_y[] = {
        -1.0 / 42.0, -3.0 / 42.0, -1.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         1.0 / 42.0,  3.0 / 42.0,  1.0 / 42.0,
        -3.0 / 42.0, -6.0 / 42.0, -3.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         3.0 / 42.0,  6.0 / 42.0,  3.0 / 42.0,
        -1.0 / 42.0, -3.0 / 42.0, -1.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         1.0 / 42.0,  3.0 / 42.0,  1.0 / 42.0,
    };

    double kernel_3_z[] = {
        -1.0 / 42.0, -3.0 / 42.0, -1.0 / 42.0,
        -3.0 / 42.0, -6.0 / 42.0, -3.0 / 42.0,
        -1.0 / 42.0, -3.0 / 42.0, -1.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         0.0 / 42.0,  0.0 / 42.0,  0.0 / 42.0,
         1.0 / 42.0,  3.0 / 42.0,  1.0 / 42.0,
         3.0 / 42.0,  6.0 / 42.0,  3.0 / 42.0,
         1.0 / 42.0,  3.0 / 42.0,  1.0 / 42.0,
    };

    template <typename T, unsigned int N>
    class CFTLETool {
    public:
        static T ComputeSigma_D2N9(const CPoint<T, N> *pts, const T *dis) {
            Eigen::MatrixXd mat;
            mat.resize(N, 2);
            mat.setZero();
            for (unsigned int i = 0; i < N; ++i) {
                for (int k = 0; k < 9; ++k) {
                    mat(i, 0) += kernel_2_x[k] * pts[k][i];
                    mat(i, 1) += kernel_2_y[k] * pts[k][i];
                }
                mat(i, 0) /= dis[0];
                mat(i, 1) /= dis[1];
            }

            Eigen::Matrix2d matCG;
            matCG = mat.transpose() * mat;
            T tmp = CEigenTool<T>::MaxEval2D(matCG);
            T sigma = std::sqrt(tmp);
            return sigma;
        }

        static T ComputeFTLE_D2N9(const CPoint<T, N> *pts, const T delta, const T *dis) {
            Eigen::MatrixXd mat;
            mat.resize(N, 2);
            mat.setZero();
            for (unsigned int i = 0; i < N; ++i) {
                for (int k = 0; k < 9; ++k) {
                    mat(i, 0) += kernel_2_x[k] * pts[k][i];
                    mat(i, 1) += kernel_2_y[k] * pts[k][i];
                }
                mat(i, 0) /= dis[0];
                mat(i, 1) /= dis[1];
            }

            Eigen::Matrix2d matCG;
            matCG = mat.transpose() * mat;
            //Eigen::EigenSolver<Eigen::Matrix2d> es(matCG);
            //T tmp = std::max(es.eigenvalues()(0).real(), es.eigenvalues()(1).real());
            T tmp = CEigenTool<T>::MaxEval2D(matCG);
            T ftle = std::log(std::sqrt(tmp)) / std::fabs(delta);
            return ftle;
        }

        static T ComputeFTLE_D2N4(const CPoint<T, N> *pts, const T delta, const T *dis) {
            Eigen::MatrixXd mat;
            mat.resize(N, 2);
            mat.setZero();
            for (unsigned int i = 0; i < N; ++i) {
                mat(i, 0) = (pts[1][i] - pts[0][i]) / dis[0];
                mat(i, 1) = (pts[3][i] - pts[2][i]) / dis[1];
            }

            Eigen::Matrix2d matCG;
            matCG = mat.transpose() * mat;
            //Eigen::EigenSolver<Eigen::Matrix2d> es(matCG);
            //T tmp = std::max(es.eigenvalues()(0).real(), es.eigenvalues()(1).real());
            T tmp = CEigenTool<T>::MaxEval2D(matCG);
            T ftle = std::log(std::sqrt(tmp)) / std::fabs(delta);
            return ftle;
        }

        static T ComputeFTLE_D3N6(const CPoint<T, N> *pts, const T delta, const T *dis) {
            Eigen::MatrixXd mat;
            mat.resize(N, 3);
            mat.setZero();
            for (unsigned int i = 0; i < N; ++i) {
                mat(i, 0) = (pts[1][i] - pts[0][i]) / dis[0];
                mat(i, 1) = (pts[3][i] - pts[2][i]) / dis[1];
                mat(i, 2) = (pts[5][i] - pts[4][i]) / dis[2];
            }

            Eigen::Matrix3d matCG;
            matCG = mat.transpose() * mat;
            Eigen::EigenSolver<Eigen::Matrix3d> es(matCG);
            T tmp = std::max(es.eigenvalues()(0).real(), std::max(es.eigenvalues()(1).real(), es.eigenvalues()(2).real()));
            T ftle = std::log(std::sqrt(tmp)) / std::fabs(delta);
            return ftle;
        }

        static T ComputeFTLE_D3N27(const CPoint<T, N> *pts, const T delta, const T *dis) {
            Eigen::MatrixXd mat;
            mat.resize(N, 3);
            mat.setZero();
            for (unsigned int i = 0; i < N; ++i) {
                for (int k = 0; k < 27; ++k) {
                    mat(i, 0) += kernel_3_x[k] * pts[k][i];
                    mat(i, 1) += kernel_3_y[k] * pts[k][i];
                    mat(i, 2) += kernel_3_z[k] * pts[k][i];
                }
                mat(i, 0) /= dis[0];
                mat(i, 1) /= dis[1];
                mat(i, 2) /= dis[2];
            }

            Eigen::Matrix3d matCG;
            matCG = mat.transpose() * mat;
            Eigen::EigenSolver<Eigen::Matrix3d> es(matCG);
            T tmp = std::max(es.eigenvalues()(0).real(), std::max(es.eigenvalues()(1).real(), es.eigenvalues()(2).real()));
            T ftle = std::log(std::sqrt(tmp)) / std::fabs(delta);
            return ftle;
        }
    };
}



#endif
