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
#ifndef _ZD_LIB_LE_TOOL_HPP_
#define _ZD_LIB_LE_TOOL_HPP_

#include "../Base/ZD_Point.hpp"

#include <Eigen/Eigenvalues> 

namespace ZD {
    extern CPoint<int, 2> neighbors_D2N4[4];
    extern CPoint<int, 2> neighbors_D2N9[9];
    extern CPoint<int, 3> neighbors_D3N6[6];
    extern CPoint<int, 3> neighbors_D3N27[27];

    enum ZD_LE_Method {
        ZD_LE_Method_Moment,        /* parameterized by moment */
        ZD_LE_Method_Fixed            /* parameterized by fixed number points */
    };

    template <typename T, unsigned int N>
    class CLETool {
    private:
        static Eigen::MatrixXd mat_27;        // 3D 
        static Eigen::MatrixXd mat_9;        // 2D

        static T dxKernel[25];
        static T dyKernel[25];

    public:
        static void Initialize(const T dx, const T dy) {
            mat_9.resize(3, 9);
            Eigen::MatrixXd temp;
            temp.resize(9, 3);
            int id = 0;
            for (int y = -1; y <= 1; ++y) {
                for (int x = -1; x <= 1; ++x) {
                    temp(id, 0) = (T)x * dx;
                    temp(id, 1) = (T)y * dy;
                    temp(id, 2) = 1.0;
                    ++id;
                }
            }

            // pinv
            Eigen::MatrixXd tmp = temp.transpose() * temp;
            mat_9 = tmp.inverse() * temp.transpose();
        }

        static void Initialize(const T dx, const T dy, const T dz) {
            mat_27.resize(4, 27);
            Eigen::MatrixXd temp;
            temp.resize(27, 4);
            int id = 0;
            for (int z = -1; z <= 1; ++z) {
                for (int y = -1; y <= 1; ++y) {
                    for (int x = -1; x <= 1; ++x) {
                        temp(id, 0) = (T)x * dx; 
                        temp(id, 1) = (T)y * dy;
                        temp(id, 2) = (T)z * dz;
                        temp(id, 3) = 1.0;
                        ++id;
                    }
                }
            }
            
            // pinv
            Eigen::MatrixXd tmp = temp.transpose() * temp;
            mat_27 = tmp.inverse() * temp.transpose();
        }

        static T LagrangianEulerian2D(ZD::CPoint<T, N> *params) {
            T fun[N][3];
            /* estimate function for each parameter */
            for (int i = 0; i < N; ++i) {
                Eigen::MatrixXd matB;
                matB.resize(9, 1);
                for (int k = 0; k < 9; ++k)
                    matB(k, 0) = params[k][i];
                Eigen::MatrixXd tmp = mat_9 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
            }

            Eigen::Matrix2d matJ;
            matJ.setZero();
            for (int i = 0; i < N; ++i) {
                matJ(0, 0) += fun[i][0] * fun[i][0];
                matJ(1, 0) += fun[i][1] * fun[i][0];
                matJ(0, 1) += fun[i][0] * fun[i][1];
                matJ(1, 1) += fun[i][1] * fun[i][1];
            }
            Eigen::EigenSolver<Eigen::Matrix2d> es(matJ);
            T tmp = std::max(es.eigenvalues()(0).real(), es.eigenvalues()(1).real());
            T le = std::sqrt(tmp);
            return le;
        }

        static T LagrangianEulerian2DSobel(ZD::CPoint<T, N> *params)
        {
            // using 5 x 5 sobel kernel to estimate the gradient
            T gradient[N][2];
            for (int i = 0; i < N; ++i) {
                gradient[i][0] = gradient[i][1] = 0.0;
                for (int k = 0; k < 25; ++k) {
                    gradient[i][0] += dxKernel[k] * params[k][i];        // dx
                    gradient[i][1] += dyKernel[k] * params[k][i];        // dy
                }
            }

            Eigen::Matrix2d matJ;
            matJ.setZero();
            for (int i = 0; i < N; ++i) {
                matJ(0, 0) += gradient[i][0] * gradient[i][0];
                matJ(1, 0) += gradient[i][1] * gradient[i][0];
                matJ(0, 1) += gradient[i][0] * gradient[i][1];
                matJ(1, 1) += gradient[i][1] * gradient[i][1];
            }
            Eigen::EigenSolver<Eigen::Matrix2d> es(matJ);
            T tmp = std::max(es.eigenvalues()(0).real(), es.eigenvalues()(1).real());
            T le = std::sqrt(tmp);
            return le;
        }

        static CPoint<T, 4> LagrangianEulerianST2D(CPoint<T, N> *params) {
            T fun[N][3];
            /* estimate function for each parameter */
            for (int i = 0; i < N; ++i) {
                Eigen::MatrixXd matB;
                matB.resize(9, 1);
                for (int k = 0; k < 9; ++k)
                    matB(k, 0) = params[k][i];
                Eigen::MatrixXd tmp = mat_9 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
            }

            CPoint<T, 4> st;
            for (int i = 0; i < N; ++i) {
                st[0] += fun[i][0] * fun[i][0];
                st[1] += fun[i][1] * fun[i][0];
                st[2] += fun[i][0] * fun[i][1];
                st[3] += fun[i][1] * fun[i][1];
            }
            return st;
        }

        static T LagrangianEulerian3D(CPoint<T, N> *params) {
            T fun[N][4];
            /* estimate function for each parameter */
            for (int i = 0; i < N; ++i) {
                Eigen::MatrixXd matB;
                matB.resize(27, 1);
                for (int k = 0; k < 27; ++k)
                    matB(k, 0) = params[k][i];
                Eigen::MatrixXd tmp = mat_27 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
                fun[i][3] = tmp(3, 0);
            }

            Eigen::Matrix3d matJ;    // not the jacobian matrix, it eqauls to J x J^{T}
            matJ.setZero();
            for (int i = 0; i < N; ++i) {
                matJ(0, 0) += fun[i][0] * fun[i][0];
                matJ(1, 0) += fun[i][1] * fun[i][0];
                matJ(2, 0) += fun[i][2] * fun[i][0];
                matJ(0, 1) += fun[i][0] * fun[i][1];
                matJ(1, 1) += fun[i][1] * fun[i][1];
                matJ(2, 1) += fun[i][2] * fun[i][1];
                matJ(0, 2) += fun[i][0] * fun[i][2];
                matJ(1, 2) += fun[i][1] * fun[i][2];
                matJ(2, 2) += fun[i][2] * fun[i][2];
            }
            Eigen::EigenSolver<Eigen::Matrix3d> es(matJ);
            T tmp = std::max(es.eigenvalues()(0).real(), std::max(es.eigenvalues()(1).real(), es.eigenvalues()(2).real()));
            T le = std::sqrt(tmp);
            return le;
        }

        static T LagrangianEulerian3DSobel(ZD::CPoint<T, N> *params)
        {
            T kernel_x[27] = {
                 -1,  0,  1,
                 -3,  0,  3, 
                 -1,  0,  1,
                 -3,  0,  3,
                 -6,  0,  6, 
                 -3,  0,  3,
                 -1,  0,  1,
                 -3,  0,  3,
                 -1,  0,  1
            };
            T kernel_y[27] = {
                -1, -3, -1,
                 0,  0,  0,
                 1,  3,  1,
                -3, -6, -3,
                 0,  0,  0,
                 3,  6,  3,
                -1, -3, -1,
                 0,  0,  0,
                 1,  3,  1
            };
            T kernel_z[27] = {
                 -1, -3, -1,
                 -3, -6, -3, 
                 -1, -3, -1,
                  0,  0,  0,
                  0,  0,  0, 
                  0,  0,  0,
                  1,  3,  1,
                  3,  6,  3,
                  1,  3,  1
            };

            // using 3 x 3 x 3 sobel kernel to estimate the gradient
            T gradient[N][3];
            for (int i = 0; i < N; ++i) {
                gradient[i][0] = gradient[i][1] = gradient[i][2] = 0.0;
                for (int k = 0; k < 27; ++k) {
                    gradient[i][0] += kernel_x[k] * params[k][i];        // dx
                    gradient[i][1] += kernel_y[k] * params[k][i];        // dy
                    gradient[i][2] += kernel_z[k] * params[k][i];        // dz
                }
            }

            Eigen::Matrix3d matJ;
            matJ.setZero();
            for (int i = 0; i < N; ++i) {
                matJ(0, 0) += gradient[i][0] * gradient[i][0];
                matJ(1, 0) += gradient[i][1] * gradient[i][0];
                matJ(2, 0) += gradient[i][2] * gradient[i][0];
                matJ(0, 1) += gradient[i][0] * gradient[i][1];
                matJ(1, 1) += gradient[i][1] * gradient[i][1];
                matJ(2, 1) += gradient[i][2] * gradient[i][1];
                matJ(0, 2) += gradient[i][0] * gradient[i][2];
                matJ(1, 2) += gradient[i][1] * gradient[i][2];
                matJ(2, 2) += gradient[i][2] * gradient[i][2];
            }
            Eigen::EigenSolver<Eigen::Matrix3d> es(matJ);
            T tmp = std::max(std::max(es.eigenvalues()(0).real(), es.eigenvalues()(1).real()), es.eigenvalues()(2).real());
            T le = std::sqrt(tmp);
            return le;
        }
    };

    template <typename T, unsigned int N> Eigen::MatrixXd CLETool<T, N>::mat_27;
    template <typename T, unsigned int N> Eigen::MatrixXd CLETool<T, N>::mat_9;
    template <typename T, unsigned int N> T CLETool<T, N>::dxKernel[25] = { 1.0,  2.0, 0.0, -2.0, -1.0,
                                                                            4.0,  8.0, 0.0, -8.0, -4.0,
                                                                            6.0, 12.0, 0.0, -12.0, -6.0,
                                                                            4.0,  8.0, 0.0, -8.0, -4.0,
                                                                            1.0,  2.0, 0.0, -2.0, -1.0 }; // 5 x 5 sobel kernel
    template <typename T, unsigned int N> T CLETool<T, N>::dyKernel[25] = { 1.0,  4.0,   6.0,  4.0,  1.0,
                                                                            2.0,  8.0,  12.0,  8.0,  2.0,
                                                                            0.0,  0.0,   0.0,  0.0,  0.0,
                                                                           -2.0, -8.0, -12.0, -8.0, -2.0,
                                                                           -1.0, -4.0,  -6.0, -4.0, -1.0 }; // 5 x 5 sobel kernel
}

#endif
