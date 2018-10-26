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
#ifndef _ZD_LIB_CS_TOOL_HPP_
#define _ZD_LIB_CS_TOOL_HPP_

#include "../Base/ZD_Point.hpp"

#include <Eigen/Eigenvalues> 

namespace ZD {
    //double cs_mat_26[4][26] = {
    //    { -0.0555555555555556,                    0,     0.0555555555555555, -0.0555555555555555,                   0,  0.0555555555555555, -0.0555555555555555,                      0,  0.0555555555555555, -0.0555555555555555,                     0,     0.0555555555555555, -0.0555555555555555, 0.0555555555555555, -0.0555555555555555,                    0, 0.0555555555555555, -0.0555555555555555,                      0,  0.0555555555555555, -0.0555555555555555,                  0,    0.0555555555555555, -0.0555555555555555,                  0, 0.0555555555555555 },
    //    { -0.0555555555555555, -0.0555555555555556, -0.0555555555555555,                      0,                   0,                   0,    0.0555555555555555,     0.0555555555555555,  0.0555555555555555, -0.0555555555555555, -0.0555555555555555, -0.0555555555555555,                   0,                  0,  0.0555555555555555, 0.0555555555555555, 0.0555555555555555, -0.0555555555555555, -0.0555555555555555, -0.0555555555555555,                    0,                  0,                   0,     0.0555555555555555, 0.0555555555555555, 0.0555555555555555 },
    //    { -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556, -0.0555555555555556,                   0,                   0,                   0,                   0,                  0,                   0,                  0,                  0,    0.0555555555555555,     0.0555555555555556,  0.0555555555555555,  0.0555555555555556, 0.0555555555555555,    0.0555555555555555,     0.0555555555555556, 0.0555555555555555, 0.0555555555555556 },
    //    {  0.0384615384615385,  0.0384615384615385,     0.0384615384615385,     0.0384615384615385,  0.0384615384615385,  0.0384615384615385,    0.0384615384615385,     0.0384615384615385,  0.0384615384615385,  0.0384615384615385,    0.0384615384615385,     0.0384615384615385,  0.0384615384615385, 0.0384615384615385,  0.0384615384615385, 0.0384615384615385, 0.0384615384615385,    0.0384615384615385,     0.0384615384615385,  0.0384615384615385,  0.0384615384615385, 0.0384615384615385,    0.0384615384615385,     0.0384615384615385, 0.0384615384615385, 0.0384615384615385 }
    //};

    template <typename T>
    class CCSTool {
    private:
        static Eigen::MatrixXd mat_ff_27;

    private:
        //static void LeastSquare(T *b, T *x)
        //{
        //    x[0] = x[1] = x[2] = x[3] = 0.0;
        //    for (int i = 0; i < 26; ++i) {
        //        x[0] += cs_mat_26[0][i] * b[i];
        //        x[1] += cs_mat_26[1][i] * b[i];
        //        x[2] += cs_mat_26[2][i] * b[i];
        //        x[3] += cs_mat_26[3][i] * b[i];
        //    }
        //}

        static T ComputeFSRValue(CPoint<T, 3> *pts, const T *dis)
        {
            Eigen::Matrix3d mat;
            for (unsigned int i = 0; i < 3; ++i) {
                mat(i, 0) = (pts[1][i] - pts[0][i]) / dis[0];
                mat(i, 1) = (pts[3][i] - pts[2][i]) / dis[1];
                mat(i, 2) = (pts[5][i] - pts[4][i]) / dis[2];
            }

            Eigen::Matrix3d matCG;
            matCG = mat.transpose() * mat;
            Eigen::EigenSolver<Eigen::Matrix3d> es(matCG);
            T value = std::max(es.eigenvalues()(0).real(), std::max(es.eigenvalues()(1).real(), es.eigenvalues()(2).real()));
            return value;
        }

        static T ComputeFSRST(CPoint<T, 3> *pts, const T *dis, CPoint<T, 7>& tensor)
        {
            Eigen::Matrix3d mat;
            for (unsigned int i = 0; i < 3; ++i) {
                mat(i, 0) = (pts[1][i] - pts[0][i]) / dis[0];
                mat(i, 1) = (pts[3][i] - pts[2][i]) / dis[1];
                mat(i, 2) = (pts[5][i] - pts[4][i]) / dis[2];
            }

            Eigen::Matrix3d matCG;
            matCG = mat.transpose() * mat;
            Eigen::EigenSolver<Eigen::Matrix3d> es(matCG);
            T value = std::max(es.eigenvalues()(0).real(), std::max(es.eigenvalues()(1).real(), es.eigenvalues()(2).real()));
            
            tensor[0] = 1.0;
            tensor[1] = matCG(0, 0);        // xx
            tensor[2] = matCG(0, 1);        // xy
            tensor[3] = matCG(0, 2);        // xz
            tensor[4] = matCG(1, 1);        // yy
            tensor[5] = matCG(1, 2);        // yz
            tensor[6] = matCG(2, 2);        // zz

            return value;
        }

    public:
        static void InitializeFiberFunction(const T dx, const T dy, const T dz)
        {
            mat_ff_27.resize(4, 27);
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
            mat_ff_27 = tmp.inverse() * temp.transpose();
        }

        //static T FiberFunction(CPoint<T, 9> **params)
        //{
        //    T fun[9][4];
        //    /* estimate function for each parameter */
        //    for (int i = 0; i < 9; ++i) {
        //        T b[26];
        //        for (int k = 0; k < 26; ++k) {
        //            b[k] = (*params[k])[i];
        //        }
        //        LeastSquare(b, fun[i]);
        //    }

        //    T dis = 0.0;
        //    for (int i = 0; i < 9; ++i) {
        //        dis += fun[i][0] * fun[i][0] + fun[i][1] * fun[i][1] + fun[i][2] * fun[i][2];
        //    }

        //    return dis;
        //}

        static T FiberFunction(CPoint<T, 3> **params)
        {
            T fun[3][4];
            Eigen::MatrixXd matB;
            matB.resize(27, 1);
            /* estimate function for each parameter */
            for (int i = 0; i < 3; ++i) {
                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = (*params[k])[i];
                }
                Eigen::MatrixXd tmp = mat_ff_27 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
                fun[i][3] = tmp(3, 0);
            }

            T dis = 0.0;
            for (int i = 0; i < 3; ++i) {
                dis += fun[i][0] * fun[i][0] + fun[i][1] * fun[i][1] + fun[i][2] * fun[i][2];
            }

            return dis;
        }

        static T FiberFunction(CPoint<T, 9> **params)
        {
            T fun[9][4];
            Eigen::MatrixXd matB;
            matB.resize(27, 1);
            /* estimate function for each parameter */
            for (int i = 0; i < 9; ++i) {
                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = (*params[k])[i];
                }
                Eigen::MatrixXd tmp = mat_ff_27 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
                fun[i][3] = tmp(3, 0);
            }

            T dis = 0.0;
            for (int i = 0; i < 9; ++i) {
                dis += fun[i][0] * fun[i][0] + fun[i][1] * fun[i][1] + fun[i][2] * fun[i][2];
            }

            return dis;
        }
        
        static T FSR(CPoint<T, 9> **params, const T *dis)
        {
            CPoint<T, 3> ref_dir = CPoint<T, 3>((*params[6])[0], (*params[6])[1], (*params[6])[2]);

            CPoint<T, 3> dir[6], f[6], b[6];
            for (int i = 0; i < 6; ++i) {
                dir[i] = CPoint<T, 3>((*params[i])[0], (*params[i])[1], (*params[i])[2]);
                if (InnerProduct(ref_dir, dir[i]) > 0.0) {
                    f[i] = CPoint<T, 3>((*params[i])[3], (*params[i])[4], (*params[i])[5]);
                    b[i] = CPoint<T, 3>((*params[i])[6], (*params[i])[7], (*params[i])[8]);
                }
                else {
                    b[i] = CPoint<T, 3>((*params[i])[3], (*params[i])[4], (*params[i])[5]);
                    f[i] = CPoint<T, 3>((*params[i])[6], (*params[i])[7], (*params[i])[8]);
                }
            }

            T fv = ComputeFSRValue(f, dis);
            T bv = ComputeFSRValue(b, dis);
            return fv > bv ? fv : bv;
        }

        //static CPoint<float, 7> FiberFunctionST(CPoint<T, 9> **params)
        //{
        //    T fun[9][4];
        //    /* estimate function for each parameter */
        //    for (int i = 0; i < 9; ++i) {
        //        T b[26];
        //        for (int k = 0; k < 26; ++k) {
        //            b[k] = (*params[k])[i];
        //        }
        //        LeastSquare(b, fun[i]);
        //    }

        //    CPoint<float, 7> tensor;
        //    tensor[0] = 1.0;
        //    for (int i = 0; i < 9; ++i) {
        //        tensor[1] += fun[i][0] * fun[i][0];        // xx
        //        tensor[2] += fun[i][0] * fun[i][1];        // xy
        //        tensor[3] += fun[i][0] * fun[i][2];        // xz
        //        tensor[4] += fun[i][1] * fun[i][1];        // yy
        //        tensor[5] += fun[i][1] * fun[i][2];        // yz
        //        tensor[6] += fun[i][2] * fun[i][2];        // zz
        //    }
        //    return tensor;
        //}

        static CPoint<T, 7> FiberFunctionST(CPoint<T, 3> **params)
        {
            T fun[3][4];
            Eigen::MatrixXd matB;
            matB.resize(27, 1);
            /* estimate function for each parameter */
            for (int i = 0; i < 3; ++i) {
                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = (*params[k])[i];
                }
                Eigen::MatrixXd tmp = mat_ff_27 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
                fun[i][3] = tmp(3, 0);
            }

            CPoint<T, 7> tensor;
            tensor[0] = 1.0;
            for (int i = 0; i < 3; ++i) {
                tensor[1] += fun[i][0] * fun[i][0];        // xx
                tensor[2] += fun[i][0] * fun[i][1];        // xy
                tensor[3] += fun[i][0] * fun[i][2];        // xz
                tensor[4] += fun[i][1] * fun[i][1];        // yy
                tensor[5] += fun[i][1] * fun[i][2];        // yz
                tensor[6] += fun[i][2] * fun[i][2];        // zz
            }
            return tensor;
        }

        static CPoint<T, 7> FiberFunctionST(CPoint<T, 9> **params)
        {
            T fun[9][4];
            Eigen::MatrixXd matB;
            matB.resize(27, 1);
            /* estimate function for each parameter */
            for (int i = 0; i < 9; ++i) {
                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = (*params[k])[i];
                }
                Eigen::MatrixXd tmp = mat_ff_27 * matB;
                fun[i][0] = tmp(0, 0);
                fun[i][1] = tmp(1, 0);
                fun[i][2] = tmp(2, 0);
                fun[i][3] = tmp(3, 0);
            }

            CPoint<T, 7> tensor;
            tensor[0] = 1.0;
            for (int i = 0; i < 9; ++i) {
                tensor[1] += fun[i][0] * fun[i][0];        // xx
                tensor[2] += fun[i][0] * fun[i][1];        // xy
                tensor[3] += fun[i][0] * fun[i][2];        // xz
                tensor[4] += fun[i][1] * fun[i][1];        // yy
                tensor[5] += fun[i][1] * fun[i][2];        // yz
                tensor[6] += fun[i][2] * fun[i][2];        // zz
            }
            return tensor;
        }

        static CPoint<T, 7> ResampleST(CPoint<T, 3> **params, const int num)
        {
            T *fun = new T[3*num*4];
            //T fun[9][4];
            Eigen::MatrixXd matB;
            matB.resize(27, 1);
            /* estimate function for each parameter */
            for (int i = 0; i < num; ++i) {
                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = params[k][i][0];
                }
                Eigen::MatrixXd tmp = mat_ff_27 * matB;
                fun[i*12+0] = tmp(0, 0);
                fun[i*12+1] = tmp(1, 0);
                fun[i*12+2] = tmp(2, 0);
                fun[i*12+3] = tmp(3, 0);

                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = params[k][i][1];
                }
                tmp = mat_ff_27 * matB;
                fun[i*12+4] = tmp(0, 0);
                fun[i*12+5] = tmp(1, 0);
                fun[i*12+6] = tmp(2, 0);
                fun[i*12+7] = tmp(3, 0);

                for (int k = 0; k < 27; ++k) {
                    matB(k, 0) = params[k][i][2];
                }
                tmp = mat_ff_27 * matB;
                fun[i*12+8] = tmp(0, 0);
                fun[i*12+9] = tmp(1, 0);
                fun[i*12+10] = tmp(2, 0);
                fun[i*12+11] = tmp(3, 0);
            }

            CPoint<T, 7> tensor;
            tensor[0] = 1.0;
            for (int i = 0; i < num*3; ++i) {
                tensor[1] += fun[i*4+0] * fun[i*4+0];        // xx
                tensor[2] += fun[i*4+0] * fun[i*4+1];        // xy
                tensor[3] += fun[i*4+0] * fun[i*4+2];        // xz
                tensor[4] += fun[i*4+1] * fun[i*4+1];        // yy
                tensor[5] += fun[i*4+1] * fun[i*4+2];        // yz
                tensor[6] += fun[i*4+2] * fun[i*4+2];        // zz
            }
            
            SafeDeleteArray(fun);
            
            return tensor;
        }

        static CPoint<T, 7> FSRST(CPoint<T, 9> **params, const T *dis)
        {
            CPoint<T, 3> ref_dir = CPoint<T, 3>((*params[6])[0], (*params[6])[1], (*params[6])[2]);

            CPoint<T, 3> dir[6], f[6], b[6];
            for (int i = 0; i < 6; ++i) {
                dir[i] = CPoint<T, 3>((*params[i])[0], (*params[i])[1], (*params[i])[2]);
                if (InnerProduct(ref_dir, dir[i]) > 0.0) {
                    f[i] = CPoint<T, 3>((*params[i])[3], (*params[i])[4], (*params[i])[5]);
                    b[i] = CPoint<T, 3>((*params[i])[6], (*params[i])[7], (*params[i])[8]);
                }
                else {
                    b[i] = CPoint<T, 3>((*params[i])[3], (*params[i])[4], (*params[i])[5]);
                    f[i] = CPoint<T, 3>((*params[i])[6], (*params[i])[7], (*params[i])[8]);
                }
            }

            T fsr[2];
            CPoint<T, 7> tensor[2];
            fsr[0] = ComputeFSRST(f, dis, tensor[0]);
            fsr[1] = ComputeFSRST(b, dis, tensor[1]);
            return fsr[0] > fsr[1] ? tensor[0] : tensor[1];
        }
    };

    template <typename T> Eigen::MatrixXd CCSTool<T>::mat_ff_27;
}

#endif
