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
#ifndef _ZD_LIB_EIGEN_TOOL_HPP_
#define _ZD_LIB_EIGEN_TOOL_HPP_

#include "../Base/ZD_Point.hpp"

#include <Eigen/Eigenvalues> 

namespace ZD {
    template <typename T>
    class CEigenTool {
    public:
        // 2D
        /* eigen-vectors and eigen-values of 2 by 2 matrix, eval1 > eval2 */
        static void Eigen2D(const CPoint<T, 4>& mat, CPoint<T, 1> *eval, CPoint<T, 2> *evec) 
        {
            Eigen::Matrix2d tmp;
            tmp(0, 0) = mat[0];
            tmp(0, 1) = mat[1];
            tmp(1, 0) = mat[2];
            tmp(1, 1) = mat[3];

            Eigen::EigenSolver<Eigen::Matrix2d> es(tmp);
            int minID, maxID;
            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
                maxID = 0;
                minID = 1;
            }
            else {
                maxID = 1;
                minID = 0;
            }

            if (eval) {
                eval[0] = CPoint<T, 1>(es.eigenvalues()(maxID).real());
                eval[1] = CPoint<T, 1>(es.eigenvalues()(minID).real());
            }
            if (evec) {
                evec[0] = CPoint<T, 2>(es.eigenvectors()(0, maxID).real(), es.eigenvectors()(1, maxID).real());
                evec[1] = CPoint<T, 2>(es.eigenvectors()(0, minID).real(), es.eigenvectors()(1, minID).real());
            }
        }
        static void Eigen2D(const Eigen::Matrix2d& mat, CPoint<T, 1> *eval, CPoint<T, 2> *evec)
        {
            Eigen::EigenSolver<Eigen::Matrix2d> es(mat);
            int minID, maxID;
            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
                maxID = 0;
                minID = 1;
            }
            else {
                maxID = 1;
                minID = 0;
            }

            if (eval) {
                eval[0] = CPoint<T, 1>(es.eigenvalues()(maxID).real());
                eval[1] = CPoint<T, 1>(es.eigenvalues()(minID).real());
            }
            if (evec) {
                evec[0] = CPoint<T, 2>(es.eigenvectors()(0, maxID).real(), es.eigenvectors()(1, maxID).real());
                evec[1] = CPoint<T, 2>(es.eigenvectors()(0, minID).real(), es.eigenvectors()(1, minID).real());
            }
        }
        //static void Eigen2D(const Eigen::Matrix2f& mat, CPoint<T, 1> *eval, CPoint<T, 2> *evec)
        //{
        //    Eigen::EigenSolver<Eigen::Matrix2d> es(mat);
        //    int minID, maxID;
        //    if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
        //        maxID = 0;
        //        minID = 1;
        //    }
        //    else {
        //        maxID = 0;
        //        minID = 1;
        //    }

        //    if (eval) {
        //        eval[0] = CPoint<T, 1>(es.eigenvalues()(maxID).real());
        //        eval[1] = CPoint<T, 1>(es.eigenvalues()(minID).real());
        //    }
        //    if (evec) {
        //        evec[0] = CPoint<T, 2>(es.eigenvectors()(0, maxID).real(), es.eigenvectors()(1, maxID).real());
        //        evec[1] = CPoint<T, 2>(es.eigenvectors()(0, minID).real(), es.eigenvectors()(1, minID).real());
        //    }
        //}
        static T MaxEval2D(const CPoint<T, 4> *mat)
        {
            Eigen::Matrix2d tmp;
            tmp(0, 0) = mat[0];
            tmp(0, 1) = mat[1];
            tmp(1, 0) = mat[2];
            tmp(1, 1) = mat[3];

            Eigen::EigenSolver<Eigen::Matrix2d> es(tmp);

            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real())
                return es.eigenvalues()(0).real();
            else
                return es.eigenvalues()(1).real();
        }
        static T MaxEval2D(const Eigen::Matrix2d& mat)
        {
            Eigen::EigenSolver<Eigen::Matrix2d> es(mat);

            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real())
                return es.eigenvalues()(0).real();
            else
                return es.eigenvalues()(1).real();
        }

        // 3D
        /* eigen-vectors and eigen-values of 3 by 3 matrix, eval1 > eval2 > eval3 */
        static void Eigen3D(const CPoint<T, 9>& mat, CPoint<T, 1> *eval, CPoint<T, 3> *evec)
        {
            Eigen::Matrix3d tmp;
            tmp(0, 0) = mat[0];
            tmp(0, 1) = mat[1];
            tmp(0, 2) = mat[2];
            tmp(1, 0) = mat[3];
            tmp(1, 1) = mat[4];
            tmp(1, 2) = mat[5];
            tmp(2, 0) = mat[6];
            tmp(2, 1) = mat[7];
            tmp(2, 2) = mat[8];

            Eigen::EigenSolver<Eigen::Matrix3d> es(tmp);
            int minID, maxID, midID;
            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
                if (es.eigenvalues()(1).real() > es.eigenvalues()(2).real()) {
                    maxID = 0;
                    midID = 1;
                    minID = 2;
                } else {
                    if (es.eigenvalues()(0).real() > es.eigenvalues()(2).real()) {
                        maxID = 0;
                        midID = 2;
                        minID = 1;
                    } else {
                        maxID = 2;
                        midID = 0;
                        minID = 1;
                    }
                }
            }
            else {
                if (es.eigenvalues()(2).real() > es.eigenvalues()(1).real()) {
                    maxID = 2;
                    midID = 1;
                    minID = 0;
                }
                else {
                    if (es.eigenvalues()(2).real() > es.eigenvalues()(0).real()) {
                        maxID = 1;
                        midID = 2;
                        minID = 0;
                    } else {
                        maxID = 1;
                        midID = 0;
                        minID = 2;
                    }
                }
                
            }

            if (eval) {
                eval[0] = CPoint<T, 1>(es.eigenvalues()(maxID).real());
                eval[1] = CPoint<T, 1>(es.eigenvalues()(midID).real());
                eval[2] = CPoint<T, 1>(es.eigenvalues()(minID).real());
            }
            if (evec) {
                evec[0] = CPoint<T, 3>(es.eigenvectors()(0, maxID).real(), es.eigenvectors()(1, maxID).real(), es.eigenvectors()(2, maxID).real());
                evec[1] = CPoint<T, 3>(es.eigenvectors()(0, midID).real(), es.eigenvectors()(1, midID).real(), es.eigenvectors()(2, midID).real());
                evec[2] = CPoint<T, 3>(es.eigenvectors()(0, minID).real(), es.eigenvectors()(1, minID).real(), es.eigenvectors()(2, minID).real());
            }
        }

        static T MaxEval3D(const CPoint<T, 9>& mat)
        {
            Eigen::Matrix3d tmp;
            tmp(0, 0) = mat[0];
            tmp(0, 1) = mat[1];
            tmp(0, 2) = mat[2];
            tmp(1, 0) = mat[3];
            tmp(1, 1) = mat[4];
            tmp(1, 2) = mat[5];
            tmp(2, 0) = mat[6];
            tmp(2, 1) = mat[7];
            tmp(2, 2) = mat[8];

            Eigen::EigenSolver<Eigen::Matrix3d> es(tmp);
            int minID, maxID, midID;
            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
                if (es.eigenvalues()(1).real() > es.eigenvalues()(2).real()) {
                    maxID = 0;
                    midID = 1;
                    minID = 2;
                }
                else {
                    if (es.eigenvalues()(0).real() > es.eigenvalues()(2).real()) {
                        maxID = 0;
                        midID = 2;
                        minID = 1;
                    }
                    else {
                        maxID = 2;
                        midID = 0;
                        minID = 1;
                    }
                }
            }
            else {
                if (es.eigenvalues()(2).real() > es.eigenvalues()(1).real()) {
                    maxID = 2;
                    midID = 1;
                    minID = 0;
                }
                else {
                    if (es.eigenvalues()(2).real() > es.eigenvalues()(0).real()) {
                        maxID = 1;
                        midID = 2;
                        minID = 0;
                    }
                    else {
                        maxID = 1;
                        midID = 0;
                        minID = 2;
                    }
                }
            }

            return es.eigenvalues()(maxID).real();
        }

        static T MaxEval3D(const Eigen::Matrix3d& mat) {
            Eigen::EigenSolver<Eigen::Matrix3d> es(mat);
            int minID, maxID, midID;
            if (es.eigenvalues()(0).real() > es.eigenvalues()(1).real()) {
                if (es.eigenvalues()(1).real() > es.eigenvalues()(2).real()) {
                    maxID = 0;
                    midID = 1;
                    minID = 2;
                }
                else {
                    if (es.eigenvalues()(0).real() > es.eigenvalues()(2).real()) {
                        maxID = 0;
                        midID = 2;
                        minID = 1;
                    }
                    else {
                        maxID = 2;
                        midID = 0;
                        minID = 1;
                    }
                }
            }
            else {
                if (es.eigenvalues()(2).real() > es.eigenvalues()(1).real()) {
                    maxID = 2;
                    midID = 1;
                    minID = 0;
                }
                else {
                    if (es.eigenvalues()(2).real() > es.eigenvalues()(0).real()) {
                        maxID = 1;
                        midID = 2;
                        minID = 0;
                    }
                    else {
                        maxID = 1;
                        midID = 0;
                        minID = 2;
                    }
                }
            }

            return es.eigenvalues()(maxID).real();
        }

    };
}

#endif

//namespace ZD {
//    template <typename T>
//    class CEigenTool {
//    public:
//        static CEigenTool* GetInstance();
//
//        /* eigen-vectors and eigen-values of 2 by 2 matrix, eval1 > eval2 */
//        void Eigen2D(const CPoint<T, 4>& mat, CPoint<T, 1> *eval, CPoint<T, 2> *evec);
//        /* eigen-vectors and eigen-values of 3 by 3 matrix, eval1 > eval2 > eval3 */
//        void Eigen3D(const CPoint<T, 9>& mat, CPoint<T, 1> *eval, CPoint<T, 3> *evec);
//
//    private:
//        CEigenTool();
//        CEigenTool(CEigenTool const&);
//        CEigenTool& operator=(CEigenTool const&);
//
//    private:
//        static CEigenTool* m_pInstance;
//    };
//}
