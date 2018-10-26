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
#ifndef _ZD_LIB_GRADIENT_TOOL_HPP_
#define _ZD_LIB_GRADIENT_TOOL_HPP_

#include "../Base/ZD_Point.hpp"
#include "../Base/ZD_Field.hpp"
#include "../define.hpp"

#include <Eigen/Eigenvalues> 

namespace ZD {
    template <typename T>
    class CGradientTool {
    private:
        // 2D
        //static void GX2D(const T *in, T *out, const int w, const int h)
        //{
        //    for (int y = 0; y < h; ++y) {
        //        // first
        //        out[y*w] = in[y*w + 1] - in[y*w];
        //        // last
        //        out[y*w + w - 1] = in[y*w + w - 1] - in[y*w + w - 2];

        //        // middle
        //        for (int x = 1; x < w - 1; ++x) {
        //            out[y*w + x] = (in[y*w + x + 1] - in[y*w + x - 1]) / 2.0;
        //        }
        //    }
        //}
        //static void GY2D(const T *in, T *out, const int w, const int h)
        //{
        //    for (int x = 0; x < w; ++x) {
        //        // first
        //        out[x] = in[w + x] - in[x];
        //        // last
        //        out[(h - 1)*w + x] = in[(h - 1)*w + x] - in[(h - 2)*w + x];

        //        // middle
        //        for (int y = 1; y < h - 1; ++y) {
        //            out[y*w + x] = (in[(y + 1)*w + x] - in[(y - 1)*w + x]) / 2.0;
        //        }
        //    }
        //}
        static void GX2D(const T *in, T *out, const int w, const int h)
        {
            memset(out, 0, sizeof(T)*w*h);
            for (int y = 1; y < h-1; ++y) {
                for (int x = 1; x < w - 1; ++x) {
                    out[y*w+x] = -3.0 * in[(y-1)*w+(x-1)] - 10.0 * in[y*w+(x-1)] - 3.0 * in[(y+1)*w+(x-1)] + 
                                  3.0 * in[(y-1)*w+(x+1)] + 10.0 * in[y*w+(x+1)] + 3.0 * in[(y+1)*w+(x+1)];
                }
            }
        }
        static void GY2D(const T *in, T *out, const int w, const int h)
        {
            memset(out, 0, sizeof(T)*w*h);
            for (int y = 1; y < h-1; ++y) {
                for (int x = 1; x < w - 1; ++x) {
                    out[y*w+x] = -3.0 * in[(y-1)*w+(x-1)] - 10.0 * in[(y-1)*w+x] - 3.0 * in[(y-1)*w+(x-1)] + 
                                  3.0 * in[(y+1)*w+(x+1)] + 10.0 * in[(y+1)*w+x] + 3.0 * in[(y+1)*w+(x+1)];
                }
            }
        }

        // 3D
        static void GX3D(const T *in, T *out, const int w, const int h, const int d)
        {
            //#pragma omp parallel for 
            for (int z = 0; z < d; ++z) {
                for (int y = 0; y < h; ++y) {
                    // first
                    out[(z*h + y)*w] = in[(z*h + y)*w + 1] - in[(z*h + y)*w];
                    // last
                    out[(z*h + y)*w + w - 1] = in[(z*h + y)*w + w - 1] - in[(z*h + y)*w + w - 2];

                    // middle
                    for (int x = 1; x < w - 1; ++x) {
                        out[(z*h + y)*w + x] = (in[(z*h + y)*w + x + 1] - in[(z*h + y)*w + x - 1]) / 2.0f;
                    }
                }
            }
        }
        static void GY3D(const T *in, T *out, const int w, const int h, const int d)
        {
            for (int z = 0; z < d; ++z) {
                for (int y = 0; y < h; ++y) {
                    for (int x = 0; x < w; ++x) {
                        if (y == 0) {
                            out[(z*h + y)*w + x] = in[(z*h + 1)*w + x] - in[(z*h + 0)*w + x];
                        }
                        else if (y == h - 1) {
                            out[(z*h + y)*w + x] = in[(z*h + h - 1)*w + x] - in[(z*h + h - 2)*w + x];
                        }
                        else {
                            out[(z*h + y)*w + x] = (in[(z*h + y + 1)*w + x] - in[(z*h + y - 1)*w + x]) / 2.0;
                        }
                    }
                }
            }
        }
        static void GZ3D(const T *in, T *out, const int w, const int h, const int d) 
        {
            for (int z = 0; z < d; ++z) {
                for (int y = 0; y < h; ++y) {
                    for (int x = 0; x < w; ++x) {
                        if (z == 0) {
                            out[(z*h + y)*w + x] = in[(1 * h + y)*w + x] - in[(0 * h + y)*w + x];
                        }
                        else if (z == d - 1) {
                            out[(z*h + y)*w + x] = in[((d - 1)*h + y)*w + x] - in[((d - 2)*h + y)*w + x];
                        }
                        else {
                            out[(z*h + y)*w + x] = (in[((z + 1)*h + y)*w + x] - in[((z - 1)*h + y)*w + x]) / 2.0;
                        }
                    }
                }
            }
        }
        static void EstimateFunction(const T *s, T *f, T *f2, const int w, const int h)
        {
            const int pts[9][2] = { {-1, -1}, {0, -1}, {1, -1}, 
                                  {-1,  0}, {0,  0}, {1,  0}, 
                                  {-1,  1}, {0,  1}, {1,  1} };
            {
                Eigen::Matrix2d matA;
                for (int i = 0; i < 9; ++i) {
                    matA(i, 0) = pts[i][0] * pts[i][0];    // x * x
                    matA(i, 1) = pts[i][0] * pts[i][1];    // x * y
                    matA(i, 2) = pts[i][1] * pts[i][1];    // y * y
                    matA(i, 3) = pts[i][0];                // x
                    matA(i, 4) = pts[i][1];              // y
                    matA(i, 5) = 1.0;                    // 1.0
                }
                // pinv
                Eigen::Matrix2d tmp = matA.transpose() * matA;
                Eigen::Matrix2d pinvA = tmp.inverse() * matA.transpose();

                for (int y = 1; y < h - 1; ++y) {
                    for (int x = 1; x < w - 1; ++x) {
                        Eigen::Matrix2d matB;
                        for (int i = 0; i < 9; ++i) {
                            matB(i, 0) = s[(y + pts[i][1])*w + (x + pts[i][0])];
                        }
                        Eigen::Matrix2d matX = pinvA * matB;

                        int index = y * w + x;
                        f[index * 6 + 0] = matX(0, 0);
                        f[index * 6 + 1] = matX(1, 0);
                        f[index * 6 + 2] = matX(2, 0);
                        f[index * 6 + 3] = matX(3, 0);
                        f[index * 6 + 4] = matX(4, 0);
                        f[index * 6 + 5] = matX(5, 0);
                    }
                }
            }

            {
                Eigen::Matrix2d matA;
                for (int i = 0; i < 9; ++i) {
                    matA(i, 0) = pts[i][0];                // x
                    matA(i, 1) = pts[i][1];              // y
                    matA(i, 2) = 1.0;                    // 1.0
                }
                // pinv
                Eigen::Matrix2d tmp = matA.transpose() * matA;
                Eigen::Matrix2d pinvA = tmp.inverse() * matA.transpose();

                for (int y = 1; y < h - 1; ++y) {
                    for (int x = 1; x < w - 1; ++x) {
                        Eigen::Matrix2d matB;
                        for (int i = 0; i < 9; ++i) {
                            matB(i, 0) = s[(y + pts[i][1])*w + (x + pts[i][0])];
                        }
                        Eigen::Matrix2d matX = pinvA * matB;

                        int index = y * w + x;
                        f2[index*3+0] = matX(0, 0);
                        f2[index*3+1] = matX(1, 0);
                        f2[index*3+2] = matX(2, 0);
                    }
                }
            }
            
        }
    public:
        // 2D
        static void ComputeGrad(const CField2<T, 1> *scalar, CField2<T, 2> *grad, const int w, const int h)
        {
            T *s = new T[w*h];
            T *gx = new T[w*h];
            T *gy = new T[w*h];

            memcpy(s, scalar->GetData(), sizeof(T)*w*h);

            GX2D(s, gx, w, h);
            GY2D(s, gy, w, h);

            SafeDeleteArray(s);

            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    grad->SetValue(x, y, CPoint<T, 2>(gx[y*w+x], gy[y*w+x]));
                }
            }

            SafeDeleteArray(gx);
            SafeDeleteArray(gy);
        }
        static void ComputeHess(const CField2<T, 1> *scalar, CField2<T, 4> *hess, const int w, const int h)
        {
            T *s = new T[w*h];
            T *gx = new T[w*h];
            T *gy = new T[w*h];

            memcpy(s, scalar->GetData(), sizeof(T)*w*h);

            GX2D(s, gx, w, h);
            GY2D(s, gy, w, h);

            SafeDeleteArray(s);

            T *hxx = new T[w*h];
            T *hxy = new T[w*h];
            T *hyx = new T[w*h];
            T *hyy = new T[w*h];

            GX2D(gx, hxx, w, h);
            GY2D(gx, hxy, w, h);
            GX2D(gy, hyx, w, h);
            GY2D(gy, hyy, w, h);

            SafeDeleteArray(gx);
            SafeDeleteArray(gy);

            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    hess->SetValue(x, y, CPoint<T, 4>(hxx[y*w+x], hxy[y*w+x], hyx[y*w+x], hyy[y*w+x]));
                }
            }

            SafeDeleteArray(hxx);
            SafeDeleteArray(hxy);
            SafeDeleteArray(hyx);
            SafeDeleteArray(hyy);
        }
        static void ComputeGradHess(const CField2<T, 1> *scalar, CField2<T, 2> *grad, CField2<T, 4> *hess, const int w, const int h)
        {
            T *s = new T[w*h];
            T *gx = new T[w*h];
            T *gy = new T[w*h];

            memcpy(s, scalar->GetData(), sizeof(T)*w*h);

            GX2D(s, gx, w, h);
            GY2D(s, gy, w, h);

            SafeDeleteArray(s);

            if (grad != nullptr) {
                for (int y = 0; y < h; ++y) {
                    for (int x = 0; x < w; ++x) {
                        grad->SetValue(x, y, CPoint<T, 2>(gx[y*w + x], gy[y*w + x]));
                    }
                }
            }

            T *hxx = new T[w*h];
            T *hxy = new T[w*h];
            T *hyx = new T[w*h];
            T *hyy = new T[w*h];

            GX2D(gx, hxx, w, h);
            GY2D(gx, hxy, w, h);
            GX2D(gy, hyx, w, h);
            GY2D(gy, hyy, w, h);

            SafeDeleteArray(gx);
            SafeDeleteArray(gy);

            if (hess != nullptr) {
                for (int y = 0; y < h; ++y) {
                    for (int x = 0; x < w; ++x) {
                        hess->SetValue(x, y, CPoint<T, 4>(hxx[y*w + x], hxy[y*w + x], hyx[y*w + x], hyy[y*w + x]));
                    }
                }
            }

            SafeDeleteArray(hxx);
            SafeDeleteArray(hxy);
            SafeDeleteArray(hyx);
            SafeDeleteArray(hyy);
        }
        //static void ComputeGradHess(const CField2<T, 1> *scalar, CField2<T, 2> *grad, CField2<T, 4> *hess, const int w, const int h)
        //{
        //    T *s = new T[w*h];
        //    T *function = new T[w*h*6];
        //    T *function2 = new T[w*h*3];

        //    memset(function, 0, sizeof(T)*w*h*6);
        //    memset(function2, 0, sizeof(T)*w*h*3);
        //    memcpy(s, scalar->GetData(), sizeof(T)*w*h);

        //    EstimateFunction(s, function, function2, w, h);

        //    SafeDeleteArray(s);

        //    if (grad != nullptr) {
        //        for (int y = 0; y < h; ++y) {
        //            for (int x = 0; x < w; ++x) {
        //                int index = y * w + x;
        //                grad->SetValue(x, y, CPoint<T, 2>(function2[index*3+0], function2[index*3+1]));
        //            }
        //        }
        //    }

        //    if (hess != nullptr) {
        //        for (int y = 0; y < h; ++y) {
        //            for (int x = 0; x < w; ++x) {
        //                int index = y * w + x;
        //                hess->SetValue(x, y, CPoint<T, 4>(2.0 * function[index*6+0], function[index*6+1], function[index*6+1], 2.0 * function[index*6+2]));
        //            }
        //        }
        //    }

        //    SafeDeleteArray(function);
        //    SafeDeleteArray(function2);
        //}
        // 3D
        static void ComputeGradHess(const CField3<T, 1> *scalar, CField3<T, 3> *grad, CField3<T, 9> *hess, const int w, const int h, const int d)
        {
            T *s = new T[w*h*d];
            T *gx = new T[w*h*d];
            T *gy = new T[w*h*d];
            T *gz = new T[w*h*d];

            memcpy(s, scalar->GetData(), sizeof(T)*w*h*d);

            GX3D(s, gx, w, h, d);
            GY3D(s, gy, w, h, d);
            GZ3D(s, gz, w, h, d);

            SafeDeleteArray(s);

            if (grad != nullptr) {
                for (int z = 0; z < d; ++z) {
                    for (int y = 0; y < h; ++y) {
                        for (int x = 0; x < w; ++x) {
                            grad->SetValue(x, y, z, CPoint<T, 3>(gx[(z*h+y)*w+x], gy[(z*h+y)*w+x], gz[(z*h+y)*w+x]));
                        }
                    }
                }
            }

            T *hxx = new T[w*h*d];
            T *hxy = new T[w*h*d];
            T *hxz = new T[w*h*d];
            T *hyy = new T[w*h*d];
            T *hyz = new T[w*h*d];
            T *hzz = new T[w*h*d];

            GX3D(gx, hxx, w, h, d);
            GY3D(gx, hxy, w, h, d);
            GZ3D(gx, hxz, w, h, d);
            GY3D(gy, hyy, w, h, d);
            GZ3D(gy, hyz, w, h, d);
            GZ3D(gz, hzz, w, h, d);

            SafeDeleteArray(gx);
            SafeDeleteArray(gy);
            SafeDeleteArray(gz);

            T H[9];
            if (hess != nullptr) {
                for (int z = 0; z < d; ++z) {
                    for (int y = 0; y < h; ++y) {
                        for (int x = 0; x < w; ++x) {
                            T H[9] = { 
                                hxx[(z*h+y)*w+x], hxy[(z*h+y)*w+x], hxz[(z*h+y)*w+x], 
                                hxy[(z*h+y)*w+x], hyy[(z*h+y)*w+x], hyz[(z*h+y)*w+x], 
                                hxz[(z*h+y)*w+x], hyz[(z*h+y)*w+x], hzz[(z*h+y)*w+x] };
                            hess->SetValue(x, y, z, CPoint<T, 9>(H));
                        }
                    }
                }
            }

            SafeDeleteArray(hxx);
            SafeDeleteArray(hxy);
            SafeDeleteArray(hxz);
            SafeDeleteArray(hyy);
            SafeDeleteArray(hyz);
            SafeDeleteArray(hzz);
        }
    };
}

#endif
