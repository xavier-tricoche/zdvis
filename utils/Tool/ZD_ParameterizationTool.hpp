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
#ifndef _ZD_LIB_PARAMETERIZATION_TOOL_HPP_
#define _ZD_LIB_PARAMETERIZATION_TOOL_HPP_

#include "../Base/ZD_Point.hpp"

#include <assert.h>

namespace ZD {
    template <typename T, unsigned int K>
    class CParameterizationTool {
    public:
        static CPoint<T, K> ParameterizeSingleScalar(T *values, const int n)
        {
            //assert(n > 1 && K >= 1);
            if (n == 1) {
                CPoint<T, K> p;
                p[0] = values[0];

                for (int k = 1; k < K; ++k) {
                    p[k] = 0.0;
                }
                return p;
            }
            else {
                CPoint<T, K> p;

                for (int i = 0; i < n; ++i)
                    p[0] += values[i];
                p[0] /= (T)n;

                for (int k = 1; k < K; ++k) {
                    for (int i = 0; i < n; ++i)
                        p[k] += std::pow((values[i] - p[0]), (k + 1));
                    p[k] /= (T)(n - 1);
                    T sign = p[k] > 0.0 ? 1.0 : -1.0;
                    p[k] = std::pow((double)std::abs(p[k]), 1.0 / (T)(k + 1)) * sign;
                }

                return p;
            }
        }

        static CPoint<T, K> ParameterizeSingleScalar(const std::vector<T>& values)
        {
            //assert(values.size() > 1 && K >= 1);
            if (values.size() == 1) {
                CPoint<T, K> p;
                p[0] = values[0];

                for (int k = 1; k < K; ++k) {
                    p[k] = 0.0;
                }
                return p;
            }
            else {
                CPoint<T, K> p;

                for (int i = 0; i < values.size(); ++i)
                    p[0] += values[i];
                p[0] /= T(values.size());

                for (int k = 1; k < K; ++k) {
                    for (int i = 0; i < values.size(); ++i)
                        p[k] += std::pow((values[i] - p[0]), (k + 1));
                    p[k] /= T(values.size() - 1);
                    T sign = p[k] > 0.0 ? 1.0 : -1.0;
                    p[k] = std::pow(std::abs(p[k]), 1.0 / (T)(k + 1)) * sign;
                }

                return p;
            }
        }
    };
}

#endif
