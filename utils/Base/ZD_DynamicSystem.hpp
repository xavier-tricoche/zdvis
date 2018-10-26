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
#ifndef _ZD_LIB_DYNAMIC_SYSTEM_HPP_
#define _ZD_LIB_DYNAMIC_SYSTEM_HPP_

#include "ZD_Point.hpp"

namespace ZD {
    template <typename T, unsigned int N>
    class CDynamicSystem {
    public:
        virtual inline const CPoint<T, N> Velocity(const CPoint<T, N>& p, 
            const T& t) const = 0;
    };
};

#endif
