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
#ifndef _ZD_LIB_RANDOM_TOOL_HPP_
#define _ZD_LIB_RANDOM_TOOL_HPP_

#include <random>

namespace ZD {
    template <typename T>
    class CRandomTool {
    private:
        std::random_device m_rd;
        std::mt19937 m_gen;
        std::uniform_real_distribution<T> m_dis;

    public:
        CRandomTool() : m_gen(m_rd()), m_dis(0, 1) {

        }

    public:
        T GetRandomNumer() {
            return m_dis(m_gen);
        }
    };
}

#endif
