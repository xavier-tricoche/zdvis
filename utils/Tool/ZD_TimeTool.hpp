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
#ifndef _ZD_LIB_TIME_TOOL_HPP_
#define _ZD_LIB_TIME_TOOL_HPP_

#include <chrono>

namespace ZD {
    template <typename TYPE>
    class CTimeTool {
    private:
        std::chrono::system_clock::time_point m_start;
        std::chrono::system_clock::time_point m_end;
        std::chrono::duration<TYPE>           m_diff;

    public:
        CTimeTool() {

        }
        ~CTimeTool() {

        }

    public:
        void StartTimer() {
            m_start = std::chrono::system_clock::now();
        }
        void StopTimer() {
            m_end = std::chrono::system_clock::now();
            m_diff = m_end - m_start;
        }
        TYPE Duration() {
            return m_diff.count();
        }
    };
}

#endif
