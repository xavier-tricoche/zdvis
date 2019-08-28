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
#ifndef _ZD_LIB_STRING_TOOL_HPP_
#define _ZD_LIB_STRING_TOOL_HPP_

#include <string>
#include <cstring>

#include "../define.hpp"

namespace ZD {
    class CStringTool {
    public:
        /* to upper */
        static std::string ToUpper(const char *str) {
            size_t len = strlen(str);
            char *tmp = new char[len + 1];
            for (int i = 0; i < len; ++i) {
                if (str[i] >= 'a' && str[i] <= 'z') {
                    tmp[i] = str[i] + 'A' - 'a';
                }
                else {
                    tmp[i] = str[i];
                }
            }
            tmp[len] = '\0';

            std::string res = tmp;
            SafeDeleteArray(tmp);

            return res;
        }

        /* to lower */
        static std::string ToLower(const char *str) {
            size_t len = strlen(str);
            char *tmp = new char[len + 1];
            for (int i = 0; i < len; ++i) {
                if (str[i] >= 'A' && str[i] <= 'Z') {
                    tmp[i] = str[i] + 'a' - 'A';
                }
                else {
                    tmp[i] = str[i];
                }
            }
            tmp[len] = '\0';

            std::string res = tmp;
            SafeDeleteArray(tmp);

            return res;
        }
    };
}

#endif
