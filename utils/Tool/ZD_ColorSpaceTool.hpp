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
#ifndef _ZD_LIB_COLOR_SPACE_TOOL_HPP_
#define _ZD_LIB_COLOR_SPACE_TOOL_HPP_

#include "../Base/ZD_Point.hpp"

namespace ZD {
    template <typename T>
    class CColorSpaceTool {
    public:
        static void RGB2HSV(CPoint<unsigned char, 3>& rgb, CPoint<T, 3>& hsv)
        {

        }

        static void RGB2HSV(CPoint<T, 3>& rgb, CPoint<T, 3>& hsv)
        {

        }

        static void HSV2RGB(CPoint<T, 3>& hsv, CPoint<unsigned char, 3>& rgb)
        {
            T      hh, p, q, t, ff;
            long   i;

            if (hsv[1] <= 0.0) {
                rgb[0] = (unsigned char)(hsv[2] * 255.0);
                rgb[1] = (unsigned char)(hsv[2] * 255.0);
                rgb[2] = (unsigned char)(hsv[2] * 255.0);
                return;
            }
            hh = hsv[0];
            if (hh >= 360.0) hh = 0.0;
            hh /= 60.0;
            i = (long)hh;
            ff = hh - i;

            p = hsv[2] * (1.0 - hsv[1]);
            q = hsv[2] * (1.0 - (hsv[1] * ff));
            t = hsv[2] * (1.0 - (hsv[1] * (1.0 - ff)));

            switch (i) {
            case 0:
                rgb[0] = (unsigned char)(hsv[2] * 255.0);
                rgb[1] = (unsigned char)(t      * 255.0);
                rgb[2] = (unsigned char)(p      * 255.0);
                break;
            case 1:
                rgb[0] = (unsigned char)(q        * 255.0);
                rgb[1] = (unsigned char)(hsv[2] * 255.0);
                rgb[2] = (unsigned char)(p        * 255.0);
                break;
            case 2:
                rgb[0] = (unsigned char)(p        * 255.0);
                rgb[1] = (unsigned char)(hsv[2] * 255.0);
                rgb[2] = (unsigned char)(t        * 255.0);
                break;

            case 3:
                rgb[0] = (unsigned char)(p        * 255.0);
                rgb[1] = (unsigned char)(q        * 255.0);
                rgb[2] = (unsigned char)(hsv[2] * 255.0);
                break;
            case 4:
                rgb[0] = (unsigned char)(t        * 255.0);
                rgb[1] = (unsigned char)(p        * 255.0);
                rgb[2] = (unsigned char)(hsv[2] * 255.0);
                break;
            case 5:
            default:
                rgb[0] = (unsigned char)(hsv[2] * 255.0);
                rgb[1] = (unsigned char)(p        * 255.0);
                rgb[2] = (unsigned char)(q        * 255.0);
                break;
            }
        }
    };
}

#endif
