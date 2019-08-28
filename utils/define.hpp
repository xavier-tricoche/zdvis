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
#ifndef _ZD_LIB_DEFINE_HPP_
#define _ZD_LIB_DEFINE_HPP_


#include <stdlib.h>
#include <vtkConfigure.h>


#ifndef ZD_CONSTANTS
#define ZD_CONSTANTS
const double ZD_PI = 3.1415926535897932384626433832795;
const double ZD_E = 2.718281828459045235;
const double ZD_EPSILON = 0.0000001;
const double ZD_EARTH_R = 6371.0;
#endif


#ifndef ZD_PATHNAME_LENGTH
#define ZD_PATHNAME_LENGTH 256
#endif

#if VTK_MAJOR_VERSION > 5
#define VTK_SET_INPUT(a, b) a->SetInputData(b)
#else
#define VTK_SET_INPUT(a, b) a->SetInput(b)
#endif

#ifndef ZD_MIN
#define ZD_MIN(a, b) (a) < (b) ? (a) : (b)
#endif

#ifndef ZD_MAX
#define ZD_MAX(a, b) (a) > (b) ? (a) : (b)
#endif


template< class T > void SafeDelete( T*& pVal )
{
    if (pVal != nullptr) {
        delete pVal;
        pVal = nullptr;
    }
}

template< class T > void SafeDeleteArray( T*& pVal )
{
    if (pVal != nullptr) {
        delete[] pVal;
        pVal = nullptr;
    }
}

template< class T > inline int ZD_Sign(T& a)
{
    if (a > 0.0)
        return 1;
    else if (a < 0.0)
        return -1;
    else
        return 0;
}

#endif
