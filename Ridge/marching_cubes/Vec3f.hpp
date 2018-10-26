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
#ifndef __VEC3F_H
#define __VEC3F_H

typedef unsigned int uint;

struct Vec3f;

#ifdef __linux
#define __forceinline
#else
#define __forceinline 
#endif

struct Vec3f
{
    //
    // Initalization
    //
    __forceinline Vec3f();
    __forceinline Vec3f(const Vec3f &V);
    __forceinline Vec3f(float _x, float _y, float _z);

    //
    // Assignment
    //
    __forceinline Vec3f& operator = (const Vec3f &V);

    //
    // Overloaded operators
    //
    __forceinline Vec3f& operator *= (float Right);
    __forceinline Vec3f& operator *= (uint Right);
    __forceinline Vec3f& operator /= (float Right);
    __forceinline Vec3f& operator /= (uint Right);
    __forceinline Vec3f& operator += (const Vec3f &Right);
    __forceinline Vec3f& operator -= (const Vec3f &Right);

    //
    // Normalization
    //
    __forceinline void SetLength(float NewLength);

    //
    // Accessors
    //
    __forceinline float Length() const;
    __forceinline float LengthSq() const;

    __forceinline bool Valid() const;

#ifdef USE_D3D
    __forceinline operator D3DXVECTOR3() const;
#endif

    __forceinline float& operator[](uint Index)
    {
        return ((float *)this)[Index];
    }
    __forceinline float operator[](uint Index) const
    {
        return ((float *)this)[Index];
    }

    //
    // Local data
    //
    float x, y, z;

    //
    // Constants
    //
    static const Vec3f Origin;
    static const Vec3f eX;
    static const Vec3f eY;
    static const Vec3f eZ;
};


__forceinline Vec3f operator * (const Vec3f &Left, float Right);
__forceinline Vec3f operator * (float Left, const Vec3f &Right);
__forceinline Vec3f operator / (const Vec3f &Left, float Right);
__forceinline Vec3f operator + (const Vec3f &Left, const Vec3f &Right);
__forceinline Vec3f operator - (const Vec3f &Left, const Vec3f &Right);

//#include "Vec3f.inl"


#endif
