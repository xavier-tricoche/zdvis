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
#ifndef _ZD_LIB_COLOR_MAP_HPP_
#define _ZD_LIB_COLOR_MAP_HPP_

#include "utils/Base/ZD_Point.hpp"
#include "utils/Tool/ZD_StringTool.hpp"
#include "utils/Tool/ZD_ColorSpaceTool.hpp"

namespace ZD {
    enum ZD_COLOR_MAP_TYPE {
        ZD_COLOR_MAP_NONE,
        ZD_COLOR_MAP_BLACK,
        ZD_COLOR_MAP_BLEU_YELLOW,
        ZD_COLOR_MAP_RAINBOW,
        ZD_COLOR_MAP_WHITE_RED,
        ZD_COLOR_MAP_CUBE_HELIX
    };

    template <typename T>
    class CColorMap {
    protected:
        T m_minValue, m_maxValue;

    public:
        CColorMap(const T min_value, const T max_value);
        ~CColorMap();

        virtual CPoint<T, 3> GetColor(const T value) = 0;

    public:
        static ZD_COLOR_MAP_TYPE ColorMapByName(const char *name) {
            std::string tmp = CStringTool::ToLower(name);

            if (tmp == "blue_yellow" || tmp == "blueyellow") {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_BLEU_YELLOW;
            } else if (tmp == "rainbow") {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_RAINBOW;
            } else if (tmp == "white_red" || tmp == "whitered") {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_WHITE_RED;
            } else if (tmp == "cube_helix" || tmp == "cubehelix") {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_CUBE_HELIX;
            } else if (tmp == "black") {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_BLACK;
            } 
            else {
                return ZD_COLOR_MAP_TYPE::ZD_COLOR_MAP_NONE;
            }
        }
    };
    
    template <typename T>
    class CColorMapBlueYellow : public CColorMap<T> {
    public:
        CColorMapBlueYellow(const T min_value, const T max_value);
        ~CColorMapBlueYellow();

        virtual CPoint<T, 3> GetColor(const T value);
    };
    
    template <typename T>
    class CColorMapRainbow : public CColorMap<T> {
    public:
        CColorMapRainbow(const T min_value, const T max_value);
        ~CColorMapRainbow();

        virtual CPoint<T, 3> GetColor(const T value);
    };
    
    template <typename T>
    class CColorMapCubeHelix : public CColorMap<T> {
    private:
        CPoint<T, 3> m_colorTable[8];
    public:
        CColorMapCubeHelix(const T min_value, const T max_value);
        ~CColorMapCubeHelix();

        virtual CPoint<T, 3> GetColor(const T value);
    };
} // namespace ZD

template<typename T>
ZD::CColorMap<T>::CColorMap(const T min_value, const T max_value)
{
    m_minValue = min_value;
    m_maxValue = max_value;
}

template<typename T>
ZD::CColorMap<T>::~CColorMap()
{
}

template<typename T>
ZD::CColorMapBlueYellow<T>::CColorMapBlueYellow(const T min_value, const T max_value) : CColorMap<T>(min_value, max_value)
{
}

template<typename T>
ZD::CColorMapBlueYellow<T>::~CColorMapBlueYellow()
{
}

template<typename T>
ZD::CPoint<T, 3> ZD::CColorMapBlueYellow<T>::GetColor(const T value)
{
    T factor;
    if (value < this->m_minValue)
        factor = 0.0;
    else if (value > this->m_maxValue)
        factor = 1.0;
    else
        factor = (value - this->m_minValue) / (this->m_maxValue - this->m_minValue);

    CPoint<unsigned char, 3> rgb;
    rgb[0] = (unsigned char)(255.0 * factor);
    rgb[1] = (unsigned char)(255.0 * factor);
    rgb[2] = (unsigned char)(255.0 * (1.0 - factor));

    CPoint<T, 3> res;
    res[0] = rgb[0] / 255.0;
    res[1] = rgb[1] / 255.0;
    res[2] = rgb[2] / 255.0;

    return res;
}

template<typename T>
ZD::CColorMapRainbow<T>::CColorMapRainbow(const T min_value, const T max_value) : CColorMap<T>(min_value, max_value)
{
}

template<typename T>
ZD::CColorMapRainbow<T>::~CColorMapRainbow()
{
}

template<typename T>
ZD::CPoint<T, 3> ZD::CColorMapRainbow<T>::GetColor(const T value)
{
    T factor;
    if (value < this->m_minValue)
        factor = 0.0;
    else if (value > this->m_maxValue)
        factor = 1.0;
    else
        factor = (value - this->m_minValue) / (this->m_maxValue - this->m_minValue);

    CPoint<T, 3> hsv;
    hsv[0] = (1.0 - factor) * 240.0;
    hsv[1] = 1.0;
    hsv[2] = 1.0;

    CPoint<unsigned char, 3> rgb;
    CColorSpaceTool<T>::HSV2RGB(hsv, rgb);

    CPoint<T, 3> res;
    res[0] = rgb[0] / 255.0;
    res[1] = rgb[1] / 255.0;
    res[2] = rgb[2] / 255.0;

    return res;
}

template<typename T>
ZD::CColorMapCubeHelix<T>::CColorMapCubeHelix(const T min_value, const T max_value) : CColorMap<T>(min_value, max_value)
{
    m_colorTable[0] = ZD::CPoint<T, 3>(0.929411765, 0.819607843, 0.796078431);
    m_colorTable[1] = ZD::CPoint<T, 3>(0.878431373, 0.694117647, 0.705882353);
    m_colorTable[2] = ZD::CPoint<T, 3>(0.811764706, 0.568627451, 0.639215686);
    m_colorTable[3] = ZD::CPoint<T, 3>(0.717647059, 0.454901961, 0.584313725);
    m_colorTable[4] = ZD::CPoint<T, 3>(0.603921569, 0.356862745, 0.533333333);
    m_colorTable[5] = ZD::CPoint<T, 3>(0.462745098, 0.266666667, 0.462745098);
    m_colorTable[6] = ZD::CPoint<T, 3>(0.317647059, 0.192156863, 0.368627451);
    m_colorTable[7] = ZD::CPoint<T, 3>(0.17254902,  0.117647059, 0.239215686);
}

template<typename T>
ZD::CColorMapCubeHelix<T>::~CColorMapCubeHelix()
{
}

template<typename T>
ZD::CPoint<T, 3> ZD::CColorMapCubeHelix<T>::GetColor(const T value)
{
    T factor;
    if (value < this->m_minValue)
        factor = 0.0;
    else if (value > this->m_maxValue)
        factor = 1.0;
    else
        factor = (value - this->m_minValue) / (this->m_maxValue - this->m_minValue);

    CPoint<T, 3> res;
    factor = factor * 7.0;
    if (factor < 1.0) {
        T f = factor;
        res = m_colorTable[0] * (1.0 - f) + m_colorTable[1] * f;
    }
    else if (factor < 2.0) {
        T f = factor - 1.0;
        res = m_colorTable[1] * (1.0 - f) + m_colorTable[2] * f;
    }
    else if (factor < 3.0) {
        T f = factor - 2.0;
        res = m_colorTable[2] * (1.0 - f) + m_colorTable[3] * f;
    }
    else if (factor < 4.0) {
        T f = factor - 3.0;
        res = m_colorTable[3] * (1.0 - f) + m_colorTable[4] * f;
    }
    else if (factor < 5.0) {
        T f = factor - 4.0;
        res = m_colorTable[4] * (1.0 - f) + m_colorTable[5] * f;
    }
    else if (factor < 6.0) {
        T f = factor - 5.0;
        res = m_colorTable[5] * (1.0 - f) + m_colorTable[6] * f;
    }
    else if (factor < 7.0) {
        T f = factor - 6.0;
        res = m_colorTable[6] * (1.0 - f) + m_colorTable[7] * f;
    }
    else {
        res = m_colorTable[7];
    }
    
    return res;
}

#endif
