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
#ifndef _ZD_LIB_LINE_ATTRIBUTE_HPP_
#define _ZD_LIB_LINE_ATTRIBUTE_HPP_

#include "ZD_Line.hpp"
#include "utils/define.hpp"

namespace ZD {
    template <typename T, unsigned int N>
    class CLineAttribute {
    public:
        CPoint<T, N> *m_pAttributes;
        int m_count;

    public:
        CLineAttribute();
        ~CLineAttribute();

    public:
        void CreateLineAttribute(const CPoint<T, N> *attributes, const int count);
        void CreateLineAttribute(const std::vector<CPoint<T, N>>& attributes);

    public:
        friend void ReadLineAttributes(const char *pathname, CLineAttribute<T, N> **lineAttributes, int &size)
        {
            FILE *fp = fopen(pathname, "rb");
            if (fp == nullptr)
                return;

            fread(&size, sizeof(int), 1, fp);
            try {
                *lineAttributes = new CLineAttribute<T, N>[size];
            }
            catch (...) {
                std::cerr << "cannot alloc memory!" << std::endl;
                fclose(fp);
                return;
            }
            for (int i = 0; i < size; ++i) {
                int count;
                fread(&count, sizeof(int), 1, fp);
                (*lineAttributes)[i].m_count = count;
                (*lineAttributes)[i].m_pAttributes = new CPoint<T, N>[count];
                fread((*lineAttributes)[i].m_pAttributes, sizeof(CPoint<T, N>), count, fp);
            }
            fclose(fp);
        }

        friend void SaveLineAttributes(const char *pathname, const CLineAttribute<T, N> *lineAttributes, const int size)
        {
            FILE *fp = fopen(pathname, "wb");
            if (fp == nullptr)
                return;

            fwrite(&size, sizeof(int), 1, fp);
            for (int i = 0; i < size; ++i) {
                fwrite(&(lineAttributes[i].m_count), sizeof(int), 1, fp);
                fwrite(lineAttributes[i].m_pAttributes, sizeof(CPoint<T, N>), lineAttributes[i].m_count, fp);
            }
            fclose(fp);
        }
    };
} // namespace ZD


template<typename T, unsigned int N>
ZD::CLineAttribute<T, N>::CLineAttribute()
{
    this->m_pAttributes = nullptr;
    this->m_count = 0;
}

template<typename T, unsigned int N>
ZD::CLineAttribute<T, N>::~CLineAttribute()
{
    SafeDeleteArray(this->m_pAttributes);
    this->m_count = 0;
}

template<typename T, unsigned int N>
void ZD::CLineAttribute<T, N>::CreateLineAttribute(const CPoint<T, N>* attributes, const int count)
{
    this->m_count = count;
    this->m_pAttributes = new CPoint<T, N>[this->m_count];
    memcpy(this->m_pAttributes, attributes, sizeof(CPoint<T, N>)*this->m_count);
}

template<typename T, unsigned int N>
void ZD::CLineAttribute<T, N>::CreateLineAttribute(const std::vector<CPoint<T, N>>& attributes)
{
    this->m_count = attributes.size();
    this->m_pAttributes = new CPoint<T, N>[this->m_count];
    int i = 0;
    for (auto iter = attributes.begin(); iter != attributes.end(); iter++) {
        this->m_pAttributes[i++] = *iter;
    }
}

#endif
