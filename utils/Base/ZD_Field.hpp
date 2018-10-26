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
#ifndef _ZD_LIB_FIELD_HPP_
#define _ZD_LIB_FIELD_HPP_

#include "ZD_Point.hpp"
#include "utils/define.hpp"
#include <teem/nrrd.h>
#include <teem/biff.h>
#include <typeinfo>

namespace ZD {
    template <typename T, unsigned int N, unsigned int DIM>
    class CField {
    protected:
        CPoint<T, N> *m_pData;
        int m_size[DIM];

    public:
        CField(void);
        ~CField(void);

        int OpenNrrdFile(const char *pathname, FILE *fp = stdout);
        int SaveNrrdFile(const char *pathname, FILE *fp = stdout) const;
        int SaveNrrdFileGzip(const char *pathname, FILE *fp = stdout) const;

        int CreateField(const int nx, ...);
        int CreateField(const int size[DIM], CPoint<T, N> *data = nullptr);
        int CreateField(const int size[DIM], T *data);

        int ClearData();        // set values to 0

        int GetDimension() const;
        int GetSize(const int k) const;
        const int *GetSize() const;

        CPoint<T, N> *GetData() const;

        CPoint<T, N> GetValue(const int x, ...) const;
        inline CPoint<T, N> GetValueFast(const int idx) const;
        void SetValue(const CPoint<T, N>& v, const int x, ...);

        void MinMax(CPoint<T, N>& min_value, CPoint<T, N>& max_value) const;
        CPoint<T, N> Mean() const;
        void Normalize();
    };

    template <typename T, unsigned int N>
    class CField2 : public CField<T, N, 2> {
        typedef CField<T, N, 2> base_type;
    public:
        CField2() : base_type() {}
        inline CPoint<T, N> GetValue(const int x, const int y) const;
        CPoint<T, N> GetValue(const float x, const float y) const;
        CPoint<T, N> GetValue(const double x, const double y) const;
        CPoint<T, N> GetValue(const CPoint<int, 2>& p) const;
        CPoint<T, N> GetValue(const CPoint<int, 2> *p) const;
        CPoint<T, N> GetValue(const CPoint<float, 2>& p) const;
        CPoint<T, N> GetValue(const CPoint<float, 2> *p) const;
        CPoint<T, N> GetValue(const CPoint<double, 2>& p) const;
        CPoint<T, N> GetValue(const CPoint<double, 2> *p) const;

        void SetValue(const CPoint<int, 2>& p, const CPoint<T, N>& v);
        void SetValue(const CPoint<int, 2> *p, const CPoint<T, N> *v);
        inline void SetValue(const int x, const int y, const CPoint<T, N>& v);
        void SetValue(const int x, const int y, const CPoint<T, N> *v);

    public:
        void GaussianFilter(const T s);
    };

    template <typename T, unsigned int N>
    class CField3 : public CField<T, N, 3> {
        typedef CField<T, N, 3> base_type;
    public:
        CField3() : base_type() {}
        inline CPoint<T, N> GetValue(const int x, const int y, const int z) const;
        CPoint<T, N> GetValue(const float x, const float y, const float z) const;
        CPoint<T, N> GetValue(const double x, const double y, const double z) const;
        CPoint<T, N> GetValue(const CPoint<int, 3>& p) const;
        CPoint<T, N> GetValue(const CPoint<int, 3> *p) const;
        CPoint<T, N> GetValue(const CPoint<float, 3>& p) const;
        CPoint<T, N> GetValue(const CPoint<float, 3> *p) const;
        CPoint<T, N> GetValue(const CPoint<double, 3>& p) const;
        CPoint<T, N> GetValue(const CPoint<double, 3> *p) const;

        void SetValue(const CPoint<int, 3>& p, const CPoint<T, N>& v);
        void SetValue(const CPoint<int, 3> *p, const CPoint<T, N> *v);
        inline void SetValue(const int x, const int y, const int z, const CPoint<T, N>& v);
        void SetValue(const int x, const int y, const int z, const CPoint<T, N> *v);

    public:
        void GaussianFilter(const T s);
        //void GaussianFilterFast(const T s);

    private:
        //void BoxFilter();
    };
} // namespace ZD

template <typename T, unsigned int N, unsigned int DIM>
ZD::CField<T, N, DIM>::CField(void)
{
    this->m_pData = nullptr;
    for (unsigned int dim = 0; dim < DIM; ++dim)
        this->m_size[dim] = 0;
}

template <typename T, unsigned int N, unsigned int DIM>
ZD::CField<T, N, DIM>::~CField(void)
{
    SafeDeleteArray(this->m_pData);
    for (unsigned int dim = 0; dim < DIM; ++dim)
        this->m_size[dim] = 0;
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::OpenNrrdFile(const char *pathname, FILE *fp)
{
    if (pathname == nullptr)
        return -1;            // no input pathname

    fprintf(fp, "Reading nrrd file %s ... ", pathname);

    Nrrd *field = nrrdNew();
    if (nrrdLoad(field, pathname, NULL)) {
        fprintf(fp, "cannot load file %s:\n%s\n", pathname, biffGetDone(NRRD));
        return -2;            // cannot load file
    }

    if (field->dim == DIM && N == 1) {
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            this->m_size[dim] = (int)(field->axis[dim].size);
        }
    }
    else if ((field->dim - 1) == DIM && N != 1) {
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            this->m_size[dim] = (int)(field->axis[dim+1].size);
        }
    }
    else {
        nrrdNuke(field);
        fprintf(fp, "Dimension doesn't match (%d != %d)!\n", field->dim, DIM);
        return -2;            // cannot load file
    }

    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim) {
        size *= this->m_size[dim];
    }

    try {
        m_pData = new CPoint<T, N>[size];
    }
    catch (...) {
        fprintf(fp, "Cannot alloc memory!\n");
        nrrdNuke(field);
        return -3;        // cannot alloc memory
    }

    if (field->type == nrrdTypeFloat) {
        float *temp = (float *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N+j];
        }
    }
    else if (field->type == nrrdTypeDouble) {
        double *temp = (double *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N + j];
        }
    }
    else if (field->type == nrrdTypeInt) {
        int *temp = (int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N + j];
        }
    }
    else if (field->type == nrrdTypeUInt) {
        unsigned int *temp = (unsigned int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N + j];
        }
    }
    else if (field->type == nrrdTypeUChar) {
        unsigned char *temp = (unsigned char *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N + j];
        }
    }
    else if (field->type == nrrdTypeUShort) {
        unsigned short *temp = (unsigned short *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j)
                this->m_pData[i][j] = temp[i*N + j];
        }
    }
    else {
        fprintf(fp, "unknown data type!\n");
        return -4;
    }

    nrrdNuke(field);

    fprintf(fp, "Success.\n");
    //fprintf(fp, "Width: %d, Height: %d\n", m_width, m_height);

    return 0;
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::SaveNrrdFile(const char *pathname, FILE *fp) const
{
    if (pathname == nullptr || this->m_pData == nullptr)
        return -1;

    std::string tname = typeid(T).name();
    
    Nrrd *field = nrrdNew();
    
    int nrrdDim     = 0;
    size_t *ss      = nullptr;
    double *spacing = nullptr;
    int *kind       = nullptr;
    char **unit     = nullptr;
    if (N == 1) {
        nrrdDim = DIM;
        ss      = new size_t[DIM];
        spacing = new double[DIM];
        kind    = new int[DIM];
        unit    = new char*[DIM];
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            ss[dim]      = this->m_size[dim];
            spacing[dim] = 1.0;
            kind[dim]    = nrrdKindSpace;
            unit[dim]    = (char *)("mm");
        }
    }
    else {
        nrrdDim = DIM + 1;
        ss      = new size_t[DIM+1];
        spacing = new double[DIM+1];
        kind    = new int[DIM+1];
        unit    = new char*[DIM+1];
        ss[0]      = N;
        spacing[0] = 1.0;
        kind[0]    = nrrdKindSpace;
        unit[0]    = (char *)("mm");
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            ss[dim+1]      = this->m_size[dim];
            spacing[dim+1] = 1.0;
            kind[dim+1]    = nrrdKindList;
            unit[dim+1]    = (char *)("");
        }
        
    }

    nrrdAxisInfoSet_nva(field, nrrdAxisInfoSpacing, spacing);
    nrrdAxisInfoSet_nva(field, nrrdAxisInfoKind, kind);
    nrrdAxisInfoSet_nva(field, nrrdAxisInfoUnits, unit);

    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim) {
        size *= this->m_size[dim];
    }

    if (strcmp(tname.c_str(), "double") == 0 || strcmp(tname.c_str(), "d") == 0) {
        // double
        nrrdAlloc_nva(field, nrrdTypeDouble, nrrdDim, ss);
        double *data = (double *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "float") == 0 || strcmp(tname.c_str(), "f") == 0) {
        // float
        nrrdAlloc_nva(field, nrrdTypeFloat, nrrdDim, ss);
        float *data = (float *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "int") == 0 || strcmp(tname.c_str(), "i") == 0) {
        // int
        nrrdAlloc_nva(field, nrrdTypeInt, nrrdDim, ss);
        int *data = (int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned char") == 0 || strcmp(tname.c_str(), "uc") == 0) {
        // unsigned char
        nrrdAlloc_nva(field, nrrdTypeUChar, nrrdDim, ss);
        unsigned char *data = (unsigned char *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned int") == 0 || strcmp(tname.c_str(), "uint") == 0) {
        // unsigned int
        nrrdAlloc_nva(field, nrrdTypeUInt, nrrdDim, ss);
        unsigned int *data = (unsigned int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned short") == 0 || strcmp(tname.c_str(), "t") == 0) {
        // unsigned short
        nrrdAlloc_nva(field, nrrdTypeUShort, nrrdDim, ss);
        unsigned short *data = (unsigned short *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else {
        ;
    }

    if (nrrdSave(pathname, field, NULL)) {
        fprintf(fp, "cannot save file %s:\n%s\n", pathname, biffGetDone(NRRD));
        return -1;
    }

    nrrdNuke(field);

    SafeDeleteArray(ss);
    SafeDeleteArray(spacing);
    SafeDeleteArray(kind);
    SafeDeleteArray(unit);

    return 0;
}

template<typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::SaveNrrdFileGzip(const char * pathname, FILE * fp) const
{
    if (pathname == nullptr || this->m_pData == nullptr)
        return -1;

    std::string tname = typeid(T).name();

    Nrrd *field = nrrdNew();
    NrrdIoState *nios = nrrdIoStateNew();

    if (nrrdEncodingGzip->available()) {
        nrrdIoStateEncodingSet(nios, nrrdEncodingGzip);
        nrrdIoStateSet(nios, nrrdIoStateZlibLevel, 9);
    }
    else {
        fprintf(fp, "WARNING: TEEM library does not support Gzip compression encoding.\n");
    }

    int nrrdDim = 0;
    size_t *ss = nullptr;
    double *spacing = nullptr;
    int *kind = nullptr;
    char **unit = nullptr;
    if (N == 1) {
        nrrdDim = DIM;
        ss = new size_t[DIM];
        spacing = new double[DIM];
        kind = new int[DIM];
        unit = new char*[DIM];
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            ss[dim] = this->m_size[dim];
            spacing[dim] = 1.0;
            kind[dim] = nrrdKindSpace;
            unit[dim] = (char *)("mm");
        }
    }
    else {
        nrrdDim = DIM + 1;
        ss = new size_t[DIM + 1];
        spacing = new double[DIM + 1];
        kind = new int[DIM + 1];
        unit = new char*[DIM + 1];
        ss[0] = N;
        spacing[0] = 1.0;
        kind[0] = nrrdKindSpace;
        unit[0] = (char *)("mm");
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            ss[dim + 1] = this->m_size[dim];
            spacing[dim + 1] = 1.0;
            kind[dim + 1] = nrrdKindList;
            unit[dim + 1] = (char *)("mm");
        }

    }

    nrrdAxisInfoSet_nva(field, nrrdAxisInfoSpacing, spacing);
    nrrdAxisInfoSet_nva(field, nrrdAxisInfoKind, kind);
    nrrdAxisInfoSet_nva(field, nrrdAxisInfoUnits, unit);

    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim) {
        size *= this->m_size[dim];
    }

    if (strcmp(tname.c_str(), "double") == 0 || strcmp(tname.c_str(), "d") == 0) {
        // double
        nrrdAlloc_nva(field, nrrdTypeDouble, nrrdDim, ss);
        double *data = (double *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "float") == 0 || strcmp(tname.c_str(), "f") == 0) {
        // float
        nrrdAlloc_nva(field, nrrdTypeFloat, nrrdDim, ss);
        float *data = (float *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "int") == 0 || strcmp(tname.c_str(), "i") == 0) {
        // int
        nrrdAlloc_nva(field, nrrdTypeInt, nrrdDim, ss);
        int *data = (int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned char") == 0 || strcmp(tname.c_str(), "uc") == 0) {
        // unsigned char
        nrrdAlloc_nva(field, nrrdTypeUChar, nrrdDim, ss);
        unsigned char *data = (unsigned char *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned int") == 0 || strcmp(tname.c_str(), "uint") == 0) {
        // unsigned int
        nrrdAlloc_nva(field, nrrdTypeUInt, nrrdDim, ss);
        unsigned int *data = (unsigned int *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else if (strcmp(tname.c_str(), "unsigned short") == 0 || strcmp(tname.c_str(), "t") == 0) {
        // unsigned short
        nrrdAlloc_nva(field, nrrdTypeUShort, nrrdDim, ss);
        unsigned short *data = (unsigned short *)(field->data);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i*N + j] = this->m_pData[i][j];
            }
        }
    }
    else {
        ;
    }

    if (nrrdSave(pathname, field, nios)) {
        fprintf(fp, "cannot save file %s:\n%s\n", pathname, biffGetDone(NRRD));
        return -1;
    }

    nrrdNuke(field);
    nrrdIoStateNix(nios);

    SafeDeleteArray(ss);
    SafeDeleteArray(spacing);
    SafeDeleteArray(kind);
    SafeDeleteArray(unit);

    return 0;
}

template<typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::CreateField(const int nx, ...)
{
    int size[DIM];
    size[0] = nx;
    va_list ap;
    va_start(ap, nx);
    for (unsigned dim = 1; dim < DIM; ++dim)
        size[dim] = va_arg(ap, int);
    CPoint<T, N> *data = (CPoint<T, N> *)va_arg(ap, void*);
    va_end(ap);

    return this->CreateField(size, data);
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::CreateField(const int size[DIM], CPoint<T, N> *data)
{
    if (this->m_pData != nullptr)
        return 0;

    int element_count = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim) {
        this->m_size[dim] = size[dim];
        element_count *= this->m_size[dim];
    }

    try {
        this->m_pData = new CPoint<T, N>[element_count];
    }
    catch (...) {
        fprintf(stderr, "Cannot alloc memory!\n");
        return -3;        // cannot alloc memory
    }

    if (data != nullptr) {
        memcpy(this->m_pData, data, sizeof(CPoint<T, N>)*element_count);
    }
    else {
        memset(this->m_pData, 0, sizeof(CPoint<T, N>)*element_count);
    }

    return 0;
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::CreateField(const int size[DIM], T *data)
{
    if (this->m_pData != nullptr)
        return 0;

    int element_count = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim) {
        this->m_size[dim] = size[dim];
        element_count *= this->m_size[dim];
    }

    try {
        this->m_pData = new CPoint<T, N>[element_count];
    }
    catch (...) {
        fprintf(stderr, "Cannot alloc memory!\n");
        return -3;        // cannot alloc memory
    }

    if (data != nullptr) {
        memcpy(this->m_pData, data, sizeof(CPoint<T, N>)*element_count);
    }
    else {
        memset(this->m_pData, 0, sizeof(CPoint<T, N>)*element_count);
    }

    return 0;
}

template<typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::ClearData()
{
    if (this->m_pData) {
        int element_count = 1;
        for (unsigned int dim = 0; dim < DIM; ++dim) {
            element_count *= this->m_size[dim];
        }
        memset(this->m_pData, 0, sizeof(CPoint<T, N>)*element_count);
    }
    return 0;
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::GetDimension() const
{
    return DIM;
}

template <typename T, unsigned int N, unsigned int DIM>
int ZD::CField<T, N, DIM>::GetSize(const int k) const
{
    return this->m_size[k];
}

template <typename T, unsigned int N, unsigned int DIM>
const int* ZD::CField<T, N, DIM>::GetSize() const
{
    return this->m_size;
}

template <typename T, unsigned int N, unsigned int DIM>
ZD::CPoint<T, N> *ZD::CField<T, N, DIM>::GetData() const
{
    return this->m_pData;
}

template <typename T, unsigned int N, unsigned int DIM>
ZD::CPoint<T, N> ZD::CField<T, N, DIM>::GetValue(const int x, ...) const
{
    int pos[DIM];
    pos[0] = x;
    va_list ap;
    va_start(ap, x);
    for (unsigned dim = 1; dim < DIM; ++dim)
        pos[dim] = va_arg(ap, int);
    va_end(ap);

    CPoint<T, N> res;
    res.SetZero();
    int index = 0;
    for (unsigned dim = 0; dim < DIM; ++dim) {
        if (pos[dim] < 0 || pos[dim] > (this->m_size[dim] - 1))
            return res;
        index = index * this->m_size[DIM - dim - 1] + pos[DIM - dim - 1];
    }

    return this->m_pData[index];
}

template <typename T, unsigned int N, unsigned int DIM>
inline ZD::CPoint<T, N> ZD::CField<T, N, DIM>::GetValueFast(const int idx) const
{
    return this->m_pData[idx];
}

template <typename T, unsigned int N, unsigned int DIM>
void ZD::CField<T, N, DIM>::SetValue(const CPoint<T, N>& v, const int x, ...)
{
    int pos[DIM];
    pos[0] = x;
    va_list ap;
    va_start(ap, x);
    for (unsigned dim = 1; dim < DIM; ++dim)
        pos[dim] = va_arg(ap, int);
    va_end(ap);

    CPoint<T, N> res;
    res.SetZero();
    int index = 0;
    for (unsigned dim = 0; dim < DIM; ++dim) {
        if (pos[dim] < 0 || pos[dim] > (this->m_size[dim] - 1))
            return;
        index = index * this->m_size[DIM - dim - 1] + pos[DIM - dim - 1];
    }

    this->m_pData[index] = v;
}

template <typename T, unsigned int N, unsigned int DIM>
void ZD::CField<T, N, DIM>::MinMax(CPoint<T, N>& min_value, CPoint<T, N>& max_value) const
{
    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim)
        size = size * this->m_size[dim];

    min_value = this->m_pData[0];
    max_value = this->m_pData[0];

    for (int i = 1; i < size; ++i) {
        for (int k = 0; k < N; ++k) {
            if (this->m_pData[i][k] < min_value[k])
                min_value[k] = this->m_pData[i][k];
            if (this->m_pData[i][k] > max_value[k])
                max_value[k] = this->m_pData[i][k];
        }
    }
}

template <typename T, unsigned int N, unsigned int DIM>
ZD::CPoint<T, N> ZD::CField<T, N, DIM>::Mean() const
{
    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim)
        size = size * this->m_size[dim];

    CPoint<T, N> mean;
    mean.SetZero();
    for (int i = 0; i < size; ++i) {
        mean += this->m_pData[i];
    }
    mean /= T(size);
    return mean;
}

template <typename T, unsigned int N, unsigned int DIM>
void ZD::CField<T, N, DIM>::Normalize()
{
    CPoint<T, N> min_value, max_value;
    this->MinMax(min_value, max_value);

    int size = 1;
    for (unsigned int dim = 0; dim < DIM; ++dim)
        size = size * this->m_size[dim];

    for (int i = 0; i < size; ++i) {
        for (int k = 0; k < N; ++k) {
            this->m_pData[i][k] = (this->m_pData[i][k] - min_value[k]) / (max_value[k] - min_value[k]);
        }
    }
}

template <typename T, unsigned int N>
inline ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const int x, const int y) const
{
    CPoint<T, N> res;
    res.SetZero();

    if (x < 0 || x > (this->m_size[0] - 1) ||
        y < 0 || y > (this->m_size[1] - 1)) {
        return res;
    }

    int index = y * this->m_size[0] + x;
    return this->m_pData[index];
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const float x, const float y) const
{
    CPoint<T, N> res;
    res.SetZero();

    if (x < -0.5 || x >(float)(this->m_size[0] - 0.5) ||
        y < -0.5 || y >(float)(this->m_size[1] - 0.5)) {
        return res;
    }

    int x1 = (int)floor(x);
    int y1 = (int)floor(y);

    int x2 = (int)ceil(x);
    int y2 = (int)ceil(y);

    float a = x - (float)x1;
    float b = y - (float)y1;

    CPoint<T, N> v1 = this->GetValue(x1, y1);
    CPoint<T, N> v2 = this->GetValue(x2, y1);
    CPoint<T, N> v3 = this->GetValue(x2, y2);
    CPoint<T, N> v4 = this->GetValue(x1, y2);

    res = v1 + a * (v2 - v1) + b * (v4 - v1) + a * b * (v1 - v2 + v3 - v4);

    return res;
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const double x, const double y) const
{
    CPoint<T, N> res;
    res.SetZero();

    if (x < -0.5 || x >(double)(this->m_size[0] - 0.5) ||
        y < -0.5 || y >(double)(this->m_size[1] - 0.5)) {
        return res;
    }

    int x1 = (int)floor(x);
    int y1 = (int)floor(y);

    int x2 = (int)ceil(x);
    int y2 = (int)ceil(y);

    double a = x - (double)x1;
    double b = y - (double)y1;

    CPoint<T, N> v1 = this->GetValue(x1, y1);
    CPoint<T, N> v2 = this->GetValue(x2, y1);
    CPoint<T, N> v3 = this->GetValue(x2, y2);
    CPoint<T, N> v4 = this->GetValue(x1, y2);

    res = v1 + a * (v2 - v1) + b * (v4 - v1) + a * b * (v1 - v2 + v3 - v4);

    return res;
}


template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<int, 2>& p) const
{
    return this->GetValue(p[0], p[1]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<int, 2> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<float, 2>& p) const
{
    return this->GetValue(p[0], p[1]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<float, 2> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<double, 2>& p) const
{
    return this->GetValue(p[0], p[1]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField2<T, N>::GetValue(const CPoint<double, 2> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1]);
}

template <typename T, unsigned int N>
void ZD::CField2<T, N>::SetValue(const int x, const int y, const CPoint<T, N>& v)
{
    if (x < 0 || x >(this->m_size[0] - 1) ||
        y < 0 || y >(this->m_size[1] - 1)) {
        return;
    }

    int index = y * this->m_size[0] + x;
    this->m_pData[index] = v;
}

template <typename T, unsigned int N>
void ZD::CField2<T, N>::SetValue(const int x, const int y, const CPoint<T, N> *v)
{
    this->SetValue(x, y, *v);
}

template <typename T, unsigned int N>
void ZD::CField2<T, N>::SetValue(const CPoint<int, 2>& p, const CPoint<T, N>& v)
{
    this->SetValue(p[0], p[1], v);
}

template <typename T, unsigned int N>
void ZD::CField2<T, N>::SetValue(const CPoint<int, 2> *p, const CPoint<T, N> *v)
{
    this->SetValue(p->m_data[0], p->m_data[1], *v);
}

template<typename T, unsigned int N>
void ZD::CField2<T, N>::GaussianFilter(const T s)
{
    int size = ceil(s * 5.0);

    T *kernel = new T[(size * 2 + 1)*(size * 2 + 1)];
    for (int y = -size; y <= size; ++y) {
        for (int x = -size; x <= size; ++x) {
            T temp = -(x * x + y * y) / (2.0 * s * s);
            kernel[(y + size)*(size * 2 + 1) + (x + size)] = 1.0 / (2.0 * ZD_PI * s * s) * std::exp((double)temp);
        }
    }

    CPoint<T, N> *results = new CPoint<T, N>[this->m_size[0]*this->m_size[1]];
    memset(results, 0, sizeof(CPoint<T, N>)*this->m_size[0]*this->m_size[1]);
    for (int y = 0; y < this->m_size[1]; ++y) {
        for (int x = 0; x < this->m_size[0]; ++x) {
            int min_x = (x - size) > 0 ? (x - size) : 0;
            int min_y = (y - size) > 0 ? (y - size) : 0;
            int max_x = (x + size) < this->m_size[0] ? (x + size) : (this->m_size[0] - 1);
            int max_y = (y + size) < this->m_size[1] ? (y + size) : (this->m_size[1] - 1);

            T weight = 0.0;
            for (int cy = min_y; cy <= max_y; ++cy) {
                for (int cx = min_x; cx <= max_x; ++cx) {
                    int ix = cx - x + size;
                    int iy = cy - y + size;
                    for (int k = 0; k < N; ++k) {
                        results[y*this->m_size[0]+x][k] += this->m_pData[cy*this->m_size[0]+cx][k] * kernel[iy*(size*2+1)+ix];
                    }
                    
                    weight += kernel[iy*(size*2+1)+ix];
                }
            }

            for (int k = 0; k < N; ++k) {
                results[y*this->m_size[0]+x] /= weight;
            }
        }
    }


    memcpy(this->m_pData, results, sizeof(CPoint<T, N>)*this->m_size[0]*this->m_size[1]);

    SafeDeleteArray(results);
    SafeDeleteArray(kernel);
}

template <typename T, unsigned int N>
inline ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const int x, const int y, const int z) const
{
    CPoint<T, N> res;
    res.SetZero();

    if (x < 0 || x >(this->m_size[0] - 1) ||
        y < 0 || y >(this->m_size[1] - 1) || 
        z < 0 || z >(this->m_size[2] - 1)) {
        return res;
    }

    int index = (z * this->m_size[1] + y) * this->m_size[0] + x;
    return this->m_pData[index];
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const float x, const float y, const float z) const
{
    CPoint<T, N> res;
    res.SetZero();

    //if (x < -0.5 || x > (float)(this->m_size[0] - 0.5) ||
    //    y < -0.5 || y > (float)(this->m_size[1] - 0.5) ||
    //    z < -0.5 || z > (float)(this->m_size[2] - 0.5)) {
    //    return res;
    //}
    if (x < 0.0 || x > (float)(this->m_size[0] - 1.0) ||
        y < 0.0 || y > (float)(this->m_size[1] - 1.0) ||
        z < 0.0 || z > (float)(this->m_size[2] - 1.0)) {
        return res;
    }

    int x1 = (int)std::floor(x);
    int y1 = (int)std::floor(y);
    int z1 = (int)std::floor(z);

    int x2 = (int)std::ceil(x);
    int y2 = (int)std::ceil(y);
    int z2 = (int)std::ceil(z);

    float a = x - (float)x1;
    float b = y - (float)y1;
    float c = z - (float)z1;

    //CPoint<T, N> v1 = this->GetValue(x1, y1, z1);
    //CPoint<T, N> v2 = this->GetValue(x2, y1, z1);
    //CPoint<T, N> v3 = this->GetValue(x2, y2, z1);
    //CPoint<T, N> v4 = this->GetValue(x1, y2, z1);
    //CPoint<T, N> v5 = this->GetValue(x1, y1, z2);
    //CPoint<T, N> v6 = this->GetValue(x2, y1, z2);
    //CPoint<T, N> v7 = this->GetValue(x2, y2, z2);
    //CPoint<T, N> v8 = this->GetValue(x1, y2, z2);

    const CPoint<T, N>& v1 = this->m_pData[(z1 * this->m_size[1] + y1) * this->m_size[0] + x1];
    const CPoint<T, N>& v2 = this->m_pData[(z1 * this->m_size[1] + y1) * this->m_size[0] + x2];
    const CPoint<T, N>& v3 = this->m_pData[(z1 * this->m_size[1] + y2) * this->m_size[0] + x2];
    const CPoint<T, N>& v4 = this->m_pData[(z1 * this->m_size[1] + y2) * this->m_size[0] + x1];
    const CPoint<T, N>& v5 = this->m_pData[(z2 * this->m_size[1] + y1) * this->m_size[0] + x1];
    const CPoint<T, N>& v6 = this->m_pData[(z2 * this->m_size[1] + y1) * this->m_size[0] + x2];
    const CPoint<T, N>& v7 = this->m_pData[(z2 * this->m_size[1] + y2) * this->m_size[0] + x2];
    const CPoint<T, N>& v8 = this->m_pData[(z2 * this->m_size[1] + y2) * this->m_size[0] + x1];

    //res = v1 + a * (v2 - v1) + b * (v4 - v1) + c * (v5 - v1) +
    //    a * b * (v1 - v2 + v3 - v4) +
    //    a * c * (v1 - v2 + v6 - v5) +
    //    b * c * (v1 - v4 + v8 - v5) +
    //    a * b * c * (-v1 + v2 - v3 + v4 + v5 - v6 + v7 - v8);

    for (unsigned int i = 0; i < N; ++i) {
        res.m_data[i] = v1.m_data[i] + a * (v2.m_data[i] - v1.m_data[i]) + b * (v4.m_data[i] - v1.m_data[i]) + c * (v5.m_data[i] - v1.m_data[i]) +
            a * b * (v1.m_data[i] - v2.m_data[i] + v3.m_data[i] - v4.m_data[i]) +
            a * c * (v1.m_data[i] - v2.m_data[i] + v6.m_data[i] - v5.m_data[i]) +
            b * c * (v1.m_data[i] - v4.m_data[i] + v8.m_data[i] - v5.m_data[i]) +
            a * b * c * (-v1.m_data[i] + v2.m_data[i] - v3.m_data[i] + v4.m_data[i] + v5.m_data[i] - v6.m_data[i] + v7.m_data[i] - v8.m_data[i]);
    }
    

    return res;
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const double x, const double y, const double z) const
{
    CPoint<T, N> res;
    res.SetZero();

    //if (x < -0.5 || x > (double)(this->m_size[0] - 0.5) ||
    //    y < -0.5 || y > (double)(this->m_size[1] - 0.5) ||
    //    z < -0.5 || z > (double)(this->m_size[2] - 0.5)) {
    //    return res;
    //}
    if (x < 0.0 || x > (double)(this->m_size[0] - 1.0) ||
        y < 0.0 || y > (double)(this->m_size[1] - 1.0) ||
        z < 0.0 || z > (double)(this->m_size[2] - 1.0)) {
        return res;
    }

    int x1 = (int)std::floor(x);
    int y1 = (int)std::floor(y);
    int z1 = (int)std::floor(z);

    int x2 = (int)std::ceil(x);
    int y2 = (int)std::ceil(y);
    int z2 = (int)std::ceil(z);

    double a = x - (double)x1;
    double b = y - (double)y1;
    double c = z - (double)z1;

    //CPoint<T, N> v1 = this->GetValue(x1, y1, z1);
    //CPoint<T, N> v2 = this->GetValue(x2, y1, z1);
    //CPoint<T, N> v3 = this->GetValue(x2, y2, z1);
    //CPoint<T, N> v4 = this->GetValue(x1, y2, z1);
    //CPoint<T, N> v5 = this->GetValue(x1, y1, z2);
    //CPoint<T, N> v6 = this->GetValue(x2, y1, z2);
    //CPoint<T, N> v7 = this->GetValue(x2, y2, z2);
    //CPoint<T, N> v8 = this->GetValue(x1, y2, z2);

    //res = v1 + a * (v2 - v1) + b * (v4 - v1) + c * (v5 - v1) +
    //    a * b * (v1 - v2 + v3 - v4) +
    //    a * c * (v1 - v2 + v6 - v5) +
    //    b * c * (v1 - v4 + v8 - v5) +
    //    a * b * c * (-v1 + v2 - v3 + v4 + v5 - v6 + v7 - v8);

    const CPoint<T, N>& v1 = this->m_pData[(z1 * this->m_size[1] + y1) * this->m_size[0] + x1];
    const CPoint<T, N>& v2 = this->m_pData[(z1 * this->m_size[1] + y1) * this->m_size[0] + x2];
    const CPoint<T, N>& v3 = this->m_pData[(z1 * this->m_size[1] + y2) * this->m_size[0] + x2];
    const CPoint<T, N>& v4 = this->m_pData[(z1 * this->m_size[1] + y2) * this->m_size[0] + x1];
    const CPoint<T, N>& v5 = this->m_pData[(z2 * this->m_size[1] + y1) * this->m_size[0] + x1];
    const CPoint<T, N>& v6 = this->m_pData[(z2 * this->m_size[1] + y1) * this->m_size[0] + x2];
    const CPoint<T, N>& v7 = this->m_pData[(z2 * this->m_size[1] + y2) * this->m_size[0] + x2];
    const CPoint<T, N>& v8 = this->m_pData[(z2 * this->m_size[1] + y2) * this->m_size[0] + x1];

    for (unsigned int i = 0; i < N; ++i) {
        res.m_data[i] = v1.m_data[i] + a * (v2.m_data[i] - v1.m_data[i]) + b * (v4.m_data[i] - v1.m_data[i]) + c * (v5.m_data[i] - v1.m_data[i]) +
            a * b * (v1.m_data[i] - v2.m_data[i] + v3.m_data[i] - v4.m_data[i]) +
            a * c * (v1.m_data[i] - v2.m_data[i] + v6.m_data[i] - v5.m_data[i]) +
            b * c * (v1.m_data[i] - v4.m_data[i] + v8.m_data[i] - v5.m_data[i]) +
            a * b * c * (-v1.m_data[i] + v2.m_data[i] - v3.m_data[i] + v4.m_data[i] + v5.m_data[i] - v6.m_data[i] + v7.m_data[i] - v8.m_data[i]);
    }


    return res;
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<int, 3>& p) const
{
    return this->GetValue(p[0], p[1], p[2]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<int, 3> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1], p->m_data[2]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<float, 3>& p) const
{
    return this->GetValue(p[0], p[1], p[2]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<float, 3> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1], p->m_data[2]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<double, 3>& p) const
{
    return this->GetValue(p[0], p[1], p[2]);
}

template <typename T, unsigned int N>
ZD::CPoint<T, N> ZD::CField3<T, N>::GetValue(const CPoint<double, 3> *p) const
{
    return this->GetValue(p->m_data[0], p->m_data[1], p->m_data[2]);
}

template <typename T, unsigned int N>
void ZD::CField3<T, N>::SetValue(const int x, const int y, const int z, const CPoint<T, N>& v)
{
    if (x < 0 || x >(this->m_size[0] - 1) ||
        y < 0 || y >(this->m_size[1] - 1) ||
        z < 0 || z >(this->m_size[2] - 1)) {
        return;
    }

    int index = (z * this->m_size[1] + y) * this->m_size[0] + x;
    this->m_pData[index] = v;
}

template <typename T, unsigned int N>
void ZD::CField3<T, N>::SetValue(const int x, const int y, const int z, const CPoint<T, N> *v)
{
    this->SetValue(x, y, z, *v);
}

template <typename T, unsigned int N>
void ZD::CField3<T, N>::SetValue(const CPoint<int, 3>& p, const CPoint<T, N>& v)
{
    this->SetValue(p[0], p[1], p[2], v);
}

template <typename T, unsigned int N>
void ZD::CField3<T, N>::SetValue(const CPoint<int, 3> *p, const CPoint<T, N> *v)
{
    this->SetValue(p->m_data[0], p->m_data[1], p->m_data[2], *v);
}

template<typename T, unsigned int N>
void ZD::CField3<T, N>::GaussianFilter(const T s)
{
    // 3D gaussian filter
    int size = std::ceil(s * 5.0);
    int ss = size * 2 + 1;
    T *kernel = new T[ss];

    for (int x = -size; x <= size; ++x) {
        T temp = -(x * x) / (2.0 * s * s);
        kernel[x + size] = 1.0 / (2.0 * ZD_PI * s * s) * std::exp((double)temp);
    }
    CPoint<T, N> *results = new CPoint<T, N>[this->m_size[0]*this->m_size[1]*this->m_size[2]];
    
    for (int z = 0; z < this->m_size[2]; ++z) {
        for (int y = 0; y < this->m_size[1]; ++y) {
            for (int x = 0; x < this->m_size[0]; ++x) {
                results[(z*this->m_size[1] + y)*this->m_size[0] + x].SetZero();
                T weight = 0.0;
                for (int ix = x - size; ix <= x + size; ++ix) {
                    int xx = ZD_MIN(this->m_size[0]-1, ZD_MAX(0, ix));
                    results[(z*this->m_size[1]+y)*this->m_size[0]+x] += this->m_pData[(z*this->m_size[1]+y)*this->m_size[0]+xx] * kernel[ix+size-x];
                    weight += kernel[ix+size-x];
                }
                results[(z*this->m_size[1]+y)*this->m_size[0]+x] /= weight;
            }
        }
    }
    memcpy(this->m_pData, results, sizeof(CPoint<T, N>)*this->m_size[0]*this->m_size[1]*this->m_size[2]);

    for (int z = 0; z < this->m_size[2]; ++z) {
        for (int y = 0; y < this->m_size[1]; ++y) {
            for (int x = 0; x < this->m_size[0]; ++x) {
                results[(z*this->m_size[1] + y)*this->m_size[0] + x].SetZero();
                T weight = 0.0;
                for (int iy = y - size; iy <= y + size; ++iy) {
                    int yy = ZD_MIN(this->m_size[1] - 1, ZD_MAX(0, iy));
                    results[(z*this->m_size[1]+y)*this->m_size[0]+x] += this->m_pData[(z*this->m_size[1]+yy)*this->m_size[0]+x] * kernel[iy + size - y];
                    weight += kernel[iy + size - y];
                }
                results[(z*this->m_size[1]+y)*this->m_size[0]+x] /= weight;
            }
        }
    }
    memcpy(this->m_pData, results, sizeof(CPoint<T, N>)*this->m_size[0]*this->m_size[1]*this->m_size[2]);

    for (int z = 0; z < this->m_size[2]; ++z) {
        for (int y = 0; y < this->m_size[1]; ++y) {
            for (int x = 0; x < this->m_size[0]; ++x) {
                results[(z*this->m_size[1] + y)*this->m_size[0] + x].SetZero();
                T weight = 0.0;
                for (int iz = z - size; iz <= z + size; ++iz) {
                    int zz = ZD_MIN(this->m_size[2] - 1, ZD_MAX(0, iz));
                    results[(z*this->m_size[1]+y)*this->m_size[0]+x] += this->m_pData[(zz*this->m_size[1]+y)*this->m_size[0]+x] * kernel[iz + size - z];
                    weight += kernel[iz + size - z];
                }
                results[(z*this->m_size[1]+y)*this->m_size[0]+x] /= weight;
            }
        }
    }
    memcpy(this->m_pData, results, sizeof(CPoint<T, N>)*this->m_size[0]*this->m_size[1]*this->m_size[2]);

    SafeDeleteArray(results);
    SafeDeleteArray(kernel);
}

#endif // _ZD_LIB_FIELD_HPP_
