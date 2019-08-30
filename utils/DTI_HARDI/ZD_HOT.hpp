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
#ifndef _ZD_LIB_HOT_HPP_
#define _ZD_LIB_HOT_HPP_

#include "utils/Base/ZD_Field.hpp"
#include "utils/define.hpp"
#include "ZD_Glyph.hpp"

#include <teem/tijk.h>


namespace ZD {
    template <typename T>
    class CHOT {
    private:
        CField3<T, 15> *m_pHOT4Volume;        // 4th order
        CField3<T, 28> *m_pHOT6Volume;        // 6th order
        unsigned int m_order;

    public:
        //static const int factor[2][4] = {
        //    {1,4,6,4},
        //    {1,6,15,20}
        //}

    public:
        CHOT(const int w, const int h, const int d, const unsigned int order);
        CHOT(const char *hotPathname);
        ~CHOT();

    public:
        const unsigned int GetOrder() const;
        const int GetSize(const int k) const;

        void SaveNrrdFile(const char *pathname) const;

    public:
        void GetHOT(const CPoint<int, 3> p, T *hot) const;
        void GetHOT(const int x, const int y, const int z, T *hot) const;
        void GetHOT(const CPoint<T, 3> p, T *hot) const;
        void GetHOT(const T x, const T y, const T z, T *hot) const;

    public:
        void SetHOT(const CPoint<int, 3> p, T *hot);
        void SetHOT(const int x, const int y, const int z, T *hot);

    public:
        static void SH2HOT(const T *sh, T *hot, const unsigned int order);
        static void HOT2SH(const T *hot, T *sh, const unsigned int order);

        static const int HOT2Directions(const T *hot, CPoint<T, 3> *dirs, const unsigned int order);

        static const int ModelSelection(const T *hot, const unsigned int order);

        static const T Norm(const T *hot, const unsigned int order);

        static T RankKDecomposition(const T *hot, const unsigned int order, 
            const unsigned int k, CPoint<T, 3> *vectors, T *weights);

        static T MaxAxis(const T *hot, const unsigned int order);

    };
}

template<typename T>
ZD::CHOT<T>::CHOT(const int w, const int h, const int d, const unsigned int order)
{
    if (order == 4) {
        m_pHOT4Volume = new CField3<T, 15>();
        m_pHOT4Volume->CreateField(w, h, d, nullptr);
        m_order = 4;
    }
    else if (order == 6) {
        m_pHOT6Volume = new CField3<T, 28>();
        m_pHOT6Volume->CreateField(w, h, d, nullptr);
        m_order = 6;
    }
    else {
        ;
    }
}

template<typename T>
ZD::CHOT<T>::CHOT(const char * hotPathname)
{
    m_pHOT4Volume = nullptr;
    m_pHOT6Volume = nullptr;
    m_order = 0;

    Nrrd *field = nrrdNew();
    if (nrrdLoad(field, hotPathname, NULL)) {
        printf("cannot load file %s:\n%s\n", hotPathname, biffGetDone(NRRD));
        exit(-1);
    }
    unsigned int n = field->axis[0].size;
    nrrdNuke(field);

    if (n == 15) {
        m_pHOT4Volume = new CField3<T, 15>();
        m_pHOT4Volume->OpenNrrdFile(hotPathname);
        m_order = 4;
    }
    else if (n == 28) {
        m_pHOT6Volume = new CField3<T, 28>();
        m_pHOT6Volume->OpenNrrdFile(hotPathname);
        m_order = 6;
    }
    else {
        printf("unknown order. (n=%d).\n", n);
        exit(-1);
    }
}

template<typename T>
ZD::CHOT<T>::~CHOT()
{
    SafeDelete(m_pHOT4Volume);
    SafeDelete(m_pHOT6Volume);
    m_order = 0;
}

template<typename T>
const unsigned int ZD::CHOT<T>::GetOrder() const
{
    return m_order;
}

template<typename T>
const int ZD::CHOT<T>::GetSize(const int k) const
{
    if (m_order == 4)
        return m_pHOT4Volume->GetSize(k);
    else if (m_order == 6)
        return m_pHOT6Volume->GetSize(k);
    else
        return -1;
}


template<typename T>
void ZD::CHOT<T>::SaveNrrdFile(const char * pathname) const
{
    if (m_order == 4) {
        m_pHOT4Volume->SaveNrrdFile(pathname);
    }
    else if (m_order == 6) {
        m_pHOT6Volume->SaveNrrdFile(pathname);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CHOT<T>::GetHOT(const CPoint<int, 3> p, T * hot) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pHOT4Volume->GetValue(p);
        memcpy(hot, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pHOT6Volume->GetValue(p);
        memcpy(hot, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CHOT<T>::GetHOT(const int x, const int y, const int z, T * hot) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pHOT4Volume->GetValue(x, y, z);
        memcpy(hot, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pHOT6Volume->GetValue(x, y, z);
        memcpy(hot, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CHOT<T>::GetHOT(const CPoint<T, 3> p, T * hot) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pHOT4Volume->GetValue(p);
        memcpy(hot, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pHOT6Volume->GetValue(p);
        memcpy(hot, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CHOT<T>::GetHOT(const T x, const T y, const T z, T * hot) const
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        tmp = m_pHOT4Volume->GetValue(x, y, z);
        memcpy(hot, tmp.m_data, sizeof(T) * 15);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        tmp = m_pHOT6Volume->GetValue(x, y, z);
        memcpy(hot, tmp.m_data, sizeof(T) * 28);
    }
    else {
        ;
    }
}


template<typename T>
void ZD::CHOT<T>::SetHOT(const CPoint<int, 3> p, T * hot)
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        memcpy(tmp.m_data, hot, sizeof(T) * 15);
        m_pHOT4Volume->SetValue(p, tmp);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        memcpy(tmp.m_data, hot, sizeof(T) * 28);
        m_pHOT6Volume->SetValue(p, tmp);
    }
    else {
        ;
    }
}

template<typename T>
void ZD::CHOT<T>::SetHOT(const int x, const int y, const int z, T * hot)
{
    if (m_order == 4) {
        CPoint<T, 15> tmp;
        memcpy(tmp.m_data, hot, sizeof(T) * 15);
        m_pHOT4Volume->SetValue(x, y, z, tmp);
    }
    else if (m_order == 6) {
        CPoint<T, 28> tmp;
        memcpy(tmp.m_data, hot, sizeof(T) * 28);
        m_pHOT6Volume->SetValue(x, y, z, tmp);
    }
    else {
        ;
    }
}

template <>
void ZD::CHOT<float>::SH2HOT(const float * sh, float * hot, const unsigned int order)
{
    if (order == 4 || order == 6) {
        tijk_esh_to_3d_sym_f(hot, sh, order);
    }
    else {
        ;
        // not implemented
    }
}

template <>
void ZD::CHOT<double>::SH2HOT(const double * sh, double * hot, const unsigned int order)
{
    if (order == 4 || order == 6) {
        tijk_esh_to_3d_sym_d(hot, sh, order);
    }
    else {
        // not implemented
    }
}

template <>
void ZD::CHOT<float>::HOT2SH(const float * hot, float * sh, const unsigned int order)
{
    if (order == 4 || order == 6) {
        if (order == 4)
            tijk_3d_sym_to_esh_f(sh, hot, tijk_4o3d_sym);
        else
            tijk_3d_sym_to_esh_f(sh, hot, tijk_6o3d_sym);
    }
    else {
        // not implemented
    }
}

template <>
void ZD::CHOT<double>::HOT2SH(const double * hot, double * sh, const unsigned int order)
{
    if (order == 4 || order == 6) {
        if (order == 4)
            tijk_3d_sym_to_esh_d(sh, hot, tijk_4o3d_sym);
        else
            tijk_3d_sym_to_esh_d(sh, hot, tijk_6o3d_sym);
    }
    else {
        // not implemented
    }
}

template <>
const int ZD::CHOT<float>::HOT2Directions(const float * hot, CPoint<float, 3>* dirs, const unsigned int order)
{
    return 0;
}

template <>
const int ZD::CHOT<double>::HOT2Directions(const double * hot, CPoint<double, 3>* dirs, const unsigned int order)
{
    return 0;
}


template <>
const int ZD::CHOT<float>::ModelSelection(const float * hot, const unsigned int order)
{
    return 0;
}

template <>
const int ZD::CHOT<double>::ModelSelection(const double * hot, const unsigned int order)
{
    const tijk_type *type = nullptr;
    if (order == 4)
        type = tijk_4o3d_sym;
    else if (order == 6)
        type = tijk_4o3d_sym;
    else
        ;

    tijk_refine_rank1_parm *rank1Param = tijk_refine_rank1_parm_new();
    rank1Param->eps_impr = 0.001;
    tijk_refine_rankk_parm *rankkParam = tijk_refine_rankk_parm_new();
    rankkParam->pos = 1;
    rankkParam->eps_impr = 0.001;
    rankkParam->rank1_parm = rank1Param;
    tijk_approx_heur_parm *parm = tijk_approx_heur_parm_new();
    parm->eps_res = 0.02;
    parm->eps_impr = 0.01;
    parm->refine_parm = rankkParam;

    double res[28], ten[28];
    /* first get the positive approximation of the tensor */
    //tijk_approx_rankk_3d_d(NULL, NULL, res, hot, type, 6, rankkParam);
    //tijk_sub_d(ten, hot, res, type);

    double ls[4], vs[12];
    int kk = tijk_approx_heur_3d_d(ls, vs, res, hot, type, 4, parm);

    tijk_approx_heur_parm_nix(parm);

    return kk;
}

template<>
const float ZD::CHOT<float>::Norm(const float * hot, const unsigned int order)
{
    const tijk_type *type = nullptr;
    if (order == 4)
        type = tijk_4o3d_sym;
    else if (order == 6)
        type = tijk_6o3d_sym;
    else
        ;

    return type->norm_f(hot);
}

template<>
const double ZD::CHOT<double>::Norm(const double * hot, const unsigned int order)
{
    const tijk_type *type = nullptr;
    if (order == 4)
        type = tijk_4o3d_sym;
    else if (order == 6)
        type = tijk_6o3d_sym;
    else
        ;

    return type->norm_d(hot);
}


template<>
float ZD::CHOT<float>::RankKDecomposition(const float * hot, const unsigned int order, const unsigned int k, CPoint<float, 3> *vectors, float *weights)
{
    const tijk_type *type = nullptr;
    if (order == 4)
        type = tijk_4o3d_sym;
    else if (order == 6)
        type = tijk_4o3d_sym;
    else
        ;

    tijk_refine_rank1_parm *rank1Param = tijk_refine_rank1_parm_new();
    rank1Param->eps_impr = 0.0001;
    tijk_refine_rankk_parm *rankkParam = tijk_refine_rankk_parm_new();
    rankkParam->pos = 1;
    rankkParam->eps_impr = 0.0001;
    rankkParam->rank1_parm = rank1Param;

    float res[28], ten[28];

    /* first get the positive approximation of the tensor */
    tijk_approx_rankk_3d_f(NULL, NULL, res, hot, type, 6, rankkParam);
    tijk_sub_f(ten, hot, res, type);

    /* run rank-k decomposition */
    float ls[3];
    float vs[9];
    tijk_approx_rankk_3d_f(ls, vs, res, ten, type, k, rankkParam);

    tijk_refine_rank1_parm_nix(rank1Param);
    tijk_refine_rankk_parm_nix(rankkParam);

    for (int i = 0; i < k; ++i) {
        vectors[i] = CPoint<float, 3>(vs[i*3+0], vs[i*3+1], vs[i*3+2]);
        weights[i] = ls[i];
    }

    return type->norm_f(res);
}

template<>
double ZD::CHOT<double>::RankKDecomposition(const double * hot, const unsigned int order, const unsigned int k,  CPoint<double, 3> *vectors, double *weights)
{
    const tijk_type *type = nullptr;
    if (order == 4)
        type = tijk_4o3d_sym;
    else if (order == 6)
        type = tijk_6o3d_sym;
    else
        ;

    tijk_refine_rank1_parm *rank1Param = tijk_refine_rank1_parm_new();
    rank1Param->eps_impr = 0.0001;
    tijk_refine_rankk_parm *rankkParam = tijk_refine_rankk_parm_new();
    rankkParam->pos = 1;
    rankkParam->eps_impr = 0.0001;
    rankkParam->rank1_parm = rank1Param;

    double res[28], ten[28];

    /* first get the positive approximation of the tensor */
    tijk_approx_rankk_3d_d(NULL, NULL, res, hot, type, 6, rankkParam);
    tijk_sub_d(ten, hot, res, type);

    /* run rank-k decomposition */
    double ls[3];
    double vs[9];
    tijk_approx_rankk_3d_d(ls, vs, res, ten, type, k, rankkParam);

    //tijk_refine_rank1_parm_nix(rank1Param);
    tijk_refine_rankk_parm_nix(rankkParam);

    for (int i = 0; i < k; ++i) {
        vectors[i] = CPoint<double, 3>(vs[i * 3 + 0], vs[i * 3 + 1], vs[i * 3 + 2]);
        weights[i] = ls[i];
    }

    return type->norm_d(res);
}

template<>
float ZD::CHOT<float>::MaxAxis(const float *hot, const unsigned int order)
{
    CPoint<float, 3> vectors[1];
    float weights[1];
    RankKDecomposition(hot, order, 1, vectors, weights);
    
    float sh[28];
    HOT2SH(hot, sh, order);

    float phi = std::atan2(vectors[0][1], vectors[0][0]);
    float tmp = std::sqrt(vectors[0][0] * vectors[0][0] + vectors[0][1] * vectors[0][1]);
    float theta = std::atan2(tmp, vectors[0][2]);
    float length = tijk_eval_esh_f(sh, order, theta, phi);
    length = length < 0.0 ? 0.0 : length;
    return length;
}

template<>
double ZD::CHOT<double>::MaxAxis(const double *hot, const unsigned int order)
{
    CPoint<double, 3> vectors[1];
    double weights[1];
    RankKDecomposition(hot, order, 1, vectors, weights);

    double sh[28];
    HOT2SH(hot, sh, order);

    double phi = std::atan2(vectors[0][1], vectors[0][0]);
    double tmp = std::sqrt(vectors[0][0] * vectors[0][0] + vectors[0][1] * vectors[0][1]);
    double theta = std::atan2(tmp, vectors[0][2]);
    double length = tijk_eval_esh_d(sh, order, theta, phi);
    length = length < 0.0 ? 0.0 : length;
    return length;
}


#endif
