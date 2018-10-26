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
//#include "ZD_DTI.hpp"
//#include "ZD_HOT.hpp"
#ifndef _ZD_LIB_PHANTOM_HPP_
#define _ZD_LIB_PHANTOM_HPP_

#include "utils/Base/ZD_Point.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "ZD_DWI.hpp"
#include "utils/define.hpp"
#include "utils/configure.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"

#include "ZD_DTI.hpp"
#include "ZD_SH.hpp"
#include "ZD_HOT.hpp"
#include "ZD_Glyph.hpp"

#include <teem/tijk.h>

#include "utils/UnitSphere/ZD_UnitSphereMedian.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"
#include "utils/configure.hpp"

#include <vector>



namespace ZD {
    template <typename T>
    class CPhantom  : public CDWI<T> {
    private:
        CField3<T, PHANTOM_DTI_N> *m_pDTI;            // dti volume
        CField3<T, PHANTOM_SH_N>  *m_pSH;            // sh volume
        CField3<T, PHANTOM_HOT_N> *m_pHOT;            // hot volume
        CField3<T, 1>             *m_pMask;            // mask volume
        CField3<T, 1>             *m_pSHMaxValue;    // sh maximum value volume

        /* ground truth */
        CPoint<T, 3> **m_fibers;        // ground truth fibers
        T *m_radius;                    // ground truth fiber radius
        CPoint<T, 3> **m_tangents;        // tangents

    public:
        CPhantom();
        CPhantom(const ZD_DWI_MODE mode);
        ~CPhantom();
    protected:
        virtual inline const bool Init(const ZD_DWI_MODE mode);

    public:
        virtual void LoadDTI();
        virtual void LoadHOT();
        virtual void LoadSH();
        void LoadMask();
        void LoadSHMaxValue();
        void LoadGroundTruth();

    public:
        virtual inline const bool HasDTI() const;
        virtual inline const bool HasHOT() const;
        virtual inline const bool HasSH() const;

    public:
        virtual const int ModelSelection(const CPoint<T, 3>& p) const;
        virtual const int ModelSelection(const T x, const T y, const T z) const;
        virtual const int ModelSelection(const CPoint<int, 3>& p) const;
        virtual const int ModelSelection(const int x, const int y, const int z) const;

        virtual inline void SetModel(const CPoint<int, 3>& p, const int model);
        virtual inline void SetModel(const int x, const int y, const int z, const int model);


    public:
        /* get local DTI value */
        virtual inline const bool GetDTI(const CPoint<int, 3>& p, T *dti) const;
        virtual inline const bool GetDTI(const int x, const int y, const int z, T *dti) const;
        virtual inline const bool GetDTI(const CPoint<T, 3>& p, T *dti) const;
        virtual inline const bool GetDTI(const T x, const T y, const T z, T *dti) const;
        virtual inline const bool GetDTIMLS(const CPoint<T, 3>& p, T *dti) const;
        virtual inline const bool GetDTIMLS(const T x, const T y, const T z, T *dti) const;

    public:
        /* get local HOT value */
        virtual inline const bool GetHOT(const CPoint<int, 3>& p, T *hot) const;
        virtual inline const bool GetHOT(const int x, const int y, const int z, T *hot) const;
        virtual inline const bool GetHOT(const CPoint<T, 3>& p, T *hot) const;
        virtual inline const bool GetHOT(const T x, const T y, const T z, T *hot) const;
        virtual inline const bool GetHOTMLS(const CPoint<T, 3>& p, T *hot) const;
        virtual inline const bool GetHOTMLS(const T x, const T y, const T z, T *hot) const;
        virtual inline const int GetHOTOrder() const;

    public:
        /* get local SH value */
        virtual inline const bool GetSH(const CPoint<int, 3>& p, T *sh) const;
        virtual inline const bool GetSH(const int x, const int y, const int z, T *sh) const;
        virtual inline const bool GetSH(const CPoint<T, 3>& p, T *sh) const;
        virtual inline const bool GetSH(const T x, const T y, const T z, T *sh) const;
        virtual inline const bool GetSHMLS(const CPoint<T, 3>& p, T *sh) const;
        virtual inline const bool GetSHMLS(const T x, const T y, const T z, T *sh) const;
        virtual inline const int GetSHOrder() const;
        virtual inline const T GetSHMaxValue(const CPoint<T, 3>& p) const;
        virtual inline const T GetSHMaxValue(const T x, const T y, const T z) const;

    public:
        virtual inline const bool IsWhiteMatter(const CPoint<T, 3>& p) const;
        virtual inline const bool IsWhiteMatter(const T x, const T y, const T z) const;

        //virtual const int GetOrder() const;
        
        //virtual void EnableMLSFilter();
        //virtual void DisableMLSFilter();

        /* world coordinate to volume grid */
        //virtual const bool W2G(const CPoint<T, 3>& wp, CPoint<T, 3>& gp) const;
        //virtual const bool W2G(const T wx, const T wy, const T wz, CPoint<T, 3>& gp) const;
        //virtual const bool W2G(const CPoint<T, 3>& wp, CPoint<int, 3>& gp) const;
        //virtual const bool W2G(const T wx, const T wy, const T wz, CPoint<int, 3>& gp) const;

        /* check if one position is inside a ground truth fiber */
        const bool IsInsideFiber(const CPoint<T, 3>& p, const int id) const;

        ///* mask */
        //virtual inline const bool LoadMask(const char *pathname);
        //virtual inline const bool SaveMask(const char *pathname) const;
        //virtual inline const int GetMask(const CPoint<T, 3>& p) const;
        //virtual inline const int GetMask(const T x, const T y, const T z) const;
        //virtual inline const int GetMask(const CPoint<int, 3>& p) const;
        //virtual inline const int GetMask(const int x, const int y, const int z) const;
        //virtual inline const int SetMask(const CPoint<int, 3>& p, const int mask);
        //virtual inline const int SetMask(const int x, const int y, const int z, const int mask);

    private:
        const bool IsInsideMask(const CPoint<T, 3>& p) const;

        /* fiber tracking */
    protected:
        virtual inline const bool SHNext(CRandomTool<T>& randomTool, CPoint<T, 3>& pos, 
            CPoint<T, 3>& dir, const T ss, const T alpha) const;
        virtual inline const bool SHDirection(CRandomTool<T>& randomTool, 
            const CPoint<T, 3>& pos, CPoint<T, 3>& dir) const;
        virtual inline const bool SHNextDirection(CRandomTool<T>& randomTool, const CPoint<T, 3>& pos, 
            const CPoint<T, 3>& in_dir, CPoint<T, 3>& out_dir, const T alpha) const;

        /* ground truth fibers */
    public:
        const bool GetFibers(std::vector<std::vector<CPoint<T, 3>>>& fibers) const;

    };
} // namespace ZD

template<typename T>
ZD::CPhantom<T>::CPhantom()
{
    CPhantom(ZD_DWI_NONE);
}

template<typename T>
ZD::CPhantom<T>::CPhantom(const ZD_DWI_MODE mode)
{
    m_fibers = nullptr;
    m_radius = nullptr;
    m_tangents = nullptr;

    m_pDTI = nullptr;
    m_pHOT = nullptr;
    m_pSH = nullptr;
    m_pMask = nullptr;
    m_pSHMaxValue = nullptr;

    if (Init(mode) == false) {
        printf("Error in func \"%s\" in file \"%s\", unknow mode.\n", __func__, __FILE__);
        exit(-1);
    }

    this->m_mode = mode;
}

template<typename T>
ZD::CPhantom<T>::~CPhantom()
{
    SafeDelete(m_pDTI);
    SafeDelete(m_pHOT);
    SafeDelete(m_pSH);
    SafeDelete(m_pMask);
    SafeDelete(m_pSHMaxValue);

    this->m_bbox[0] = this->m_bbox[1] = CPoint<T, 3>(0.0, 0.0, 0.0);

    if (m_fibers) {
        for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
            SafeDeleteArray(m_fibers[i]);
        }
        SafeDeleteArray(m_fibers);
    }

    SafeDeleteArray(m_radius);

    if (m_tangents) {
        for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
            SafeDeleteArray(m_tangents[i]);
        }
        SafeDeleteArray(m_tangents);
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::Init(const ZD_DWI_MODE mode)
{
    if (mode != ZD_DWI_NONE && mode != ZD_DWI_DTI && mode != ZD_DWI_HOT && mode != ZD_DWI_SH &&
        mode != ZD_DWI_DTI_MLS && mode != ZD_DWI_HOT_MLS && mode != ZD_DWI_SH_MLS && mode != ZD_DWI_ALL)
        return false;

    this->m_bbox[0] = this->m_bbox[1] = CPoint<T, 3>(0.0, 0.0, 0.0);
    this->m_size = CPoint<int, 3>(0, 0, 0);

    if (mode & 0x0001 || mode & 0x1000) {
        LoadDTI();
        this->m_size[0] = m_pDTI->GetSize(0);
        this->m_size[1] = m_pDTI->GetSize(1);
        this->m_size[2] = m_pDTI->GetSize(2);
    }
    if (mode & 0x0010) {
        LoadHOT();
        this->m_size[0] = m_pHOT->GetSize(0);
        this->m_size[1] = m_pHOT->GetSize(1);
        this->m_size[2] = m_pHOT->GetSize(2);
    }
    if (mode & 0x0100) {
        LoadSH();
        LoadSHMaxValue();
        this->m_size[0] = m_pSH->GetSize(0);
        this->m_size[1] = m_pSH->GetSize(1);
        this->m_size[2] = m_pSH->GetSize(2);
    }
    
    LoadMask();
    
    LoadGroundTruth();

    if (mode & 0x1000) {
        this->EnableMLSFilter();
    }

    return true;
}

template<typename T>
void ZD::CPhantom<T>::LoadDTI()
{
    SafeDelete(m_pDTI);
    m_pDTI = new CField3<T, PHANTOM_DTI_N>();
    m_pDTI->OpenNrrdFile(PHANTOM_DTI_PATHNAME);
    //m_bbox[0] = ZD::CPoint<T, 3>(-24.5, -24.5, -24.5);
    //m_bbox[1] = ZD::CPoint<T, 3>(24.5, 24.5, 24.5);
    this->m_bbox[0] = ZD::CPoint<T, 3>(-24.75, -24.75, -24.75);
    this->m_bbox[1] = ZD::CPoint<T, 3>(24.75, 24.75, 24.75);
    this->m_psize = ZD::CPoint<T, 3>(50.0, 50.0, 50.0);
}

template<typename T>
void ZD::CPhantom<T>::LoadHOT()
{
    SafeDelete(m_pHOT);
    m_pHOT = new CField3<T, PHANTOM_HOT_N>();
    m_pHOT->OpenNrrdFile(PHANTOM_HOT_PATHNAME);
    //m_bbox[0] = ZD::CPoint<T, 3>(-24.5, -24.5, -24.5);
    //m_bbox[1] = ZD::CPoint<T, 3>(24.5, 24.5, 24.5);
    this->m_bbox[0] = ZD::CPoint<T, 3>(-24.75, -24.75, -24.75);
    this->m_bbox[1] = ZD::CPoint<T, 3>(24.75, 24.75, 24.75);
    this->m_psize = ZD::CPoint<T, 3>(50.0, 50.0, 50.0);
}

template<typename T>
void ZD::CPhantom<T>::LoadSH()
{
    SafeDelete(m_pSH);
    m_pSH = new CField3<T, PHANTOM_SH_N>();
    m_pSH->OpenNrrdFile(PHANTOM_SH_PATHNAME);
    //m_bbox[0] = ZD::CPoint<T, 3>(-24.5, -24.5, -24.5);
    //m_bbox[1] = ZD::CPoint<T, 3>(24.5, 24.5, 24.5);
    this->m_bbox[0] = ZD::CPoint<T, 3>(-24.75, -24.75, -24.75);
    this->m_bbox[1] = ZD::CPoint<T, 3>(24.75, 24.75, 24.75);
    this->m_psize = ZD::CPoint<T, 3>(50.0, 50.0, 50.0);
}

template<typename T>
void ZD::CPhantom<T>::LoadMask()
{
    SafeDelete(m_pMask);
    m_pMask = new CField3<T, 1>();
    m_pMask->OpenNrrdFile(PHANTOM_MASK_PATHNAME);
}

template<typename T>
void ZD::CPhantom<T>::LoadSHMaxValue()
{
    SafeDelete(m_pSHMaxValue);
    m_pSHMaxValue = new CField3<T, 1>();
    m_pSHMaxValue->OpenNrrdFile(PHANTOM_SH_MAX_VALUE_PATHNAME);
}

template<typename T>
void ZD::CPhantom<T>::LoadGroundTruth()
{
    m_fibers = new CPoint<T, 3>*[PHANTOM_GROUND_TRUTH_COUNT];

    for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
        m_fibers[i] = new CPoint<T, 3>[PHANTOM_GROUND_TRUTH_VERTEX_COUNT];
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%sfiber-%02d.txt", PHANTOM_GROUND_TRUTH_PATH, i + 1);
        FILE *fp = fopen(pathname, "r");
        if (fp) {
            for (int k = 0; k < PHANTOM_GROUND_TRUTH_VERTEX_COUNT; ++k) {
                double x, y, z;
                fscanf(fp, "%lf %lf %lf", &x, &y, &z);
                m_fibers[i][k] = CPoint<T, 3>(x, y, z);
            }
            fclose(fp);
        }
    }

    m_radius = new T[PHANTOM_GROUND_TRUTH_COUNT];
    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sfibers-radii.txt", PHANTOM_GROUND_TRUTH_PATH);
    FILE *fp = fopen(pathname, "r");
    if (fp) {
        for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
            double r;
            int tmp;
            fscanf(fp, "%d %lf", &tmp, &r);
            m_radius[i] = r;
        }
        fclose(fp);
    }

    m_tangents = new CPoint<T, 3>*[PHANTOM_GROUND_TRUTH_COUNT];
    for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
        m_tangents[i] = new CPoint<T, 3>[PHANTOM_GROUND_TRUTH_VERTEX_COUNT];
        for (int k = 1; k < PHANTOM_GROUND_TRUTH_VERTEX_COUNT-1; ++k) {
            m_tangents[i][k] = m_fibers[i][k + 1] - m_fibers[i][k - 1];
            m_tangents[i][k].Normalize();
        }
        m_tangents[i][0] = m_fibers[i][1] - m_fibers[i][0];
        m_tangents[i][0].Normalize();
        m_tangents[i][PHANTOM_GROUND_TRUTH_VERTEX_COUNT-1] = m_fibers[i][PHANTOM_GROUND_TRUTH_VERTEX_COUNT-1] - m_fibers[i][PHANTOM_GROUND_TRUTH_VERTEX_COUNT-2];
        m_tangents[i][PHANTOM_GROUND_TRUTH_VERTEX_COUNT-1].Normalize();
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::HasDTI() const
{
    return this->m_pDTI != nullptr;
}

template<typename T>
inline const bool ZD::CPhantom<T>::HasHOT() const
{
    return this->m_pHOT != nullptr;
}


template<typename T>
inline const bool ZD::CPhantom<T>::HasSH() const
{
    return this->m_pSH != nullptr;
}



/* DTI section */

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTI(const CPoint<int, 3>& p, T *dti) const
{
    CPoint<T, PHANTOM_DTI_N> tmp = m_pDTI->GetValue(p);
    memcpy(dti, tmp.m_data, sizeof(T) * PHANTOM_DTI_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTI(const int x, const int y, const int z, T *dti) const
{
    CPoint<T, PHANTOM_DTI_N> tmp = m_pDTI->GetValue(x, y, z);
    memcpy(dti, tmp.m_data, sizeof(T) * PHANTOM_DTI_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTI(const CPoint<T, 3>& p, T *dti) const
{
    CPoint<T, 3> gp;
    if (this->W2G(p, gp)) {
        /* inside the domain */
        CPoint<T, PHANTOM_DTI_N> tmp = m_pDTI->GetValue(gp);
        memcpy(dti, tmp.m_data, sizeof(T) * PHANTOM_DTI_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTI(const T x, const T y, const T z, T *dti) const
{
    CPoint<T, 3> gp;
    if (this->W2G(x, y, z, gp)) {
        /* inside the domain */
        CPoint<T, PHANTOM_DTI_N> tmp = m_pDTI->GetValue(gp);
        memcpy(dti, tmp.m_data, sizeof(T) * PHANTOM_DTI_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTIMLS(const CPoint<T, 3>& p, T *dti) const
{
    return false;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetDTIMLS(const T x, const T y, const T z, T *dti) const
{
    return false;
}

/* end of DTI section */


/* HOT section */

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOT(const CPoint<int, 3>& p, T *hot) const
{
    CPoint<T, PHANTOM_HOT_N> tmp = m_pHOT->GetValue(p);
    memcpy(hot, tmp.m_data, sizeof(T) * PHANTOM_HOT_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOT(const int x, const int y, const int z, T *hot) const
{
    CPoint<T, PHANTOM_HOT_N> tmp = m_pHOT->GetValue(x, y, z);
    memcpy(hot, tmp.m_data, sizeof(T) * PHANTOM_HOT_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOT(const CPoint<T, 3>& p, T *hot) const
{
    CPoint<T, 3> gp;
    if (this->W2G(p, gp)) {
        CPoint<T, PHANTOM_HOT_N> tmp = m_pHOT->GetValue(gp);
        memcpy(hot, tmp.m_data, sizeof(T) * PHANTOM_HOT_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOT(const T x, const T y, const T z, T *hot) const
{
    CPoint<T, 3> gp;
    if (this->W2G(x, y, z, gp)) {
        CPoint<T, PHANTOM_HOT_N> tmp = m_pHOT->GetValue(gp);
        memcpy(hot, tmp.m_data, sizeof(T) * PHANTOM_HOT_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOTMLS(const CPoint<T, 3>& p, T *hot) const
{
    return false;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetHOTMLS(const T x, const T y, const T z, T *hot) const
{
    return false;
}

template<typename T>
inline const int ZD::CPhantom<T>::GetHOTOrder() const
{
    return PHANTOM_HOT_ORDER;
}

/* end of HOT section */


/* SH section */

template<typename T>
inline const bool ZD::CPhantom<T>::GetSH(const CPoint<int, 3>& p, T *sh) const
{
    CPoint<T, PHANTOM_SH_N> tmp = m_pSH->GetValue(p);
    memcpy(sh, tmp.m_data, sizeof(T)*PHANTOM_SH_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetSH(const int x, const int y, const int z, T *sh) const
{
    CPoint<T, PHANTOM_SH_N> tmp = m_pSH->GetValue(x, y, z);
    memcpy(sh, tmp.m_data, sizeof(T)*PHANTOM_SH_N);
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetSH(const CPoint<T, 3>& p, T *sh) const
{
    CPoint<T, 3> gp;
    if (this->W2G(p, gp)) {
        /* inside the domain */
        CPoint<T, PHANTOM_SH_N> tmp = m_pSH->GetValue(gp);
        memcpy(sh, tmp.m_data, sizeof(T)*PHANTOM_SH_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetSH(const T x, const T y, const T z, T *sh) const
{
    CPoint<T, 3> gp;
    if (this->W2G(x, y, z, gp)) {
        /* inside the domain */
        CPoint<T, PHANTOM_SH_N> tmp = m_pSH->GetValue(gp);
        memcpy(sh, tmp.m_data, sizeof(T)*PHANTOM_SH_N);
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetSHMLS(const CPoint<T, 3>& p, T *sh) const
{
    return false;
}

template<typename T>
inline const bool ZD::CPhantom<T>::GetSHMLS(const T x, const T y, const T z, T *sh) const
{
    return false;
}


template<typename T>
inline const int ZD::CPhantom<T>::GetSHOrder() const
{
    return PHANTOM_SH_ORDER;
}


template<typename T>
inline const T ZD::CPhantom<T>::GetSHMaxValue(const CPoint<T, 3>& p) const
{
    CPoint<T, 3> gp;
    if (this->W2G(p, gp)) {
        /* inside the domain */
        return m_pSHMaxValue->GetValue(gp)[0];
    }
    else {
        return 0.0;
    }
}

template<typename T>
inline const T ZD::CPhantom<T>::GetSHMaxValue(const T x, const T y, const T z) const
{
    CPoint<T, 3> gp;
    if (this->W2G(x, y, z, gp)) {
        /* inside the domain */
        return m_pSHMaxValue->GetValue(gp)[0];
    }
    else {
        return 0.0;
    }
}

/* end of SH section */


template<typename T>
const int ZD::CPhantom<T>::ModelSelection(const CPoint<T, 3>& p) const
{
    const T threshold = 0.05;

    CPoint<T, PHANTOM_HOT_N> hot;
    GetHOT(p, hot.m_data);
    T residuals[4];
    residuals[0] = CHOT<T>::Norm(hot.m_data, this->GetHOTOrder());

    for (int k = 1; k <= 3; ++k) {
        ZD::CPoint<T, 3> directions[3];
        T weights[3];
        residuals[k] = CHOT<T>::RankKDecomposition(hot.m_data, this->GetHOTOrder(), k, directions, weights);
    }
    
    if (residuals[3] / residuals[0] > threshold)
        return 0;
    if (residuals[1] / residuals[0] < threshold)
        return 1;
    if (residuals[2] / residuals[0] < threshold)
        return 2;
    if (residuals[3] / residuals[0] < threshold)
        return 3;
    return 0;
}

template<typename T>
const int ZD::CPhantom<T>::ModelSelection(const T x, const T y, const T z) const
{
    const T threshold = 0.05;

    CPoint<T, 3> p = CPoint<T, 3>(x, y, z);
    CPoint<T, PHANTOM_HOT_N> hot;
    GetHOT(p, hot.m_data);
    T residuals[4];
    residuals[0] = CHOT<T>::Norm(hot.m_data, this->GetHOTOrder());

    for (int k = 1; k <= 3; ++k) {
        ZD::CPoint<T, 3> directions[3];
        T weights[3];
        residuals[k] = CHOT<T>::RankKDecomposition(hot.m_data, this->GetHOTOrder(), k, directions, weights);
    }

    if (residuals[3] / residuals[0] > threshold)
        return 0;
    if (residuals[1] / residuals[0] < threshold)
        return 1;
    if (residuals[2] / residuals[0] < threshold)
        return 2;
    if (residuals[3] / residuals[0] < threshold)
        return 3;
}

template<typename T>
const int ZD::CPhantom<T>::ModelSelection(const CPoint<int, 3>& p) const
{
    return 0;
}

template<typename T>
const int ZD::CPhantom<T>::ModelSelection(const int x, const int y, const int z) const
{
    return 0;
}

template<typename T>
inline const bool ZD::CPhantom<T>::IsWhiteMatter(const CPoint<T, 3>& p) const
{
    CPoint<T, 3> gp;
    if (this->W2G(p, gp) && m_pMask->GetValue(gp)[0] > 0.5)
        return true;
    else
        return false;
}

template<typename T>
inline const bool ZD::CPhantom<T>::IsWhiteMatter(const T x, const T y, const T z) const
{
    CPoint<T, 3> p = CPoint<T, 3>(x, y, z);
    CPoint<T, 3> gp;
    if (this->W2G(p, gp) && m_pMask->GetValue(gp)[0] > 0.5)
        return true;
    else
        return false;
}

template<typename T>
inline void ZD::CPhantom<T>::SetModel(const CPoint<int, 3>& p, const int model)
{
}

template<typename T>
inline void ZD::CPhantom<T>::SetModel(const int x, const int y, const int z, const int model)
{
}


template<typename T>
const bool ZD::CPhantom<T>::IsInsideFiber(const CPoint<T, 3>& p, const int fiberID) const
{
    T minDis = FLT_MAX;
    CPoint<T, 3> t;
    CPoint<T, 3> q;

    for (int i = 0; i < PHANTOM_GROUND_TRUTH_VERTEX_COUNT-1; ++i) {
        CPoint<T, 3> q0 = m_fibers[fiberID][i];
        CPoint<T, 3> q1 = m_fibers[fiberID][i + 1];

        CPoint<T, 3> v0 = p - q0;
        CPoint<T, 3> v1 = q1 - q0;

        T tmp = InnerProduct(v0, v1) / v1.Length();
        T d = std::sqrt(v0.Length2() - tmp * tmp);
        if (tmp >= 0.0 && tmp <= v1.Length() && d < minDis) {
            t = v1;
            t.Normalize();
            q = tmp / v1.Length() * q1 + (1.0 - tmp / v1.Length()) * q0;
            minDis = d;
        }
    }

    for (int i = 0; i < PHANTOM_GROUND_TRUTH_VERTEX_COUNT; ++i) {
        CPoint<T, 3> v = m_fibers[fiberID][i] - p;
        if (v.Length() < minDis) {
            t = m_tangents[fiberID][i];
            q = m_fibers[fiberID][i];
            minDis = v.Length();
        }
    }

    if (minDis < m_radius[fiberID]) {
        CPoint<T, 3> d = p - q;
        d.Normalize();
        T tmp = InnerProduct(d, t);
        if (std::fabs(tmp) < 0.1)
            return true;
        else
            return false;
    }
    else {
        return false;
    }
}

template<typename T>
const bool ZD::CPhantom<T>::IsInsideMask(const CPoint<T, 3>& p) const
{
    const T radius = 5.0;
    const CPoint<T, 3> center(-0.5, -7.5, -5.5);
    CPoint<T, 3> d = p - center;
    if (d.Length() < radius)
        return true;
    else
        return false;
}


template<typename T>
inline const bool ZD::CPhantom<T>::SHNext(CRandomTool<T>& randomTool, CPoint<T, 3>& pos, CPoint<T, 3>& dir, const T ss, const T alpha) const
{
    CPoint<T, 3> next_dir;
    if (SHNextDirection(randomTool, pos, dir, next_dir, alpha) == false)
        return false;

    dir = next_dir;
    pos = pos + ss * dir;
    return true;
}

template<typename T>
inline const bool ZD::CPhantom<T>::SHDirection(CRandomTool<T>& randomTool, const CPoint<T, 3>& pos, CPoint<T, 3>& dir) const
{
    if (IsWhiteMatter(pos) == false)
        return false;

    CPoint<T, PHANTOM_SH_N> sh;
    this->GetSH(pos, sh.m_data);
    
    //T maxValue = CSH<T>::MaxValue(sh.m_data, PHANTOM_SH_ORDER);
    const T maxValue = this->GetSHMaxValue(pos);

    const T threshold = 0.1;
    bool success = false;
    int tryCount = 0;
    const int maxTry = 1000;
    do {
        T theta = randomTool.GetRandomNumer() * 2.0 * ZD_PI;
        T phi = randomTool.GetRandomNumer() * ZD_PI;
        T value = CSH<T>::Evaluate(sh.m_data, PHANTOM_SH_ORDER, theta, phi);
        if (value > threshold && value / maxValue > randomTool.GetRandomNumer()) {
            dir[0] = std::sin(theta) * std::cos(phi);
            dir[1] = std::sin(theta) * std::sin(phi);
            dir[2] = std::cos(theta);
            //T tmp_value = CSH<T>::Evaluate(sh.m_data, PHANTOM_SH_ORDER, dir);
            //tmp_value = CSH<T>::MaxValue(sh.m_data, PHANTOM_SH_ORDER, &dir);
            success = true;
        }
        tryCount++;
    } while (success == false && tryCount <= maxTry);
    
    return success;
}

template<typename T>
inline const bool ZD::CPhantom<T>::SHNextDirection(CRandomTool<T>& randomTool, const CPoint<T, 3>& pos, 
    const CPoint<T, 3>& in_dir, CPoint<T, 3>& out_dir, const T alpha) const
{
    if (IsWhiteMatter(pos) == false)
        return false;

    CPoint<T, PHANTOM_SH_N> sh;
    this->GetSH(pos, sh.m_data);

    //T maxValue = CSH<T>::MaxValue(sh.m_data, PHANTOM_SH_ORDER);
    const T maxValue = this->GetSHMaxValue(pos);

    const T threshold = 0.1;
    bool success = false;
    int tryCount = 0;
    const int maxTry = 1000;
    CPoint<T, 3> axisZ(0.0, 0.0, 1.0);
    do {
        T theta = randomTool.GetRandomNumer() * alpha;
        T phi = randomTool.GetRandomNumer() * 2.0 * ZD_PI;
        CPoint<T, 3> v;
        v[0] = std::tan(theta) * std::cos(phi);
        v[1] = std::tan(theta) * std::sin(phi);
        v[2] = 1.0;
        v.Normalize();
        CPoint<T, 3> vv;
        vv = this->RotateVector(axisZ, in_dir, v);
        vv.Normalize();
        T value = CSH<T>::Evaluate(sh.m_data, PHANTOM_SH_ORDER, vv);
        if (value > threshold && value / maxValue > randomTool.GetRandomNumer()) {
            out_dir = vv;
            success = true;
        }
        tryCount++;
    } while (success == false && tryCount <= maxTry);

    return success;
}


template<typename T>
const bool ZD::CPhantom<T>::GetFibers(std::vector<std::vector<CPoint<T, 3>>>& fibers) const
{
    fibers.resize(PHANTOM_GROUND_TRUTH_COUNT);
    for (int i = 0; i < PHANTOM_GROUND_TRUTH_COUNT; ++i) {
        fibers[i].resize(PHANTOM_GROUND_TRUTH_VERTEX_COUNT);
        for (int j = 0; j < PHANTOM_GROUND_TRUTH_VERTEX_COUNT; ++j) {
            fibers[i][j] = m_fibers[i][j];
        }
    }
    return true;
}


#endif
