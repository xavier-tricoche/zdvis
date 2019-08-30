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
#ifndef _ZD_LIB_DWI_HPP_
#define _ZD_LIB_DWI_HPP_

#include "utils/Base/ZD_Point.hpp"
#include "utils/Base/ZD_Field.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"
#include "ZD_Fiber.hpp"
#include "utils/define.hpp"
#include "utils/Tool/ZD_RandomTool.hpp"


namespace ZD {
    enum ZD_DWI_MODE {
        ZD_DWI_NONE    = 0x0000,        // default mode, nothing is loaded
        ZD_DWI_DTI     = 0x0001,        // dti mode, dti data is loaded
        ZD_DWI_HOT     = 0x0010,        // hot mode, hot data is loaded
        ZD_DWI_SH      = 0x0100,        // sh mode, sh data is loaded
        ZD_DWI_DTI_MLS = 0x1001,        // dti with MLS mode
        ZD_DWI_HOT_MLS = 0x1010,        // hot with MLS mode
        ZD_DWI_SH_MLS  = 0x1100,        // sh with MLS mode
        ZD_DWI_ALL     = 0x1111            // dti, hot, and sh data are loaded
    };

    template <typename T>
    class CDWI {
    protected:
        ZD_DWI_MODE m_mode;
    public:
        CDWI();
        ~CDWI();
    protected:
        virtual inline const bool Init(const ZD_DWI_MODE mode) = 0;

    public:
        virtual void LoadDTI() = 0;
        virtual void LoadHOT() = 0;
        virtual void LoadSH() = 0;

    public:
        virtual inline const bool HasDTI() const = 0;
        virtual inline const bool HasHOT() const = 0;
        virtual inline const bool HasSH() const = 0;

    public:
        /* get local DTI value */
        virtual inline const bool GetDTI(const CPoint<int, 3>& p, T *dti) const = 0;
        virtual inline const bool GetDTI(const int x, const int y, const int z, T *dti) const = 0;
        virtual inline const bool GetDTI(const CPoint<T, 3>& p, T *dti) const = 0;
        virtual inline const bool GetDTI(const T x, const T y, const T z, T *dti) const = 0;
        virtual inline const bool GetDTIMLS(const CPoint<T, 3>& p, T *dti) const = 0;
        virtual inline const bool GetDTIMLS(const T x, const T y, const T z, T *dti) const = 0;
    
    public:
        /* get local HOT value */
        virtual inline const bool GetHOT(const CPoint<int, 3>& p, T *hot) const = 0;
        virtual inline const bool GetHOT(const int x, const int y, const int z, T *hot) const = 0;
        virtual inline const bool GetHOT(const CPoint<T, 3>& p, T *hot) const = 0;
        virtual inline const bool GetHOT(const T x, const T y, const T z, T *hot) const = 0;
        virtual inline const bool GetHOTMLS(const CPoint<T, 3>& p, T *hot) const = 0;
        virtual inline const bool GetHOTMLS(const T x, const T y, const T z, T *hot) const = 0;
        virtual inline const int GetHOTOrder() const = 0;

    public:
        /* get local SH value */
        virtual inline const bool GetSH(const CPoint<int, 3>& p, T *sh) const = 0;
        virtual inline const bool GetSH(const int x, const int y, const int z, T *sh) const = 0;
        virtual inline const bool GetSH(const CPoint<T, 3>& p, T *sh) const = 0;
        virtual inline const bool GetSH(const T x, const T y, const T z, T *sh) const = 0;
        virtual inline const bool GetSHMLS(const CPoint<T, 3>& p, T *sh) const = 0;
        virtual inline const bool GetSHMLS(const T x, const T y, const T z, T *sh) const = 0;
        virtual inline const int GetSHOrder() const = 0;
        virtual inline const T GetSHMaxValue(const CPoint<T, 3>& p) const = 0;
        virtual inline const T GetSHMaxValue(const T x, const T y, const T z) const = 0;

        /* DWI model selection, 0 -- not white matter, 1 or grater -- white matter */
        virtual inline const int ModelSelection(const CPoint<T, 3>& p) const = 0;
        virtual inline const int ModelSelection(const T x, const T y, const T z) const = 0;
        virtual inline const int ModelSelection(const CPoint<int, 3>& p) const = 0;
        virtual inline const int ModelSelection(const int x, const int y, const int z) const = 0;

        /* check if inside white matter */
        virtual inline const bool IsWhiteMatter(const CPoint<T, 3>& p) const = 0;
        virtual inline const bool IsWhiteMatter(const T x, const T y, const T z) const = 0;

        virtual inline void SetModel(const CPoint<int, 3>& p, const int model) = 0;
        virtual inline void SetModel(const int x, const int y, const int z, const int model) = 0;
        
        /* define the size of the dataset */
    protected:
        CPoint<T, 3> m_bbox[2];        // bounding box
        CPoint<int, 3> m_size;        // grid size 
        CPoint<T, 3> m_psize;        // physical size
    public:
        const int GetGridSize(const int dim) const;
        const T GetPhysicalSize(const int dim) const;
        const void GetBBox(CPoint<T, 3>& min, CPoint<T, 3>& max) const;
        /* world coordinate to volume grid */
        const bool W2G(const CPoint<T, 3>& wp, CPoint<T, 3>& gp) const;
        const bool W2G(const T wx, const T wy, const T wz, CPoint<T, 3>& gp) const;
        const bool W2G(const CPoint<T, 3>& wp, CPoint<int, 3>& gp) const;
        const bool W2G(const T wx, const T wy, const T wz, CPoint<int, 3>& gp) const;

        /* MLS filter */
    protected:
        CPoint<T, 3> *m_pMLSPoints;
        T *m_pMLSWeights;
        int m_MLSPointCount;
        T m_MLSSize;
    public:
        virtual void EnableMLSFilter();
        virtual void DisableMLSFilter();

        /* fiber tracking (DTI) */
    public:
        const T FiberTrackingDTIStream(const CPoint<T, 3>& seed, CFiber<T>& fiber,
            const T tl, const T ss, const T at) const;
        const T FiberTrackingDTIStreamMLS(const CPoint<T, 3>& seed, CFiber<T>& fiber,
            const T length, const T ss, const T at) const;
    private:
        inline const bool DTINext(CPoint<T, 3>& pos, CPoint<T, 3>& dir, const T ss, const T at) const;
        inline const bool DTINextDirection(const CPoint<T, 3>& pos, const CPoint<T, 3>& in_dir,
                CPoint<T, 3>& out_dir, const T at) const;
        
        /* fiber tracking (SH) */
    public:
        const T FiberTrackingSHProb(const CPoint<T, 3>& seed, CFiber<T>& fiber,
            const T tl, const T ss, const T alpha) const;
        const T FiberTrackingSHProbMLS(const CPoint<T, 3>& seed, CFiber<T>& fiber,
            const T tl, const T ss, const T alpha) const;
    protected:
        virtual inline const bool SHNext(CRandomTool<T>& randomTool, CPoint<T, 3>& pos, 
            CPoint<T, 3>& dir, const T ss, const T alpha) const = 0;
        virtual inline const bool SHDirection(CRandomTool<T>& randomTool, 
            const CPoint<T, 3>& pos, CPoint<T, 3>& dir) const = 0;
        virtual inline const bool SHNextDirection(CRandomTool<T>& randomTool, const CPoint<T, 3>& pos, 
            const CPoint<T, 3>& in_dir,    CPoint<T, 3>& out_dir, const T alpha) const = 0;

        //inline CPoint<T, 3>& RandomDirection(CRandomTool<T>& randomTool) const;
        //inline CPoint<T, 3>& RandomDirection(CRandomTool<T>& randomTool, const CPoint<T, 3>& in, const T theta) const;
        inline CPoint<T, 3> RotateVector(const CPoint<T, 3>& fromV, 
            const CPoint<T, 3>& toV, const CPoint<T, 3>& v) const;

    };
}

template<typename T>
ZD::CDWI<T>::CDWI()
{
    m_mode = ZD_DWI_MODE::ZD_DWI_NONE;
    m_pMLSPoints = nullptr;
    m_pMLSWeights = nullptr;
}

template<typename T>
ZD::CDWI<T>::~CDWI()
{
    SafeDeleteArray(m_pMLSPoints);
    SafeDeleteArray(m_pMLSWeights);
}

template<typename T>
const int ZD::CDWI<T>::GetGridSize(const int dim) const
{
    if (dim >= 0 && dim < 3)
        return m_size[dim];
    else
        return 0;
}

template<typename T>
const T ZD::CDWI<T>::GetPhysicalSize(const int dim) const
{
    if (dim >= 0 && dim < 3)
        return m_psize[dim];
    else
        return 0;
}

template<typename T>
const void ZD::CDWI<T>::GetBBox(CPoint<T, 3>& min, CPoint<T, 3>& max) const
{
    min = m_bbox[0];
    max = m_bbox[1];
}


template<typename T>
const bool ZD::CDWI<T>::W2G(const CPoint<T, 3>& wp, CPoint<T, 3>& gp) const
{
    if (wp[0] < m_bbox[0][0] || wp[1] < m_bbox[0][1] || wp[2] < m_bbox[0][2] ||
        wp[0] > m_bbox[1][0] || wp[1] > m_bbox[1][1] || wp[2] > m_bbox[1][2]) {
        return false;
    }
    else {
        CPoint<T, 3> tmp = wp - m_bbox[0];
        tmp[0] /= m_psize[0];
        tmp[1] /= m_psize[1];
        tmp[2] /= m_psize[2];
        tmp[0] = tmp[0] < 0.0 ? 0.0 : tmp[0];
        tmp[1] = tmp[1] < 0.0 ? 0.0 : tmp[1];
        tmp[2] = tmp[2] < 0.0 ? 0.0 : tmp[2];
        tmp[0] = tmp[0] > 1.0 ? 1.0 : tmp[0];
        tmp[1] = tmp[1] > 1.0 ? 1.0 : tmp[1];
        tmp[2] = tmp[2] > 1.0 ? 1.0 : tmp[2];
        gp[0] = tmp[0] * T(m_size[0] - 1) + 0.5;
        gp[1] = tmp[1] * T(m_size[1] - 1) + 0.5;
        gp[2] = tmp[2] * T(m_size[2] - 1) + 0.5;
        return true;
    }
}

template<typename T>
const bool ZD::CDWI<T>::W2G(const T wx, const T wy, const T wz, CPoint<T, 3>& gp) const
{
    if (wx < m_bbox[0][0] || wy < m_bbox[0][1] || wz < m_bbox[0][2] ||
        wx > m_bbox[1][0] || wy > m_bbox[1][1] || wz > m_bbox[1][2]) {
        return false;
    }
    else {
        CPoint<T, 3> tmp;
        tmp[0] = (wx - m_bbox[0][0]) / m_psize[0];
        tmp[1] = (wy - m_bbox[0][1]) / m_psize[1];
        tmp[2] = (wz - m_bbox[0][2]) / m_psize[2];
        tmp[0] = tmp[0] < 0.0 ? 0.0 : tmp[0];
        tmp[1] = tmp[1] < 0.0 ? 0.0 : tmp[1];
        tmp[2] = tmp[2] < 0.0 ? 0.0 : tmp[2];
        tmp[0] = tmp[0] > 1.0 ? 1.0 : tmp[0];
        tmp[1] = tmp[1] > 1.0 ? 1.0 : tmp[1];
        tmp[2] = tmp[2] > 1.0 ? 1.0 : tmp[2];
        gp[0] = tmp[0] * T(m_size[0] - 1) + 0.5;
        gp[1] = tmp[1] * T(m_size[1] - 1) + 0.5;
        gp[2] = tmp[2] * T(m_size[2] - 1) + 0.5;
        return true;
    }
}

template<typename T>
const bool ZD::CDWI<T>::W2G(const CPoint<T, 3>& wp, CPoint<int, 3>& gp) const
{
    if (wp[0] < m_bbox[0][0] || wp[1] < m_bbox[0][1] || wp[2] < m_bbox[0][2] ||
        wp[0] > m_bbox[1][0] || wp[1] > m_bbox[1][1] || wp[2] > m_bbox[1][2]) {
        return false;
    }
    else {
        CPoint<T, 3> tmp = wp - m_bbox[0];
        tmp[0] /= m_psize[0];
        tmp[1] /= m_psize[1];
        tmp[2] /= m_psize[2];
        tmp[0] = tmp[0] < 0.0 ? 0.0 : tmp[0];
        tmp[1] = tmp[1] < 0.0 ? 0.0 : tmp[1];
        tmp[2] = tmp[2] < 0.0 ? 0.0 : tmp[2];
        tmp[0] = tmp[0] > 1.0 ? 1.0 : tmp[0];
        tmp[1] = tmp[1] > 1.0 ? 1.0 : tmp[1];
        tmp[2] = tmp[2] > 1.0 ? 1.0 : tmp[2];
        gp[0] = int(tmp[0] * T(m_size[0] - 1) + 0.5);
        gp[1] = int(tmp[1] * T(m_size[1] - 1) + 0.5);
        gp[2] = int(tmp[2] * T(m_size[2] - 1) + 0.5);
        return true;
    }
}

template<typename T>
const bool ZD::CDWI<T>::W2G(const T wx, const T wy, const T wz, CPoint<int, 3>& gp) const
{
    if (wx < m_bbox[0][0] || wy < m_bbox[0][1] || wz < m_bbox[0][2] ||
        wx > m_bbox[1][0] || wy > m_bbox[1][1] || wz > m_bbox[1][2]) {
        return false;
    }
    else {
        CPoint<T, 3> tmp;
        tmp[0] = (wx - m_bbox[0][0]) / m_psize[0];
        tmp[1] = (wy - m_bbox[0][1]) / m_psize[1];
        tmp[2] = (wz - m_bbox[0][2]) / m_psize[2];
        tmp[0] = tmp[0] < 0.0 ? 0.0 : tmp[0];
        tmp[1] = tmp[1] < 0.0 ? 0.0 : tmp[1];
        tmp[2] = tmp[2] < 0.0 ? 0.0 : tmp[2];
        tmp[0] = tmp[0] > 1.0 ? 1.0 : tmp[0];
        tmp[1] = tmp[1] > 1.0 ? 1.0 : tmp[1];
        tmp[2] = tmp[2] > 1.0 ? 1.0 : tmp[2];
        gp[0] = int(tmp[0] * T(m_size[0] - 1) + 0.5);
        gp[1] = int(tmp[1] * T(m_size[1] - 1) + 0.5);
        gp[2] = int(tmp[2] * T(m_size[2] - 1) + 0.5);
        return true;
    }
}

template<typename T>
void ZD::CDWI<T>::EnableMLSFilter()
{
    // need DTI volume to perform MLS filter
    if (this->HasDTI() == false) {
        LoadDTI();
    }

    m_MLSPointCount = 50;
    m_MLSSize = 1.0;

    m_pMLSPoints = new CPoint<T, 3>[m_MLSPointCount];
    m_pMLSWeights = new T[m_MLSPointCount];

    m_pMLSPoints[0] = CPoint<T, 3>(0.0, 0.0, 0.0);
    m_pMLSWeights[0] = 1.0;

    CRandomTool<T> rt;

    for (int i = 1; i < m_MLSPointCount; ++i) {
        bool flag = false;
        CPoint<T, 3> pt;
        while (!flag) {
            T a = rt.GetRandomNumer() * 2.0 * ZD_PI;
            T b = rt.GetRandomNumer() * ZD_PI;
            T r = rt.GetRandomNumer();

            pt[0] = r * sin(a) * sin(b);
            pt[1] = r * sin(a) * cos(b);
            pt[2] = r * cos(a);

            flag = true;
            for (int j = 0; j < i; ++j) {
                T dis = Distance(m_pMLSPoints[j], pt);
                if (dis < 0.3) {
                    flag = false;
                    break;
                }
            }
        }

        m_pMLSPoints[i] = pt;
        m_pMLSWeights[i] = InnerProduct(m_pMLSPoints[i], m_pMLSPoints[i]) * 4.0;
        m_pMLSWeights[i] = exp(m_pMLSWeights[i] / -2.5);
    }
}

template<typename T>
void ZD::CDWI<T>::DisableMLSFilter()
{
    SafeDeleteArray(m_pMLSPoints);
    SafeDeleteArray(m_pMLSWeights);
    m_MLSPointCount = 0;
    m_MLSSize = 0.0;
}



template<typename T>
const T ZD::CDWI<T>::FiberTrackingDTIStream(const CPoint<T, 3>& seed, CFiber<T>& fiber, 
    const T tl, const T ss, const T at) const
{
    CPoint<T, 7> dti;
    this->GetDTI(seed, dti.m_data);
    CPoint<T, 3> dir = CDTI<T>::DTI2Direction(dti.m_data);

    CPoint<T, 3> p, last_dir;

    //const T deltaLength = 1.0;

    // forward
    std::vector<CPoint<T, 3>> forward;
    p = seed;
    last_dir = dir;
    forward.push_back(p);
    T forward_length = 0.0;
    for (; forward_length < tl * (1.0 - ZD_EPSILON); forward_length += ss) {
        if (DTINext(p, last_dir, ss, at) == false)
            break;
        forward.push_back(p);
    }

    // backward
    std::vector<CPoint<T, 3>> backward;
    p = seed;
    last_dir = -dir;
    backward.push_back(p);
    T backward_length = 0.0;
    for (;  backward_length < tl * (1.0 - ZD_EPSILON); backward_length += ss) {
        if (DTINext(p, last_dir, ss, at) == false)
            break;
        backward.push_back(p);
    }

    fiber.CreateFiber(forward, backward, seed, dir);

    return (forward_length + backward_length);
}


template<typename T>
const T ZD::CDWI<T>::FiberTrackingDTIStreamMLS(const CPoint<T, 3>& seed, CFiber<T>& fiber, 
    const T length, const T ss, const T at) const
{
    return 0.0;
}


template<typename T>
const T ZD::CDWI<T>::FiberTrackingSHProb(const CPoint<T, 3>& seed, CFiber<T>& fiber, 
    const T tl, const T ss, const T alpha) const
{
    CRandomTool<T> rt;
    CPoint<T, 3> dir;
    if (SHDirection(rt, seed, dir) == false)
        return 0.0;

    CPoint<T, 3> p, last_dir;

    // forward
    std::vector<CPoint<T, 3>> forward;
    p = seed;
    last_dir = dir;
    forward.push_back(p);
    T forward_length = 0.0;
    for (; forward_length < tl * (1.0 - ZD_EPSILON); forward_length += ss) {
        if (SHNext(rt, p, last_dir, ss, alpha) == false)
            break;
        forward.push_back(p);
    }

    // backward
    std::vector<CPoint<T, 3>> backward;
    p = seed;
    last_dir = -dir;
    backward.push_back(p);
    T backward_length = 0.0;
    for (; backward_length < tl * (1.0 - ZD_EPSILON); backward_length += ss) {
        if (SHNext(rt, p, last_dir, ss, alpha) == false)
            break;
        backward.push_back(p);
    }

    fiber.CreateFiber(forward, backward, seed, dir);

    return (forward_length + backward_length);
}

template<typename T>
const T ZD::CDWI<T>::FiberTrackingSHProbMLS(const CPoint<T, 3>& seed, CFiber<T>& fiber, 
    const T tl, const T ss, const T at) const
{
    return 0.0;
}


template<typename T>
inline const bool ZD::CDWI<T>::DTINext(CPoint<T, 3>& pos, CPoint<T, 3>& dir, const T ss, const T at) const
{
    CPoint<T, 3> k1, k2, temp;
    
    if (DTINextDirection(pos, dir, k1, at) == false)
        return false;

    temp = pos + k1 * (ss * 0.5);
    if (DTINextDirection(temp, dir, k2, at) == false)
        return false;

    dir = k2;                    // update last_dir
    pos = pos + ss * dir;        // update pos
    return true;
}

template<typename T>
inline const bool ZD::CDWI<T>::DTINextDirection(const CPoint<T, 3>& pos, const CPoint<T, 3>& in_dir, 
    CPoint<T, 3>& out_dir, const T at) const
{
    if (IsWhiteMatter(pos) == false)
        return false;

    CPoint<T, 7> dti;
    this->GetDTI(pos, dti.m_data);

    out_dir = CDTI<T>::DTI2Direction(dti.m_data);
    T ipv = InnerProduct(in_dir, out_dir);
    if (std::abs(ipv) < at)
        return false;
    if (ipv < 0.0)
        out_dir = -out_dir;
    return true;
}

template<typename T>
inline ZD::CPoint<T, 3> ZD::CDWI<T>::RotateVector(const CPoint<T, 3>& fromV, 
    const CPoint<T, 3>& toV, const CPoint<T, 3>& v) const
{
    T dotValue = InnerProduct(fromV, toV);
    CPoint<T, 3> vv;
    if (std::abs(dotValue) > 0.9999) {
        if (dotValue > 0.0)
            vv = v;
        else
            vv = -v;
    }
    else {
        CPoint<T, 3> u = CrossProduct(fromV, toV);        // rotation axis
        T C = dotValue;                                    // cos(theta)
        T S = u.Length();                                // sin(theta)
        u.Normalize();

        T m[3][3];
        m[0][0] = C + u[0] * u[0] * (1.0 - C);        m[1][0] = u[0] * u[1] * (1.0 - C) + u[2] * S; m[2][0] = u[0] * u[2] * (1.0 - C) - u[1] * S;
        m[0][1] = u[0] * u[1] * (1.0 - C) - u[2] * S; m[1][1] = C + u[1] * u[1] * (1.0 - C);        m[2][1] = u[1] * u[2] * (1.0 - C) + u[0] * S;
        m[0][2] = u[0] * u[2] * (1.0 - C) + u[1] * S; m[1][2] = u[1] * u[2] * (1.0 - C) - u[0] * S; m[2][2] = C + u[2] * u[2] * (1.0 - C);

        vv[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
        vv[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
        vv[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
    }

    return vv;
}


#endif
