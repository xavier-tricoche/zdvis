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
#include "utils.hpp"

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

#include <vector>

char *pScalarPathname = nullptr;
char *pStrePathname = nullptr;
char *pGradPathname = nullptr;
char *pEvecPathname = nullptr;
char *pResultPathname = nullptr;

DATA_TYPE threshold = 0.0;
DATA_TYPE gaussianKernel = 0.0;


void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpGK, tmpT;

    hestOptAdd(&hopt, "i",    "filename",    airTypeString, 1, 1, &pScalarPathname, NULL,   "input 2D scalar field file.");
    hestOptAdd(&hopt, "o",    "filename",    airTypeString, 1, 1, &pResultPathname, NULL,   "output ridge lines file.");
    hestOptAdd(&hopt, "stre", "filename",    airTypeString, 1, 1, &pStrePathname,   "NULL", "output ridge strength file.");
    hestOptAdd(&hopt, "grad", "filename",    airTypeString, 1, 1, &pGradPathname,   "NULL", "output gradient file.");
    hestOptAdd(&hopt, "evec", "filename",    airTypeString, 1, 1, &pEvecPathname,   "NULL", "output eigenvector file.");
    hestOptAdd(&hopt, "g",    "kernel size", airTypeDouble, 1, 1, &tmpGK,           "0.0",  "Gaussian filter kernel size.");
    hestOptAdd(&hopt, "t",    "threshold",   airTypeDouble, 1, 1, &tmpT,            "0.0",  "ridge strength threshold.");
    
    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Extract ridge lines from a 2D scalar field",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    if (strcmp(pStrePathname, "NULL") == 0)
        pStrePathname = nullptr;
    if (strcmp(pGradPathname, "NULL") == 0)
        pGradPathname = nullptr;
    if (strcmp(pEvecPathname, "NULL") == 0)
        pEvecPathname = nullptr;

    threshold = tmpT;
    gaussianKernel = tmpGK;
}


int MarchingSquare(DATA_TYPE *value, ZD::CPoint<DATA_TYPE, 2> *pos, const int code, ZD::CPoint<DATA_TYPE, 2> *lines)
{
    ZD::CPoint<DATA_TYPE, 2> pt1, pt2, pt3, pt4;
    DATA_TYPE a, b;
    int count = 0;

    switch (code) {
    case 1:
    case 14:
        a = fabs(value[0]) / (fabs(value[1]) + fabs(value[0]));
        b = fabs(value[0]) / (fabs(value[3]) + fabs(value[0]));
        pt1 = pos[0] * (1.0f - a) + pos[1] * a;
        pt2 = pos[0] * (1.0f - b) + pos[3] * b;
        count = 1;
        break;
    case 2:
    case 13:
        a = fabs(value[1]) / (fabs(value[0]) + fabs(value[1]));
        b = fabs(value[1]) / (fabs(value[2]) + fabs(value[1]));
        pt1 = pos[1] * (1.0f - a) + pos[0] * a;
        pt2 = pos[1] * (1.0f - b) + pos[2] * b;
        count = 1;
        break;
    case 4:
    case 11:
        a = fabs(value[2]) / (fabs(value[1]) + fabs(value[2]));
        b = fabs(value[2]) / (fabs(value[3]) + fabs(value[2]));
        pt1 = pos[2] * (1.0f - a) + pos[1] * a;
        pt2 = pos[2] * (1.0f - b) + pos[3] * b;
        count = 1;
        break;
    case 8:
    case 7:
        a = fabs(value[3]) / (fabs(value[0]) + fabs(value[3]));
        b = fabs(value[3]) / (fabs(value[2]) + fabs(value[3]));
        pt1 = pos[3] * (1.0f - a) + pos[0] * a;
        pt2 = pos[3] * (1.0f - b) + pos[2] * b;
        count = 1;
        break;
    case 12:
    case 3:
        a = fabs(value[0]) / (fabs(value[0]) + fabs(value[3]));
        b = fabs(value[1]) / (fabs(value[1]) + fabs(value[2]));
        pt1 = pos[0] * (1.0f - a) + pos[3] * a;
        pt2 = pos[1] * (1.0f - b) + pos[2] * b;
        count = 1;
        break;
    case 9:
    case 6:
        a = fabs(value[0]) / (fabs(value[0]) + fabs(value[1]));
        b = fabs(value[3]) / (fabs(value[3]) + fabs(value[2]));
        pt1 = pos[0] * (1.0f - a) + pos[1] * a;
        pt2 = pos[3] * (1.0f - b) + pos[2] * b;
        count = 1;
        break;
    case 10:
    case 5:
    {
        DATA_TYPE v = 0.25 * (value[0] + value[1] + value[2] + value[3]);
        if (v * value[0] > 0.0f) {
            a = fabs(value[1]) / (fabs(value[0]) + fabs(value[1]));
            b = fabs(value[1]) / (fabs(value[2]) + fabs(value[1]));
            pt1 = pos[1] * (1.0f - a) + pos[0] * a;
            pt2 = pos[1] * (1.0f - b) + pos[2] * b;
            a = fabs(value[3]) / (fabs(value[0]) + fabs(value[3]));
            b = fabs(value[3]) / (fabs(value[2]) + fabs(value[3]));
            pt3 = pos[3] * (1.0f - a) + pos[0] * a;
            pt4 = pos[3] * (1.0f - b) + pos[2] * b;
        }
        else {
            a = fabs(value[0]) / (fabs(value[1]) + fabs(value[0]));
            b = fabs(value[0]) / (fabs(value[3]) + fabs(value[0]));
            pt1 = pos[0] * (1.0f - a) + pos[1] * a;
            pt2 = pos[0] * (1.0f - b) + pos[3] * b;
            a = fabs(value[2]) / (fabs(value[1]) + fabs(value[2]));
            b = fabs(value[2]) / (fabs(value[3]) + fabs(value[2]));
            pt3 = pos[2] * (1.0f - a) + pos[1] * a;
            pt4 = pos[2] * (1.0f - b) + pos[3] * b;
        }
        count = 2;
    }
    break;
    default:
        count = 0;
        break;
    }

    if (count == 1) {
        lines[0] = pt1;
        lines[1] = pt2;
    }
    else if (count == 2) {
        lines[0] = pt1;
        lines[1] = pt2;
        lines[2] = pt3;
        lines[3] = pt4;
    }

    return count;
}


int MarchingRidges(ZD::CPoint<DATA_TYPE, 2> *grad, ZD::CPoint<DATA_TYPE, 2> *evec, ZD::CPoint<DATA_TYPE, 2> *pos, ZD::CPoint<DATA_TYPE, 2>  *lines)
{
    // re-orientation eigen vector
    ZD::CPoint<DATA_TYPE, 4> mat;
    for (int i = 0; i < 4; ++i) {
        mat[0] += evec[i][0] * evec[i][0];
        mat[1] += evec[i][0] * evec[i][1];
        mat[2] += evec[i][1] * evec[i][0];
        mat[3] += evec[i][1] * evec[i][1];
    }

    ZD::CPoint<DATA_TYPE, 2> ec[2]; 
    ZD::CPoint<DATA_TYPE, 1> el[2];
    ZD::CEigenTool<DATA_TYPE>::Eigen2D(mat, el, ec);

    ZD::CPoint<DATA_TYPE, 2>  oevec[4];
    for (int i = 0; i < 4; ++i) {
        if (InnerProduct(evec[i], ec[0]) < 0.0f) {
            oevec[i][0] = -evec[i][0];
            oevec[i][1] = -evec[i][1];
        }
        else {
            oevec[i][0] = evec[i][0];
            oevec[i][1] = evec[i][1];
        }
    }

    // calculate zero-crossing
    DATA_TYPE value[4];
    int flag[4];
    int sum = 0;
    for (int i = 0; i < 4; ++i) {
        value[i] = InnerProduct(oevec[i], grad[i]);
        if (value[i] < 0.0)
            flag[i] = 0;
        else
            flag[i] = 1;

        sum += flag[i] * 1 << i;
    }

    // marching square
    int count = MarchingSquare(value, pos, sum, lines);
    return count;
}


int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    ZD::CField2<DATA_TYPE, 1> *pScalar = new ZD::CField2<DATA_TYPE, 1>();
    pScalar->OpenNrrdFile(pScalarPathname);
    if (gaussianKernel > 0.1)
        pScalar->GaussianFilter(gaussianKernel);

    int w = pScalar->GetSize(0);
    int h = pScalar->GetSize(1);

    ZD::CField2<DATA_TYPE, 4> *pHess = new ZD::CField2<DATA_TYPE, 4>();
    pHess->CreateField(w, h, nullptr);
    ZD::CField2<DATA_TYPE, 2> *pGrad = new ZD::CField2<DATA_TYPE, 2>();
    pGrad->CreateField(w, h, nullptr);
    ZD::CGradientTool<DATA_TYPE>::ComputeGradHess(pScalar, pGrad, pHess, w, h);

    ZD::CField2<DATA_TYPE, 2> *pEvec = new ZD::CField2<DATA_TYPE, 2>();
    pEvec->CreateField(w, h, nullptr);
    ZD::CField2<DATA_TYPE, 1> *pStre = new ZD::CField2<DATA_TYPE, 1>();
    pStre->CreateField(w, h, nullptr);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            ZD::CPoint<DATA_TYPE, 4> hess = pHess->GetValue(x, y);
            ZD::CPoint<DATA_TYPE, 1> eval[2];
            ZD::CPoint<DATA_TYPE, 2> evec[2];
            ZD::CEigenTool<DATA_TYPE>::Eigen2D(hess, eval, evec);
            pEvec->SetValue(x, y, evec[1]);
            if (eval[1][0] > 0.0)
                pStre->SetValue(x, y, ZD::CPoint<DATA_TYPE, 1>(0.0));
            else
                pStre->SetValue(x, y, ZD::CPoint<DATA_TYPE, 1>(std::abs(eval[1][0])));
        }
    }

    SafeDelete(pHess);

    ZD::CPoint<DATA_TYPE, 1> max_s, min_s;
    pStre->MinMax(min_s, max_s);
    DATA_TYPE s_thre = max_s[0] * threshold;

    // marching ridges
    std::vector<ZD::CPoint<DATA_TYPE, 2> > ridge_lines;
    for (int y = 6; y < h - 5; ++y) {
        for (int x = 6; x < w - 5; ++x) {
            ZD::CPoint<DATA_TYPE, 2>  pos[4];
            ZD::CPoint<DATA_TYPE, 2> grad[4];
            ZD::CPoint<DATA_TYPE, 2> evec[4];

            // position 
            pos[0][0] = x - 1; pos[0][1] = y - 1;
            pos[1][0] = x;     pos[1][1] = y - 1;
            pos[2][0] = x;     pos[2][1] = y;
            pos[3][0] = x - 1; pos[3][1] = y;
            // gradient
            grad[0] = pGrad->GetValue(x - 1, y - 1);
            grad[1] = pGrad->GetValue(x, y - 1);
            grad[2] = pGrad->GetValue(x, y);
            grad[3] = pGrad->GetValue(x - 1, y);
            // eigenvector
            evec[0] = pEvec->GetValue(x - 1, y - 1);
            evec[1] = pEvec->GetValue(x, y - 1);
            evec[2] = pEvec->GetValue(x, y);
            evec[3] = pEvec->GetValue(x - 1, y);

            ZD::CPoint<DATA_TYPE, 2> buffer[4];
            int count = MarchingRidges(grad, evec, pos, buffer);
            for (int i = 0; i < count; ++i) {
                if (pStre->GetValue(buffer[i*2+0])[0] > s_thre &&
                    pStre->GetValue(buffer[i*2+1])[0] > s_thre) {
                    ridge_lines.push_back(buffer[i*2+0]);
                    ridge_lines.push_back(buffer[i*2+1]);
                }
            }
        }
    }

    // save results
    size_t line_count = ridge_lines.size() / 2;
    ZD::CLine<DATA_TYPE, 2> *lines = new ZD::CLine<DATA_TYPE, 2>[line_count];
    for (size_t i = 0; i < line_count; ++i) {
        ZD::CPoint<DATA_TYPE, 2> pts[2];
        pts[0] = ridge_lines.at(i * 2 + 0);
        pts[1] = ridge_lines.at(i * 2 + 1);
        lines[i].CreateLine(pts, 2, pts[0]);
    }

    std::string ext = ZD::CFileTool::GetFileExtension(pResultPathname);
    if (strcmp(ext.c_str(), "vtk") == 0 || strcmp(ext.c_str(), "VTK") == 0) {
        SaveLinesVTK(pResultPathname, lines, line_count);
    }
    else {
        SaveLines(pResultPathname, lines, line_count);
    }

    if (pStrePathname != nullptr)
        pStre->SaveNrrdFile(pStrePathname);
    if (pEvecPathname != nullptr)
        pEvec->SaveNrrdFile(pEvecPathname);
    if (pGradPathname != nullptr)
        pGrad->SaveNrrdFile(pGradPathname);

    SafeDeleteArray(lines);
    SafeDelete(pScalar);
    SafeDelete(pGrad);
    SafeDelete(pEvec);
    SafeDelete(pStre);

    return EXIT_SUCCESS;
}
