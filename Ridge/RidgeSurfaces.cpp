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

#include "marching_cubes/MarchingCubes.h"

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

    hestOptAdd(&hopt, "i",    "filename",    airTypeString, 1, 1, &pScalarPathname, NULL,   "input 3D scalar field file.");
    hestOptAdd(&hopt, "o",    "filename",    airTypeString, 1, 1, &pResultPathname, NULL,   "output ridge surfaces file.");
    hestOptAdd(&hopt, "stre", "filename",    airTypeString, 1, 1, &pStrePathname,   "NULL", "output ridge strength file.");
    hestOptAdd(&hopt, "grad", "filename",    airTypeString, 1, 1, &pGradPathname,   "NULL", "output gradient file.");
    hestOptAdd(&hopt, "evec", "filename",    airTypeString, 1, 1, &pEvecPathname,   "NULL", "output eigenvector file.");
    hestOptAdd(&hopt, "g",    "kernel size", airTypeDouble, 1, 1, &tmpGK,           "0.0",  "Gaussian filter kernel size.");
    hestOptAdd(&hopt, "t",    "threshold",   airTypeDouble, 1, 1, &tmpT,            "0.0",  "ridge strength threshold.");
    
    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Extract ridge surfaces from a 3D scalar field",
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


int MarchingRidges(ZD::CPoint<DATA_TYPE, 3> *grad, ZD::CPoint<DATA_TYPE, 3> *evec,
    ZD::CPoint<DATA_TYPE, 3> *pos, ZD::CPoint<DATA_TYPE, 3> *triangles)
{
    // re-orientation eigen vector
    ZD::CPoint<DATA_TYPE, 9> mat;
    mat.SetZero();
    for (int i = 0; i < 8; ++i) {
        mat[0] += evec[i][0] * evec[i][0];
        mat[1] += evec[i][0] * evec[i][1];
        mat[2] += evec[i][0] * evec[i][2];
        mat[4] += evec[i][1] * evec[i][1];
        mat[5] += evec[i][1] * evec[i][2];
        mat[8] += evec[i][2] * evec[i][2];
    }
    mat[3] = mat[1];
    mat[6] = mat[2];
    mat[7] = mat[5];

    ZD::CPoint<DATA_TYPE, 1> tmp_eval[3];
    ZD::CPoint<DATA_TYPE, 3> tmp_evec[3];
    ZD::CEigenTool<DATA_TYPE>::Eigen3D(mat, tmp_eval, tmp_evec);
    ZD::CPoint<DATA_TYPE, 3> ec = tmp_evec[0];

    ZD::CPoint<DATA_TYPE, 3> oevec[8];
    for (int i = 0; i < 8; ++i) {
        if (InnerProduct(evec[i], ec) < 0.0f) {
            oevec[i] = -evec[i];
        }
        else {
            oevec[i] = evec[i];
        }
    }

    // calculate zero-crossing
    float sum = 0.0f;
    float value[8];
    for (int i = 0; i < 8; ++i) {
        //grad[i].Normalize();
        value[i] = InnerProduct(oevec[i], grad[i]);
        sum += value[i];
    }

    GRIDCELL grid;
    for (int i = 0; i < 8; ++i) {
        grid.p[i] = Vec3f(pos[i][0], pos[i][1], pos[i][2]);
        if (sum < 0.0f)
            grid.val[i] = -value[i];
        else
            grid.val[i] = value[i];
    }

    // marching cube
    TriMeshFace _FaceStorage[10];
    Vec3f _VertexStorage[15];
    int NewVertexCount = 0;
    int count = Polygonise(grid, _FaceStorage, NewVertexCount, _VertexStorage);
    for (int i = 0; i < count; ++i) {
        triangles[i * 3 + 0][0] = _VertexStorage[_FaceStorage[i].I[0]].x;
        triangles[i * 3 + 0][1] = _VertexStorage[_FaceStorage[i].I[0]].y;
        triangles[i * 3 + 0][2] = _VertexStorage[_FaceStorage[i].I[0]].z;

        triangles[i * 3 + 1][0] = _VertexStorage[_FaceStorage[i].I[1]].x;
        triangles[i * 3 + 1][1] = _VertexStorage[_FaceStorage[i].I[1]].y;
        triangles[i * 3 + 1][2] = _VertexStorage[_FaceStorage[i].I[1]].z;

        triangles[i * 3 + 2][0] = _VertexStorage[_FaceStorage[i].I[2]].x;
        triangles[i * 3 + 2][1] = _VertexStorage[_FaceStorage[i].I[2]].y;
        triangles[i * 3 + 2][2] = _VertexStorage[_FaceStorage[i].I[2]].z;
    }
    return count;
}

void RidgeSurfaces(ZD::CField3<DATA_TYPE, 3> *pGrad, ZD::CField3<DATA_TYPE, 3> *pEvec, 
    ZD::CField3<DATA_TYPE, 1> *pStre, std::vector<ZD::CPoint<DATA_TYPE, 3>>& triangles)
{
    const int w = pStre->GetSize(0);
    const int h = pStre->GetSize(1);
    const int d = pStre->GetSize(2);

    ZD::CPoint<DATA_TYPE, 1> max_s, min_s;
    pStre->MinMax(min_s, max_s);
    DATA_TYPE s_thre = max_s[0] * threshold;

    // marching ridges
    int skip = 0;

    for (int z = skip; z < d - skip; ++z) {
        for (int y = skip; y < h - skip; ++y) {
            for (int x = skip; x < w - skip; ++x) {
                ZD::CPoint<DATA_TYPE, 3>  pos[8];
                ZD::CPoint<DATA_TYPE, 3> grad[8];
                ZD::CPoint<DATA_TYPE, 3> evec[8];

                // position 
                pos[0] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x - 1), DATA_TYPE(y - 1), DATA_TYPE(z - 1));
                pos[1] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x),     DATA_TYPE(y - 1), DATA_TYPE(z - 1));
                pos[2] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x),     DATA_TYPE(y),     DATA_TYPE(z - 1));
                pos[3] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x - 1), DATA_TYPE(y),     DATA_TYPE(z - 1));
                pos[4] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x - 1), DATA_TYPE(y - 1), DATA_TYPE(z));
                pos[5] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x),     DATA_TYPE(y - 1), DATA_TYPE(z));
                pos[6] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x),     DATA_TYPE(y),     DATA_TYPE(z));
                pos[7] = ZD::CPoint<DATA_TYPE, 3>(DATA_TYPE(x - 1), DATA_TYPE(y),     DATA_TYPE(z));

                // gradient
                grad[0] = pGrad->GetValue(x - 1, y - 1, z - 1);
                grad[1] = pGrad->GetValue(x, y - 1, z - 1);
                grad[2] = pGrad->GetValue(x, y, z - 1);
                grad[3] = pGrad->GetValue(x - 1, y, z - 1);
                grad[4] = pGrad->GetValue(x - 1, y - 1, z);
                grad[5] = pGrad->GetValue(x, y - 1, z);
                grad[6] = pGrad->GetValue(x, y, z);
                grad[7] = pGrad->GetValue(x - 1, y, z);

                // eigenvector
                evec[0] = pEvec->GetValue(x - 1, y - 1, z - 1);
                evec[1] = pEvec->GetValue(x, y - 1, z - 1);
                evec[2] = pEvec->GetValue(x, y, z - 1);
                evec[3] = pEvec->GetValue(x - 1, y, z - 1);
                evec[4] = pEvec->GetValue(x - 1, y - 1, z);
                evec[5] = pEvec->GetValue(x, y - 1, z);
                evec[6] = pEvec->GetValue(x, y, z);
                evec[7] = pEvec->GetValue(x - 1, y, z);

                // marching ridges
                ZD::CPoint<DATA_TYPE, 3> buffer[100];
                int count = MarchingRidges(grad, evec, pos, buffer);
                for (int i = 0; i < count; ++i) {
                    if (pStre->GetValue(buffer[i * 3 + 0])[0] > s_thre &&
                        pStre->GetValue(buffer[i * 3 + 1])[0] > s_thre &&
                        pStre->GetValue(buffer[i * 3 + 2])[0] > s_thre) {
                        triangles.push_back(buffer[i * 3 + 0]);
                        triangles.push_back(buffer[i * 3 + 1]);
                        triangles.push_back(buffer[i * 3 + 2]);
                    }
                }
            }
        }
    }
}


int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    ZD::CField3<DATA_TYPE, 1> *pScalar = new ZD::CField3<DATA_TYPE, 1>();
    pScalar->OpenNrrdFile(pScalarPathname);
    if (gaussianKernel > 0.1)
        pScalar->GaussianFilter(gaussianKernel);

    const int w = pScalar->GetSize(0);
    const int h = pScalar->GetSize(1);
    const int d = pScalar->GetSize(2);

    ZD::CField3<DATA_TYPE, 9> *pHess = new ZD::CField3<DATA_TYPE, 9>();
    pHess->CreateField(w, h, d, nullptr);
    ZD::CField3<DATA_TYPE, 3> *pGrad = new ZD::CField3<DATA_TYPE, 3>();
    pGrad->CreateField(w, h, d, nullptr);
    ZD::CGradientTool<DATA_TYPE>::ComputeGradHess(pScalar, pGrad, pHess, w, h, d);

    ZD::CField3<DATA_TYPE, 3> *pEvec = new ZD::CField3<DATA_TYPE, 3>();
    pEvec->CreateField(w, h, d, nullptr);
    ZD::CField3<DATA_TYPE, 1> *pStre = new ZD::CField3<DATA_TYPE, 1>();
    pStre->CreateField(w, h, d, nullptr);
    for (int z = 0; z < d; ++z) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                ZD::CPoint<DATA_TYPE, 9> hess = pHess->GetValue(x, y, z);
                ZD::CPoint<DATA_TYPE, 1> eval[3];
                ZD::CPoint<DATA_TYPE, 3> evec[3];
                ZD::CEigenTool<DATA_TYPE>::Eigen3D(hess, eval, evec);
                pEvec->SetValue(x, y, z, evec[2]);
                if (eval[2][0] > 0.0)
                    pStre->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(0.0));
                else
                    pStre->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(std::abs(eval[2][0])));
            }
        }
    }

    SafeDelete(pHess);

    std::vector<ZD::CPoint<DATA_TYPE, 3> > triangles;
    RidgeSurfaces(pGrad, pEvec, pStre, triangles);

    if (pGradPathname != nullptr)
        pGrad->SaveNrrdFile(pGradPathname);
    if (pEvecPathname != nullptr)
        pEvec->SaveNrrdFile(pEvecPathname);
    if (pStrePathname != nullptr)
        pStre->SaveNrrdFile(pStrePathname);

    SafeDelete(pScalar);
    SafeDelete(pGrad);
    SafeDelete(pEvec);
    SafeDelete(pStre);

    ZD::CMesh<DATA_TYPE> *pMesh = new ZD::CMesh<DATA_TYPE>(triangles);
    pMesh->SaveMeshVTK(pResultPathname);
    SafeDelete(pMesh);

    return EXIT_SUCCESS;
}
