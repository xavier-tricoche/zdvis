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
#include "configure.hpp"

#include <Eigen/Eigen>

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

char *pFiberConfPathname = nullptr;
char *pResultPathname = nullptr;

int nx, ny, nz;
DATA_TYPE dx, dy, dz;
char **fiberFilePathnames;
int order;

void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    hestOptAdd(&hopt, "i",  "input",            airTypeString, 1, 1, &pFiberConfPathname, NULL, "input HARDI fiber configuration file.");
    hestOptAdd(&hopt, "o",  "result",           airTypeString, 1, 1, &pResultPathname,    NULL, "output result coherent structure file (structure tensor).");
    hestOptAdd(&hopt, "k",  "moment order",     airTypeInt,    1, 1, &order,              "1",  "moment order used in fiber_function method, 1 or 2.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "HOT fiber tracking",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    if (order > 2) {
        order = 2;
    } else if (order < 1) {
        order = 1;
    } else {
        ;
    }
}

void ReadFiberConfigureFile()
{
    FILE *fp = fopen(pFiberConfPathname, "r");
    if (fp == NULL) {
        printf("Error, cannot open fiber configure file %s\n", pFiberConfPathname);
        exit(-1);
    }

    ZD::CFileTool::SkipLines(fp, 3);

    char temp[256];
    double tempDX, tempDY, tempDZ;
    fscanf(fp, "%s %s %s %s", temp, temp, temp, temp);
    fscanf(fp, "%d %s %s %s %d %s %s %s %d", &nx, temp, temp, temp, &ny, temp, temp, temp, &nz);
    fscanf(fp, "%s %s %s %lf %s %s %s %lf %s %s %s %lf", temp, temp, temp, &tempDX, temp, temp, temp, &tempDY, temp, temp, temp, &tempDZ);
    dx = tempDX;
    dy = tempDY;
    dz = tempDZ;

    fiberFilePathnames = new char*[nz];
    for (int z = 0; z < nz; ++z) {
        fiberFilePathnames[z] = new char[ZD_PATHNAME_LENGTH];
        fscanf(fp, "%s\n", fiberFilePathnames[z]);
    }

    fclose(fp);
}

unsigned int FindCloestFiber(const ZD::CPoint<DATA_TYPE, 3> &srcDir, const ZD::CPoint<DATA_TYPE, 3> *dirs, const unsigned int count)
{
    unsigned int id = 0;
    DATA_TYPE maxV = 0.0;
    for (unsigned int i = 0; i < count; ++i) {
        DATA_TYPE tmpV = std::abs(ZD::InnerProduct(srcDir, dirs[i]));
        if (tmpV > maxV) {
            maxV = tmpV;
            id = i;
        }
    }
    return id;
}

void HARDI_CS_FiberFunction_ST_1(ZD::CField3<DATA_TYPE, 7> *pST)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> **fibers = nullptr;
    unsigned int size;
    unsigned int *fiberCounts[3];
    ZD::CPoint<DATA_TYPE, 3> **fiber_directions[3];
    ZD::CPoint<DATA_TYPE, 3> **fiber_parameters[3];

    for (int z = 0; z < 2; ++z) {
        OpenHARDIFibers(fiberFilePathnames[z], fibers, fiberCounts[z], size);
        fiber_parameters[z] = new ZD::CPoint<DATA_TYPE, 3>*[size];
        fiber_directions[z] = new ZD::CPoint<DATA_TYPE, 3>*[size];
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fiber_parameters[z][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[z][i]];
            fiber_directions[z][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[z][i]];
            for (unsigned int j = 0; j < fiberCounts[z][i]; ++j) {
                fibers[i][j].ParameterizeFF(fiber_parameters[z][i][j], -1);
                fiber_directions[z][i][j] = fibers[i][j].m_seedDir;
            }
            SafeDeleteArray(fibers[i]);
        }
        SafeDeleteArray(fibers);
    }

    for (int z = 1; z < nz - 1; ++z) {
        OpenHARDIFibers(fiberFilePathnames[z + 1], fibers, fiberCounts[2], size);
        fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 3>*[size];
        fiber_directions[2] = new ZD::CPoint<DATA_TYPE, 3>*[size];
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fiber_parameters[2][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[2][i]];
            fiber_directions[2][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[2][i]];
            for (unsigned int j = 0; j < fiberCounts[2][i]; ++j) {
                fibers[i][j].ParameterizeFF(fiber_parameters[2][i][j], -1);
                fiber_directions[2][i][j] = fibers[i][j].m_seedDir;
            }
            SafeDeleteArray(fibers[i]);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 7> st;
                st.SetZero();
                st[0] = 1.0;
                for (unsigned int j = 0; j < fiberCounts[1][y*nx + x]; ++j) {
                    // for each fiber at (x, y, z), find the closest fibers from its neighbors, and compute the strcuture tensor from them
                    ZD::CPoint<DATA_TYPE, 3> srcDir = fiber_directions[1][y*nx + x][j];
                    ZD::CPoint<DATA_TYPE, 3> *buffer[27];
                    int index = 0;
                    for (int kk = 0; kk < 3; ++kk) {
                        for (int yy = y - 1; yy <= y + 1; ++yy) {
                            for (int xx = x - 1; xx <= x + 1; ++xx) {
                                unsigned int dstIndex = FindCloestFiber(srcDir, fiber_directions[kk][yy*nx + xx], fiberCounts[kk][yy*nx + xx]);
                                buffer[index] = &(fiber_parameters[kk][yy*nx + xx][dstIndex]);
                                index++;
                            }
                        }
                    }
                    ZD::CPoint<DATA_TYPE, 7> tmpST = ZD::CCSTool<DATA_TYPE>::FiberFunctionST(buffer);

                    // get the sum of all structure tensors at (x, y, z)
                    for (unsigned int j = 1; j < 7; ++j) {
                        st[j] += tmpST[j];
                    }
                }

                pST->SetValue(x, y, z, st);
            }
        }

        // free all memory allocated for the previous layer
        for (unsigned int i = 0; i < size; ++i) {
            SafeDeleteArray(fiber_parameters[0][i]);
            SafeDeleteArray(fiber_directions[0][i]);
        }
        SafeDeleteArray(fiber_parameters[0]);
        SafeDeleteArray(fiber_directions[0]);
        SafeDeleteArray(fiberCounts[0]);


        fiber_parameters[0] = fiber_parameters[1];
        fiber_directions[0] = fiber_directions[1];
        fiberCounts[0] = fiberCounts[1];

        fiber_parameters[1] = fiber_parameters[2];
        fiber_directions[1] = fiber_directions[2];
        fiberCounts[1] = fiberCounts[2];
    }

    for (unsigned int i = 0; i < size; ++i) {
        SafeDeleteArray(fiber_parameters[0][i]);
        SafeDeleteArray(fiber_parameters[1][i]);
        SafeDeleteArray(fiber_directions[0][i]);
        SafeDeleteArray(fiber_directions[1][i]);
    }
    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_directions[0]);
    SafeDeleteArray(fiber_directions[1]);
    SafeDeleteArray(fiberCounts[0]);
    SafeDeleteArray(fiberCounts[1]);
}

void HARDI_CS_FiberFunction_ST_2(ZD::CField3<DATA_TYPE, 7> *pST)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> **fibers = nullptr;
    unsigned int size;
    unsigned int *fiberCounts[3];
    ZD::CPoint<DATA_TYPE, 3> **fiber_directions[3];
    ZD::CPoint<DATA_TYPE, 9> **fiber_parameters[3];

    for (int z = 0; z < 2; ++z) {
        OpenHARDIFibers(fiberFilePathnames[z], fibers, fiberCounts[z], size);
        fiber_parameters[z] = new ZD::CPoint<DATA_TYPE, 9>*[size];
        fiber_directions[z] = new ZD::CPoint<DATA_TYPE, 3>*[size];
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fiber_parameters[z][i] = new ZD::CPoint<DATA_TYPE, 9>[fiberCounts[z][i]];
            fiber_directions[z][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[z][i]];
            for (unsigned int j = 0; j < fiberCounts[z][i]; ++j) {
                fibers[i][j].ParameterizeFF(fiber_parameters[z][i][j], -1);
                fiber_directions[z][i][j] = fibers[i][j].m_seedDir;
            }
            SafeDeleteArray(fibers[i]);
        }
        SafeDeleteArray(fibers);
    }

    for (int z = 1; z < nz - 1; ++z) {
        OpenHARDIFibers(fiberFilePathnames[z + 1], fibers, fiberCounts[2], size);
        fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 9>*[size];
        fiber_directions[2] = new ZD::CPoint<DATA_TYPE, 3>*[size];
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fiber_parameters[2][i] = new ZD::CPoint<DATA_TYPE, 9>[fiberCounts[2][i]];
            fiber_directions[2][i] = new ZD::CPoint<DATA_TYPE, 3>[fiberCounts[2][i]];
            for (unsigned int j = 0; j < fiberCounts[2][i]; ++j) {
                fibers[i][j].ParameterizeFF(fiber_parameters[2][i][j], -1);
                fiber_directions[2][i][j] = fibers[i][j].m_seedDir;
            }
            SafeDeleteArray(fibers[i]);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 7> st;
                st.SetZero();
                st[0] = 1.0;
                for (unsigned int j = 0; j < fiberCounts[1][y*nx + x]; ++j) {
                    // for each fiber at (x, y, z), find the closest fibers from its neighbors, and compute the strcuture tensor from them
                    ZD::CPoint<DATA_TYPE, 3> srcDir = fiber_directions[1][y*nx + x][j];
                    ZD::CPoint<DATA_TYPE, 9> *buffer[27];
                    int index = 0;
                    for (int kk = 0; kk < 3; ++kk) {
                        for (int yy = y - 1; yy <= y + 1; ++yy) {
                            for (int xx = x - 1; xx <= x + 1; ++xx) {
                                unsigned int dstIndex = FindCloestFiber(srcDir, fiber_directions[kk][yy*nx + xx], fiberCounts[kk][yy*nx + xx]);
                                buffer[index] = &(fiber_parameters[kk][yy*nx + xx][dstIndex]);
                                index++;
                            }
                        }
                    }
                    ZD::CPoint<DATA_TYPE, 7> tmpST = ZD::CCSTool<DATA_TYPE>::FiberFunctionST(buffer);

                    // get the sum of all structure tensors at (x, y, z)
                    for (unsigned int j = 1; j < 7; ++j) {
                        st[j] += tmpST[j];
                    }
                }

                pST->SetValue(x, y, z, st);
            }
        }

        // free all memory allocated for the previous layer
        for (unsigned int i = 0; i < size; ++i) {
            SafeDeleteArray(fiber_parameters[0][i]);
            SafeDeleteArray(fiber_directions[0][i]);
        }
        SafeDeleteArray(fiber_parameters[0]);
        SafeDeleteArray(fiber_directions[0]);
        SafeDeleteArray(fiberCounts[0]);


        fiber_parameters[0] = fiber_parameters[1];
        fiber_directions[0] = fiber_directions[1];
        fiberCounts[0] = fiberCounts[1];

        fiber_parameters[1] = fiber_parameters[2];
        fiber_directions[1] = fiber_directions[2];
        fiberCounts[1] = fiberCounts[2];
    }

    for (unsigned int i = 0; i < size; ++i) {
        SafeDeleteArray(fiber_parameters[0][i]);
        SafeDeleteArray(fiber_parameters[1][i]);
        SafeDeleteArray(fiber_directions[0][i]);
        SafeDeleteArray(fiber_directions[1][i]);
    }
    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_directions[0]);
    SafeDeleteArray(fiber_directions[1]);
    SafeDeleteArray(fiberCounts[0]);
    SafeDeleteArray(fiberCounts[1]);
}



int main(int argc, char *argv[])
{
    Initialize(argc, argv);

#if _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    ReadFiberConfigureFile();

    ZD::CField3<DATA_TYPE, 7> *pST = new ZD::CField3<DATA_TYPE, 7>();
    pST->CreateField(nx, ny, nz, nullptr);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x) {
                DATA_TYPE tmp[7] = {
                    (DATA_TYPE)1.0, (DATA_TYPE)1.0, (DATA_TYPE)0.0,
                    (DATA_TYPE)0.0, (DATA_TYPE)1.0, (DATA_TYPE)0.0,
                    (DATA_TYPE)1.0
                };
                pST->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 7>(tmp));
            }

    if (order == 1) {
        HARDI_CS_FiberFunction_ST_1(pST);
    }
    else if (order == 2) {
        HARDI_CS_FiberFunction_ST_2(pST);
    }
    else {
        ;
    }

    pST->SaveNrrdFile(pResultPathname);
    SafeDelete(pST);

    for (int z = 0; z < nz; ++z) {
        SafeDeleteArray(fiberFilePathnames[z]);
    }
    SafeDeleteArray(fiberFilePathnames);

    return EXIT_SUCCESS;
}
