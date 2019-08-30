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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <teem/hest.h>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

enum ZD_DTI_CS_METHOD {
    ZD_DTI_CS_NONE,
    ZD_DTI_CS_FSR,
    ZD_DTI_CS_FIBER_FUNCTION
};

char *pFiberConfPathname = NULL;
char *pResultPathname = NULL;
bool STFlag = false;            // structure tensor
ZD_DTI_CS_METHOD method = ZD_DTI_CS_METHOD::ZD_DTI_CS_NONE;

int nx, ny, nz;
DATA_TYPE dx, dy, dz;
int length;
char **fiberFilePathnames;
int order;

void Initialize(const int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    char *tmpMethod;

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    hestOptAdd(&hopt, "i",  "filename",         airTypeString, 1, 1, &pFiberConfPathname, NULL, "input fiber configure file.");
    hestOptAdd(&hopt, "o",  "filename",         airTypeString, 1, 1, &pResultPathname,    NULL, "output coherent structure file.");
    hestOptAdd(&hopt, "m",  "method",           airTypeString, 1, 1, &tmpMethod,          NULL, "method: FSR and fiber_function.");
    hestOptAdd(&hopt, "k",  "moment order",     airTypeInt,    1, 1, &order,              "1",  "moment order used in fiber_function method, 1 or 2.");
    hestOptAdd(&hopt, "l",  "length",           airTypeInt,    1, 1, &length,             "-1", "fiber length used in computing coherent structures, -1 means using the whole fiber.");
    hestOptAdd(&hopt, "st", "structure tensor", airTypeBool,   0, 0, &STFlag,             NULL, "structure tensor flag.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute DTI coherent structure",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    if (strcmp(tmpMethod, "FSR") == 0 || strcmp(tmpMethod, "fsr") == 0)
        method = ZD_DTI_CS_METHOD::ZD_DTI_CS_FSR;
    else if (strcmp(tmpMethod, "fiber_function") == 0)
        method = ZD_DTI_CS_METHOD::ZD_DTI_CS_FIBER_FUNCTION;
    else
        method = ZD_DTI_CS_METHOD::ZD_DTI_CS_FIBER_FUNCTION;

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


void DTI_CS_FSR_ST(ZD::CField3<DATA_TYPE, 7> *pST)
{
    DATA_TYPE dis[3] = { dx, dy, dz };

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 9> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
    for (unsigned int i = 0; i < size; ++i) {
        fibers[i].ParameterizeFSR(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
    for (unsigned int i = 0; i < size; ++i) {
        fibers[i].ParameterizeFSR(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z + 1], fibers, size);
        for (unsigned int i = 0; i < size; ++i) {
            fibers[i].ParameterizeFSR(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 9> *buffer[7];
                buffer[0] = &(fiber_parameters[1][y*nx + (x - 1)]);
                buffer[1] = &(fiber_parameters[1][y*nx + (x + 1)]);
                buffer[2] = &(fiber_parameters[1][(y - 1)*nx + x]);
                buffer[3] = &(fiber_parameters[1][(y + 1)*nx + x]);
                buffer[4] = &(fiber_parameters[0][y*nx + x]);
                buffer[5] = &(fiber_parameters[2][y*nx + x]);
                buffer[6] = &(fiber_parameters[1][y*nx + x]);

                pST->SetValue(x, y, z, ZD::CCSTool<DATA_TYPE>::FSRST(buffer, dis));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}

void DTI_CS_FiberFunction_ST_1(ZD::CField3<DATA_TYPE, 7> *pST)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 3> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z + 1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fibers[i].ParameterizeFF(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 3> *buffer[27];
                int index = 0;
                for (int kk = 0; kk < 3; ++kk) {
                    for (int yy = y - 1; yy <= y + 1; ++yy) {
                        for (int xx = x - 1; xx <= x + 1; ++xx) {
                            buffer[index] = &(fiber_parameters[kk][yy*nx + xx]);
                            index++;
                        }
                    }
                }

                pST->SetValue(x, y, z, ZD::CCSTool<DATA_TYPE>::FiberFunctionST(buffer));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}

void DTI_CS_FiberFunction_ST_2(ZD::CField3<DATA_TYPE, 7> *pST)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 9> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z + 1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fibers[i].ParameterizeFF(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 9> *buffer[27];
                int index = 0;
                for (int kk = 0; kk < 3; ++kk) {
                    for (int yy = y - 1; yy <= y + 1; ++yy) {
                        for (int xx = x - 1; xx <= x + 1; ++xx) {
                            buffer[index] = &(fiber_parameters[kk][yy*nx + xx]);
                            index++;
                        }
                    }
                }

                pST->SetValue(x, y, z, ZD::CCSTool<DATA_TYPE>::FiberFunctionST(buffer));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}


void DTI_CS_FSR(ZD::CField3<DATA_TYPE, 1> *pCS)
{
    DATA_TYPE dis[3] = { dx, dy, dz };

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 9> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
    for (unsigned int i = 0; i < size; ++i) {
        fibers[i].ParameterizeFSR(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
    for (unsigned int i = 0; i < size; ++i) {
        fibers[i].ParameterizeFSR(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z + 1], fibers, size);
        for (unsigned int i = 0; i < size; ++i) {
            fibers[i].ParameterizeFSR(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 9> *buffer[7];
                buffer[0] = &(fiber_parameters[1][y*nx+(x-1)]);
                buffer[1] = &(fiber_parameters[1][y*nx+(x+1)]);
                buffer[2] = &(fiber_parameters[1][(y-1)*nx+x]);
                buffer[3] = &(fiber_parameters[1][(y+1)*nx+x]);
                buffer[4] = &(fiber_parameters[0][y*nx+x]);
                buffer[5] = &(fiber_parameters[2][y*nx+x]);
                buffer[6] = &(fiber_parameters[1][y*nx+x]);

                pCS->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(ZD::CCSTool<DATA_TYPE>::FSR(buffer, dis)));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}

void DTI_CS_FiberFunction_1(ZD::CField3<DATA_TYPE, 1> *pCS)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 3> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 3>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z+1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fibers[i].ParameterizeFF(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 3> *buffer[27];
                int index = 0;
                for (int kk = 0; kk < 3; ++kk) {
                    for (int yy = y - 1; yy <= y + 1; ++yy) {
                        for (int xx = x - 1; xx <= x + 1; ++xx) {
                            buffer[index] = &(fiber_parameters[kk][yy*nx + xx]);
                            index++;
                        }
                    }
                }

                pCS->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(ZD::CCSTool<DATA_TYPE>::FiberFunction(buffer)));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 3>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}

void DTI_CS_FiberFunction_2(ZD::CField3<DATA_TYPE, 1> *pCS)
{
    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    ZD::CFiber<DATA_TYPE> *fibers = nullptr;
    unsigned int size;
    ZD::CPoint<DATA_TYPE, 9> *fiber_parameters[3];
    fiber_parameters[0] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[1] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];
    fiber_parameters[2] = new ZD::CPoint<DATA_TYPE, 9>[nx*ny];

    OpenDTIFibers(fiberFilePathnames[0], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[0][i], length);
    }
    SafeDeleteArray(fibers);

    OpenDTIFibers(fiberFilePathnames[1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i = 0; i < (int)size; ++i) {
        fibers[i].ParameterizeFF(fiber_parameters[1][i], length);
    }
    SafeDeleteArray(fibers);

    for (int z = 1; z < nz - 1; ++z) {
        OpenDTIFibers(fiberFilePathnames[z + 1], fibers, size);
#if _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (int i = 0; i < (int)size; ++i) {
            fibers[i].ParameterizeFF(fiber_parameters[2][i], length);
        }
        SafeDeleteArray(fibers);

        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 9> *buffer[27];
                int index = 0;
                for (int kk = 0; kk < 3; ++kk) {
                    for (int yy = y - 1; yy <= y + 1; ++yy) {
                        for (int xx = x - 1; xx <= x + 1; ++xx) {
                            buffer[index] = &(fiber_parameters[kk][yy*nx + xx]);
                            index++;
                        }
                    }
                }

                pCS->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(ZD::CCSTool<DATA_TYPE>::FiberFunction(buffer)));
            }
        }

        memcpy(fiber_parameters[0], fiber_parameters[1], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memcpy(fiber_parameters[1], fiber_parameters[2], sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
        memset(fiber_parameters[2], 0, sizeof(ZD::CPoint<DATA_TYPE, 9>)*nx*ny);
    }

    SafeDeleteArray(fiber_parameters[0]);
    SafeDeleteArray(fiber_parameters[1]);
    SafeDeleteArray(fiber_parameters[2]);
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

#if _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    ReadFiberConfigureFile();

    if (STFlag) {
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
        if (method == ZD_DTI_CS_METHOD::ZD_DTI_CS_FSR) {
            DTI_CS_FSR_ST(pST);
        }
        else if (method == ZD_DTI_CS_METHOD::ZD_DTI_CS_FIBER_FUNCTION) {
            if (order == 1) {
                DTI_CS_FiberFunction_ST_1(pST);
            } else if (order == 2) {
                DTI_CS_FiberFunction_ST_2(pST);
            }
            else {
                ;
            }
        }
        else {
            ;
        }

        pST->SaveNrrdFile(pResultPathname);
        SafeDelete(pST);
    }
    else {
        ZD::CField3<DATA_TYPE, 1> *pCS = new ZD::CField3<DATA_TYPE, 1>();
        pCS->CreateField(nx, ny, nz, nullptr);

        if (method == ZD_DTI_CS_METHOD::ZD_DTI_CS_FSR) {
            DTI_CS_FSR(pCS);
        }
        else if (method == ZD_DTI_CS_METHOD::ZD_DTI_CS_FIBER_FUNCTION){
            if (order == 1) {
                DTI_CS_FiberFunction_1(pCS);
            }
            else if (order == 2) {
                DTI_CS_FiberFunction_2(pCS);
            }
            else {
                ;
            }
        }
        else {
            ;
        }

        pCS->SaveNrrdFile(pResultPathname);
        SafeDelete(pCS);
    }

    for (int z = 0; z < nz; ++z) {
        SafeDeleteArray(fiberFilePathnames[z]);
    }
    SafeDeleteArray(fiberFilePathnames);

    return EXIT_SUCCESS;
}
