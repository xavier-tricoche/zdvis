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

#include <vector>

#include <teem/hest.h>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif


char *pResultPathname = nullptr;
char *pFlowName = nullptr;
char **ppScalarNames = nullptr;
int scalarsCount = 0;
int nx, ny;
DATA_TYPE minx, maxx, miny, maxy;
DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;


void Initialize(const int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpMinX, tmpMaxX, tmpMinY, tmpMaxY;
    double tmpST, tmpET, tmpSS;

    hestOptAdd(&hopt, "o",    "filename",     airTypeString, 1, 1,  &pResultPathname,   NULL,  "output Lagrangian/Eulerian file.");
    hestOptAdd(&hopt, "f",    "flow name",    airTypeString, 1, 1,  &pFlowName,         NULL,  "2D flow name, support: double_gyre, convection, and boussinesq.");
    hestOptAdd(&hopt, "s",    "scalar name",  airTypeString, 1, -1, &ppScalarNames,     NULL,  "scalar quantities names", &scalarsCount);
    hestOptAdd(&hopt, "nx",   "resolution x", airTypeInt,    1, 1,  &nx,                NULL,  "resolution x");
    hestOptAdd(&hopt, "ny",   "resolution y", airTypeInt,    1, 1,  &ny,                NULL,  "resolution y");
    hestOptAdd(&hopt, "minx", "min x coord",  airTypeDouble, 1, 1,  &tmpMinX,           "nan", "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord",  airTypeDouble, 1, 1,  &tmpMaxX,           "nan", "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord",  airTypeDouble, 1, 1,  &tmpMinY,           "nan", "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord",  airTypeDouble, 1, 1,  &tmpMaxY,           "nan", "max y coord. of bounding box.");
    hestOptAdd(&hopt, "st",   "start time",   airTypeDouble, 1, 1,  &tmpST,             NULL,  "pathline integration starting time.");
    hestOptAdd(&hopt, "et",   "end time",     airTypeDouble, 1, 1,  &tmpET,             NULL,  "pathline integration ending time.");
    hestOptAdd(&hopt, "ss",   "step size",    airTypeDouble, 1, 1,  &tmpSS,             NULL,  "pathline integration step size.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute 2D Lagrangian/Eulerian",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    minx = tmpMinX;
    maxx = tmpMaxX;
    miny = tmpMinY;
    maxy = tmpMaxY;
    startTime = tmpST;
    endTime = tmpET;
    stepSize = tmpSS;
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    const int K = 2;

    DATA_TYPE dx, dy;
    dx = (maxx - minx) / DATA_TYPE(nx - 1);
    dy = (maxy - miny) / DATA_TYPE(ny - 1);
    ZD::CLETool<DATA_TYPE, K>::Initialize(dx, dy);

#if _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    for (int i = 0; i < scalarsCount; ++i) {
        ZD::CFlow<DATA_TYPE, 2> *pFlow = nullptr;
        if (strcmp(pFlowName, "double_gyre") == 0) {
            pFlow = new ZD::CFlowDoubleGyre<DATA_TYPE>();
        }
        else if (strcmp(pFlowName, "convection") == 0) {
            DATA_TYPE minTime = std::min(startTime, endTime) - stepSize;
            DATA_TYPE maxTime = std::max(startTime, endTime) + stepSize;
            pFlow = new ZD::CFlowConvection<DATA_TYPE>(minTime, maxTime, false, ppScalarNames[i]);
        }
        else if (strcmp(pFlowName, "boussinesq") == 0) {
            DATA_TYPE minTime = std::min(startTime, endTime) - stepSize;
            DATA_TYPE maxTime = std::max(startTime, endTime) + stepSize;
            pFlow = new ZD::CFlowBoussinesq<DATA_TYPE>(minTime, maxTime, ppScalarNames[i]);
        }
        else {
            std::cerr << "Unknown 2D flow name: " << pFlowName << std::endl;
            return EXIT_FAILURE;
        }

        pFlow->SetStepSize(stepSize);

        // tracking pathlines & parameterization
        DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
        DATA_TYPE totalTime = std::abs(endTime - startTime);
        DATA_TYPE dt = (endTime - startTime) / DATA_TYPE(PATHLINE_MAX_COUNT);

        ZD::CField2<DATA_TYPE, K> *parameters = new ZD::CField2<DATA_TYPE, K>();
        int size[2] = { nx, ny };
        parameters->CreateField(size);


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                ZD::CPoint<DATA_TYPE, 2> seed;
                seed[0] = minx + (maxx - minx) * DATA_TYPE(x) / DATA_TYPE(nx - 1);
                seed[1] = miny + (maxy - miny) * DATA_TYPE(y) / DATA_TYPE(ny - 1);
                DATA_TYPE integratedTime = 0.0;

                std::vector<DATA_TYPE> buffer;
                buffer.push_back(pFlow->Scalar(seed, startTime + integratedTime * direction));

                for (int j = 0; j < PATHLINE_MAX_COUNT; ++j) {
                    DATA_TYPE st = startTime + dt * DATA_TYPE(j);
                    DATA_TYPE tt = pFlow->NextTime(seed, st, std::fabs(dt), direction);
                    integratedTime += tt;
                    buffer.push_back(pFlow->Scalar(seed, startTime + integratedTime * direction));
                }
                parameters->SetValue(x, y, ZD::CParameterizationTool<DATA_TYPE, K>::ParameterizeSingleScalar(buffer));
            }
        }

        // compute Lagrangian-Eulerian
        ZD::CField2<DATA_TYPE, 1> *le = new ZD::CField2<DATA_TYPE, 1>();
        le->CreateField(nx, ny, nullptr);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, K> buffer[9];
                for (int k = 0; k < 9; ++k) {
                    int xx = x + ZD::neighbors_D2N9[k][0];
                    int yy = y + ZD::neighbors_D2N9[k][1];
                    buffer[k] = parameters->GetValue(xx, yy);
                }

                le->SetValue(x, y, ZD::CPoint<DATA_TYPE, 1>(ZD::CLETool<DATA_TYPE, K>::LagrangianEulerian2D(buffer)));
            }
        }

        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%sle_%s.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), ppScalarNames[i]);
        le->SaveNrrdFile(pathname);

        sprintf(pathname, "%sle_%s_params.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), ppScalarNames[i]);
        parameters->SaveNrrdFile(pathname);

        SafeDelete(parameters);
        SafeDelete(le);
        SafeDelete(pFlow);
    }

    return EXIT_SUCCESS;
}
