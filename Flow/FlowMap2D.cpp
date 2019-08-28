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

#include <limits>

#include <teem/hest.h>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

typedef ZD::CIntegratorRK45<DATA_TYPE, 2> integ_type;

ZD::CFlow<DATA_TYPE, 2> *pFlow = nullptr;

char *pFlowName = nullptr;
char *pResultPathname = nullptr;
char *pDataPath = nullptr;

int nx, ny;

DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;

DATA_TYPE minx, maxx, miny, maxy;

DATA_TYPE deltaTime = 0.0;

void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpST, tmpET, tmpSS;
    double tmpDT;
    double tmpMinX, tmpMaxX, tmpMinY, tmpMaxY;

    hestOptAdd(&hopt, "o",  "pathname",           airTypeString, 1, 1, &pResultPathname,    NULL,   "output configure file pathname.");
    hestOptAdd(&hopt, "nx", "resolution x",       airTypeInt,    1, 1, &nx,                 NULL,   "sampling resolution along x.");
    hestOptAdd(&hopt, "ny", "resolution y",       airTypeInt,    1, 1, &ny,                 NULL,   "sampling resolution along y.");
    hestOptAdd(&hopt, "f",  "2D unsteady flow",   airTypeString, 1, 1, &pFlowName,          NULL,   "2D unsteady flow name, support: double_gyre, steady_double_gyre, meandering_jet, boussinesq, convection, and gaussian_vortices.");

    hestOptAdd(&hopt, "st", "start time", airTypeDouble, 1, 1, &tmpST, NULL, "integration starting time.");
    hestOptAdd(&hopt, "et", "end time",   airTypeDouble, 1, 1, &tmpET, NULL, "integration ending time.");
    hestOptAdd(&hopt, "ss", "step size",  airTypeDouble, 1, 1, &tmpSS, NULL, "integration step size.");

    hestOptAdd(&hopt, "dt", "output delta time", airTypeDouble, 1, 1, &tmpDT, NULL, "Output Flow Map delta time.");

    hestOptAdd(&hopt, "minx", "min x coord", airTypeDouble, 1, 1, &tmpMinX, "nan", "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord", airTypeDouble, 1, 1, &tmpMaxX, "nan", "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord", airTypeDouble, 1, 1, &tmpMinY, "nan", "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord", airTypeDouble, 1, 1, &tmpMaxY, "nan", "max y coord. of bounding box.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute Flow Map from 2D unsteady flow",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    startTime = tmpST;
    endTime   = tmpET;
    stepSize  = tmpSS;

    deltaTime = tmpDT;

    minx = tmpMinX;
    maxx = tmpMaxX;
    miny = tmpMinY;
    maxy = tmpMaxY;
}

void SaveConfigureFile()
{
    DATA_TYPE d[2];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);

    FILE *fp = fopen(pResultPathname, "w");
    if (fp == NULL) {
        printf("Cannot open file %s to write.\n", pResultPathname);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "Computing Flow Map from %s\n", pFlowName);
    fprintf(fp, "\tBounding box (%f, %f, %f, %f)\n", minx, maxx, miny, maxy);
    fprintf(fp, "\tSampling resolution %d x %d\n", nx, ny);
    fprintf(fp, "\tdx = %.10lf, dy = %.10lf\n", d[0], d[1]);
    fprintf(fp, "\tStarting time %f, ending time %f, time delta %f\n", startTime, endTime, deltaTime);
    fprintf(fp, "End\n");
    fclose(fp);
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    DATA_TYPE minTime = std::min(startTime, endTime);
    DATA_TYPE maxTime = std::max(startTime, endTime);
    if (strcmp(pFlowName, "double_gyre") == 0) {
        pFlow = new ZD::CFlowDoubleGyre<DATA_TYPE>();
    }
    else if (strcmp(pFlowName, "steady_double_gyre") == 0) {
        pFlow = new ZD::CFlowSteadyDoubleGyre<DATA_TYPE>();
    }
    else if (strcmp(pFlowName, "meandering_jet") == 0) {
        pFlow = new ZD::CFlowMeanderingJet<DATA_TYPE>();
    }
    else if (strcmp(pFlowName, "convection") == 0) {
        pFlow = new ZD::CFlowConvection<DATA_TYPE>(minTime, maxTime, true);
    }
    else if (strcmp(pFlowName, "boussinesq") == 0) {
        pFlow = new ZD::CFlowBoussinesq<DATA_TYPE>(minTime, maxTime);
    }
    else if (strcmp(pFlowName, "gaussian_vortices") == 0) {
        pFlow = new ZD::CFlowGaussianVortices<DATA_TYPE>(minTime, maxTime);
    }
    else {
        std::cerr << "Unknown 2D flow name: " << pFlowName << std::endl;
        return EXIT_FAILURE;
    }

    if (nx <= 3 || ny <= 3) {
        std::cerr << "Error! sampling resolution ( " << nx << " x " << ny << " ) is too small." << std::endl;
        SafeDelete(pFlow);
        return EXIT_FAILURE;
    }

    ZD::CPoint<DATA_TYPE, 2> bbox[2];
    pFlow->GetBBox(bbox[0], bbox[1]);
    minx = std::isnan(minx) ? bbox[0][0] : minx;
    miny = std::isnan(miny) ? bbox[0][1] : miny;
    maxx = std::isnan(maxx) ? bbox[1][0] : maxx;
    maxy = std::isnan(maxy) ? bbox[1][1] : maxy;


    SaveConfigureFile();

#ifdef _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    ZD::CPoint<DATA_TYPE, 2> *particles = new ZD::CPoint<DATA_TYPE, 2>[nx*ny];

    // initialize seeding position
    DATA_TYPE d[2];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            particles[y*nx+x] = ZD::CPoint<DATA_TYPE, 2>(minx + d[0] * x, miny + d[1] * y);
        }
    }

    // compute flow map
    pFlow->SetStepSize(stepSize);
    DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    DATA_TYPE remainTime = std::abs(endTime - startTime);
    DATA_TYPE integratedTime = 0.0;
    while (remainTime > 0.0) {
        DATA_TYPE currentStartTime = startTime + direction * integratedTime;

#ifdef _OPENMP
        omp_set_num_threads(OMP_MAX_THREADS);
        int thread_count[OMP_MAX_THREADS];
        memset(thread_count, 0, sizeof(int)*OMP_MAX_THREADS);
        int last_prec = 0;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int i = 0; i < nx*ny; ++i) {

            // std::cout << "integrating particle #" << i << '\n';

            pFlow->NextTime(particles[i], currentStartTime, deltaTime, direction,
                            new integ_type(stepSize, deltaTime/100.));

#ifdef _OPENMP
            int threadID = omp_get_thread_num();
            thread_count[threadID]++;
            if (threadID == 0) {
                int total = 0;
                for (int id = 0; id < OMP_MAX_THREADS; ++id)
                    total += thread_count[id];
                int prec = int(10000.0 * (double)total / double(nx*ny));
                if (prec > last_prec) {
                    printf("\r(%d / %d) pathlines finished.", total, nx*ny);
                    last_prec = prec;
                }
            }
#endif
        }

        integratedTime += deltaTime;
        remainTime -= deltaTime;

        // save flow map
        std::string temp = ZD::CFileTool::GetFilePath(pResultPathname);
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%sflowmap_dt=%09.04f.nrrd", temp.c_str(), integratedTime*direction);
        ZD::CField2<DATA_TYPE, 2> *pFlowMap = new ZD::CField2<DATA_TYPE, 2>();
        pFlowMap->CreateField(nx, ny, particles);
        pFlowMap->SaveNrrdFile(pathname);
        SafeDelete(pFlowMap);
    }


    SafeDelete(pFlow);
    SafeDeleteArray(particles);
    return EXIT_SUCCESS;
}
