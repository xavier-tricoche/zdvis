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


typedef ZD::CIntegratorRK45<DATA_TYPE, 3> integ_type;

ZD::CFlow<DATA_TYPE, 3> *pFlow = nullptr;

char *pFlowName = nullptr;
char *pResultPathname = nullptr;
char *pDataPath = nullptr;
char *pCtName = nullptr;

int nx, ny, nz;

DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;

DATA_TYPE minx, maxx, miny, maxy, minz, maxz;

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
    double tmpMinX, tmpMaxX, tmpMinY, tmpMaxY, tmpMinZ, tmpMaxZ;

    hestOptAdd(&hopt, "o",  "pathname",           airTypeString, 1, 1, &pResultPathname,    NULL,   "output configure file pathname.");
    hestOptAdd(&hopt, "nx", "resolution x",       airTypeInt,    1, 1, &nx,                 NULL,   "sampling resolution along x.");
    hestOptAdd(&hopt, "ny", "resolution y",       airTypeInt,    1, 1, &ny,                 NULL,   "sampling resolution along y.");
    hestOptAdd(&hopt, "nz", "resolution z",       airTypeInt,    1, 1, &nz,                 NULL,   "sampling resolution along z.");
    hestOptAdd(&hopt, "f",  "3D unsteady flow",   airTypeString, 1, 1, &pFlowName,          NULL,   "3D unsteady flow name, support: abc, delta_wing.");

    hestOptAdd(&hopt, "st", "start time", airTypeDouble, 1, 1, &tmpST, NULL, "integration starting time.");
    hestOptAdd(&hopt, "et", "end time",   airTypeDouble, 1, 1, &tmpET, NULL, "integration ending time.");
    hestOptAdd(&hopt, "ss", "step size",  airTypeDouble, 1, 1, &tmpSS, NULL, "integration step size.");

    hestOptAdd(&hopt, "dt", "output delta time", airTypeDouble, 1, 1, &tmpDT, NULL, "Output Flow Map delta time.");

    hestOptAdd(&hopt, "minx", "min x coord", airTypeDouble, 1, 1, &tmpMinX, "nan", "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord", airTypeDouble, 1, 1, &tmpMaxX, "nan", "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord", airTypeDouble, 1, 1, &tmpMinY, "nan", "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord", airTypeDouble, 1, 1, &tmpMaxY, "nan", "max y coord. of bounding box.");
    hestOptAdd(&hopt, "minz", "min z coord", airTypeDouble, 1, 1, &tmpMinZ, "nan", "min z coord. of bounding box.");
    hestOptAdd(&hopt, "maxz", "max z coord", airTypeDouble, 1, 1, &tmpMaxZ, "nan", "max z coord. of bounding box.");

    hestOptAdd(&hopt, "ctname", "celltree filename", airTypeString, 0, 1, &pCtName, "(none)", "Name of file in which to store celltree after construction or from which to import it if available.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute Flow Map from 3D unsteady flow",
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
    minz = tmpMinZ;
    maxz = tmpMaxZ;
}

void SaveConfigureFile()
{
    DATA_TYPE d[3];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);
    d[2] = (maxz - minz) / DATA_TYPE(nz - 1);

    FILE *fp = fopen(pResultPathname, "w");
    if (fp == NULL) {
        printf("Cannot open file %s to write.\n", pResultPathname);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "Computing Flow Map from %s\n", pFlowName);
    fprintf(fp, "\tBounding box (%f, %f, %f, %f, %f, %f)\n", minx, maxx, miny, maxy, minz, maxz);
    fprintf(fp, "\tSampling resolution %d x %d x %d\n", nx, ny, nz);
    fprintf(fp, "\tdx = %.10lf, dy = %.10lf, dz = %.10lf\n", d[0], d[1], d[2]);
    fprintf(fp, "\tStarting time %f, ending time %f, time delta %f\n", startTime, endTime, deltaTime);
    fprintf(fp, "End\n");
    fclose(fp);
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    std::string flow_name(pFlowName);

    bool failed = false;

    if (flow_name == "ABC" || flow_name == "abc") {
        pFlow = new ZD::CFlowABC<DATA_TYPE>(1.0, 1.0, 1.0);
    }
    else if (flow_name == "delta_wing" || flow_name == "delta") {
        DATA_TYPE minTime = std::min(startTime, endTime);
        DATA_TYPE maxTime = std::max(startTime, endTime);
        std::string ctname = pCtName;
        if (ctname == "(none)" || ctname.empty())
            pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, true);
        else
            pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, true, nullptr, ctname);
    }
    else {
        std::cerr << "Unknown 3D flow name: " << pFlowName << std::endl;
        return EXIT_FAILURE;
    }

    if (nx <= 3 || ny <= 3 || nz <= 3) {
        std::cerr << "Error! sampling resolution ( " << nx << " x " << ny << "x" << nz << " ) is too small." << std::endl;
        SafeDelete(pFlow);
        return EXIT_FAILURE;
    }

    ZD::CPoint<DATA_TYPE, 3> bbox[2];
    pFlow->GetBBox(bbox[0], bbox[1]);
    minx = std::isnan(minx) ? bbox[0][0] : minx;
    miny = std::isnan(miny) ? bbox[0][1] : miny;
    minz = std::isnan(minz) ? bbox[0][2] : minz;
    maxx = std::isnan(maxx) ? bbox[1][0] : maxx;
    maxy = std::isnan(maxy) ? bbox[1][1] : maxy;
    maxz = std::isnan(maxz) ? bbox[1][2] : maxz;

    SaveConfigureFile();

// #ifdef _OPENMP
//     omp_set_num_threads(OMP_MAX_THREADS);
// #endif

    ZD::CPoint<DATA_TYPE, 3> *particles = new ZD::CPoint<DATA_TYPE, 3>[nx*ny*nz];

    // initialize seeding position
    DATA_TYPE d[3];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);
    d[2] = (maxz - minz) / DATA_TYPE(nz - 1);
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles[(z*ny+y)*nx+x] = ZD::CPoint<DATA_TYPE, 3>(minx + d[0] * x, miny + d[1] * y, minz + d[2] * z);
            }
        }
    }

    // compute flow map
    pFlow->SetStepSize(stepSize);
    DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    DATA_TYPE remainTime = std::abs(endTime - startTime);
    DATA_TYPE integratedTime = 0.0;
    while (remainTime > ZD_EPSILON) {
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
        for (int i = 0; i < nx*ny*nz; ++i) {
            pFlow->NextTime(particles[i], currentStartTime, deltaTime, direction,
                    new integ_type(stepSize, deltaTime/100.));

#ifdef _OPENMP
            int threadID = omp_get_thread_num();
            thread_count[threadID]++;
            if (threadID == 0) {
                int total = 0;
                for (int id = 0; id < OMP_MAX_THREADS; ++id)
                    total += thread_count[id];
                int prec = (double)total / double(nx*ny*nz) * 1000.0;
                if (prec > last_prec) {
                    //printf("\r(%d / %d) pathlines finished.", total, nx*ny*nz);
                    printf("\r%.1f%% pathlines finished.", prec/10.0);
                    fflush(stdout);
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
        if (direction > 0.0)
            sprintf(pathname, "%sflowmap_dt=%011.05f.nrrd", temp.c_str(), integratedTime*direction);
        else
            sprintf(pathname, "%sflowmap_dt=%012.05f.nrrd", temp.c_str(), integratedTime*direction);

        ZD::CField3<DATA_TYPE, 3> *pFlowMap = new ZD::CField3<DATA_TYPE, 3>();

        pFlowMap->CreateField(nx, ny, nz, particles);
        pFlowMap->SaveNrrdFile(pathname);
        SafeDelete(pFlowMap);
    }


    SafeDelete(pFlow);
    SafeDeleteArray(particles);

    return EXIT_SUCCESS;

}
