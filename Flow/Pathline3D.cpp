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

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <string>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

ZD::CFlow<DATA_TYPE, 3> *pFlow = nullptr;

DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;

char *pFlowName = nullptr;
char *pPathlinePathname = nullptr;
char *pDataPath = nullptr;

ZD::CPoint<DATA_TYPE, 3> *pSeeds = nullptr;
int seedCount = 0;
int timeFlag = 0;
int uniformSeedFlag = 0;

int nx, ny, nz;
DATA_TYPE minx, maxx, miny, maxy, minz, maxz;


void Initialize(int argc, char *argv[])
{
hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double seed_pos[3];
    char *pSeedPathname = nullptr;
    double tmpMinX, tmpMaxX, tmpMinY, tmpMaxY, tmpMinZ, tmpMaxZ;
    double tmpST, tmpET, tmpSS;

    hestOptAdd(&hopt, "o",    "pathline",         airTypeString, 1, 1, &pPathlinePathname, NULL, "output pathline name.");
    hestOptAdd(&hopt, "f",    "3D unsteady flow", airTypeString, 1, 1, &pFlowName, NULL, "3D unsteady flow name, support: abc, delta_wing, test.");
    hestOptAdd(&hopt, "sp",   "seed point",       airTypeDouble, 3, 3, &seed_pos, "nan nan nan", "input seed point.");
    hestOptAdd(&hopt, "su",   "seed uniform",     airTypeBool,   0, 0, &uniformSeedFlag, NULL, "uniform seed, using nx, ny, minx, maxx, miny, maxy.");
    hestOptAdd(&hopt, "nx",   "resolution x",     airTypeInt,    1, 1, &nx, "0", "sampling resolution along x.");
    hestOptAdd(&hopt, "ny",   "resolution y",     airTypeInt,    1, 1, &ny, "0", "sampling resolution along y.");
    hestOptAdd(&hopt, "nz",   "resolution z",     airTypeInt,    1, 1, &nz, "0", "sampling resolution along z.");
    hestOptAdd(&hopt, "minx", "min x coord",      airTypeDouble, 1, 1, &tmpMinX, "nan", "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord",      airTypeDouble, 1, 1, &tmpMaxX, "nan", "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord",      airTypeDouble, 1, 1, &tmpMinY, "nan", "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord",      airTypeDouble, 1, 1, &tmpMaxY, "nan", "max y coord. of bounding box.");
    hestOptAdd(&hopt, "minz", "min z coord",      airTypeDouble, 1, 1, &tmpMinZ, "nan", "min z coord. of bounding box.");
    hestOptAdd(&hopt, "maxz", "max z coord",      airTypeDouble, 1, 1, &tmpMaxZ, "nan", "max z coord. of bounding box.");
    hestOptAdd(&hopt, "time", "output time file", airTypeBool,   0, 0, &timeFlag, NULL, "output the time information for each vertice in pathline");
    hestOptAdd(&hopt, "st",   "start time",       airTypeDouble, 1, 1, &tmpST, NULL, "pathline integration starting time.");
    hestOptAdd(&hopt, "et",   "end time",         airTypeDouble, 1, 1, &tmpET, NULL, "pathline integration ending time.");
    hestOptAdd(&hopt, "ss",   "step size",        airTypeDouble, 1, 1, &tmpSS, NULL, "pathline integration step size.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute Pathline from 3D unsteady flow",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    minx = tmpMinX;
    maxx = tmpMaxX;
    miny = tmpMinY;
    maxy = tmpMaxY;
    minz = tmpMinZ;
    maxz = tmpMaxZ;

    startTime = tmpST;
    endTime   = tmpET;
    stepSize  = tmpSS;

    if (std::isnan(seed_pos[0]) != true) {
        seedCount = 1;
        pSeeds = new ZD::CPoint<DATA_TYPE, 3>[seedCount];
        pSeeds[0][0] = seed_pos[0];
        pSeeds[0][1] = seed_pos[1];
        pSeeds[0][2] = seed_pos[2];
    }
    else if (uniformSeedFlag == 1) {
        seedCount = nx * ny * nz;
        pSeeds = new ZD::CPoint<DATA_TYPE, 3>[seedCount];
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    pSeeds[(z*ny+y)*nx+x][0] = minx + (maxx - minx) * DATA_TYPE(x) / DATA_TYPE(nx - 1);
                    pSeeds[(z*ny+y)*nx+x][1] = miny + (maxy - miny) * DATA_TYPE(y) / DATA_TYPE(ny - 1);
                    pSeeds[(z*ny+y)*nx+x][2] = minz + (maxz - minz) * DATA_TYPE(z) / DATA_TYPE(nz - 1);
                }
            }
        }
    }
    else {
        printf("Please specify a seed method\n");
        exit(EXIT_FAILURE);
    }
}


void SaveConfigureFile()
{
    DATA_TYPE d[3];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);
    d[2] = (maxz - minz) / DATA_TYPE(nz - 1);

    std::string tmp = ZD::CFileTool::RemoveFileExtension(pPathlinePathname);
    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%s_conf.txt", tmp.c_str());
    FILE *fp = fopen(pathname, "w");
    if (fp == NULL) {
        printf("Cannot open file %s to write.\n", pathname);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "Computing Pathlines from %s\n", pFlowName);
    fprintf(fp, "\tBounding box (%f, %f, %f, %f, %f, %f)\n", minx, maxx, miny, maxy, minz, maxz);
    fprintf(fp, "\tSampling resolution %d x %d x %d\n", nx, ny, nz);
    fprintf(fp, "\tdx = %.10lf, dy = %.10lf, dz = %.10lf\n", d[0], d[1], d[2]);
    fprintf(fp, "\tStarting time %f, ending time %f, step size %f\n", startTime, endTime, stepSize);
    fprintf(fp, "End\n");

    std::string path = ZD::CFileTool::GetFilePath(pPathlinePathname);
    std::string ext = ZD::CFileTool::GetFileExtension(pPathlinePathname);
    for (int z = 0; z < nz; ++z) {
        DATA_TYPE zz = minz + (maxz - minz) * DATA_TYPE(z) / DATA_TYPE(nz - 1);
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%spathline_z=%09.04f.%s", path.c_str(), zz, ext.c_str());
        fprintf(fp, "%s", pathname);

        if (timeFlag == 1) {
            sprintf(pathname, "%stime_z=%09.04f.la", path.c_str(), zz);
            fprintf(fp, " %s", pathname);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void Pathline3DSingleSeeding()
{
    const DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    const DATA_TYPE totalTime = std::abs(endTime - startTime);
    const DATA_TYPE dt = (endTime - startTime) / DATA_TYPE(PATHLINE_MAX_COUNT);

    ZD::CLine<DATA_TYPE, 3> *pLines = nullptr;
    ZD::CLineAttribute<DATA_TYPE, 1> *pTimes = nullptr;

    pLines = new ZD::CLine<DATA_TYPE, 3>[seedCount];
    if (timeFlag)
        pTimes = new ZD::CLineAttribute<DATA_TYPE, 1>[seedCount];

    for (int i = 0; i < seedCount; ++i) {
        std::vector<ZD::CPoint<DATA_TYPE, 3>> buffer;
        std::vector<ZD::CPoint<DATA_TYPE, 1>> timeBuffer;
        buffer.reserve(PATHLINE_MAX_COUNT + 5);
        timeBuffer.reserve(PATHLINE_MAX_COUNT + 5);

        DATA_TYPE integratedTime = 0.0;

        buffer.push_back(pSeeds[i]);
        timeBuffer.push_back(startTime + integratedTime * direction);

        for (int j = 0; j < PATHLINE_MAX_COUNT; ++j) {
            DATA_TYPE st = startTime + dt * DATA_TYPE(j);
            DATA_TYPE tt = pFlow->NextTime(pSeeds[i], st, std::fabs(dt), direction);
            buffer.push_back(pSeeds[i]);
            integratedTime += tt;
            timeBuffer.push_back(startTime + integratedTime * direction);
        }
        pLines[i].CreateLine(buffer);
        pLines[i].ComputeColor();
        if (timeFlag)
            pTimes[i].CreateLineAttribute(timeBuffer);
    }

    std::string ext = ZD::CFileTool::GetFileExtension(pPathlinePathname);
    if (ext == "vtk" || ext == "VTK")
        SaveLinesVTK(pPathlinePathname, pLines, seedCount);
    else
        SaveLines(pPathlinePathname, pLines, seedCount);

    if (timeFlag) {
        std::string tmp = ZD::CFileTool::RemoveFileExtension(pPathlinePathname);
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%s_time.la", tmp.c_str());
        SaveLineAttributes(pathname, pTimes, seedCount);
    }

    SafeDeleteArray(pLines);
    SafeDeleteArray(pTimes);
}

void Pathline3DUniformSeeding()
{
    const DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    const DATA_TYPE totalTime = std::abs(endTime - startTime);

    DATA_TYPE dt = (endTime - startTime) / DATA_TYPE(PATHLINE_MAX_COUNT);
    int maxCount = PATHLINE_MAX_COUNT;
    if (std::fabs(dt) < stepSize) {
        dt = dt > 0.0 ? stepSize : -stepSize;
        maxCount = int(std::fabs((endTime - startTime) / stepSize) + ZD_EPSILON);
    }

    ZD::CTimeTool<DATA_TYPE> timer;
    timer.StartTimer();

#ifdef _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    for (int z = 0; z < nz; ++z) {
        ZD::CLine<DATA_TYPE, 3> *pLines = nullptr;
        ZD::CLineAttribute<DATA_TYPE, 1> *pTimes = nullptr;

        pLines = new ZD::CLine<DATA_TYPE, 3>[nx*ny];
        if (timeFlag)
            pTimes = new ZD::CLineAttribute<DATA_TYPE, 1>[nx*ny];

#ifdef _OPENMP
        int thread_count[OMP_MAX_THREADS];
        memset(thread_count, 0, sizeof(int)*OMP_MAX_THREADS);
        int last_prec = 0;
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int i = 0; i < nx*ny; ++i) {
            std::vector<ZD::CPoint<DATA_TYPE, 3>> buffer;
            std::vector<ZD::CPoint<DATA_TYPE, 1>> timeBuffer;
            buffer.reserve(maxCount + 5);
            timeBuffer.reserve(maxCount + 5);

            DATA_TYPE integratedTime = 0.0;

            ZD::CPoint<DATA_TYPE, 3> seed = pSeeds[z*nx*ny + i];

            buffer.push_back(seed);
            timeBuffer.push_back(startTime + integratedTime * direction);

            for (int j = 0; j < maxCount; ++j) {
                DATA_TYPE st = startTime + dt * DATA_TYPE(j);
                DATA_TYPE tt = pFlow->NextTime(seed, st, std::fabs(dt), direction);
                if (pFlow->CheckPosition(seed) == false)
                    break;
                buffer.push_back(seed);
                integratedTime += tt;
                timeBuffer.push_back(startTime + integratedTime * direction);
            }
            pLines[i].CreateLine(buffer);
            pLines[i].ComputeColor();
            if (timeFlag)
                pTimes[i].CreateLineAttribute(timeBuffer);

#ifdef _OPENMP
            int threadID = omp_get_thread_num();
            thread_count[threadID]++;
            if (threadID == 0) {
                int total = 0;
                for (int id = 0; id < OMP_MAX_THREADS; ++id)
                    total += thread_count[id];
                int prec = (double)total / double(nx*ny) * 100.0;
                if (prec > last_prec) {
                    printf("\r(%d / %d) pathlines finished.", total, nx*ny);
                    last_prec = prec;
                }
            }
#endif
        }

#ifdef _OPENMP
        printf("\n");
#endif

        std::string ext = ZD::CFileTool::GetFileExtension(pPathlinePathname);
        std::string path = ZD::CFileTool::GetFilePath(pPathlinePathname);
        DATA_TYPE zz = minz + (maxz - minz) * DATA_TYPE(z) / DATA_TYPE(nz - 1);
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%spathline_z=%09.04f.%s", path.c_str(), zz, ext.c_str());
        if (ext == "vtk" || ext == "VTK")
            SaveLinesVTK(pathname, pLines, nx*ny);
        else
            SaveLines(pathname, pLines, nx*ny);

        if (timeFlag) {
            sprintf(pathname, "%stime_z=%09.04f.la", path.c_str(), zz);
            SaveLineAttributes(pathname, pTimes, nx*ny);
        }

        SafeDeleteArray(pLines);
        SafeDeleteArray(pTimes);
    }

    timer.StopTimer();
    std::cout << "Running time " << timer.Duration() << " seconds." << std::endl;
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    if (strcmp(pFlowName, "abc") == 0) {
        pFlow = new ZD::CFlowABC<DATA_TYPE>();
    }
    else if (strcmp(pFlowName, "delta_wing") == 0) {
        DATA_TYPE minTime = std::min(startTime, endTime);
        DATA_TYPE maxTime = std::max(startTime, endTime);
        pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, true);
    }
    else if (strcmp(pFlowName, "test") == 0) {
        pFlow = new ZD::CFlowTest3D<DATA_TYPE>();
    }
    else {
        std::cerr << "Unknown 3D flow name: " << pFlowName << std::endl;
        return EXIT_FAILURE;
    }

    if (uniformSeedFlag == 1)
        SaveConfigureFile();

    pFlow->SetStepSize(stepSize);

    if (seedCount == 1) {
        /* single seeding */
        Pathline3DSingleSeeding();
    }
    else {
        /* uniform seeding */
        Pathline3DUniformSeeding();
    }

    SafeDelete(pFlow);
    SafeDeleteArray(pSeeds);

    return EXIT_SUCCESS;
}
