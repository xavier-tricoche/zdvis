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

const DATA_TYPE VALUE_THRESHOLD = 0.00001;
const DATA_TYPE TIME_THRESHOLD = 0.00001;

ZD::CFlow<DATA_TYPE, 2> *pFlow = nullptr;

char *pFlowName = nullptr;
char *pResultPathname = nullptr;
char *pDataPath = nullptr;

int nx, ny;

DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;

DATA_TYPE minx, maxx, miny, maxy;

DATA_TYPE sigma;

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
    double tmpS;

    hestOptAdd(&hopt, "o",    "pathname",         airTypeString, 1, 1, &pResultPathname, NULL,  "output FSLE result pathname.");
    hestOptAdd(&hopt, "nx",   "resolution x",     airTypeInt,    1, 1, &nx,              NULL,  "sampling resolution along x.");
    hestOptAdd(&hopt, "ny",   "resolution y",     airTypeInt,    1, 1, &ny,              NULL,  "sampling resolution along y.");
    hestOptAdd(&hopt, "f",    "2D unsteady flow", airTypeString, 1, 1, &pFlowName,       NULL,  "2D unsteady flow name, support: double_gyre, boussinesq, meandering_jet, convection.");
    hestOptAdd(&hopt, "st",   "start time",       airTypeDouble, 1, 1, &tmpST,           NULL,  "integration starting time.");
    hestOptAdd(&hopt, "et",   "end time",         airTypeDouble, 1, 1, &tmpET,           NULL,  "integration ending time.");
    hestOptAdd(&hopt, "ss",   "step size",        airTypeDouble, 1, 1, &tmpSS,           NULL,  "integration step size.");
    hestOptAdd(&hopt, "s",    "sigma",            airTypeDouble, 1, 1, &tmpS,            NULL,  "sigma (FSLE parameter).");
    hestOptAdd(&hopt, "minx", "min x coord",      airTypeDouble, 1, 1, &tmpMinX,         "nan", "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord",      airTypeDouble, 1, 1, &tmpMaxX,         "nan", "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord",      airTypeDouble, 1, 1, &tmpMinY,         "nan", "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord",      airTypeDouble, 1, 1, &tmpMaxY,         "nan", "max y coord. of bounding box.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute FSLE from 2D unsteady flow",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    startTime = tmpST;
    endTime = tmpET;
    stepSize  = tmpSS;

    minx = tmpMinX;
    maxx = tmpMaxX;
    miny = tmpMinY;
    maxy = tmpMaxY;

    sigma = tmpS;
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    DATA_TYPE minTime = std::min(startTime, endTime);
    DATA_TYPE maxTime = std::max(startTime, endTime);
    if (strcmp(pFlowName, "double_gyre") == 0) {
        pFlow = new ZD::CFlowDoubleGyre<DATA_TYPE>();
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

    pFlow->SetStepSize(stepSize);

    ZD::CField2<DATA_TYPE, 1> *pFSLE = new ZD::CField2<DATA_TYPE, 1>();
    pFSLE->CreateField(nx - 2, ny - 2, nullptr);

    DATA_TYPE d[2];
    d[0] = (maxx - minx) / DATA_TYPE(nx - 1);
    d[1] = (maxy - miny) / DATA_TYPE(ny - 1);

    DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    const int TOTAL_STEP = 100;
    DATA_TYPE deltaTime = std::abs(endTime - startTime) / (DATA_TYPE)TOTAL_STEP;
    //deltaTime = deltaTime > std::abs(stepSize) ? deltaTime : std::abs(stepSize);

    DATA_TYPE dis[2];
    dis[0] = 2.0 * d[0];
    dis[1] = 2.0 * d[1];

#ifdef _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    for (int y = 1; y < ny-1; ++y) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int x = 1; x < nx-1; ++x) {
            // initialize seeds
            ZD::CPoint<DATA_TYPE, 2> seeds[9];
            for (int i = 0; i < 9; ++i) {
                seeds[i] = ZD::CPoint<DATA_TYPE, 2>(minx + d[0] * (x + ZD::neighbors_D2N9[i][0]),
                    miny + d[1] * (y + ZD::neighbors_D2N9[i][1]));
            }

            // find the section
            DATA_TYPE lastValue = 0.0;
            DATA_TYPE currValue = 0.0;
            DATA_TYPE lastTime = startTime;
            DATA_TYPE currTime = 0.0;
            ZD::CPoint<DATA_TYPE, 2> lastPos[9];
            ZD::CPoint<DATA_TYPE, 2> currPos[9];
            memcpy(lastPos, seeds, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
            DATA_TYPE totalTime = 0.0;
            int step = 0;
            for (; step < TOTAL_STEP; ++step) {
                ///* debug */
                //ZD::CPoint<DATA_TYPE, 2> tmpPos[9];
                //memcpy(tmpPos, seeds, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
                //for (int i = 0; i < 9; ++i) {
                //    pFlow->NextTime(tmpPos[i], startTime, DATA_TYPE(step + 1) * deltaTime, direction);
                //}
                //DATA_TYPE tmpValue = ZD::CFTLETool<DATA_TYPE, 2>::ComputeSigma_D2N9(tmpPos, dis);
                ///* end of debug */

                currTime = startTime + direction * (DATA_TYPE)step * deltaTime;
                memcpy(currPos, lastPos, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
                for (int i = 0; i < 9; ++i) {
                    pFlow->NextTime(currPos[i], currTime, deltaTime, direction);
                }
                currValue = ZD::CFTLETool<DATA_TYPE, 2>::ComputeSigma_D2N9(currPos, dis);
                if (currValue > sigma) {
                    lastTime = currTime;
                    break;
                } else {
                    memcpy(lastPos, currPos, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
                    lastValue = currValue;
                    lastTime = currTime + deltaTime;
                    totalTime += deltaTime;
                }
            }

            if (step < TOTAL_STEP) {
                // newton's method
                DATA_TYPE localDeltaTime = deltaTime / 2.0;
                while (true) {
                    memcpy(currPos, lastPos, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
                    for (int i = 0; i < 9; ++i) {
                        pFlow->NextTime(currPos[i], lastTime, localDeltaTime, direction);
                    }
                    currValue = ZD::CFTLETool<DATA_TYPE, 2>::ComputeSigma_D2N9(currPos, dis);
                    if (currValue > sigma) {
                        localDeltaTime /= 2.0;
                    }
                    else {
                        lastTime += direction * localDeltaTime;
                        lastValue = currValue;
                        memcpy(lastPos, currPos, sizeof(ZD::CPoint<DATA_TYPE, 2>) * 9);
                        localDeltaTime /= 2.0;
                    }
                    if (localDeltaTime < TIME_THRESHOLD || std::abs(currValue - sigma) < VALUE_THRESHOLD) {
                        break;
                    }
                }
            }

            //for (int i = 0; i < 9; ++i) {
            //    pFlow->NextTime(seeds[i], startTime, std::abs(lastTime - startTime), direction);
            //}
            //DATA_TYPE tmpValue = ZD::CFTLETool<DATA_TYPE, 2>::ComputeSigma_D2N9(seeds, dis);

            DATA_TYPE tmpValue = ZD::CFTLETool<DATA_TYPE, 2>::ComputeFTLE_D2N9(lastPos, std::abs(lastTime - startTime), dis);
            if (std::isnan(tmpValue)) {
                std::cout << x << "," << y << std::endl;
            } else {
                pFSLE->SetValue(x - 1, y - 1, ZD::CPoint<DATA_TYPE, 1>(tmpValue));
            }

            //pFSLE->SetValue(x - 1, y - 1,
            //    ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 2>::ComputeFTLE_D2N9(lastPos, std::abs(lastTime-startTime), dis)));
        }
        printf("y == %d\n", y);
    }

    pFSLE->SaveNrrdFile(pResultPathname);

    SafeDelete(pFSLE);
    SafeDelete(pFlow);

    return EXIT_SUCCESS;

    SafeDelete(pFSLE);
    SafeDelete(pFlow);

    return EXIT_SUCCESS;

}
