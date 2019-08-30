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
#elif defined _TBB
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#endif

#include <limits>
#include <ctime>

#include <teem/hest.h>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

typedef DATA_TYPE                          value_type;
typedef ZD::CIntegratorRK45<value_type, 2> integrator_type;
typedef ZD::CPoint<value_type, 2>          point_type;
typedef ZD::CFlow<value_type, 2>           flow_type;


flow_type *pFlow = nullptr;

char *pFlowName = nullptr;
char *pResultPathname = nullptr;
char *pDataPath = nullptr;

int nx, ny;

value_type startTime = 0.0;
value_type endTime = 0.0;
value_type stepSize = 0.0;

value_type minx, maxx, miny, maxy;

value_type deltaTime = 0.0;

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
    value_type d[2];
    d[0] = (maxx - minx) / value_type(nx - 1);
    d[1] = (maxy - miny) / value_type(ny - 1);

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

void flowmap_iteration(point_type& p, value_type t0, value_type dt, 
                       value_type step_size, value_type dir) {
    pFlow->NextTime(p, t0, dt, dir, 
                    new integrator_type(step_size, dt/100.));
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    value_type minTime = std::min(startTime, endTime);
    value_type maxTime = std::max(startTime, endTime);
    if (strcmp(pFlowName, "double_gyre") == 0) {
        pFlow = new ZD::CFlowDoubleGyre<value_type>();
    }
    else if (strcmp(pFlowName, "steady_double_gyre") == 0) {
        pFlow = new ZD::CFlowSteadyDoubleGyre<value_type>();
    }
    else if (strcmp(pFlowName, "meandering_jet") == 0) {
        pFlow = new ZD::CFlowMeanderingJet<value_type>();
    }
    else if (strcmp(pFlowName, "convection") == 0) {
        pFlow = new ZD::CFlowConvection<value_type>(minTime, maxTime, true);
    }
    else if (strcmp(pFlowName, "boussinesq") == 0) {
        pFlow = new ZD::CFlowBoussinesq<value_type>(minTime, maxTime);
    }
    else if (strcmp(pFlowName, "gaussian_vortices") == 0) {
        pFlow = new ZD::CFlowGaussianVortices<value_type>(minTime, maxTime);
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

    point_type bbox[2];
    pFlow->GetBBox(bbox[0], bbox[1]);
    minx = std::isnan(minx) ? bbox[0][0] : minx;
    miny = std::isnan(miny) ? bbox[0][1] : miny;
    maxx = std::isnan(maxx) ? bbox[1][0] : maxx;
    maxy = std::isnan(maxy) ? bbox[1][1] : maxy;


    SaveConfigureFile();

    point_type *particles = new point_type[nx*ny];

    // initialize seeding position
    value_type d[2];
    d[0] = (maxx - minx) / value_type(nx - 1);
    d[1] = (maxy - miny) / value_type(ny - 1);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            particles[y*nx+x] = point_type(minx + d[0] * x, miny + d[1] * y);
        }
    }

    // compute flow map
    pFlow->SetStepSize(stepSize);
    value_type direction = endTime > startTime ? 1.0 : -1.0;
    value_type remainTime = std::abs(endTime - startTime);
    value_type integratedTime = 0.0;
    size_t step=0;
    while (remainTime > 0.0) {
        ++step;
        std::cout << "\ncurrently computing step #" << step << ". " << remainTime << " s. remaining\n";
        
        value_type currentStartTime = startTime + direction * integratedTime;

        std::clock_t start_time = std::clock();
#ifdef _OPENMP
        int last_prec = 0;
        omp_set_num_threads(OMP_MAX_THREADS);
        int thread_count[OMP_MAX_THREADS];
        memset(thread_count, 0, sizeof(int)*OMP_MAX_THREADS);
        
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i=0; i<nx*ny; ++i) {
            int threadID = omp_get_thread_num();
            thread_count[threadID]++;
            flowmap_iteration(particles[i], currentStartTime, deltaTime, 
                              stepSize, direction);
            
            if (threadID == 0) {
                int total = 0;
                for (int id = 0; id < OMP_MAX_THREADS; ++id)
                    total += thread_count[id];
                ++nb_done;
                int total = nb_done;
                int prec = int(10000.0 * (double)total / double(nx*ny));
                
                if (prec > last_prec) {
                    printf("\r(%d / %d) pathlines completed.", total, nx*ny);
                    last_prec = prec;
                }
            }
        }
#elif defined _TBB
        tbb::atomic<int> tbb_progress_counter = 0;
        tbb::atomic<int> tbb_last_done=0;
        int stride = nx*ny/100;
        
        std::cout << "Running flow map computation with TBB\n";
    
        tbb::parallel_for(tbb::blocked_range<int>(0,nx*ny),
                           [&](tbb::blocked_range<int> r) {
            for (int i=r.begin(); i!=r.end(); ++i) {
                int ndone = tbb_progress_counter++;
                flowmap_iteration(particles[i], currentStartTime, deltaTime, 
                                  stepSize, direction);
                int last_done = tbb_last_done;
                if (true || ndone == last_done + stride) {
                    printf("\r(%d / %d) pathlines completed.", ndone, nx*ny);
                    tbb_last_done.fetch_and_store(ndone);
                }
            }
        });
#else
        int ndone = 0;
        int last_prec = 0;
        for (int i=0; i<nx*ny; ++i) {
            flowmap_iteration(particles[i], currentStartTime, deltaTime, 
                              stepSize, direction);
            int prec = int(10000.0 * (double)++ndone / double(nx*ny));
                
            if (prec > last_prec) {
                printf("\r(%d / %d) pathlines completed.", ndone, nx*ny);
                last_prec = prec;
            }
        }
#endif
        std::clock_t end_time = std::clock();
        auto elapsed = (end_time-start_time)/CLOCKS_PER_SEC;
        std::cout << "time: " << elapsed
            << " (" << (double)(nx*ny)/(double)elapsed << " Hz)\n";

        integratedTime += deltaTime;
        remainTime -= deltaTime;

        // save flow map
        std::string temp = ZD::CFileTool::GetFilePath(pResultPathname);
        char pathname[ZD_PATHNAME_LENGTH];
        sprintf(pathname, "%sflowmap_dt=%09.04f.nrrd", temp.c_str(), integratedTime*direction);
        ZD::CField2<value_type, 2> *pFlowMap = new ZD::CField2<value_type, 2>();
        pFlowMap->CreateField(nx, ny, particles);
        pFlowMap->SaveNrrdFile(pathname);
        SafeDelete(pFlowMap);
    }


    SafeDelete(pFlow);
    SafeDeleteArray(particles);
    return EXIT_SUCCESS;
}
