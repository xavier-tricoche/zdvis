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

#include <teem/hest.h>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif


char *pFlowMapPathname = nullptr;
char *pResultPathname = nullptr;
DATA_TYPE dt = 0.0;
DATA_TYPE dx, dy, dz;


void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpDX, tmpDY, tmpDZ, tmpDT;

    hestOptAdd(&hopt, "i",  "flow map",              airTypeString, 1, 1, &pFlowMapPathname, NULL, "input flow map name.");
    hestOptAdd(&hopt, "o",  "FTLE",                  airTypeString, 1, 1, &pResultPathname,  NULL, "output FTLE name.");
    hestOptAdd(&hopt, "dx", "sampling resolution x", airTypeDouble, 1, 1, &tmpDX,            NULL, "sampling resolution x.");
    hestOptAdd(&hopt, "dy", "sampling resolution y", airTypeDouble, 1, 1, &tmpDY,            NULL, "sampling resolution y.");
    hestOptAdd(&hopt, "dz", "sampling resolution z", airTypeDouble, 1, 1, &tmpDZ,            NULL, "sampling resolution z.");
    hestOptAdd(&hopt, "dt", "integration time",      airTypeDouble, 1, 1, &tmpDT,            NULL, "integration time.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute FTLE from 3D unsteady flow",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    dt = tmpDT;
    dx = tmpDX;
    dy = tmpDY;
    dz = tmpDZ;
}


int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    ZD::CField3<DATA_TYPE, 3> *pFlowMap = new ZD::CField3<DATA_TYPE, 3>();
    pFlowMap->OpenNrrdFile(pFlowMapPathname);
    const int *size = pFlowMap->GetSize();

    ZD::CField3<DATA_TYPE, 1> *ftle = new ZD::CField3<DATA_TYPE, 1>();
    int ftleSize[3] = { size[0] - 2, size[1] - 2, size[2] - 2 };
    ftle->CreateField(ftleSize);

    DATA_TYPE dis[3];
    dis[0] = 2.0 * dx;
    dis[1] = 2.0 * dy;
    dis[2] = 2.0 * dz;

#ifdef _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int z = 1; z < size[2] - 1; ++z) {
        for (int y = 1; y < size[1] - 1; ++y) {
            for (int x = 1; x < size[0] - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 3> buffer[27];
                for (int i = 0; i < 27; ++i) {
                    buffer[i] = pFlowMap->GetValue(x + ZD::neighbors_D3N27[i][0], y + ZD::neighbors_D3N27[i][1], z + ZD::neighbors_D3N27[i][2]);
                }
                ftle->SetValue(x - 1, y - 1, z-1, ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 3>::ComputeFTLE_D3N27(buffer, dt, dis)));
            }
        }
    }

    ftle->SaveNrrdFile(pResultPathname);

    SafeDelete(pFlowMap);
    SafeDelete(ftle);

    return EXIT_SUCCESS;
}
