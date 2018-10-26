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
DATA_TYPE dx, dy;

void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpDX, tmpDY, tmpDT;

    hestOptAdd(&hopt, "i",  "flow map",              airTypeString, 1, 1, &pFlowMapPathname, NULL, "input flow map name.");
    hestOptAdd(&hopt, "o",  "FTLE",                  airTypeString, 1, 1, &pResultPathname,  NULL, "output FTLE name.");
    hestOptAdd(&hopt, "dx", "sampling resolution x", airTypeDouble, 1, 1, &tmpDX,            NULL, "sampling resolution x.");
    hestOptAdd(&hopt, "dy", "sampling resolution y", airTypeDouble, 1, 1, &tmpDY,            NULL, "sampling resolution y.");
    hestOptAdd(&hopt, "dt", "integration time",      airTypeDouble, 1, 1, &tmpDT,            NULL, "integration time.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute FTLE from 2D unsteady flow",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    dt = tmpDT;
    dx = tmpDX;
    dy = tmpDY;
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    ZD::CField2<DATA_TYPE, 2> *pFlowMap = new ZD::CField2<DATA_TYPE, 2>();
    pFlowMap->OpenNrrdFile(pFlowMapPathname);
    const int *size = pFlowMap->GetSize();

    ZD::CField2<DATA_TYPE, 1> *ftle = new ZD::CField2<DATA_TYPE, 1>();
    int ftleSize[2] = { size[0], size[1]};
    ftle->CreateField(ftleSize);

    DATA_TYPE dis[2];
    dis[0] = 2.0 * dx;
    dis[1] = 2.0 * dy;

    for (int y = 1; y < size[1] - 1; ++y) {
        for (int x = 1; x < size[0] - 1; ++x) {
            ZD::CPoint<DATA_TYPE, 2> buffer[9];
            for (int i = 0; i < 9; ++i) {
                buffer[i] = pFlowMap->GetValue(x + ZD::neighbors_D2N9[i][0], y + ZD::neighbors_D2N9[i][1]);
            }
            ftle->SetValue(x, y, ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 2>::ComputeFTLE_D2N9(buffer, dt, dis)));

            //ZD::CPoint<DATA_TYPE, 2> buffer[4];
            //for (int i = 0; i < 4; ++i) {
            //    buffer[i] = pFlowMap->GetValue(x + ZD::neighbors_D2N4[i][0], y + ZD::neighbors_D2N4[i][1]);
            //}
            //ftle->SetValue(x-1, y-1, ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 2>::ComputeFTLE_D2N4(buffer, dt, dis)));
        }
    }

    ftle->SaveNrrdFile(pResultPathname);

    SafeDelete(pFlowMap);
    SafeDelete(ftle);

    return EXIT_SUCCESS;
}
