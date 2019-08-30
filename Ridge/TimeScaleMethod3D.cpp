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

#include <cmath>

#ifndef DATA_TYPE
#define DATA_TYPE double
#endif

// change me!
// #define TIMESCALE3D_OUTPUT_PATH "C:\\Ziang\\project\\dataset\\double_point_load\\non_normalized.xls"

char **pFTLEPathnames = nullptr;
char *pResultPathname = nullptr;
char *pTimePathname = nullptr;
DATA_TYPE *pTimes = nullptr;
int count = 0;
DATA_TYPE gaussianKernel = 0.0;
DATA_TYPE threshold = 0.0;
int w, h, d;

enum TS_METHOD {
    TE_METHOD_TS,
    TE_METHOD_MF
};
TS_METHOD method;

void Initialize(int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    char *pConfPathname;
    double tmpGK, tmpT;
    char *tmpMethod;

    hestOptAdd(&hopt, "i",    "filename",    airTypeString, 1, 1, &pConfPathname,   NULL,  "input FTLE configure file.");
    hestOptAdd(&hopt, "o",    "filename",    airTypeString, 1, 1, &pResultPathname, NULL,  "output temporal-scale feature field filename.");
    hestOptAdd(&hopt, "time", "filename",    airTypeString, 1, 1, &pTimePathname,   NULL,  "output temporal-scale time field filename.");
    hestOptAdd(&hopt, "g",    "kernel size", airTypeDouble, 1, 1, &tmpGK,           "0.0", "Gaussian filter kernel size.");
    hestOptAdd(&hopt, "m",    "method",      airTypeString, 1, 1, &tmpMethod,       "ts",  "method: time scale (ts) or maximum FTLE (mf).");
    //hestOptAdd(&hopt, "t",    "threshold",   airTypeDouble, 1, 1, &tmpT,            "0.0", "ridge strength threshold.");

    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "3D Temporal-Scale method",
        AIR_TRUE, AIR_TRUE, AIR_TRUE);

    hestParmFree(hparm);

    //threshold = tmpT;
    gaussianKernel = tmpGK;

    std::ifstream file_ifs;
    file_ifs.open(pConfPathname, std::ifstream::in);
    if (!file_ifs.is_open()) {
        std::cerr << "Error! cannot open configure file" << pConfPathname << std::endl;
        exit(EXIT_FAILURE);
    }
    file_ifs >> count;
    if (count < 1) {
        std::cerr << "Error! no file field in configure file" << std::endl;
        file_ifs.close();
        exit(EXIT_FAILURE);
    }
    pFTLEPathnames = new char*[count];
    pTimes = new DATA_TYPE[count];
    for (int i = 0; i < count; ++i) {
        pFTLEPathnames[i] = new char[ZD_PATHNAME_LENGTH];
        file_ifs >> pTimes[i] >> pFTLEPathnames[i];
    }
    file_ifs.close();

    if (ZD::CStringTool::ToLower(tmpMethod) == "mf") {
        method = TS_METHOD::TE_METHOD_MF;
    }
    else {
        method = TS_METHOD::TE_METHOD_TS;
    }
}

void ComputeRidgeStrength(ZD::CField3<DATA_TYPE, 9> *pHess,
    ZD::CField3<DATA_TYPE, 1> *pStre)
{
#if _OPENMP
#pragma omp parallel for
#endif
    for (int z = 0; z < d; ++z) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                ZD::CPoint<DATA_TYPE, 3> evec[3];
                ZD::CPoint<DATA_TYPE, 1> eval[3];
                ZD::CPoint<DATA_TYPE, 9> tensor;
                tensor = pHess->GetValue(x, y, z);
                ZD::CEigenTool<DATA_TYPE>::Eigen3D(tensor, eval, evec);
                if (eval[2][0] > 0.0)
                    pStre->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(0.0));
                else
                    pStre->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(std::abs(eval[2][0])));

            }
        }
    }
}

void TimeScaleMethod3D()
{
    //
    ZD::CPoint<DATA_TYPE, 3> long_time_point  = ZD::CPoint<DATA_TYPE, 3>(165.0, 104.0, 90.0);
    ZD::CPoint<DATA_TYPE, 3> short_time_point = ZD::CPoint<DATA_TYPE, 3>(153.0, 149.0, 90.0);
    double ftle_values[2][100];
    double stre_values[2][100];
    //

    ZD::CField3<DATA_TYPE, 1> *pTime = new ZD::CField3<DATA_TYPE, 1>();
    pTime->CreateField(w, h, d, nullptr);
    ZD::CField3<DATA_TYPE, 1> *pFeature = new ZD::CField3<DATA_TYPE, 1>();
    pFeature->CreateField(w, h, d, nullptr);

    for (int z = 0; z < d; ++z) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                pFeature->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(0.0));
                pTime->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(pTimes[0]));
            }
        }
    }

    for (int i = 0; i < count; ++i) {
        ZD::CField3<DATA_TYPE, 1> *pTmpFTLE = new ZD::CField3<DATA_TYPE, 1>();
        pTmpFTLE->OpenNrrdFile(pFTLEPathnames[i]);

        if (gaussianKernel > 0.1)
            pTmpFTLE->GaussianFilter(gaussianKernel);

        ZD::CField3<DATA_TYPE, 9> *pTmpHess = new ZD::CField3<DATA_TYPE, 9>();
        ZD::CField3<DATA_TYPE, 1> *pTmpStre = new ZD::CField3<DATA_TYPE, 1>();
        pTmpHess->CreateField(w, h, d, nullptr);
        pTmpStre->CreateField(w, h, d, nullptr);

        ZD::CGradientTool<DATA_TYPE>::ComputeGradHess(pTmpFTLE, nullptr, pTmpHess, w, h, d);
        ComputeRidgeStrength(pTmpHess, pTmpStre);

        for (int z = 0; z < d; ++z) {
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    DATA_TYPE stre_value = pTmpStre->GetValue(x, y, z)[0];
                    DATA_TYPE ftle_value = pTmpFTLE->GetValue(x, y, z)[0];
                    DATA_TYPE feature_value = stre_value * ftle_value;
                    if (feature_value > pFeature->GetValue(x, y, z)[0]) {
                        pFeature->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(feature_value));
                        pTime->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(pTimes[i]));
                    }
                }
            }
        }

        //
        ftle_values[0][i] = pTmpFTLE->GetValue(short_time_point)[0];
        ftle_values[1][i] = pTmpFTLE->GetValue(long_time_point)[0];
        stre_values[0][i] = pTmpStre->GetValue(short_time_point)[0];
        stre_values[1][i] = pTmpStre->GetValue(long_time_point)[0];
        //

        SafeDelete(pTmpFTLE);
        SafeDelete(pTmpStre);
        SafeDelete(pTmpHess);
    }

    pFeature->SaveNrrdFile(pResultPathname);
    pTime->SaveNrrdFile(pTimePathname);

    SafeDelete(pTime);
    SafeDelete(pFeature);

#ifdef TIMESCALE3D_OUTPUT_PATH
    //
    std::fstream fp(TIMESCALE3D_OUTPUT_PATH, std::ios::out);
    for (int i = 0; i < count; ++i) {
        fp << ftle_values[0][i] << '\t'
           << stre_values[0][i] << '\t'
           << ftle_values[1][i] << '\t'
           << stre_values[1][i] << '\n';
    }
    fp.close();
    //
#endif
}

void MaximumFTLEMethod3D()
{
    ZD::CField3<DATA_TYPE, 1> *pTime = new ZD::CField3<DATA_TYPE, 1>();
    pTime->CreateField(w, h, d, nullptr);
    ZD::CField3<DATA_TYPE, 1> *pFeature = new ZD::CField3<DATA_TYPE, 1>();
    pFeature->CreateField(w, h, d, nullptr);

    for (int z = 0; z < d; ++z) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                pFeature->SetValue(x, y, d, ZD::CPoint<DATA_TYPE, 1>(0.0));
                pTime->SetValue(x, y, d, ZD::CPoint<DATA_TYPE, 1>(pTimes[0]));
            }
        }
    }


    for (int i = 0; i < count; ++i) {
        ZD::CField3<DATA_TYPE, 1> *pTmpFTLE = new ZD::CField3<DATA_TYPE, 1>();
        pTmpFTLE->OpenNrrdFile(pFTLEPathnames[i]);

        if (gaussianKernel > 0.1)
            pTmpFTLE->GaussianFilter(gaussianKernel);

        for (int z = 0; z < d; ++z) {
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    DATA_TYPE ftle_value = pTmpFTLE->GetValue(x, y, z)[0];
                    DATA_TYPE feature_value = ftle_value;
                    if (feature_value > pFeature->GetValue(x, y, z)[0]) {
                        pFeature->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(feature_value));
                        pTime->SetValue(x, y, z, ZD::CPoint<DATA_TYPE, 1>(pTimes[i]));
                    }
                }
            }
        }

        SafeDelete(pTmpFTLE);
    }

    pFeature->SaveNrrdFile(pResultPathname);
    pTime->SaveNrrdFile(pTimePathname);

    SafeDelete(pTime);
    SafeDelete(pFeature);
}

int main(int argc, char *argv[])
{
    Initialize(argc, argv);

    ZD::CTimeTool<DATA_TYPE> timer;
    timer.StartTimer();

#if _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    ZD::CField3<DATA_TYPE, 1> *tmp = new ZD::CField3<DATA_TYPE, 1>();
    tmp->OpenNrrdFile(pFTLEPathnames[0]);
    const int *size = tmp->GetSize();
    w = size[0];
    h = size[1];
    d = size[2];
    SafeDelete(tmp);

    if (method == TS_METHOD::TE_METHOD_TS) {
        TimeScaleMethod3D();
    } else if (method == TS_METHOD::TE_METHOD_MF) {
        MaximumFTLEMethod3D();
    } else {
        ;
    }

    for (int i = 0; i < count; ++i) {
        SafeDeleteArray(pFTLEPathnames[i]);
    }
    SafeDeleteArray(pFTLEPathnames);
    SafeDeleteArray(pTimes);

    timer.StopTimer();
    std::cout << "\nTime = " << timer.Duration() << " seconds." << std::endl;


    return EXIT_SUCCESS;
}
