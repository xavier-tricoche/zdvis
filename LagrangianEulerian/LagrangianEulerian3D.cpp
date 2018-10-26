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

char *pConfigurePathname = nullptr;
char *pResultPathname = nullptr;
char *pFlowName = nullptr;
char **ppScalarNames = nullptr;
int scalarsCount = 0;

char **pPathlinePathnames = nullptr;
char **pTimePathnames = nullptr;

int stFlag = 0;

int nx, ny, nz;
DATA_TYPE dx, dy, dz;
DATA_TYPE minx, maxx, miny, maxy, minz, maxz;

DATA_TYPE startTime = 0.0;
DATA_TYPE endTime = 0.0;
DATA_TYPE stepSize = 0.0;

const int K = 2;
const int num = 21;

ZD::ZD_LE_Method method = ZD::ZD_LE_Method::ZD_LE_Method_Moment;

void Initialize(const int argc, char *argv[])
{
    hestOpt *hopt = nullptr;
    hestParm *hparm = hestParmNew();
    airArray *mop = airMopNew();
    char *me = argv[0];

    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;

    double tmpMinX, tmpMaxX, tmpMinY, tmpMaxY, tmpMinZ, tmpMaxZ;
    double tmpST, tmpET, tmpSS;
    char *tmpMethod;

    hestOptAdd(&hopt, "i",  "configure filename", airTypeString, 1, 1,  &pConfigurePathname, "NULL", "input configure file.");
    hestOptAdd(&hopt, "o",  "filename",           airTypeString, 1, 1,  &pResultPathname,    NULL,   "output Lagrangian/Eulerian file.");
    hestOptAdd(&hopt, "s",  "scalar name",        airTypeString, 1, -1, &ppScalarNames,      NULL,   "scalar quantities names", &scalarsCount);
    //hestOptAdd(&hopt, "st", "structure tensor",   airTypeInt,    0, 0,  &stFlag,             "0",    "output structure tensor");

    hestOptAdd(&hopt, "f",    "3D unsteady flow", airTypeString, 1, 1, &pFlowName, "NULL",   "3D unsteady flow name, support: delta_wing");
    hestOptAdd(&hopt, "nx",   "resolution x",     airTypeInt,    1, 1, &nx,        "0",      "sampling resolution along x.");
    hestOptAdd(&hopt, "ny",   "resolution y",     airTypeInt,    1, 1, &ny,        "0",      "sampling resolution along y.");
    hestOptAdd(&hopt, "nz",   "resolution z",     airTypeInt,    1, 1, &nz,        "0",      "sampling resolution along z.");
    hestOptAdd(&hopt, "minx", "min x coord",      airTypeDouble, 1, 1, &tmpMinX,   "nan",    "min x coord. of bounding box.");
    hestOptAdd(&hopt, "maxx", "max x coord",      airTypeDouble, 1, 1, &tmpMaxX,   "nan",    "max x coord. of bounding box.");
    hestOptAdd(&hopt, "miny", "min y coord",      airTypeDouble, 1, 1, &tmpMinY,   "nan",    "min y coord. of bounding box.");
    hestOptAdd(&hopt, "maxy", "max y coord",      airTypeDouble, 1, 1, &tmpMaxY,   "nan",    "max y coord. of bounding box.");
    hestOptAdd(&hopt, "minz", "min z coord",      airTypeDouble, 1, 1, &tmpMinZ,   "nan",    "min z coord. of bounding box.");
    hestOptAdd(&hopt, "maxz", "max z coord",      airTypeDouble, 1, 1, &tmpMaxZ,   "nan",    "max z coord. of bounding box.");
    hestOptAdd(&hopt, "st",   "start time",       airTypeDouble, 1, 1, &tmpST,     "nan",    "pathline integration starting time.");
    hestOptAdd(&hopt, "et",   "end time",         airTypeDouble, 1, 1, &tmpET,     "nan",    "pathline integration ending time.");
    hestOptAdd(&hopt, "ss",   "step size",        airTypeDouble, 1, 1, &tmpSS,     "nan",    "pathline integration step size.");
    hestOptAdd(&hopt, "m",    "method",           airTypeString, 1, 1, &tmpMethod, "moment", "Lagrangian-Eulerian method, moment or fixed number points");


    hestParseOrDie(hopt, argc - 1, (const char **)argv + 1, hparm,
        (const char *)me, "Compute 3D Lagrangian/Eulerian",
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

    if (strcmp(tmpMethod, "moment") == 0) {
        method = ZD::ZD_LE_Method::ZD_LE_Method_Moment;
    }
    else if (strcmp(tmpMethod, "fixed") == 0) {
        method = ZD::ZD_LE_Method::ZD_LE_Method_Fixed;
    }
    else {
        method = ZD::ZD_LE_Method::ZD_LE_Method_Moment;
    }
}

void ComputeLagrangianEulerian(const ZD::CPoint<DATA_TYPE, K> *parameters, const char *scalarName)
{
    ZD::CLETool<DATA_TYPE, K>::Initialize(std::abs(dx), std::abs(dy), std::abs(dz));

    ZD::CField3<DATA_TYPE, 1> *le = new ZD::CField3<DATA_TYPE, 1>();
    le->CreateField(nx-2, ny-2, nz-2, nullptr);

    for (int z = 1; z < nz - 1; ++z) {
        for (int y = 1; y < ny - 1; ++y) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, K> buffer[27];
                for (int k = 0; k < 27; ++k) {
                    int cx = x + ZD::neighbors_D3N27[k][0];
                    int cy = y + ZD::neighbors_D3N27[k][1];
                    int cz = z + ZD::neighbors_D3N27[k][2];
                    buffer[k] = parameters[(cz*ny+cy)*nx+cx];
                }

                DATA_TYPE value = ZD::CLETool<DATA_TYPE, K>::LagrangianEulerian3D(buffer);
                le->SetValue(x-1, y-1, z-1, ZD::CPoint<DATA_TYPE, 1>(value));
            }
        }
    }

    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sle_%s.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), scalarName);
    le->SaveNrrdFile(pathname);

    SafeDelete(le);
}

void ComputeLagrangianEulerianFixed(const ZD::CPoint<DATA_TYPE, num> *parameters, const char *scalarName)
{
    ZD::CLETool<DATA_TYPE, num>::Initialize(std::abs(dx), std::abs(dy), std::abs(dz));

    ZD::CField3<DATA_TYPE, 1> *le = new ZD::CField3<DATA_TYPE, 1>();
    le->CreateField(nx-2, ny-2, nz-2, nullptr);

    for (int z = 1; z < nz - 1; ++z) {
        for (int y = 1; y < ny - 1; ++y) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, num> buffer[27];
                for (int k = 0; k < 27; ++k) {
                    int cx = x + ZD::neighbors_D3N27[k][0];
                    int cy = y + ZD::neighbors_D3N27[k][1];
                    int cz = z + ZD::neighbors_D3N27[k][2];
                    buffer[k] = parameters[(cz*ny + cy)*nx + cx];
                }

                DATA_TYPE value = ZD::CLETool<DATA_TYPE, num>::LagrangianEulerian3D(buffer);
                le->SetValue(x-1, y-1, z-1, ZD::CPoint<DATA_TYPE, 1>(value));
            }
        }
    }

    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sle_%s.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), scalarName);
    le->SaveNrrdFile(pathname);

    SafeDelete(le);
}

void ComputeFTLE(const ZD::CField3<DATA_TYPE, 3> *flowMap)
{
    ZD::CField3<DATA_TYPE, 1> *ftle = new ZD::CField3<DATA_TYPE, 1>();
    ftle->CreateField(nx, ny, nz, nullptr);

    DATA_TYPE dis[3];
    dis[0] = 2.0 * dx;
    dis[1] = 2.0 * dy;
    dis[2] = 2.0 * dz;

    DATA_TYPE dt = std::fabs(startTime - endTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int z = 1; z < nz - 1; ++z) {
        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 3> buffer[27];
                for (int i = 0; i < 27; ++i) {
                    buffer[i] = flowMap->GetValue(x + ZD::neighbors_D3N27[i][0], y + ZD::neighbors_D3N27[i][1], z + ZD::neighbors_D3N27[i][2]);
                }
                ftle->SetValue(x - 1, y - 1, z - 1, ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 3>::ComputeFTLE_D3N27(buffer, dt, dis)));
            }
        }
    }

    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sle_%s.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), "ftle");
    ftle->SaveNrrdFile(pathname);

    SafeDelete(ftle);
}

void ComputeSpatialDistance(const ZD::CField3<DATA_TYPE, 9> *parameters)
{
    ZD::CField3<DATA_TYPE, 1> *cs = new ZD::CField3<DATA_TYPE, 1>();
    cs->CreateField(nx, ny, nz, nullptr);


    ZD::CCSTool<DATA_TYPE>::InitializeFiberFunction(dx, dy, dz);

    DATA_TYPE dt = std::fabs(startTime - endTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int z = 1; z < nz - 1; ++z) {
        for (int y = 1; y < ny - 1; ++y) {
            for (int x = 1; x < nx - 1; ++x) {
                ZD::CPoint<DATA_TYPE, 9> *buffer = new ZD::CPoint<DATA_TYPE, 9>[27];
                for (int i = 0; i < 27; ++i) {
                    buffer[i] = parameters->GetValue(x + ZD::neighbors_D3N27[i][0], y + ZD::neighbors_D3N27[i][1], z + ZD::neighbors_D3N27[i][2]);
                }

                cs->SetValue(x - 1, y - 1, z - 1, ZD::CCSTool<DATA_TYPE>::FiberFunction(&buffer));
                SafeDeleteArray(buffer);
                //ftle->SetValue(x - 1, y - 1, z - 1, ZD::CPoint<DATA_TYPE, 1>(ZD::CFTLETool<DATA_TYPE, 3>::ComputeFTLE_D3N27(buffer, dt, dis)));
            }
        }
    }

    char pathname[ZD_PATHNAME_LENGTH];
    sprintf(pathname, "%sle_%s.nrrd", ZD::CFileTool::GetFilePath(pResultPathname).c_str(), "spatial_distance");
    cs->SaveNrrdFile(pathname);

    SafeDelete(cs);
}

void ComputeLagrangianEulerianByFlow()
{
    dx = (maxx - minx) / DATA_TYPE(nx - 1);
    dy = (maxy - miny) / DATA_TYPE(ny - 1);
    dz = (maxz - minz) / DATA_TYPE(nz - 1);

    const DATA_TYPE direction = endTime > startTime ? 1.0 : -1.0;
    const DATA_TYPE totalTime = std::abs(endTime - startTime);

    DATA_TYPE dt = (endTime - startTime) / DATA_TYPE(PATHLINE_MAX_COUNT);
    int maxCount = PATHLINE_MAX_COUNT;
    if (std::fabs(dt) < stepSize) {
        dt = dt > 0.0 ? stepSize : -stepSize;
        maxCount = int(std::fabs((endTime - startTime) / stepSize) + ZD_EPSILON);
    }

    for (int i = 0; i < scalarsCount; ++i) {
        ZD::CFlow<DATA_TYPE, 3> *pFlow = nullptr;

        if (strcmp(pFlowName, "delta_wing") == 0) {
            DATA_TYPE minTime = std::min(startTime, endTime);
            DATA_TYPE maxTime = std::max(startTime, endTime);
            pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, true, ppScalarNames[i]);
        }
        else {
            std::cerr << "Unknown 3D flow name: " << pFlowName << std::endl;
            return;
        }

        pFlow->SetStepSize(stepSize);

        ZD::CPoint<DATA_TYPE, K> *parameters = new ZD::CPoint<DATA_TYPE, K>[nx*ny*nz];

        printf("\n");
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
                for (int x = 0; x < nx; ++x) {
                    ZD::CPoint<DATA_TYPE, 3> seed;
                    seed[0] = minx + dx * (DATA_TYPE)x;
                    seed[1] = miny + dy * (DATA_TYPE)y;
                    seed[2] = minz + dz * (DATA_TYPE)z;

                    std::vector<DATA_TYPE> buffer;
                    buffer.push_back(pFlow->Scalar(seed, startTime));
                    DATA_TYPE integratedTime = 0.0;
                    for (int j = 0; j < maxCount; ++j) {
                        DATA_TYPE st = startTime + dt * DATA_TYPE(j);
                        DATA_TYPE tt = pFlow->NextTime(seed, st, std::fabs(dt), direction);
                        integratedTime += tt;
                        buffer.push_back(pFlow->Scalar(seed, startTime + integratedTime * direction));
                    }
                    parameters[(z*ny+y)*nx+x] = ZD::CParameterizationTool<DATA_TYPE, K>::ParameterizeSingleScalar(buffer);
                }
                printf("\r%.02f%% finished", float(z*ny + y) / float(ny*nz) * 100.0);
                fflush(stdout);
            }
        }

        ComputeLagrangianEulerian(parameters, ppScalarNames[i]);

        SafeDelete(pFlow);
        SafeDeleteArray(parameters);
    }
}

void ReadConfigureFile()
{
    std::fstream fs;
    fs.open(pConfigurePathname);
    char temp[128];

    if (fs.is_open() == true) {
        /* flow name */
        pFlowName = new char[ZD_PATHNAME_LENGTH];
        fs >> temp >> temp >> temp;
        fs >> pFlowName;

        /* skip bounding box */
        fs >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;

        /* nx, ny, nz */
        fs >> temp >> temp >> nx >> temp >> ny >> temp >> nz;

        /* dx, dy, dz */
        fs >> temp >> temp >> dx >> temp;
        fs >> temp >> temp >> dy >> temp;
        fs >> temp >> temp >> dz;

        /* start time, end time, time step */
        fs >> temp >> temp >> startTime >> temp;
        fs >> temp >> temp >> endTime >> temp;
        fs >> temp >> temp >> stepSize;

        fs >> temp;

        /* pathline and time file pathnames */
        pPathlinePathnames = new char*[nz];
        pTimePathnames = new char*[nz];
        for (int z = 0; z < nz; ++z) {
            pPathlinePathnames[z] = new char[ZD_PATHNAME_LENGTH];
            pTimePathnames[z] = new char[ZD_PATHNAME_LENGTH];
            fs >> pPathlinePathnames[z] >> pTimePathnames[z];
        }

        fs.close();
    }
}



void ParameterizePathlinesMoment(const ZD::CFlow<DATA_TYPE, 3> *flow, const ZD::CLine<DATA_TYPE, 3> *lines,
    const ZD::CLineAttribute<DATA_TYPE, 1> *times, const int size, ZD::CPoint<DATA_TYPE, K> *params)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < size; ++i) {
        std::vector<DATA_TYPE> buffer;
        for (int v = 0; v < lines[i].m_count; ++v) {
            buffer.push_back(flow->Scalar(lines[i].m_pVertices[v], times[i].m_pAttributes[v][0]));
        }
        params[i] = ZD::CParameterizationTool<DATA_TYPE, K>::ParameterizeSingleScalar(buffer);
    }
}

void ParameterizePathlinesFixed(const ZD::CFlow<DATA_TYPE, 3> *flow, const ZD::CLine<DATA_TYPE, 3> *lines,
    const ZD::CLineAttribute<DATA_TYPE, 1> *times, const int size, ZD::CPoint<DATA_TYPE, num> *params)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < size; ++i) {
        ZD::CPoint<DATA_TYPE, 3> *vertexBuffer = new ZD::CPoint<DATA_TYPE, 3>[num];
        ZD::CPoint<DATA_TYPE, 1> *timeBuffer = new ZD::CPoint<DATA_TYPE, 1>[num];

        DATA_TYPE totalLength = 0.0;
        for (int j = 1; j < lines[i].m_count; ++j) {
            totalLength += (lines[i].m_pVertices[j] - lines[i].m_pVertices[j - 1]).Length();
        }

        DATA_TYPE ss = totalLength / (DATA_TYPE)(num - 1);
        vertexBuffer[0] = lines[i].m_pVertices[0];
        timeBuffer[0] = times[i].m_pAttributes[0];
        int count = 1;
        DATA_TYPE lastLength = 0.0;
        ZD::CPoint<DATA_TYPE, 3> lastPoint = vertexBuffer[0];
        ZD::CPoint<DATA_TYPE, 1> lastTime = timeBuffer[0];
        for (int j = 1; j < lines[i].m_count;) {
            DATA_TYPE length = (lines[i].m_pVertices[j] - lastPoint).Length();
            if (length > (stepSize - lastLength)) {
                DATA_TYPE factor = stepSize / length;
                vertexBuffer[count] = factor * lines[i].m_pVertices[j] + (1.0 - factor) * lastPoint;
                timeBuffer[count] = factor * times[i].m_pAttributes[j] + (1.0 - factor) * lastTime;
                lastLength = 0;
                count++;
            }
            else {
                lastLength += length;
                lastPoint = lines[i].m_pVertices[j];

            }
        }

        //lines[i].ParameterizeFixed(vertexBuffer, num);
        for (int v = 0; v < num; ++v) {
            params[i][v] = flow->Scalar(vertexBuffer[v], timeBuffer[v][0]);
        }
        SafeDeleteArray(vertexBuffer);
        SafeDeleteArray(timeBuffer);
    }
}

void ComputeLagrangianEulerianByPathlinesMoment()
{
    ReadConfigureFile();

    for (int i = 0; i < scalarsCount; ++i) {
        if (strcmp(ppScalarNames[i], "ftle") == 0 || strcmp(ppScalarNames[i], "FTLE") == 0) {
            ZD::CField3<DATA_TYPE, 3> *flowMap = new ZD::CField3<DATA_TYPE, 3>();
            flowMap->CreateField(nx, ny, nz, nullptr);
            for (int z = 0; z < nz; ++z) {
                int size = 0;
                ZD::CLine<DATA_TYPE, 3> *lines = nullptr;
                ReadLines(pPathlinePathnames[z], &lines, size);
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        flowMap->SetValue(x, y, z, lines[y*nx+x].m_pVertices[lines[y*nx+x].m_count-1]);
                    }
                }
                SafeDeleteArray(lines);
            }
            ComputeFTLE(flowMap);

            SafeDelete(flowMap);
        } else {
            ZD::CFlow<DATA_TYPE, 3> *pFlow = nullptr;
            if (strcmp(pFlowName, "delta_wing") == 0) {
                DATA_TYPE minTime = std::min(startTime, endTime);
                DATA_TYPE maxTime = std::max(startTime, endTime);
                pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, false, ppScalarNames[i]);
            }
            else {
                std::cerr << "Unknown 3D flow name: " << pFlowName << std::endl;
                return;
            }

            ZD::CPoint<DATA_TYPE, K> *parameters = new ZD::CPoint<DATA_TYPE, K>[nx*ny*nz];
            for (int z = 0; z < nz; ++z) {
                int size = 0;
                ZD::CLine<DATA_TYPE, 3> *lines = nullptr;
                ZD::CLineAttribute<DATA_TYPE, 1> *times = nullptr;
                ReadLines(pPathlinePathnames[z], &lines, size);
                ReadLineAttributes(pTimePathnames[z], &times, size);
                ParameterizePathlinesMoment(pFlow, lines, times, size, &(parameters[z*nx*ny]));
                SafeDeleteArray(lines);
                SafeDeleteArray(times);
            }

            ComputeLagrangianEulerian(parameters, ppScalarNames[i]);

            SafeDelete(pFlow);
            SafeDeleteArray(parameters);
        }
    }

    SafeDeleteArray(pFlowName);
    for (int z = 0; z < nz; ++z) {
        SafeDeleteArray(pPathlinePathnames[z]);
        SafeDeleteArray(pTimePathnames[z]);
    }
    SafeDeleteArray(pPathlinePathnames);
    SafeDeleteArray(pTimePathnames);
}


void ComputeLagrangianEulerianByPathlinesFixed()
{
    ReadConfigureFile();

    for (int i = 0; i < scalarsCount; ++i) {
        if (strcmp(ppScalarNames[i], "ftle") == 0 || strcmp(ppScalarNames[i], "FTLE") == 0) {
            ZD::CField3<DATA_TYPE, 3> *flowMap = new ZD::CField3<DATA_TYPE, 3>();
            flowMap->CreateField(nx, ny, nz, nullptr);
            for (int z = 0; z < nz; ++z) {
                int size = 0;
                ZD::CLine<DATA_TYPE, 3> *lines = nullptr;
                ReadLines(pPathlinePathnames[z], &lines, size);
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        flowMap->SetValue(x, y, z, lines[y*nx + x].m_pVertices[lines[y*nx + x].m_count - 1]);
                    }
                }
                SafeDeleteArray(lines);
            }
            ComputeFTLE(flowMap);

            SafeDelete(flowMap);
        }
        else if (strcmp(ZD::CStringTool::ToLower(ppScalarNames[i]).c_str(), "spatial_distance") == 0 ) {
            ZD::CField3<DATA_TYPE, 9> *parameters = new ZD::CField3<DATA_TYPE, 9>();
            parameters->CreateField(nx, ny, nz, nullptr);
            for (int z = 0; z < nz; ++z) {
                int size = 0;
                ZD::CLine<DATA_TYPE, 3> *lines = nullptr;
                ReadLines(pPathlinePathnames[z], &lines, size);
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        parameters->SetValue(x, y, z, lines[y*nx+x].Parameterization3());
                    }
                }
                SafeDeleteArray(lines);
            }
            ComputeSpatialDistance(parameters);

            SafeDelete(parameters);
        }
        else            {
            ZD::CFlow<DATA_TYPE, 3> *pFlow = nullptr;
            if (strcmp(pFlowName, "delta_wing") == 0) {
                DATA_TYPE minTime = std::min(startTime, endTime);
                DATA_TYPE maxTime = std::max(startTime, endTime);
                pFlow = new ZD::CFlowDeltaWing<DATA_TYPE>(minTime, maxTime, false, ppScalarNames[i]);
            }
            else {
                std::cerr << "Unknown 3D flow name: " << pFlowName << std::endl;
                return;
            }

            ZD::CPoint<DATA_TYPE, num> *parameters = new ZD::CPoint<DATA_TYPE, num>[nx*ny*nz];
            for (int z = 0; z < nz; ++z) {
                int size = 0;
                ZD::CLine<DATA_TYPE, 3> *lines = nullptr;
                ZD::CLineAttribute<DATA_TYPE, 1> *times = nullptr;
                ReadLines(pPathlinePathnames[z], &lines, size);
                ReadLineAttributes(pTimePathnames[z], &times, size);
                ParameterizePathlinesFixed(pFlow, lines, times, size, &(parameters[z*nx*ny]));
                SafeDeleteArray(lines);
                SafeDeleteArray(times);
            }

            ComputeLagrangianEulerianFixed(parameters, ppScalarNames[i]);

            SafeDelete(pFlow);
            SafeDeleteArray(parameters);
        }
    }

    SafeDeleteArray(pFlowName);
    for (int z = 0; z < nz; ++z) {
        SafeDeleteArray(pPathlinePathnames[z]);
        SafeDeleteArray(pTimePathnames[z]);
    }
    SafeDeleteArray(pPathlinePathnames);
    SafeDeleteArray(pTimePathnames);
}


int main(int argc, char *argv[])
{
    Initialize(argc, argv);

#if _OPENMP
    omp_set_num_threads(OMP_MAX_THREADS);
#endif

    if (strcmp(pConfigurePathname, "NULL") == 0) {
        ComputeLagrangianEulerianByFlow();
    }
    else {
        if (method == ZD::ZD_LE_Method::ZD_LE_Method_Moment)
            ComputeLagrangianEulerianByPathlinesMoment();
        else if (method == ZD::ZD_LE_Method::ZD_LE_Method_Fixed)
            ComputeLagrangianEulerianByPathlinesFixed();
    }

    return EXIT_SUCCESS;
}
