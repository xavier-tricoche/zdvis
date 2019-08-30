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


// Note: This file is automatically edited by CMake based on the provided
//       value of the variable PROJECT_DATASET_PATH

#ifndef _CMAKE_CONFIGURE_H_
#define _CMAKE_CONFIGURE_H_


#ifndef ON
#define ON 1
#endif

#ifndef OFF
#define OFF 0
#endif

#define ZD_ENABLE_CUDA   
#define ZD_ENABLE_OGLUI  

/* maximum openmp threads */
#define OMP_MAX_THREADS            8

/* pre-allocate pathline vertices */
#define PATHLINE_MAX_COUNT         500


/* dataset pathnames */
#define PROJECT_DATASET_PATH                         "/Users/xmt/data/zdvis"
/* boussinesq path */
#define BOUSSINESQ_TIME_PATHNAME                     "/Users/xmt/data/zdvis/boussinesq/times.txt"
#define BOUSSINESQ_DATA_PATH                         "/Users/xmt/data/zdvis/boussinesq/data/" /* vector fields */
#define BOUSSINESQ_TOTAL_COUNT                       201
/* convection path */
#define CONVECTION_TIME_PATHNAME                     "/Users/xmt/data/zdvis/convection/list.txt"
#define CONVECTION_VECTOR_PATH                       "/Users/xmt/data/zdvis/convection/velocity"
#define CONVECTION_TEMPERATURE_PATH                  "/Users/xmt/data/zdvis/convection/temperature/"
#define CONVECTION_TOTAL_COUNT                       201
/* delta wing path */
#define DELTA_WING_DLR_FILE                          "/Users/xmt/data/zdvis/delta_wing/dlr/deltawing.dlr"
#define DELTA_WING_TOTAL_COUNT                       6
/* shedder tensor */
#define SHEDDER_TENSOR_PATHNAME                      "/Users/xmt/data/zdvis/shedder/tensor/strain_1108.nrrd"
/* phantom */
#define PHANTOM_DTI_PATHNAME                         "/Users/xmt/data/zdvis/phantom/dti_30/dti.nrrd"
#define PHANTOM_HOT_PATHNAME                         "/Users/xmt/data/zdvis/phantom/hardi_30/nrrd_denoise/hot6.nrrd"
#define PHANTOM_SH_PATHNAME                          "/Users/xmt/data/zdvis/phantom/hardi_30/nrrd_denoise/sh8.nrrd"
#define PHANTOM_SH_MAX_VALUE_PATHNAME                 "/Users/xmt/data/zdvis/phantom/hardi_30/nrrd_denoise/SHMaxValue.nrrd"
#define PHANTOM_MASK_PATHNAME                         "/Users/xmt/data/zdvis/phantom/hardi_30/nrrd_denoise/mask.nrrd"
#define PHANTOM_DTI_N                                7
#define PHANTOM_HOT_N                                28
#define PHANTOM_SH_N                                 45
#define PHANTOM_HOT_ORDER                            6
#define PHANTOM_SH_ORDER                             8
#define PHANTOM_GROUND_TRUTH_PATH                    "/Users/xmt/data/zdvis/phantom/ground_truth"
#define PHANTOM_GROUND_TRUTH_COUNT                   20
#define PHANTOM_GROUND_TRUTH_VERTEX_COUNT            200
/* gaussian vortices */
#define GAUSSIAN_VORTICES_TOTAL_COUNT                1001
#define GAUSSIAN_VORTICES_DATA_PATH                  "/Users/xmt/data/zdvis/gaussian_vortices/data/"
#define GAUSSIAN_VORTICES_TIME_PATHNAME              "/Users/xmt/data/zdvis/gaussian_vortices/time/times.txt"

#endif
