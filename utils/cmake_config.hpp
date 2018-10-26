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
#define OMP_MAX_THREADS            12

/* pre-allocate pathline vertices */
#define PATHLINE_MAX_COUNT         500


/* dataset pathnames */
#define PROJECT_DATASET_PATH                         "/homes/xmt/visdata/Public/zdvis/datasets"
/* boussinesq path */
#define BOUSSINESQ_TIME_PATHNAME                     "/homes/xmt/visdata/Public/zdvis/datasets/boussinesq/times.txt"
#define BOUSSINESQ_DATA_PATH                         "/homes/xmt/visdata/Public/zdvis/datasets/boussinesq/data/" /* vector fields */
#define BOUSSINESQ_TOTAL_COUNT                       201
/* convection path */
#define CONVECTION_TIME_PATHNAME                     "/homes/xmt/visdata/Public/zdvis/datasets/convection/list.txt"
#define CONVECTION_VECTOR_PATH                       "/homes/xmt/visdata/Public/zdvis/datasets/convection/velocity"
#define CONVECTION_TEMPERATURE_PATH                  "/homes/xmt/visdata/Public/zdvis/datasets/convection/temperature/"
#define CONVECTION_TOTAL_COUNT                       201
/* delta wing path */
#define DELTA_WING_DLR_FILE                          "/homes/xmt/visdata/Public/zdvis/datasets/delta_wing/dlr/deltawing.dlr"
#define DELTA_WING_TOTAL_COUNT                       6
/* shedder tensor */
#define SHEDDER_TENSOR_PATHNAME                      "/homes/xmt/visdata/Public/zdvis/datasets/shedder/tensor/strain_1108.nrrd"
/* phantom */
#define PHANTOM_DTI_PATHNAME                         "/homes/xmt/visdata/Public/zdvis/datasets/phantom/dti_30/dti.nrrd"
#define PHANTOM_HOT_PATHNAME                         "/homes/xmt/visdata/Public/zdvis/datasets/phantom/hardi_30/nrrd_denoise/hot6.nrrd"
#define PHANTOM_SH_PATHNAME                          "/homes/xmt/visdata/Public/zdvis/datasets/phantom/hardi_30/nrrd_denoise/sh8.nrrd"
#define PHANTOM_SH_MAX_VALUE_PATHNAME                 "/homes/xmt/visdata/Public/zdvis/datasets/phantom/hardi_30/nrrd_denoise/SHMaxValue.nrrd"
#define PHANTOM_MASK_PATHNAME                         "/homes/xmt/visdata/Public/zdvis/datasets/phantom/hardi_30/nrrd_denoise/mask.nrrd"
#define PHANTOM_DTI_N                                7
#define PHANTOM_HOT_N                                28
#define PHANTOM_SH_N                                 45
#define PHANTOM_HOT_ORDER                            6
#define PHANTOM_SH_ORDER                             8
#define PHANTOM_GROUND_TRUTH_PATH                    "/homes/xmt/visdata/Public/zdvis/datasets/phantom/ground_truth"
#define PHANTOM_GROUND_TRUTH_COUNT                   20
#define PHANTOM_GROUND_TRUTH_VERTEX_COUNT            200
/* gaussian vortices */
#define GAUSSIAN_VORTICES_TOTAL_COUNT                1001
#define GAUSSIAN_VORTICES_DATA_PATH                  "/homes/xmt/visdata/Public/zdvis/datasets/gaussian_vortices/data/"
#define GAUSSIAN_VORTICES_TIME_PATHNAME              "/homes/xmt/visdata/Public/zdvis/datasets/gaussian_vortices/time/times.txt"

#endif
