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
#define PROJECT_DATASET_PATH        "homes/xmt/github/zdvis/data"
/* boussinesq path */
#define BOUSSINESQ_TIME_PATHNAME                     "E:/Ziang/project/dataset/boussinesq/times.txt"
#define BOUSSINESQ_DATA_PATH                         "E:/Ziang/project/dataset/boussinesq/data/"
#define BOUSSINESQ_VECTOR_PATH                         "E:/Ziang/project/dataset/boussinesq/resample/"        /* uniform resampled vecotr fields */
#define BOUSSINESQ_TOTAL_COUNT                         1001
/* convection path */
#define CONVECTION_TIME_PATHNAME                     "E:/Ziang/project/dataset/convection/time.txt"
#define CONVECTION_VECTOR_PATH                       "E:/Ziang/project/dataset/convection/velocity/"
#define CONVECTION_TEMPERATURE_PATH                  "E:/Ziang/project/dataset/convection/temperature/"
#define CONVECTION_TOTAL_COUNT                       3001
/* cylinder path */
#define CYLINDER_TOTAL_COUNT                         1201
#define CYLINDER_TIME_PATHNAME                         "E:/Ziang/project/dataset/cylinder/times.txt"
#define CYLINDER_DATA_PATH                             "E:/Ziang/project/dataset/cylinder/data/"
/* cylinder 2D */
#define CYLINDER2D_TOTAL_COUNT                       1001
#define CYLINDER2D_VECTOR_PATH                       "E:/Ziang/project/dataset/cylinder_2D/vector/"
/* delta wing path */
#define DELTA_WING_DLR_PATH                             "E:/Ziang/project/dataset/delta_wing/dlr/"
#define DELTA_WING_TOTAL_COUNT                         86
/* hurricane Isabel path */
#define HURRICANE_ISABEL_HEIGHT_MAP_PATHNAME         "E:/Ziang/project/dataset/vis2004/height_map/height_map.nrrd"
#define HURRICANE_ISABEL_VECTOR_PATH                 "E:/Ziang/project/dataset/vis2004/velocity/"
#define HURRICANE_ISABEL_PRESSURE_PATH               "E:/Ziang/project/dataset/vis2004/scalar/p/"
#define HURRICANE_ISABEL_TEMPERATURE_PATH            "E:/Ziang/project/dataset/Vis2004/TC/"
#define HURRICANE_ISABEL_CLOUD_PATH                  "E:/Ziang/project/dataset/Vis2004/CLOUD/"
#define HURRICANE_ISABEL_PRECIP_PATH                 "E:/Ziang/project/dataset/Vis2004/PRECIP/"
#define HURRICANE_ISABEL_SNOW_PATH                   "E:/Ziang/project/dataset/Vis2004/QSNOW/"
#define HURRICANE_ISABEL_TOTAL_COUNT                 48
/* shedder tensor */
#define SHEDDER_TENSOR_PATHNAME                         "E:/Ziang/project/dataset/shedder/tensor/strain_1108.nrrd"
/* phantom */
#define PHANTOM_DTI_PATHNAME                         "E:/Ziang/project/dataset/phantom/dti_30/dti.nrrd"
#define PHANTOM_HOT_PATHNAME                         "E:/Ziang/project/dataset/phantom/hardi_30/nrrd_denoise/hot6.nrrd"
#define PHANTOM_SH_PATHNAME                             "E:/Ziang/project/dataset/phantom/hardi_30/nrrd_denoise/sh8.nrrd"
#define PHANTOM_SH_MAX_VALUE_PATHNAME                 "E:/Ziang/project/dataset/phantom/hardi_30/nrrd_denoise/SHMaxValue.nrrd"
#define PHANTOM_MASK_PATHNAME                         "E:/Ziang/project/dataset/phantom/hardi_30/nrrd_denoise/mask.nrrd"
#define PHANTOM_DTI_N                                7
#define PHANTOM_HOT_N                                28
#define PHANTOM_SH_N                                 45
#define PHANTOM_HOT_ORDER                            6
#define PHANTOM_SH_ORDER                             8
#define PHANTOM_GROUND_TRUTH_PATH                     "E:/Ziang/project/dataset/phantom/ground_truth"
#define    PHANTOM_GROUND_TRUTH_COUNT                     20
#define PHANTOM_GROUND_TRUTH_VERTEX_COUNT             200
/* ionization front path */
#define IONIZATION_FRONT_TIME_PATHNAME               "E:/Ziang/project/dataset/vis2008/time/time.txt"
#define IONIZATION_FRONT_VECTOR_PATH                 "E:/Ziang/project/dataset/vis2008/velocity/"
#define IONIZATION_FRONT_MULTIFIELD_PATH             "E:/Ziang/project/dataset/vis2008/multifield/"
#define IONIZATION_FRONT_TOTAL_PARTICLE_DENSITY_PATH "E:/Ziang/project/dataset/vis2008/multifield/total_particle_density/"
#define IONIZATION_FRONT_GAS_TEMPERATURE_PATH        "E:/Ziang/project/dataset/vis2008/multifield/gas_temperature/"
#define IONIZATION_FRONT_H_2_MASS_PATH               "E:/Ziang/project/dataset/vis2008/multifield/H_2_mass/"
#define IONIZATION_FRONT_H_2P_MASS_PATH              "E:/Ziang/project/dataset/vis2008/multifield/H_2+_mass/"
#define IONIZATION_FRONT_H_MASS_PATH                 "E:/Ziang/project/dataset/vis2008/multifield/H_mass/"
#define IONIZATION_FRONT_HM_MASS_PATH                "E:/Ziang/project/dataset/vis2008/multifield/H-_mass/"
#define IONIZATION_FRONT_HP_MASS_PATH                "E:/Ziang/project/dataset/vis2008/multifield/H+_mass/"
#define IONIZATION_FRONT_HE_MASS_PATH                "E:/Ziang/project/dataset/vis2008/multifield/He_mass/"
#define IONIZATION_FRONT_HEP_MASS_PATH               "E:/Ziang/project/dataset/vis2008/multifield/He+_mass/"
#define IONIZATION_FRONT_HEPP_MASS_PATH              "E:/Ziang/project/dataset/vis2008/multifield/He++_mass/"
#define IONIZATION_FRONT_TOTAL_COUNT                 200
/* gaussian vortices */
#define GAUSSIAN_VORTICES_TOTAL_COUNT                 1001
#define GAUSSIAN_VORTICES_DATA_PATH                  "E:/Ziang/project/dataset/gaussian_vortices/data/"
#define GAUSSIAN_VORTICES_TIME_PATHNAME                 "E:/Ziang/project/dataset/gaussian_vortices/time/times.txt"
/* ocean */
#define OCEAN_TOTAL_COUNT                             1224
#define OCEAN_VECTOR_PATH                             "E:/Ziang/project/dataset/NCOM/nrrd/"

#endif
