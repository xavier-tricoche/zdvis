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
/*******************************************************/
/* header file for all zdvis utility files             */
/* By: Zi'ang Ding                                     */
/* 2016-4-13                                           */
/*******************************************************/

#ifndef _ZD_LIB_HPP_
#define _ZD_LIB_HPP_

#define ZD_LIB_VERSION 0.3

#include "define.hpp"

/* base */
#include "Base/ZD_DynamicSystem.hpp"
#include "Base/ZD_Field.hpp"
#include "Base/ZD_Line.hpp"
#include "Base/ZD_LineAttribute.hpp"
#include "Base/ZD_Mesh.hpp"
#include "Base/ZD_Point.hpp"

/* integrator */
#include "Integrator/ZD_Integrator.hpp"
#include "Integrator/ZD_Integrator_Euler.hpp"
#include "Integrator/ZD_Integrator_RK4.hpp"
#include "Integrator/ZD_Integrator_RK45.hpp"
#endif

/* HamiltonianSystem */
#include "HamiltonianSystem/ZD_CRTBP.hpp"
#include "HamiltonianSystem/ZD_HamiltonianSystem.hpp"
#include "HamiltonianSystem/ZD_KinematicChaosModel.hpp"
#include "HamiltonianSystem/ZD_StandardMap.hpp"

/* DTI & HARDI */
#include "DTI_HARDI/ZD_DTI.hpp"
#include "DTI_HARDI/ZD_DWI.hpp"
#include "DTI_HARDI/ZD_Fiber.hpp"
#include "DTI_HARDI/ZD_Glyph.hpp"
#include "DTI_HARDI/ZD_HOT.hpp"
#include "DTI_HARDI/ZD_Phantom.hpp"
#include "DTI_HARDI/ZD_SH.hpp"

/* Tensor */
#include "Tensor/ZD_Tensor.hpp"
#include "Tensor/ZD_DoublePointLoad.hpp"

/* Flow */
#include "Flow/ZD_Flow.hpp"
#include "Flow/ZD_Flow_ABC.hpp"
#include "Flow/ZD_Flow_Boussinesq.hpp"
#include "Flow/ZD_Flow_Convection.hpp"
#include "Flow/ZD_Flow_Delta_Wing.hpp"
#include "Flow/ZD_Flow_Double_Gyre.hpp"
#include "Flow/ZD_Flow_Gaussian_Vortices.hpp"
#include "Flow/ZD_Flow_Meandering_Jet.hpp"
#include "Flow/ZD_Flow_Steady_Double_Gyre.hpp"
#include "Flow/ZD_Flow_Test.hpp"

/* Tool */
#include "Tool/ZD_ColorSpaceTool.hpp"
#include "Tool/ZD_CSTool.hpp"
#include "Tool/ZD_EigenTool.hpp"
#include "Tool/ZD_FileTool.hpp"
#include "Tool/ZD_FTLETool.hpp"
#include "Tool/ZD_GradientTool.hpp"
#include "Tool/ZD_ImageTool.hpp"
#include "Tool/ZD_LETool.hpp"
#include "Tool/ZD_ParameterizationTool.hpp"
#include "Tool/ZD_RandomTool.hpp"
#include "Tool/ZD_StringTool.hpp"
#include "Tool/ZD_TimeTool.hpp"

/* color map */
#include "ColorMap/ZD_ColorMap.hpp"
