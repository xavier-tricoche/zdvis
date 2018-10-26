#*************************************************************************
#zdvis: Lagrangian Visualization for Vector, Tensor, and Multifield Data.
#
#Author: Zi'ang Ding
#
#Copyright (c) 2016-2018, Purdue University
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.
#**************************************************************************
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

add_project(FlowMap2D  FlowMap2D.cpp)
add_project(FlowMap3D  FlowMap3D.cpp)
add_project(FTLE2D     FTLE2D.cpp)
add_project(FTLE3D     FTLE3D.cpp)
add_project(Pathline3D Pathline3D.cpp)
add_project(Pathline2D Pathline2D.cpp)
add_project(FSLE2D     FSLE2D.cpp)

