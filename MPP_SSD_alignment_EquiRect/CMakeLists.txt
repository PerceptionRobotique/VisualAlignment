#############################################################################
#
# This file is part of the libPR software.
# Copyright (C) 2021 by MIS lab (UPJV) & JRL (CNRS-AIST). All rights reserved.
#
# See http://mis.u-picardie.fr/~g-caron/fr/index.php?page=7 for more information.
#
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# configuration file. Freely inspired from the CMakeList tree of the 
# ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################

cmake_minimum_required(VERSION 4.0)

project(MPP_SSD_alignment_EquiRect)

find_package(PER REQUIRED per_core per_core_extended per_io per_features per_sensor_pose_estimation per_estimation)
if(PER_FOUND)
	include(${PER_USE_FILE})
endif(PER_FOUND)

# ViSP (to do: list modules only)
find_package(VISP REQUIRED)
if(VISP_FOUND)
	include(${VISP_USE_FILE})
endif(VISP_FOUND)


# Boost
FIND_PACKAGE(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIRS}) # /opt/local/include/ might be needed as well under MacOS

link_directories(${Boost_LIBRARY_DIRS}) # /opt/local/lib/ might be needed as well under MacOS # similar /Users/guillaume/Developpement/librairies/visp-3.0.1/build/lib/Release/ might be needed as well under MacOS 

set(MPPSSDalignment_EquiRect_cpp
  MPPSSDalignment_EquiRect.cpp
)

foreach(cpp ${MPPSSDalignment_EquiRect_cpp})
  per_add_target(${cpp})
  if(COMMAND pr_add_dependency)
    pr_add_dependency(${cpp} "MPPSSDalignment_EquiRect")
  endif()
endforeach()

#target_link_libraries(MPPSSDgyroEstim libboost_system-mt.dylib libboost_filesystem-mt.dylib libboost_regex-mt.dylib)

#target_link_libraries(MPPSSDgyroEstim libboost_system-mt.so libboost_filesystem-mt.so libboost_regex-mt.so)

target_link_libraries(MPPSSDalignment_EquiRect libboost_system.so libboost_filesystem.so libboost_regex.so)
