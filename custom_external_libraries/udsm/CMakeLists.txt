# CMakeLists.txt for inclusion of Plaxis soil models to Kratos
# Author: hbui
# Date: Apr 27, 2015
# Remarks:
#   +   to add new soil model, simply append to PLAXIS_SOIL_MODEL_SRC
#   +   to link with new model in static/shared library, use target_link_libraries

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set( PLAXIS_SOIL_MODEL_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/usrlib.for
    ${CMAKE_CURRENT_SOURCE_DIR}/mc_tens.for
    ${CMAKE_CURRENT_SOURCE_DIR}/PlDP.for
    ${CMAKE_CURRENT_SOURCE_DIR}/PlModMC.for
    ${CMAKE_CURRENT_SOURCE_DIR}/UdsmDP.for
    ${CMAKE_CURRENT_SOURCE_DIR}/UdsmMC.for
)

###############################################################
add_library(PlaxisUDSM STATIC ${PLAXIS_SOIL_MODEL_SRC} )

