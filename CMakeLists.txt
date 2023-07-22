set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring ConstitutiveLawsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR}/applications/structural_application )

## generate variables with the sources
set( CONSTITUTIVE_LAWS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/constitutive_laws_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/umat.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/mises_umat.f
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/umat_hypo.f
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/umat_fortran_wrapper.f
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/xit.f
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/udsm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/udsme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/umat2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/umat3.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/umat3e.cpp
)

if(DEFINED UMAT_LIBRARY)
    add_definitions(-DKRATOS_UMAT_LIBRARY_IS_PROVIDED)
endif()
if(DEFINED UDSM_LIBRARY)
    add_definitions(-DKRATOS_UDSM_LIBRARY_IS_PROVIDED)
endif()

if(DEFINED OPENSEES_DIR)
    include_directories( ${OPENSEES_DIR} )
    include_directories( ${OPENSEES_DIR}/SRC )
    include_directories( ${OPENSEES_DIR}/SRC/material )
    include_directories( ${OPENSEES_DIR}/SRC/material/nD )
    include_directories( ${OPENSEES_DIR}/SRC/material/nD/UWmaterials )
    include_directories( ${OPENSEES_DIR}/SRC/matrix )
    include_directories( ${OPENSEES_DIR}/SRC/domain/component )
    include_directories( ${OPENSEES_DIR}/SRC/tagged )
    include_directories( ${OPENSEES_DIR}/SRC/handler )
    include_directories( ${OPENSEES_DIR}/SRC/actor/actor )
    include_directories( ${OPENSEES_DIR}/SRC/recorder/response )
    include_directories( ${OPENSEES_DIR}/SRC/element )
    find_library( OPENSEES_LIBRARY NAMES OpenSees PATHS ${OPENSEES_DIR}/lib NO_DEFAULT_PATH )
    find_library( OPENSEES_BLAS_LIBRARY NAMES Blas PATHS ${OPENSEES_DIR}/lib NO_DEFAULT_PATH )
    find_library( OPENSEES_LAPACK_LIBRARY NAMES Lapack PATHS ${OPENSEES_DIR}/lib NO_DEFAULT_PATH )
    set(OPENSEES_LIBRARIES ${OPENSEES_LIBRARY} ${OPENSEES_LAPACK_LIBRARY} ${OPENSEES_BLAS_LIBRARY})
    add_definitions( -DKRATOS_USE_OPENSEES )
    message("OpenSees directory:" ${OPENSEES_DIR})
    message("OpenSees libraries:" ${OPENSEES_LIBRARIES})
    set(CONSTITUTIVE_LAWS_APPLICATION_SOURCES
        ${CONSTITUTIVE_LAWS_APPLICATION_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/opensees_mat.cpp
    )
endif()

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosConstitutiveLawsApplication SHARED ${CONSTITUTIVE_LAWS_APPLICATION_SOURCES} )
target_link_libraries( KratosConstitutiveLawsApplication KratosCore )
target_link_libraries( KratosConstitutiveLawsApplication KratosStructuralApplication )
target_link_libraries( KratosConstitutiveLawsApplication gfortran )
target_link_libraries( KratosConstitutiveLawsApplication dl )
target_link_libraries( KratosConstitutiveLawsApplication ${OPENSEES_LIBRARIES} )
if(DEFINED UMAT_LIBRARY)
    target_link_libraries( KratosConstitutiveLawsApplication ${UMAT_LIBRARY} )
endif()
if(DEFINED UDSM_LIBRARY)
    target_link_libraries( KratosConstitutiveLawsApplication ${UDSM_LIBRARY} )
endif()
set_target_properties(KratosConstitutiveLawsApplication PROPERTIES PREFIX "")
install(TARGETS KratosConstitutiveLawsApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveLawsApplication.py" DESTINATION KratosMultiphysics )
