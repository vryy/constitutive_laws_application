message("")
message("**** configuring KratosConstitutiveLawsApplication ****")

if(${CMAKE_SYSTEM_NAME} MATCHES "Android")
    message("KratosConstitutiveLawsApplication is compiled on Android")
    include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists2.txt)
else()
    if((${PYTHON_VERSION_MAJOR} MATCHES "2") OR (${KRATOS_USE_BOOST_PYTHON} MATCHES ON))
        if((${PYTHON_VERSION_MAJOR} MATCHES "3") AND (${KRATOS_USE_BOOST_PYTHON} MATCHES ON))
            message("KratosConstitutiveLawsApplication is compiled with Python 3 (Boost.Python)")
        else()
            message("KratosConstitutiveLawsApplication is compiled with Python 2")
        endif()
        include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists2.txt)
    elseif(${PYTHON_VERSION_MAJOR} MATCHES "3")
        message("KratosConstitutiveLawsApplication is not compilable with Python 3")
    else()
        message("KratosConstitutiveLawsApplication does not detect the Python version")
    endif()
endif()

message("**** finish configuring KratosConstitutiveLawsApplication ****")
message("")
