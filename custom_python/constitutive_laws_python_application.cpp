//
//   Project Name:        Kratos
//   Last modified by:    $Author: stasch $
//   Date:                $Date: 2007-09-26 13:57:49 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"
#include "custom_python/add_constitutive_laws_to_python.h"

namespace Kratos
{

namespace Python
{

BOOST_PYTHON_MODULE( KratosConstitutiveLawsApplication )
{

    using namespace boost::python;

    class_ < KratosConstitutiveLawsApplication,
           KratosConstitutiveLawsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable > ( "KratosConstitutiveLawsApplication" )
           ;

    AddConstitutiveLawsToPython();

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SOIL_MODEL_NUMBER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_UNDRAINED )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BULK_W )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLAXIS_LIBRARY_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( USERMOD_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ABAQUS_LIBRARY_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_NDI )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_NSHR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_NSTATV )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_CMNAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UMAT_STATEV )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( OPENSEES_MATERIAL_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UDSM_RELAXATION_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_AXISYMMETRIC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEED_DETERMINE_INTERNAL_PARAMS )

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
