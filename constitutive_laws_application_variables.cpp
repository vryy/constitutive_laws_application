//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:54:46 $
//   Revision:            $Revision: 1.3 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(std::string, PLAXIS_LIBRARY_NAME)
KRATOS_CREATE_VARIABLE(std::string, USERMOD_NAME)
KRATOS_CREATE_VARIABLE(int, SOIL_MODEL_NUMBER)
KRATOS_CREATE_VARIABLE(bool, IS_UNDRAINED)
KRATOS_CREATE_VARIABLE(std::string, ABAQUS_LIBRARY_NAME)
KRATOS_CREATE_VARIABLE(std::string, UMAT_NAME)
KRATOS_CREATE_VARIABLE(int, UMAT_NDI)
KRATOS_CREATE_VARIABLE(int, UMAT_NSHR)
KRATOS_CREATE_VARIABLE(int, UMAT_NSTATV)
KRATOS_CREATE_VARIABLE(std::string, UMAT_CMNAME)
KRATOS_CREATE_VARIABLE(Vector, UMAT_STATEV)
KRATOS_CREATE_VARIABLE(std::string, OPENSEES_MATERIAL_NAME)
KRATOS_CREATE_VARIABLE(double, UDSM_RELAXATION_FACTOR)

}  // namespace Kratos.


