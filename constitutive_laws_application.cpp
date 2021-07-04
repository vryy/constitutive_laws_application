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
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

KratosConstitutiveLawsApplication::KratosConstitutiveLawsApplication():
    mUmat3e()
{}

void KratosConstitutiveLawsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosConstitutiveLawsApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE(PLAXIS_LIBRARY_NAME)
    KRATOS_REGISTER_VARIABLE(USERMOD_NAME)
    KRATOS_REGISTER_VARIABLE(SOIL_MODEL_NUMBER)
    KRATOS_REGISTER_VARIABLE(IS_UNDRAINED)
    KRATOS_REGISTER_VARIABLE(ABAQUS_LIBRARY_NAME)
    KRATOS_REGISTER_VARIABLE(UMAT_NAME)
    KRATOS_REGISTER_VARIABLE(UMAT_NDI)
    KRATOS_REGISTER_VARIABLE(UMAT_NSHR)
    KRATOS_REGISTER_VARIABLE(UMAT_NSTATV)
    KRATOS_REGISTER_VARIABLE(UMAT_CMNAME)
    KRATOS_REGISTER_VARIABLE(UMAT_STATEV)
    KRATOS_REGISTER_VARIABLE(OPENSEES_MATERIAL_NAME)

    Serializer::Register( "Umat3e", mUmat3e );
}


}  // namespace Kratos.


