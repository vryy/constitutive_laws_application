//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:54:44 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_CONSTITUTIVE_LAWS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_LAWS_APPLICATION_VARIABLES_H_INCLUDED



// System includes
#include <string>


// External includes


// Project includes
#include "includes/define.h"

namespace Kratos
{

KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, bool,   IS_UNDRAINED )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, double, BULK_W ) // we should use BULK_WATER, but it requires SoilMechanicsApplication. To make this application independent, BULK_W is used.
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, PLAXIS_LIBRARY_NAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, ABAQUS_LIBRARY_NAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, USERMOD_NAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, UMAT_NAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, int, SOIL_MODEL_NUMBER )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, int, UMAT_NDI )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, int, UMAT_NSHR )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, int, UMAT_NSTATV )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, UMAT_CMNAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, Vector, UMAT_STATEV )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, std::string, OPENSEES_MATERIAL_NAME )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, double, UDSM_RELAXATION_FACTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, bool,   IS_AXISYMMETRIC )
KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_LAWS_APPLICATION, bool,   NEED_DETERMINE_INTERNAL_PARAMS )

/**
 * type definition for Umat
 * @param STRESS ......... the vector of stresses
 * @param STATEV ......... the vector of state variables
 * @param DDSDDE ......... the material tangent
 * @param SSE ............
 * @param SPD ............
 * @param SCD ............
 * @param RPL ............
 * @param DDSDDT .........
 * @param DRPLDE .........
 * @param DRPLDT .........
 * @param STRAN .......... the vector of total strains
 * @param DSTRAN ......... the vector of incremental strains
 * @param TIME ........... current time
 * @param DTIME .......... current time increment
 * @param TEMP ........... current temperature
 * @param DTEMP .......... current increment of temperature
 * @param PREDEF .........
 * @param DPRED ..........
 * @param CMNAME ......... material name
 * @param NDI ............ number of direct strain components (3 in 3D)
 * @param NSHR ........... number if shear strain components (3 in 3D)
 * @param NTENS .......... number of stress components (6 in 3D)
 * @param NSTATV ......... number of state variables (size of STATEV)
 * @param PROPS .......... material parameters
 * @param NPROPS ......... number of material paramters (size of PROPS)
 * @param COORDS ......... coordinates of the integration point
 * @param DROT ........... rotation increment (3, 3)
 * @param PNEWDT .........
 * @param CELENT ......... characteristic element length
 * @param DFGRD0 .........
 * @param DFGRD1 .........
 * @param NOEL ........... element number
 * @param NPT ............ integration point number
 * @param KSLAY ..........
 * @param KSPT ...........
 * @param KSTEP .......... step number
 * @param KINC ........... increment number
 */
typedef void (*umat_t)(double* STRESS, double* STATEV, double* DDSDDE, double* SSE, double* SPD, double* SCD,
        double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
        double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
        char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
        double* COORDS, double* DROT, double* PNEWDT, double* CELENT, double* DFGRD0,
        double* DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC);

/**
 * type definition for UDSM
 * */
typedef void (*udsm_t)(int* IDTask, int* iMod, int* IsUndr, int* iStep, int* iTer, int* iEl,
        int* Int, double* X, double* Y, double* Z, double* Time0, double* dTime,
        double* Props, double* Sig0, double* Swp0, double* StVar0, double* dEps,
        double* D, double* BulkW, double* Sig, double* Swp, double* StVar,
        int* ipl, int* nStat, int* NonSym, int* iStrsDep, int* iTimeDep,
        int* iTang, int* iPrjDir, int* iPrjLen, int* iAbort);

}  // namespace Kratos.

#endif // KRATOS_EXTERNAL_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED  defined
