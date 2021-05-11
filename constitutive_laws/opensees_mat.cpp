#include <iostream>
#include "opensees_mat.h"
#include "constitutive_laws_application.h"
#include "structural_application/structural_application_variables.h"

#include "SRC/material/nD/UWmaterials/DruckerPrager3D.h"
#include "SRC/material/nD/UWmaterials/BoundingCamClay3D.h"

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

extern "C" int OPS_GetNumRemainingInputArgs()
{
    return 0;
}

extern "C" int ops_getintinput_(int *numData, int* data)
{
  return 0;
}

extern "C" int ops_getdoubleinput_(int *numData, double* data)
{
  return 0;
}

extern "C" const char* ops_getstring(void)
{
    const char *res = 0;
    return res;
}

void OpenSees_setTrialStrain(boost::shared_ptr<NDMaterial> pMat, const Kratos::Vector& rValue)
{
    Vector strain_from_element(6);
    for(int i = 0; i < 6; ++i)
        strain_from_element(i) = rValue(i);

    pMat->setTrialStrain(strain_from_element);
}

void OpenSees_getStress(boost::shared_ptr<NDMaterial> pMat, Kratos::Vector& rValue)
{
    const Vector& stress = pMat->getStress();

    for(int i = 0; i < 6; ++i)
        rValue(i) = stress(i);
}

void OpenSees_getTangent(boost::shared_ptr<NDMaterial> pMat, Kratos::Matrix& rValue)
{
    const Matrix& tangent = pMat->getTangent();

    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            rValue(i, j) = tangent(i, j);
}

void OpenSees_getInitialTangent(boost::shared_ptr<NDMaterial> pMat, Kratos::Matrix& rValue)
{
    const Matrix& tangent = pMat->getInitialTangent();

    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            rValue(i, j) = tangent(i, j);
}

namespace Kratos
{

OpenSeesMat::OpenSeesMat()
{
}

OpenSeesMat::~OpenSeesMat()
{
}

int OpenSeesMat::Check( const Kratos::Properties& props, const Kratos::ConstitutiveLaw::GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo )
{
    if(!props.Has(MATERIAL_PARAMETERS))
        KRATOS_THROW_ERROR(std::logic_error, "MATERIAL_PARAMETERS is not given to the properties", "")
    return 0;
}

bool OpenSeesMat::Has( const Kratos::Variable<int>& rThisVariable )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        return true;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        return true;
    return false;
}

bool OpenSeesMat::Has( const Kratos::Variable<double>& rThisVariable )
{
    return false;
}

bool OpenSeesMat::Has( const Kratos::Variable<Vector>& rThisVariable )
{
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
        return true;
    if(rThisVariable == STRAIN)
        return true;
    return false;
}

bool OpenSeesMat::Has( const Kratos::Variable<Matrix>& rThisVariable )
{
    return false;
}

int& OpenSeesMat::GetValue( const Kratos::Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& OpenSeesMat::GetValue( const Kratos::Variable<double>& rThisVariable, double& rValue )
{
    return rValue;
}

Vector& OpenSeesMat::GetValue( const Kratos::Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable == INTERNAL_VARIABLES)
    {
        if(rValue.size() != mCurrentStateVariables.size())
            rValue.resize(mCurrentStateVariables.size(), false);
        noalias(rValue) = mCurrentStateVariables;
    }
    if(rThisVariable == STRESSES)
    {
        if(rValue.size() != mCurrentStress.size())
            rValue.resize(mCurrentStress.size(), false);
        noalias(rValue) = mCurrentStress;
    }
    if(rThisVariable == STRAIN)
    {
        if(rValue.size() != mCurrentStrain.size())
            rValue.resize(mCurrentStrain.size(), false);
        noalias(rValue) = mCurrentStrain;
    }

    return rValue;
}

Matrix& OpenSeesMat::GetValue( const Kratos::Variable<Matrix>& rThisVariable, Kratos::Matrix& rValue )
{
    return rValue;
}

void OpenSeesMat::SetValue( const Kratos::Variable<int>& rVariable, const int& rValue, const Kratos::ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PARENT_ELEMENT_ID)
        mParentElementId = rValue;
    if(rVariable == INTEGRATION_POINT_INDEX)
        mIntegrationPointId = rValue;
}

void OpenSeesMat::SetValue( const Kratos::Variable<double>& rVariable, const double& rValue, const Kratos::ProcessInfo& rCurrentProcessInfo)
{
}

void OpenSeesMat::SetValue( const Kratos::Variable<Vector>& rVariable, const Kratos::Vector& rValue, const Kratos::ProcessInfo& rCurrentProcessInfo)
{
}

void OpenSeesMat::SetValue( const Kratos::Variable<Matrix>& rVariable, const Kratos::Matrix& rValue, const Kratos::ProcessInfo& rCurrentProcessInfo)
{
}

void OpenSeesMat::InitializeMaterial( const Kratos::Properties& props,
                               const Kratos::ConstitutiveLaw::GeometryType& geom,
                               const Kratos::Vector& ShapeFunctionsValues )
{
    const Kratos::Vector& mat_params = props[MATERIAL_PARAMETERS];

    // be careful here, this interface now only supports the 3D constitutive law
    if(props[OPENSEES_MATERIAL_NAME] == std::string("DruckerPrager3D"))
    {
        // REF: http://opensees.berkeley.edu/wiki/index.php/Drucker_Prager
        // REMARKS:
        //  +   the yield function in the wiki page is wrong, the sign of q_iso is minus (code checked)
        //  +   the f_2 yield surface equation may be wrong, it should be f2 = -Invariant_1 - T(mAlpha2_n1) <= 0
        double bulk = mat_params(0);
        double shear = mat_params(1);
		double s_y = mat_params(2);
		double r = mat_params(3);
		double r_bar = mat_params(4);
		double Kinfinity = mat_params(5);
		double Kinit = mat_params(6); 
		double d1 = mat_params(7);
		double d2 = mat_params(8);
		double H = mat_params(9);
		double t = mat_params(10);
		double massDen = mat_params(11);
		double atm = mat_params(12);
        mpConstitutiveLaw = boost::shared_ptr<NDMaterial>(new DruckerPrager3D(0, bulk, shear, s_y, r, r_bar, Kinfinity, Kinit, d1, d2, H, t, massDen, atm));
    }
    else if(props[OPENSEES_MATERIAL_NAME] == std::string("BoundingCamClay3D"))
    {
        double mDen = mat_params(0);
        double C = mat_params(1);
        double bulk = mat_params(2);
        double OCR = mat_params(3);
        double mu_o = mat_params(4);
        double alpha = mat_params(5);
        double lambda = mat_params(6);
        double h = mat_params(7);
        double m = mat_params(8);
        mpConstitutiveLaw = boost::shared_ptr<NDMaterial>(new BoundingCamClay3D(0, mDen, C, bulk, OCR, mu_o, alpha, lambda, h, m));
    }
    else
        KRATOS_THROW_ERROR(std::runtime_error, "The material is not supported:", props[OPENSEES_MATERIAL_NAME])

    mpChannel = boost::shared_ptr<MyOpenSeesChannel>(new MyOpenSeesChannel());

    mCurrentStress.resize(6);
    mCurrentStrain.resize(6);
}

void OpenSeesMat::ResetMaterial( const Kratos::Properties& props,
                           const Kratos::ConstitutiveLaw::GeometryType& geom,
                           const Kratos::Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = ZeroVector(6);
    noalias(mCurrentStrain) = ZeroVector(6);
    mpConstitutiveLaw->revertToStart();
}

void OpenSeesMat::InitializeSolutionStep( const Kratos::Properties& props,
                                    const Kratos::ConstitutiveLaw::GeometryType& geom,
                                    const Kratos::Vector& ShapeFunctionsValues ,
                                    const Kratos::ProcessInfo& CurrentProcessInfo )
{
}

void OpenSeesMat::InitializeNonLinearIteration ( const Kratos::Properties& props,
                                           const Kratos::ConstitutiveLaw::GeometryType& geom,
                                           const Kratos::Vector& ShapeFunctionsValues,
                                           const Kratos::ProcessInfo& CurrentProcessInfo )
{
}

void OpenSeesMat::CalculateMaterialResponse( const Kratos::Vector& StrainVector,
                                       const Kratos::Matrix& DeformationGradient,
                                       Kratos::Vector& StressVector,
                                       Kratos::Matrix& AlgorithmicTangent,
                                       const Kratos::ProcessInfo& CurrentProcessInfo,
                                       const Kratos::Properties& props,
                                       const Kratos::ConstitutiveLaw::GeometryType& geom,
                                       const Kratos::Vector& ShapeFunctionsValues,
                                       bool CalculateStresses,
                                       int CalculateTangent,
                                       bool SaveInternalVariables )
{
    OpenSees_setTrialStrain(mpConstitutiveLaw, StrainVector);
    OpenSees_getStress(mpConstitutiveLaw, StressVector);
    OpenSees_getTangent(mpConstitutiveLaw, AlgorithmicTangent);
//    OpenSees_getInitialTangent(mpConstitutiveLaw, AlgorithmicTangent);

    if(mParentElementId == 1)
    {
        KRATOS_WATCH(StressVector)
        KRATOS_WATCH(AlgorithmicTangent)
    }

    mpConstitutiveLaw->sendSelf(0, *mpChannel);
    mpChannel->copyVector(mCurrentStateVariables);
}

void OpenSeesMat::FinalizeNonLinearIteration ( const Kratos::Properties& props,
                                         const Kratos::ConstitutiveLaw::GeometryType& geom,
                                         const Kratos::Vector& ShapeFunctionsValues,
                                         const Kratos::ProcessInfo& CurrentProcessInfo )
{
}

void OpenSeesMat::FinalizeSolutionStep( const Kratos::Properties& props,
                                 const Kratos::ConstitutiveLaw::GeometryType& geom,
                                 const Kratos::Vector& ShapeFunctionsValues ,
                                 const Kratos::ProcessInfo& CurrentProcessInfo )
{
    mpConstitutiveLaw->commitState();
}

} // Namespace Kratos
