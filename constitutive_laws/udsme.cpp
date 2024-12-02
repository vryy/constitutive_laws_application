#include <iomanip>

#ifdef _MSC_VER
#include <Windows.h>
#else
#include <dlfcn.h>
#endif
#include "constitutive_laws/udsme.h"
#include "structural_application/structural_application_variables.h"

namespace Kratos
{

bool UDSMe::Has ( const Variable<int>& rThisVariable )
{
    if(rThisVariable == SOIL_MODEL_NUMBER)
        return true;
    return false;
}

bool UDSMe::Has ( const Variable<bool>& rThisVariable )
{
    if(rThisVariable == IS_UNDRAINED)
        return true;
    return false;
}

bool UDSMe::Has ( const Variable<double>& rThisVariable )
{
    if(rThisVariable == BULK_W)
        return true;
    return BaseType::Has( rThisVariable );
}

bool UDSMe::Has ( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == MATERIAL_PARAMETERS)
        return true;
    return BaseType::Has( rThisVariable );
}

// bool UDSMe::Has ( const Variable<std::string>& rThisVariable )
// {
//     if(rThisVariable == PLAXIS_LIBRARY_NAME)
//         return true;
//     if(rThisVariable == USERMOD_NAME)
//         return true;
//     return false;
// }

int& UDSMe::GetValue ( const Variable<int>& rThisVariable, int& rValue )
{
    if(rThisVariable == SOIL_MODEL_NUMBER)
        rValue = mModelNumber;
    return rValue;
}

bool& UDSMe::GetValue ( const Variable<bool>& rThisVariable, bool& rValue )
{
    if(rThisVariable == IS_UNDRAINED)
        rValue = mIsUndr;
    return rValue;
}

double& UDSMe::GetValue ( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == BULK_W)
        rValue = mBulkW;
    return BaseType::GetValue( rThisVariable, rValue );
}

Vector& UDSMe::GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable == MATERIAL_PARAMETERS)
    {
        if(rValue.size() != mProps.size())
            rValue.resize(mProps.size(), false);
        noalias(rValue) = mProps;
    }

    return BaseType::GetValue( rThisVariable, rValue );
}

std::string& UDSMe::GetValue ( const Variable<std::string>& rThisVariable, std::string& rValue )
{
    if(rThisVariable == PLAXIS_LIBRARY_NAME)
        rValue = mLibName;
    if(rThisVariable == USERMOD_NAME)
        rValue = mName;
    return rValue;
}

void UDSMe::SetValue ( const Variable<int>& rThisVariable,
                       const int& rValue,
                       const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == SOIL_MODEL_NUMBER)
        mModelNumber = rValue;
    else
        BaseType::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void UDSMe::SetValue ( const Variable<bool>& rThisVariable,
                       const bool& rValue,
                       const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == IS_UNDRAINED)
        mIsUndr = static_cast<int>(rValue);
}

void UDSMe::SetValue ( const Variable<double>& rThisVariable,
                       const double& rValue,
                       const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == BULK_W)
        mBulkW = rValue;
    else
        BaseType::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void UDSMe::SetValue ( const Variable<Vector>& rThisVariable,
                       const Vector& rValue,
                       const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == MATERIAL_PARAMETERS )
    {
        if(mProps.size() != rValue.size())
            mProps.resize(rValue.size(), false);
        noalias(mProps) = rValue;
    }
    else
        BaseType::SetValue( rThisVariable, rValue, rCurrentProcessInfo );
}

void UDSMe::SetValue ( const Variable<std::string>& rThisVariable,
                       const std::string& rValue,
                       const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PLAXIS_LIBRARY_NAME )
    {
        mLibName = rValue;
    }
    else if ( rThisVariable == USERMOD_NAME )
    {
        mName = rValue;
    }
}

void UDSMe::InitializeMaterial ( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues )
{
    // retrieve soil model number
    if (props.Has(SOIL_MODEL_NUMBER))
        mModelNumber = props[SOIL_MODEL_NUMBER];

    // retrieve undrained case
    if (props.Has(IS_UNDRAINED))
        mIsUndr = static_cast<int>(props[IS_UNDRAINED]);

    // retrieve material parameters
    if (props.Has(MATERIAL_PARAMETERS))
        mProps = props[MATERIAL_PARAMETERS];

    // retrieve library name
    if (props.Has(PLAXIS_LIBRARY_NAME))
        mLibName = props[PLAXIS_LIBRARY_NAME];

    // retrieve subroutine name
    if (props.Has(USERMOD_NAME))
        mName = props[USERMOD_NAME];

    // retrieve Bulk modulus of water
    if (props.Has(BULK_W))
        mBulkW = props[BULK_W];

#ifdef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    UserMod = udsm_;
#else
    if (minstances == 0)
    {
#ifdef _MSC_VER
        // TODO
#else
        mp_udsm_handle = dlopen(mLibName.c_str(), RTLD_NOW | RTLD_GLOBAL);
#endif
        if(mp_udsm_handle == 0)
        {
            KRATOS_ERROR << "Error loading Plaxis material library " << mLibName;
        }
#ifdef _MSC_VER
        // TODO
#else
        char* error;
        UserMod = (udsm_t) dlsym(mp_udsm_handle, mName.c_str());
        error = dlerror();
        if(error != NULL)
        {
            KRATOS_ERROR << "Error loading subroutine " << mName << " in the " << mLibName << " library, error message = " << error;
        }
        else
        {
            std::cout << "Loading subroutine " << mName << " in the " << mLibName << " library successfully" << std::endl;
        }
#endif
    }
    #pragma omp atomic
    ++minstances;
#endif

    // initialize the material
    BaseType::InitializeMaterial( mModelNumber, mIsUndr, mProps );
}

int UDSMe::Check ( const Properties& props,
                   const GeometryType& geom,
                   const ProcessInfo& CurrentProcessInfo ) const
{
    // DO NOTHING
    return 0;
}

/*****************************************************************************/
/*********** UDSMe Implicit-Explicit ******************************************/
/*****************************************************************************/

void UDSMeImplex::CalculateMaterialResponse ( const Vector& StrainVector,
                                              const Matrix& DeformationGradient,
                                              Vector& StressVector,
                                              Matrix& AlgorithmicTangent,
                                              const ProcessInfo& CurrentProcessInfo,
                                              const Properties& props,
                                              const GeometryType& geom,
                                              const Vector& ShapeFunctionsValues,
                                              bool CalculateStresses,
                                              int CalculateTangent,
                                              bool SaveInternalVariables )
{
    const double omega = 1.0;

    UDSM::DebugInfo proc_info;
    proc_info.time = CurrentProcessInfo[TIME];
    proc_info.delta_time = CurrentProcessInfo[DELTA_TIME];
    proc_info.step = mStep;
    proc_info.iteration = mIter;
    proc_info.element_id = mElemId;
    proc_info.point_id = mGaussId;

    UDSM::Input input( mLastStrain, mLastStress, mLastStateVariables,
            mLastExcessPorePressure);

    UDSM::Output output( mCurrentStrain, mCurrentStress, mCurrentStateVariables,
            mCurrentExcessPorePressure, mPlasticState, AlgorithmicTangent );

    UDSMImplex::StressIntegration( StrainVector, StressVector, input, output,
            UserMod, omega, proc_info, mModelNumber, mIsUndr, mBulkW, mProps );
}

/*****************************************************************************/
/*********** UDSMe Implicit ***************************************************/
/*****************************************************************************/

void UDSMeImplicit::CalculateMaterialResponse ( const Vector& StrainVector,
                                                const Matrix& DeformationGradient,
                                                Vector& StressVector,
                                                Matrix& AlgorithmicTangent,
                                                const ProcessInfo& CurrentProcessInfo,
                                                const Properties& props,
                                                const GeometryType& geom,
                                                const Vector& ShapeFunctionsValues,
                                                bool CalculateStresses,
                                                int CalculateTangent,
                                                bool SaveInternalVariables )
{
    if (CurrentProcessInfo[SET_CALCULATE_REACTION])
    {
        if (CalculateStresses)
        {
            noalias(StressVector) = mCurrentStress;
            return;
        }
    }

    UDSM::DebugInfo proc_info;
    proc_info.time = CurrentProcessInfo[TIME];
    proc_info.delta_time = CurrentProcessInfo[DELTA_TIME];
    proc_info.step = mStep;
    proc_info.iteration = mIter;
    proc_info.element_id = mElemId;
    proc_info.point_id = mGaussId;

    // save the current strain
    // KRATOS_WATCH(StrainVector)
    UDSM::VectorTo3DVector(StrainVector, mCurrentStrain);

    UDSM::Input input( mLastStrain, mLastStress, mLastStateVariables,
            mLastExcessPorePressure);

    UDSM::Output output( mCurrentStrain, mCurrentStress, mCurrentStateVariables,
            mCurrentExcessPorePressure, mPlasticState, AlgorithmicTangent );

    UDSMImplicit::StressIntegration( input, output, UserMod, proc_info,
            mModelNumber, mIsUndr, mBulkW, mProps );

    if (CalculateStresses)
    {
        // export the stress
        UDSM::Vector3DToVector(mCurrentStress, StressVector);
    }
}

} // Namespace Kratos
