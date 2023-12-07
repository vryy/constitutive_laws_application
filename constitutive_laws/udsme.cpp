#include <iomanip>

#include <dlfcn.h>
#include "constitutive_laws/umat3.h"
#include "constitutive_laws/udsme.h"
#include "structural_application/structural_application_variables.h"

namespace Kratos
{

UDSMe::UDSMe() : BaseType()
{}

UDSMe::~UDSMe()
{}

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
//    if(rThisVariable == PLASTICITY_INDICATOR)
//        return true;
    if(rThisVariable == BULK_W)
        return true;
    return BaseType::Has( rThisVariable );
}

bool UDSMe::Has ( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == MATERIAL_PARAMETERS)
        return true;
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
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
//    if(rThisVariable == PLASTICITY_INDICATOR)
//        rValue = mPlasticState;
    if(rThisVariable == PRESTRESS_FACTOR)
        rValue = mPrestressFactor;
    if(rThisVariable == BULK_W)
        rValue = mBulkW;
    return BaseType::GetValue( rThisVariable, rValue );
}

Vector& UDSMe::GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable == INTERNAL_VARIABLES)
    {
        if(rValue.size() != mCurrentStateVariables.size())
            rValue.resize(mCurrentStateVariables.size(), false);
        noalias(rValue) = mCurrentStateVariables;
    }

    if(rThisVariable == MATERIAL_PARAMETERS)
    {
        if(rValue.size() != mProps.size())
            rValue.resize(mProps.size(), false);
        noalias(rValue) = mProps;
    }

    if(rThisVariable == STRESSES)
    {
        if(rValue.size() != mCurrentStress.size())
            rValue.resize(mCurrentStress.size(), false);
        noalias(rValue) = mCurrentStress;
    }

    if(rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS)
    {
        if(rValue.size() != mPrestress.size())
            rValue.resize(mPrestress.size(), false);
        noalias(rValue) = mPrestress;
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
        mp_udsm_handle = dlopen(mLibName.c_str(), RTLD_NOW | RTLD_GLOBAL);
        if(mp_udsm_handle == 0)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Error loading Plaxis material library", mLibName)
        }

        char* error;
        UserMod = (udsm_t) dlsym(mp_udsm_handle, mName.c_str());
        error = dlerror();
        if(error != NULL)
        {
            std::stringstream ss;
            ss << "Error loading subroutine " << mName << " in the " << mLibName << " library, error message = " << error;
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }
        else
        {
            std::cout << "Loading subroutine " << mName << " in the " << mLibName << " library successfully" << std::endl;
        }
    }
    #pragma omp atomic
    ++minstances;
#endif

    // initialize the material
    BaseType::InitializeMaterial( mModelNumber, mIsUndr, mProps );
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
    int IDTask;
    Vector D(36);
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen; // can use this for debugging purpose (TODO)
    int nStat = mCurrentStateVariables.size();

    // compute the vector of previous effective stress component
    Vector Sig0(20);
    noalias(Sig0) = ZeroVector(20);
    Sig0[0] = mLastStress[0];
    Sig0[1] = mLastStress[1];
    Sig0[2] = mLastStress[2];
    Sig0[3] = mLastStress[3];
    Sig0[4] = mLastStress[4];
    Sig0[5] = mLastStress[5];

    // compute the previous excess pore pressure
    double Swp0 = mLastExcessPorePressure;

    // compute the previous value of state variables
    Vector StVar0 = mLastStateVariables;

    // compute the vector of incremental strain
    Vector dEps(12);
    dEps[0] = StrainVector[0] - mLastStrain[0];
    dEps[1] = StrainVector[1] - mLastStrain[1];
    dEps[2] = StrainVector[2] - mLastStrain[2];
    dEps[3] = StrainVector[3] - mLastStrain[3];
    dEps[4] = StrainVector[4] - mLastStrain[4];
    dEps[5] = StrainVector[5] - mLastStrain[5];
    dEps[6] = mLastStrain[0];
    dEps[7] = mLastStrain[1];
    dEps[8] = mLastStrain[2];
    dEps[9] = mLastStrain[3];
    dEps[10] = mLastStrain[4];
    dEps[11] = mLastStrain[5];

    // initialize D with elastic matrix
    IDTask = 6;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // export the stiffness matrix
    IDTask = 5;
    UserMod(&IDTask, &mModelNumber, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            &NonSym,
            &iStrsDep,
            &iTimeDep,
            &iTang,
            &iPrjDir,
            &iPrjLen,
            &iAbort );

    if(NonSym != 0)
    {
        for(std::size_t i = 0; i < 6; ++i)
            for(std::size_t j = 0; j < 6; ++j)
                AlgorithmicTangent(i, j) = D[i + 6*j];
    }
    else
    {
        for(std::size_t i = 0; i < 6; ++i)
            for(std::size_t j = i; j < 6; ++j)
                AlgorithmicTangent(i, j) = D[i + 6*j];
        for(std::size_t i = 0; i < 6; ++i)
            for(std::size_t j = 0; j < i; ++j)
                AlgorithmicTangent(i, j) = AlgorithmicTangent(j, i);
    }

    // compute delta_strain
    Vector delta_strain = StrainVector - mCurrentStrain;
//    KRATOS_WATCH(delta_strain)

    // compute equilibrium stress
    // double omega = Props[UDSM_RELAXATION_FACTOR]; // relaxation factor
    // noalias(StressVector) = mCurrentStress + omega*prod(AlgorithmicTangent, delta_strain);
    noalias(StressVector) = mCurrentStress + prod(AlgorithmicTangent, delta_strain);

    // calculate the constitutive stress
    IDTask = 2;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

//    KRATOS_WATCH(mCurrentStress)
//    KRATOS_WATCH(iPl)

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating stress, iAbort =", iAbort)

    // export the plastic state
    mPlasticState = iPl;

    // compute for the undrained state
    if(mIsUndr != 0)
    {
        // TODO add the bulk stiffness of water to the material stiffness matrix
    }

    // save the current strain
    noalias(mCurrentStrain) = StrainVector;

//    if(mPlasticState)
//    {
//        KRATOS_WATCH(StressVector)
//        KRATOS_WATCH(AlgorithmicTangent)
//        KRATOS_WATCH(mCurrentStateVariables)
//    }
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

    // save the current strain
    // KRATOS_WATCH(StrainVector)
    this->VectorTo3DVector(StrainVector, mCurrentStrain);

    int IDTask;
    Vector D(36);
    double time = CurrentProcessInfo[TIME];
    double delta_time = CurrentProcessInfo[DELTA_TIME];
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen;
    int nStat = mCurrentStateVariables.size();
    double dummyX, dummyY, dummyZ;
    int dummy_iStep=0, dummy_iTer=0, dummy_Iel=0, dummy_Int=0;

    // compute the vector of previous effective stress component
    Vector Sig0(20);
    noalias(Sig0) = ZeroVector(20);
    Sig0[0] = mLastStress[0];
    Sig0[1] = mLastStress[1];
    Sig0[2] = mLastStress[2];
    Sig0[3] = mLastStress[3];
    Sig0[4] = mLastStress[4];
    Sig0[5] = mLastStress[5];

    // compute the previous excess pore pressure
    double Swp0 = mLastExcessPorePressure;

    // compute the previous value of state variables
    Vector StVar0 = mLastStateVariables;

    // KRATOS_WATCH(mCurrentStrain)

    // compute the vector of incremental strain
    Vector dEps(12);
    dEps[0] = mCurrentStrain[0] - mLastStrain[0];
    dEps[1] = mCurrentStrain[1] - mLastStrain[1];
    dEps[2] = mCurrentStrain[2] - mLastStrain[2];
    dEps[3] = mCurrentStrain[3] - mLastStrain[3];
    dEps[4] = mCurrentStrain[4] - mLastStrain[4];
    dEps[5] = mCurrentStrain[5] - mLastStrain[5];
    dEps[6] = mLastStrain[0];
    dEps[7] = mLastStrain[1];
    dEps[8] = mLastStrain[2];
    dEps[9] = mLastStrain[3];
    dEps[10] = mLastStrain[4];
    dEps[11] = mLastStrain[5];

    // initialize D with elastic matrix
    IDTask = 6;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             &dummy_iStep, // iStep
             &dummy_iTer, // iter
             &dummy_Iel, // Iel
             &dummy_Int, // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time, // Time0
             &delta_time, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // calculate the consitutive stress
    IDTask = 2;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             &dummy_iStep, // iStep
             &dummy_iTer, // iter
             &dummy_Iel, // Iel
             &dummy_Int, // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time, // Time0
             &delta_time, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

//    KRATOS_WATCH(mCurrentStress)
//    KRATOS_WATCH(iPl)

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating stress, iAbort =", iAbort)

    // export the plastic state
    mPlasticState = iPl;

    // calculate the material stiffness
    IDTask = 3;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             &dummy_iStep, // iStep
             &dummy_iTer, // iter
             &dummy_Iel, // Iel
             &dummy_Int, // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time, // Time0
             &delta_time, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );
//    KRATOS_WATCH(D)
//    iAbort=1;

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating material stiffness, iAbort =", iAbort)

    // compute for the undrained state
    if(mIsUndr != 0)
    {
        // TODO add the bulk stiffness of water to the material stiffness matrix
    }

    if (CalculateStresses)
    {
        // export the stress
        this->Vector3DToVector(mCurrentStress, StressVector);
    }

    // obtain the stiffness matrix properties
    IDTask = 5;
    UserMod( &IDTask,
             &mModelNumber,
             &mIsUndr,
             &dummy_iStep, // iStep
             &dummy_iTer, // iter
             &dummy_Iel, // Iel
             &dummy_Int, // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time, // Time0
             &delta_time, // dTime
             &mProps[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &mBulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    if (CalculateTangent)
    {
        this->Vector1DToMatrix(D, AlgorithmicTangent, NonSym);
    }

//    if(mPlasticState)
//    {
       // KRATOS_WATCH(StressVector)
       // KRATOS_WATCH(AlgorithmicTangent)
//        KRATOS_WATCH(mCurrentStateVariables)
//    }
}


} // Namespace Kratos

#ifdef USRMOD1
#undef USRMOD1
#endif

#ifdef USRMOD2
#undef USRMOD2
#endif

#ifdef DEFINE_USRMOD_PLAXIS
#undef DEFINE_USRMOD_PLAXIS
#endif
