#include <dlfcn.h>
#include "udsm.h"
#include "constitutive_laws_application_variables.h"
#include "structural_application/structural_application_variables.h"

#ifdef KRATOS_UDSM_LIBRARY_IS_PROVIDED
extern "C" void udsm_(int* IDTask, int* iMod, int* IsUndr, int* iStep, int* iTer, int* iEl,
                    int* Int, double* X, double* Y, double* Z, double* Time0, double* dTime,
                    double* Props, double* Sig0, double* Swp0, double* StVar0, double* dEps,
                    double* D, double* BulkW, double* Sig, double* Swp, double* StVar,
                    int* ipl, int* nStat, int* NonSym, int* iStrsDep, int* iTimeDep,
                    int* iTang, int* iAbort);
#endif

namespace Kratos
{

UDSM::UDSM()
{
    mModelNumber = 1; // default model number
}

UDSM::~UDSM()
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    if(mp_udsm_handle != 0)
        dlclose(mp_udsm_handle);
    #endif
}

bool UDSM::Has ( const Variable<double>& rThisVariable )
{
//    if(rThisVariable == PLASTICITY_INDICATOR)
//        return true;
    return false;
}

bool UDSM::Has ( const Variable<Vector>& rThisVariable )
{
//    if(rThisVariable == MATERIAL_PARAMETERS)
//        return true;
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
//    if(rThisVariable == STRESSES)
//        return true;
    return false;
}

double& UDSM::GetValue ( const Variable<double>& rThisVariable, double& rValue )
{
//    if(rThisVariable == PLASTICITY_INDICATOR)
//        rValue = mPlasticState;
    return rValue;
}

Vector& UDSM::GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue )
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

    return rValue;
}

void UDSM::SetValue ( const Variable<int>& rThisVariable,
                      const int& rValue,
                      const ProcessInfo& rCurrentProcessInfo )
{}

void UDSM::SetValue ( const Variable<Vector>& rThisVariable,
                      const Vector& rValue,
                      const ProcessInfo& rCurrentProcessInfo )
{}

void UDSM::ResetMaterial ( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{}

int UDSM::Check ( const Properties& props,
                  const GeometryType& geom,
                  const ProcessInfo& CurrentProcessInfo )
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    if(props.Has( PLAXIS_LIBRARY_NAME ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define PLAXIS_LIBRARY_NAME", "")
    }
    if(props.Has( USERMOD_NAME ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define USERMOD_NAME", "")
    }
    #endif
    if(props.Has( SOIL_MODEL_NUMBER ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define SOIL_MODEL_NUMBER", "")
    }
    if(props.Has( IS_UNDRAINED ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define IS_UNDRAINED", "")
    }
    if(props.Has( MATERIAL_PARAMETERS ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define MATERIAL_PARAMETERS", "")
    }
    return 0;
}

void UDSM::InitializeMaterial ( const Properties& props,
                                const GeometryType& geom,
                                const Vector& ShapeFunctionsValues )
{
    // retrieve soil model number
    mModelNumber = props[SOIL_MODEL_NUMBER];

    // retrieve undrained case
    int IsUndr = static_cast<int>(props[IS_UNDRAINED]);

    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    // get the library name and load the udsm subroutine
    std::string lib_name = props[PLAXIS_LIBRARY_NAME];
//    mp_udsm_handle = dlopen(lib_name.c_str(), RTLD_LAZY);
    mp_udsm_handle = dlopen(lib_name.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if(mp_udsm_handle == 0)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "Error loading Plaxis material library", lib_name)
    }

    std::string udsm_name = props[USERMOD_NAME];
    char* error;
    UserMod = (void (*)(int* IDTask, int* iMod, int* IsUndr, int* iStep, int* iTer, int* iEl,
                    int* Int, double* X, double* Y, double* Z, double* Time0, double* dTime,
                    double* Props, double* Sig0, double* Swp0, double* StVar0, double* dEps,
                    double* D, double* BulkW, double* Sig, double* Swp, double* StVar,
                    int* ipl, int* nStat, int* NonSym, int* iStrsDep, int* iTimeDep,
                    int* iTang, int* iAbort)) dlsym(mp_udsm_handle, udsm_name.c_str());
    error = dlerror();
    if(error != NULL)
    {
        std::stringstream ss;
        ss << "Error loading subroutine " << udsm_name << " in the " << lib_name << " library, error message = " << error;
        KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
    }
    #else
    UserMod = udsm_;
    #endif

    // retrieve number of state variables
    int nStat;
    int IDTask = 4;
    UserMod(&IDTask, &mModelNumber, &IsUndr, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, &nStat, NULL, NULL, NULL, NULL, NULL);

    if(mCurrentStateVariables.size() != nStat)
        mCurrentStateVariables.resize(nStat);
    mCurrentStateVariables = ZeroVector(nStat);

    if(mLastStateVariables.size() != nStat)
        mLastStateVariables.resize(nStat);
    mLastStateVariables = ZeroVector(nStat);

    if(mCurrentStress.size() != 6)
        mCurrentStress.resize(6);
    noalias(mCurrentStress) = ZeroVector(6);

    if(mLastStress.size() != 6)
        mLastStress.resize(6);
    noalias(mLastStress) = ZeroVector(6);

    if(mCurrentStrain.size() != 6)
        mCurrentStrain.resize(6);
    noalias(mCurrentStrain) = ZeroVector(6);

    if(mLastStrain.size() != 6)
        mLastStrain.resize(6);
    noalias(mLastStrain) = ZeroVector(6);

    mCurrentExcessPorePressure = 0.0;
    mLastExcessPorePressure = 0.0;
    mPlasticState = 0;
}

void UDSM::InitializeSolutionStep ( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{}

void UDSM::InitializeNonLinearIteration ( const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues,
                                          const ProcessInfo& CurrentProcessInfo )
{}

void UDSM::CalculateMaterialResponse ( const Vector& StrainVector,
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
    KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
}

void UDSM::FinalizeNonLinearIteration ( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues,
                                        const ProcessInfo& CurrentProcessInfo )
{}

void UDSM::FinalizeSolutionStep ( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ,
                                  const ProcessInfo& CurrentProcessInfo )
{
    noalias(mLastStrain) = mCurrentStrain;
    noalias(mLastStress) = mCurrentStress;
    noalias(mLastStateVariables) = mCurrentStateVariables;
    mLastExcessPorePressure = mCurrentExcessPorePressure;
}




/*****************************************************************************/
/*********** UDSM Implicit-Explicit ******************************************/
/*****************************************************************************/




void UDSMImplex::CalculateMaterialResponse ( const Vector& StrainVector,
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
    int IsUndr = static_cast<int>(props[IS_UNDRAINED]);
    Vector Props = props[MATERIAL_PARAMETERS];
    Vector D(36);
    double BulkW;
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
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
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &Props[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iAbort );

    // export the stiffness matrix
    IDTask = 5;
    UserMod(&IDTask, &mModelNumber, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            &NonSym,
            &iStrsDep,
            &iTimeDep,
            &iTang,
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
    double omega = Props[49]; // relaxation factor
    noalias(StressVector) = mCurrentStress + omega*prod(AlgorithmicTangent, delta_strain);

    // calculate the constitutive stress
    IDTask = 2;
    UserMod( &IDTask,
             &mModelNumber,
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &Props[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iAbort );

//    KRATOS_WATCH(mCurrentStress)
//    KRATOS_WATCH(iPl)

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating stress, iAbort =", iAbort)

    // export the plastic state
    mPlasticState = iPl;

    // compute for the undrained state
    if(IsUndr != 0)
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
/*********** UDSM Implicit ***************************************************/
/*****************************************************************************/





void UDSMImplicit::CalculateMaterialResponse ( const Vector& StrainVector,
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
    int IsUndr = static_cast<int>(props[IS_UNDRAINED]);
    Vector Props = props[MATERIAL_PARAMETERS];
    Vector D(36);
    double BulkW;
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
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
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &Props[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iAbort );

    // calculate the consitutive stress
    IDTask = 2;
    UserMod( &IDTask,
             &mModelNumber,
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &Props[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
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
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             &Props[0],
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &mCurrentStress[0],
             &mCurrentExcessPorePressure,
             &mCurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iAbort );
//    KRATOS_WATCH(D)
//    iAbort=1;

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating material stiffness, iAbort =", iAbort)

    // compute for the undrained state
    if(IsUndr != 0)
    {
        // TODO add the bulk stiffness of water to the material stiffness matrix
    }

    // save the current strain
    noalias(mCurrentStrain) = StrainVector;

    // export the stress
    noalias(StressVector) = mCurrentStress;

    // export the stiffness matrix
    IDTask = 5;
    UserMod(&IDTask, &mModelNumber, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            &NonSym,
            &iStrsDep,
            &iTimeDep,
            &iTang,
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

//    if(mPlasticState)
//    {
//        KRATOS_WATCH(StressVector)
//        KRATOS_WATCH(AlgorithmicTangent)
//        KRATOS_WATCH(mCurrentStateVariables)
//    }
}


} // Namespace Kratos

#undef USRMOD1
#undef USRMOD2
#undef DEFINE_USRMOD_PLAXIS

