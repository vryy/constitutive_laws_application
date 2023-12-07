#include <iomanip>

#include <dlfcn.h>
#include "umat3.h"
#include "udsm.h"
#include "constitutive_laws_application_variables.h"
#include "structural_application/structural_application_variables.h"

namespace Kratos
{

#ifdef KRATOS_UDSM_LIBRARY_IS_PROVIDED
extern "C" void udsm_(int* IDTask, int* iMod, int* IsUndr, int* iStep, int* iTer, int* iEl,
                    int* Int, double* X, double* Y, double* Z, double* Time0, double* dTime,
                    double* Props, double* Sig0, double* Swp0, double* StVar0, double* dEps,
                    double* D, double* BulkW, double* Sig, double* Swp, double* StVar,
                    int* ipl, int* nStat, int* NonSym, int* iStrsDep, int* iTimeDep,
                    int* iTang, int* iPrjDir, int* iPrjLen, int* iAbort);
#else
unsigned long long UDSM::minstances = 0;
void* UDSM::mp_udsm_handle = 0;
udsm_t UDSM::UserMod = 0;
#endif

UDSM::UDSM()
{
    mModelNumber = 1; // default model number
}

UDSM::~UDSM()
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    --minstances;
    if(minstances == 0)
    {
        dlclose(mp_udsm_handle);
        std::cout << "Successfully unload the UDSM shared library/DLL" << std::endl;
    }
    #endif
}

bool UDSM::Has ( const Variable<double>& rThisVariable )
{
    return false;
}

bool UDSM::Has ( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
        return true;
    return false;
}

double& UDSM::GetValue ( const Variable<double>& rThisVariable, double& rValue )
{
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
{
    if( rThisVariable == PARENT_ELEMENT_ID )
        mElemId = rValue;
    if( rThisVariable == INTEGRATION_POINT_INDEX )
        mGaussId = rValue;
}

void UDSM::SetValue ( const Variable<double>& rThisVariable,
                      const double& rValue,
                      const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
}

void UDSM::SetValue ( const Variable<Vector>& rThisVariable,
                      const Vector& rValue,
                      const ProcessInfo& rCurrentProcessInfo )
{
    if(rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS)
    {
        if (mPrestress.size() != rValue.size())
            mPrestress.resize(rValue.size());
        noalias(mPrestress) = rValue;
    }
}

void UDSM::ResetMaterial ( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    mCurrentStrain.clear();
    mLastStrain.clear();
    mCurrentStateVariables.clear();
    mLastStateVariables.clear();
    mCurrentExcessPorePressure = 0.0;
    mLastExcessPorePressure = 0.0;

    noalias(mLastStress) = -mPrestressFactor*mPrestress;
    noalias(mCurrentStress) = mLastStress;
}

int UDSM::Check ( const Properties& props,
                  const GeometryType& geom,
                  const ProcessInfo& CurrentProcessInfo ) const
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
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
    int ModelNumber = static_cast<int>(props[SOIL_MODEL_NUMBER]);

    // retrieve undrained case
    int IsUndrained = static_cast<int>(props[IS_UNDRAINED]);

    // retrieve material properties
    const Vector& Props = props[MATERIAL_PARAMETERS];

#ifdef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    UserMod = udsm_;
#else
    if (minstances == 0)
    {
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
        UserMod = (udsm_t) dlsym(mp_udsm_handle, udsm_name.c_str());
        error = dlerror();
        if(error != NULL)
        {
            std::stringstream ss;
            ss << "Error loading subroutine " << udsm_name << " in the " << lib_name << " library, error message = " << error;
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }
        else
        {
            std::cout << "Loading subroutine " << udsm_name << " in the " << lib_name << " library successfully" << std::endl;
        }
    }
    #pragma omp atomic
    ++minstances;
#endif

    // initialize the material
    this->InitializeMaterial( ModelNumber, IsUndrained, Props );
}

void UDSM::InitializeMaterial ( int ModelNumber,
                                int IsUndrained,
                                const Vector& Props )
{
// KRATOS_WATCH(__LINE__)
// KRATOS_WATCH(ModelNumber)
    // retrieve number of state variables
    int nStat;
    int IDTask = 4;
    // UserMod(&IDTask, &ModelNumber, &IsUndrained, NULL, NULL, NULL, NULL, NULL,
    //         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    //         NULL, NULL, NULL, NULL, &nStat, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    int iStep=0, iTer=0, Iel=mElemId, Int=mGaussId, dummy_ipl, dummy_nonsym,
        dummy_iStrsDep, dummy_iTimeDep, dummy_iTang, dummy_iPrjDir, dummy_iPrjLen;
    int iAbort;
    double dummyX, dummyY, dummyZ;
    double dummy_T0, dummy_dT;
    double dummy_Sig0, dummy_Swp0, dummy_Stvar0, dummy_dEps, dummy_D, dummy_BulkW,
        dummy_Sig, dummy_Swp, dummy_StVar;
    double* props = const_cast<double*>(&Props[0]); // cast out the const'ness to use the pointer value and let hope that the UDSM does not change any input parameter.
    UserMod(&IDTask, &ModelNumber, &IsUndrained,
            &iStep, &iTer, &Iel, &Int,
            &dummyX, &dummyY, &dummyZ, &dummy_T0, &dummy_dT,
            props,
            &dummy_Sig0, &dummy_Swp0, &dummy_Stvar0, &dummy_dEps, &dummy_D, &dummy_BulkW,
            &dummy_Sig, &dummy_Swp, &dummy_StVar, &dummy_ipl,
            &nStat,
            &dummy_nonsym, &dummy_iStrsDep, &dummy_iTimeDep, &dummy_iTang, &dummy_iPrjDir, &dummy_iPrjLen,
            &iAbort);

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

    if(mPrestress.size() != 6)
        mPrestress.resize(6);
    noalias(mPrestress) = ZeroVector(6);

    if(mCurrentStrain.size() != 6)
        mCurrentStrain.resize(6);
    noalias(mCurrentStrain) = ZeroVector(6);

    if(mLastStrain.size() != 6)
        mLastStrain.resize(6);
    noalias(mLastStrain) = ZeroVector(6);

    mCurrentExcessPorePressure = 0.0;
    mLastExcessPorePressure = 0.0;
    mPrestressFactor = 1.0;
    mPlasticState = 0;

    mStep = 0;
    mIter = 0;

    // TODO initialize the internal variables by calling IDTask = 1
}

void UDSM::InitializeSolutionStep ( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{
    ++mStep;
    mIter = 0;
}

void UDSM::InitializeNonLinearIteration ( const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues,
                                          const ProcessInfo& CurrentProcessInfo )
{
    ++mIter;
}

void UDSM::CalculateMaterialResponseCauchy( Parameters& parameters )
{
    this->CalculateMaterialResponse( parameters.GetStrainVector()
        , parameters.GetDeformationGradientF()
        , parameters.GetStressVector()
        , parameters.GetConstitutiveMatrix()
        , parameters.GetProcessInfo()
        , parameters.GetMaterialProperties()
        , parameters.GetElementGeometry()
        , parameters.GetShapeFunctionsValues()
        , parameters.IsSetStressVector()
        , parameters.IsSetConstitutiveMatrix()
        , true
    );
}

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

void UDSM::VectorTo3DVector(const Vector& vector, Vector& vector_3d)
{
    if (vector.size() == 3)
    {
        vector_3d[0] = vector[0];
        vector_3d[1] = vector[1];
        vector_3d[2] = 0.0;
        vector_3d[3] = vector[2];
        vector_3d[4] = 0.0;
        vector_3d[5] = 0.0;
    }
    else if (vector.size() == 6)
    {
        noalias(vector_3d) = vector;
    }
}

void UDSM::Vector3DToVector(const Vector& vector_3d, Vector& vector)
{
    if (vector.size() == 3)
    {
        vector[0] = vector_3d[0];
        vector[1] = vector_3d[1];
        vector[2] = vector_3d[3];
    }
    else if (vector.size() == 6)
    {
        noalias(vector) = vector_3d;
    }
}

void UDSM::Vector1DToMatrix(const Vector& D, Matrix& A, int non_sym)
{
    if (A.size1() == 3)
    {
        if(non_sym != 0)
        {
            for(unsigned int i = 0; i < 3; ++i)
                for(unsigned int j = 0; j < 3; ++j)
                    A(i, j) = D[Umat3::PS[i] + 6*Umat3::PS[j]];
        }
        else
        {
            for(unsigned int i = 0; i < 3; ++i)
                for(unsigned int j = i; j < 3; ++j)
                    A(i, j) = D[Umat3::PS[i] + 6*Umat3::PS[j]];
            for(unsigned int i = 0; i < 3; ++i)
                for(unsigned int j = 0; j < i; ++j)
                    A(i, j) = A(j, i);
        }
    }
    else if (A.size1() == 6)
    {
        if(non_sym != 0)
        {
            for(unsigned int i = 0; i < 6; ++i)
                for(unsigned int j = 0; j < 6; ++j)
                    A(i, j) = D[i + 6*j];
        }
        else
        {
            for(unsigned int i = 0; i < 6; ++i)
                for(unsigned int j = i; j < 6; ++j)
                    A(i, j) = D[i + 6*j];
            for(unsigned int i = 0; i < 6; ++i)
                for(unsigned int j = 0; j < i; ++j)
                    A(i, j) = A(j, i);
        }
    }
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
    int ModelNumber = static_cast<int>(props[SOIL_MODEL_NUMBER]);
    int IsUndr = static_cast<int>(props[IS_UNDRAINED]);
    Vector Props = props[MATERIAL_PARAMETERS];
    Vector D(36);
    double BulkW;
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen; // can use this for debugging purpose (TODO)
    int nStat = mCurrentStateVariables.size();

    if (IsUndr)
    {
        BulkW = props[BULK_W];
    }

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
             &ModelNumber,
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
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // export the stiffness matrix
    IDTask = 5;
    UserMod(&IDTask, &ModelNumber, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
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
    double omega = Props[UDSM_RELAXATION_FACTOR]; // relaxation factor
    noalias(StressVector) = mCurrentStress + omega*prod(AlgorithmicTangent, delta_strain);

    // calculate the stress
    IDTask = 2;
    UserMod( &IDTask,
             &ModelNumber,
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
    int ModelNumber = static_cast<int>(props[SOIL_MODEL_NUMBER]);
    int IsUndr = static_cast<int>(props[IS_UNDRAINED]);
    Vector Props = props[MATERIAL_PARAMETERS];
    Vector D(36);
    double BulkW;
    double time = CurrentProcessInfo[TIME];
    double delta_time = CurrentProcessInfo[DELTA_TIME];
    int iPl;
    int NonSym;
    int iStrsDep, iTang = 0, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen;
    int nStat = mCurrentStateVariables.size();
    double dummyX, dummyY, dummyZ;
    int iStep=mStep, iTer=mIter, Iel=mElemId, Int=mGaussId;

    if (IsUndr)
    {
        BulkW = props[BULK_W];
    }

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
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time,         // Time0
             &delta_time,   // dTime
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
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // calculate the stress
    IDTask = 2;
    UserMod( &IDTask,
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time,         // Time0
             &delta_time,   // dTime
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
             &iPrjDir,
             &iPrjLen,
             &iAbort );

//    KRATOS_WATCH(mCurrentStress)
//    KRATOS_WATCH(iPl)

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating stress, iAbort =", iAbort)

    // export the plastic state
    mPlasticState = iPl;

    // obtain the stiffness matrix properties
    IDTask = 5;
    UserMod( &IDTask,
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time,         // Time0
             &delta_time,   // dTime
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
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // if (iTang == 1)
    // {
        // calculate the material stiffness
        IDTask = 3;
        UserMod( &IDTask,
                 &ModelNumber,
                 &IsUndr,
                 &iStep,        // iStep
                 &iTer,         // iter
                 &Iel,          // Iel
                 &Int,          // Int
                 &dummyX, &dummyY, &dummyZ, // X, Y, Z
                 &time,         // Time0
                 &delta_time,   // dTime
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
                 &iPrjDir,
                 &iPrjLen,
                 &iAbort );
    // }
    // else if (iTang == 0)
    // {
    //     // TODO approximate the tangent using finite difference scheme
    // }

//    KRATOS_WATCH(D)
//    iAbort=1;

    if(iAbort != 0)
        KRATOS_THROW_ERROR(std::logic_error, "Force calculation stop after calculating material stiffness, iAbort =", iAbort)

    // compute for the undrained state
    if(IsUndr != 0)
    {
        // TODO add the bulk stiffness of water to the material stiffness matrix
    }

    if (CalculateStresses)
    {
        // export the stress
        this->Vector3DToVector(mCurrentStress, StressVector);
    }

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
