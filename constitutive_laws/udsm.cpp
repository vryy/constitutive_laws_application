#include <iomanip>

#include "umat3.h"
#include "udsm.h"
#include "includes/c2c_variables.h"
#include "custom_utilities/shared_library_handle_shield.h"
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
#endif

UDSM::UDSM()
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    mp_udsm_handle = nullptr;
    #endif
    UserMod = nullptr;
    mModelNumber = 1; // default model number
}

UDSM::~UDSM()
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    mp_udsm_handle = nullptr;
    #endif
    UserMod = nullptr;
}

bool UDSM::Has ( const Variable<double>& rThisVariable )
{
    if(rThisVariable == PLASTICITY_INDICATOR)
        return true;
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
    if(rThisVariable == PLASTICITY_INDICATOR)
        rValue = static_cast<double>(mPlasticState); // since the precised location of accumulated plastic strain
                        // in general unknown in the internal state vector, the plasticity indicator
                        // is taken as the plastic flag returned by the constitutive law
    else if(rThisVariable == PRESTRESS_FACTOR)
        rValue = mPrestressFactor;
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
    else if(rThisVariable == STRESSES)
    {
        if(rValue.size() == 6)              // 3D
        {
            noalias(rValue) = mCurrentStress;
        }
        else if(rValue.size() == 3)         // plane strain
        {
            rValue[0] = mCurrentStress[0];
            rValue[1] = mCurrentStress[1];
            rValue[2] = mCurrentStress[3];
        }
        else
            KRATOS_ERROR << "Invalid size " << rValue.size();
    }
    else if(rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS)
    {
        if(rValue.size() != mPrestress.size())
            rValue.resize(mPrestress.size(), false);
        noalias(rValue) = mPrestress;
    }

    return rValue;
}

Matrix& UDSM::GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if(rThisVariable == ELASTIC_TANGENT)
    {
        int NonSym;
        int iStrsDep, iTang = 0, iTimeDep, iAbort = 0;
        int iPrjDir, iPrjLen;
        double* props = const_cast<double*>(&mProps[0]);

        Vector Sig0(20);
        noalias(Sig0) = ZeroVector(20);
        Sig0[0] = mLastStress[0];
        Sig0[1] = mLastStress[1];
        Sig0[2] = mLastStress[2];
        Sig0[3] = mLastStress[3];
        Sig0[4] = mLastStress[4];
        Sig0[5] = mLastStress[5];

        double Swp0 = mLastExcessPorePressure;

        int nStat = mLastStateVariables.size();
        Vector StVar0(nStat), StVar(nStat);
        StVar0.clear();
        StVar.clear();

        // initialize D with elastic matrix
        int IDTask = 6;
        Vector D(36);
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
                 props,
                 &Sig0[0],
                 &Swp0,
                 &StVar0[0],
                 NULL,
                 &D[0],
                 &mBulkW,
                 &mCurrentStress[0],
                 &mCurrentExcessPorePressure,
                 &mCurrentStateVariables[0],
                 &mPlasticState,
                 &nStat,
                 &NonSym,
                 &iStrsDep,
                 &iTimeDep,
                 &iTang,
                 &iPrjDir,
                 &iPrjLen,
                 &iAbort );

        UDSM::Vector1DToMatrix(D, rValue, NonSym);
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
        if (rValue.size() == 6)         // 3D
        {
            noalias(mPrestress) = rValue;
        }
        else if (rValue.size() == 3)    // plane strain
        {
            mPrestress[0] = rValue[0];
            mPrestress[1] = rValue[1];
            mPrestress[2] = 0.0;
            mPrestress[3] = rValue[2];
            mPrestress[4] = 0.0;
            mPrestress[5] = 0.0;
        }
        else
            KRATOS_ERROR << "Invalid size " << rValue.size();
    }
}

int UDSM::Check ( const Properties& props,
                  const GeometryType& geom,
                  const ProcessInfo& CurrentProcessInfo ) const
{
    #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    if(props.Has( PLAXIS_LIBRARY_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define PLAXIS_LIBRARY_NAME";
    }
    if(props.Has( USERMOD_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define USERMOD_NAME";
    }
    #endif
    if(props.Has( SOIL_MODEL_NUMBER ) == false)
    {
        KRATOS_ERROR << "Properties must define SOIL_MODEL_NUMBER";
    }
    if(props.Has( IS_UNDRAINED ) == false)
    {
        KRATOS_ERROR << "Properties must define IS_UNDRAINED";
    }
    if(props.Has( MATERIAL_PARAMETERS ) == false)
    {
        KRATOS_ERROR << "Properties must define MATERIAL_PARAMETERS";
    }
    return 0;
}

void UDSM::InitializeMaterial ( const Properties& props,
                                const GeometryType& geom,
                                const Vector& ShapeFunctionsValues )
{
    // retrieve soil model number
    mModelNumber = static_cast<int>(props[SOIL_MODEL_NUMBER]);

    // retrieve undrained case
    mIsUndr = static_cast<int>(props[IS_UNDRAINED]);

    // retrieve material properties
    mProps = props[MATERIAL_PARAMETERS];

    // retrieve bulk modulus of water
    mBulkW = 0.0;
    if (mIsUndr) mBulkW = props[BULK_W];

#ifdef KRATOS_UDSM_LIBRARY_IS_PROVIDED
    UserMod = udsm_;
#else
    {
        // get the library name and load the udsm subroutine
        std::string lib_name = props[PLAXIS_LIBRARY_NAME];
        mp_udsm_handle = DLL::GetSharedLibraryHandle(lib_name);
        if(mp_udsm_handle == nullptr)
        {
            KRATOS_ERROR << "Error loading Plaxis material library " << lib_name;
        }

        std::string udsm_name = props[USERMOD_NAME];
        UserMod = (udsm_t) DLL::GetSymbol(mp_udsm_handle, udsm_name);
        const char* error = DLL::GetError();
        if(error != nullptr)
        {
            KRATOS_ERROR << "Error loading subroutine " << udsm_name << " in the " << lib_name << " library"
                         << ", error message = " << DLL::GetErrorMessage(error);
        }
        else
        {
            std::cout << "Loading subroutine " << udsm_name << " in the " << lib_name << " library successfully" << std::endl;
        }
    }
#endif

    bool need_determine_internal_params = false;
    if (props.Has(NEED_DETERMINE_INTERNAL_PARAMS))
        need_determine_internal_params = props[NEED_DETERMINE_INTERNAL_PARAMS];

    // initialize the material
    this->InitializeMaterial( mModelNumber, mIsUndr, mProps, need_determine_internal_params );
}

void UDSM::InitializeMaterial ( int ModelNumber,
                                int IsUndr,
                                const Vector& Props,
                                bool need_determine_internal_params )
{
    // retrieve number of state variables
    // This is obtained by using IDTask = 4
    int nStat;
    int IDTask = 4;

    int iStep=0, iTer=0, Iel=mElemId, Int=mGaussId, dummy_ipl, dummy_nonsym,
        dummy_iStrsDep, dummy_iTimeDep, dummy_iTang, dummy_iPrjDir, dummy_iPrjLen;
    int iAbort;
    double dummyX, dummyY, dummyZ;
    double dummy_T0, dummy_dT;
    double dummy_Sig0, dummy_Swp0, dummy_Stvar0, dummy_dEps, dummy_D, dummy_BulkW,
        dummy_Sig, dummy_Swp, dummy_StVar;
    double* props = const_cast<double*>(&Props[0]); // cast out the const'ness to use the pointer value and let hope that the UDSM does not change any input parameter.
    UserMod(&IDTask, &ModelNumber, &IsUndr,
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

    if (need_determine_internal_params)
    {
        IDTask = 7; // this task number is not officially supported by Plaxis. User has to implement it in UDSM specifically.
                        // Otherwise, it does nothing

        Vector NewProps = Props;

        UserMod(&IDTask, &ModelNumber, &IsUndr,
                &iStep, &iTer, &Iel, &Int,
                &dummyX, &dummyY, &dummyZ, &dummy_T0, &dummy_dT,
                &NewProps[0],
                &dummy_Sig0, &dummy_Swp0, &dummy_Stvar0, &dummy_dEps, &dummy_D, &dummy_BulkW,
                &dummy_Sig, &dummy_Swp, &dummy_StVar, &dummy_ipl,
                &nStat,
                &dummy_nonsym, &dummy_iStrsDep, &dummy_iTimeDep, &dummy_iTang, &dummy_iPrjDir, &dummy_iPrjLen,
                &iAbort);

        const double diff = norm_2(NewProps - Props);
        if (diff > 1e-10)
        {
            KRATOS_WATCH(NewProps)
            KRATOS_WATCH(Props)

            KRATOS_ERROR << "The calculated material parameters are not the same as the given one, diff = " << diff;
        }
    }
}

void UDSM::ResetMaterial ( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    // reset the internal state
    mLastStrain.clear();
    mCurrentStrain.clear();
    mLastStateVariables.clear();
    mCurrentStateVariables.clear();
    mLastExcessPorePressure = 0.0;
    mCurrentExcessPorePressure = 0.0;

    noalias(mLastStress) = -mPrestressFactor*mPrestress;
    noalias(mCurrentStress) = mLastStress;

    /* initialize the internal variables */
    // This is done by calling IDTask = 1

    int IDTask = 1;

    mModelNumber = static_cast<int>(props[SOIL_MODEL_NUMBER]);
    mIsUndr = static_cast<int>(props[IS_UNDRAINED]);
    mProps = props[MATERIAL_PARAMETERS];

    mBulkW = 0.0;
    if (mIsUndr) mBulkW = props[BULK_W];

    Vector Sig0(20);
    noalias(Sig0) = ZeroVector(20);
    Sig0[0] = mLastStress[0];
    Sig0[1] = mLastStress[1];
    Sig0[2] = mLastStress[2];
    Sig0[3] = mLastStress[3];
    Sig0[4] = mLastStress[4];
    Sig0[5] = mLastStress[5];

    double Swp0 = mLastExcessPorePressure;

    int nStat = static_cast<int>(mLastStateVariables.size());
    Vector StVar0(nStat), StVar(nStat);
    StVar0.clear();
    StVar.clear();

    int iStep=0, iTer=0, Iel=mElemId, Int=mGaussId, dummy_ipl, dummy_nonsym,
        dummy_iStrsDep, dummy_iTimeDep, dummy_iTang, dummy_iPrjDir, dummy_iPrjLen;
    int iAbort;
    double dummyX, dummyY, dummyZ;
    double dummy_T0, dummy_dT;
    double dummy_dEps, dummy_D, dummy_BulkW,
        dummy_Sig, dummy_Swp, dummy_StVar;
    double* props_pointer = const_cast<double*>(&mProps[0]); // cast out the const'ness to use the pointer value and let hope that the UDSM does not change any input parameter.
    UserMod(&IDTask, &mModelNumber, &mIsUndr,
            &iStep, &iTer, &Iel, &Int,
            &dummyX, &dummyY, &dummyZ, &dummy_T0, &dummy_dT,
            props_pointer,
            &Sig0[0],
            &Swp0,
            &StVar0[0],
            &dummy_dEps, &dummy_D,
            &mBulkW,
            &dummy_Sig, &dummy_Swp,
            &StVar[0],
            &dummy_ipl,
            &nStat,
            &dummy_nonsym, &dummy_iStrsDep, &dummy_iTimeDep, &dummy_iTang, &dummy_iPrjDir, &dummy_iPrjLen,
            &iAbort);

    noalias(mLastStateVariables) = StVar0;
    noalias(mCurrentStateVariables) = StVar;
    // KRATOS_WATCH(mLastStateVariables)
    // KRATOS_WATCH(mCurrentStateVariables)
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
    KRATOS_ERROR << "Calling base class function";
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
    double omega = 1.0;
    if (props.Has(UDSM_RELAXATION_FACTOR))
        omega = props[UDSM_RELAXATION_FACTOR]; // relaxation factor

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

void UDSMImplex::StressIntegration ( const Vector& StrainVector, Vector& StressVector,
                                     const UDSM::Input& rInput, UDSM::Output& rOutput,
                                     const udsm_t UserFun, double omega, const UDSM::DebugInfo& pinfo,
                                     int ModelNumber, int IsUndr, double BulkW, const Vector& Props )
{
    int IDTask;
    Vector D(36);
    int iPl;
    int NonSym;
    int iStrsDep, iTang, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen; // can use this for debugging purpose (TODO)
    int nStat = rInput.LastStateVariables.size();
    double* props = const_cast<double*>(&Props[0]);

    // compute the vector of previous effective stress component
    Vector Sig0(20);
    noalias(Sig0) = ZeroVector(20);
    Sig0[0] = rInput.LastStress[0];
    Sig0[1] = rInput.LastStress[1];
    Sig0[2] = rInput.LastStress[2];
    Sig0[3] = rInput.LastStress[3];
    Sig0[4] = rInput.LastStress[4];
    Sig0[5] = rInput.LastStress[5];

    // compute the previous excess pore pressure
    double Swp0 = rInput.LastExcessPorePressure;

    // compute the previous value of state variables
    Vector StVar0 = rInput.LastStateVariables;

    // compute the vector of incremental strain
    Vector dEps(12);
    dEps[0] = StrainVector[0] - rInput.LastStrain[0];
    dEps[1] = StrainVector[1] - rInput.LastStrain[1];
    dEps[2] = StrainVector[2] - rInput.LastStrain[2];
    dEps[3] = StrainVector[3] - rInput.LastStrain[3];
    dEps[4] = StrainVector[4] - rInput.LastStrain[4];
    dEps[5] = StrainVector[5] - rInput.LastStrain[5];
    dEps[6] = rInput.LastStrain[0];
    dEps[7] = rInput.LastStrain[1];
    dEps[8] = rInput.LastStrain[2];
    dEps[9] = rInput.LastStrain[3];
    dEps[10] = rInput.LastStrain[4];
    dEps[11] = rInput.LastStrain[5];

    // initialize D with elastic matrix
    IDTask = 6;
    UserFun( &IDTask,
             &ModelNumber,
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             props,
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &rOutput.CurrentStress[0],
             &rOutput.CurrentExcessPorePressure,
             &rOutput.CurrentStateVariables[0],
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
    UserFun( &IDTask, &ModelNumber, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    UDSM::Vector1DToMatrix(D, rOutput.AlgorithmicTangent, NonSym);

    // compute delta_strain
    Vector delta_strain = StrainVector - rOutput.CurrentStrain;
//    KRATOS_WATCH(delta_strain)

    // compute equilibrium stress
    noalias(StressVector) = rOutput.CurrentStress + omega*prod(rOutput.AlgorithmicTangent, delta_strain);

    // calculate the stress
    IDTask = 2;
    UserFun( &IDTask,
             &ModelNumber,
             &IsUndr,
             NULL, // iStep
             NULL, // iter
             NULL, // Iel
             NULL, // Int
             NULL, NULL, NULL, // X, Y, Z
             NULL, // Time0
             NULL, // dTime
             props,
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &rOutput.CurrentStress[0],
             &rOutput.CurrentExcessPorePressure,
             &rOutput.CurrentStateVariables[0],
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
        KRATOS_ERROR << "Force calculation stop after calculating stress, iAbort = " << iAbort;

    // export the plastic state
    rOutput.PlasticState = iPl;

    // compute for the undrained state
    if(IsUndr != 0)
    {
        // TODO add the bulk stiffness of water to the material stiffness matrix
    }

    // save the current strain
    noalias(rOutput.CurrentStrain) = StrainVector;

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

    if (props.Has(USE_NUMERICAL_TANGENT))
    {
        if (props[USE_NUMERICAL_TANGENT])
        {
            if (props.Has(EPSILON))
                CalculateTangent = static_cast<int>(-std::ceil(std::log10(props[EPSILON])));
            if (CalculateTangent < 2)
                CalculateTangent = 2;
            if (CalculateTangent > 10)
                std::cout << "WARNING!!!The delta epsilon to compute numerical is too small" << std::endl;
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
            mModelNumber, mIsUndr, mBulkW, mProps, CalculateTangent );

    if (CalculateStresses)
    {
        // export the stress
        UDSM::Vector3DToVector(mCurrentStress, StressVector);
    }

    // KRATOS_WATCH(AlgorithmicTangent)
}

void UDSMImplicit::StressIntegration ( const UDSM::Input& rInput, UDSM::Output& rOutput,
                                       const udsm_t UserFun, const UDSM::DebugInfo& pinfo,
                                       int ModelNumber, int IsUndr, double BulkW, const Vector& Props,
                                       int CalculateTangent )
{
    int IDTask;
    Vector D(36);
    int iPl;
    int NonSym;
    int iStrsDep, iTang = 0, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen;
    int nStat = rInput.LastStateVariables.size();
    double dummyX, dummyY, dummyZ;
    double time0 = pinfo.time, dtime = pinfo.delta_time;
    int iStep = pinfo.step, iTer = pinfo.iteration,
        Iel = pinfo.element_id, Int = pinfo.point_id;
    double* props = const_cast<double*>(&Props[0]);

    // initialize the vector of previous effective stress component
    Vector Sig0(20);
    noalias(Sig0) = ZeroVector(20);
    noalias(subrange(Sig0, 0, 6)) = rInput.LastStress;

    // initialize the previous excess pore pressure
    double Swp0 = rInput.LastExcessPorePressure;

    // initialize the previous value of state variables
    Vector StVar0 = rInput.LastStateVariables;

    // KRATOS_WATCH(CurrentStrain)

    // initialize the vector of incremental strain
    Vector dEps(12);
    noalias(subrange(dEps, 0, 6)) = rOutput.CurrentStrain - rInput.LastStrain;
    noalias(subrange(dEps, 6, 12)) = rInput.LastStrain;

    // initialize D with elastic matrix
    IDTask = 6;
    UserFun( &IDTask,
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time0,        // Time0
             &dtime,        // dTime
             props,
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &rOutput.CurrentStress[0],
             &rOutput.CurrentExcessPorePressure,
             &rOutput.CurrentStateVariables[0],
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
    UserFun( &IDTask,
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time0,        // Time0
             &dtime,        // dTime
             props,
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &rOutput.CurrentStress[0],
             &rOutput.CurrentExcessPorePressure,
             &rOutput.CurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // KRATOS_WATCH(rOutput.CurrentStress)
    // KRATOS_WATCH(iPl)

    if (std::isnan(norm_2(rOutput.CurrentStress)))
    {
        KRATOS_ERROR << "NaN is detected at element " << Iel
                     << ", point " << Int
                     << ", initial stress: " << Sig0
                     << ", initial int vars: " << StVar0
                     << ", strain vector: " << dEps
                     << std::endl;
    }

    if(iAbort != 0)
        KRATOS_ERROR << "Force calculation stop after calculating stress, iAbort = " << iAbort;

    // export the plastic state
    rOutput.PlasticState = iPl;

    // obtain the stiffness matrix properties
    IDTask = 5;
    UserFun( &IDTask,
             &ModelNumber,
             &IsUndr,
             &iStep,        // iStep
             &iTer,         // iter
             &Iel,          // Iel
             &Int,          // Int
             &dummyX, &dummyY, &dummyZ, // X, Y, Z
             &time0,        // Time0
             &dtime,        // dTime
             props,
             &Sig0[0],
             &Swp0,
             &StVar0[0],
             &dEps[0],
             &D[0],
             &BulkW,
             &rOutput.CurrentStress[0],
             &rOutput.CurrentExcessPorePressure,
             &rOutput.CurrentStateVariables[0],
             &iPl,
             &nStat,
             &NonSym,
             &iStrsDep,
             &iTimeDep,
             &iTang,
             &iPrjDir,
             &iPrjLen,
             &iAbort );

    // KRATOS_WATCH(iTang)

    if (CalculateTangent)
    {
        if (CalculateTangent > 2)
        {
            // estimate the tangent using finite difference scheme
            const double epsilon = std::pow(10, -CalculateTangent);

            int aux = ComputeNumericalTangent(rInput, rOutput, rOutput.AlgorithmicTangent,
                    UserFun, pinfo, ModelNumber, IsUndr, BulkW, Props, epsilon);

            if (aux > 0)
            {
                std::cout << "WARNING!!!The numerical tangent is not valid at element " << Iel << ", point " << Int << std::endl;
            }
        }
        else
        {
            // calculate the material stiffness
            IDTask = 3;
            UserFun( &IDTask,
                     &ModelNumber,
                     &IsUndr,
                     &iStep,        // iStep
                     &iTer,         // iter
                     &Iel,          // Iel
                     &Int,          // Int
                     &dummyX, &dummyY, &dummyZ, // X, Y, Z
                     &time0,        // Time0
                     &dtime,        // dTime
                     props,
                     &Sig0[0],
                     &Swp0,
                     &StVar0[0],
                     &dEps[0],
                     &D[0],
                     &BulkW,
                     &rOutput.CurrentStress[0],
                     &rOutput.CurrentExcessPorePressure,
                     &rOutput.CurrentStateVariables[0],
                     &iPl,
                     &nStat,
                     &NonSym,
                     &iStrsDep,
                     &iTimeDep,
                     &iTang,
                     &iPrjDir,
                     &iPrjLen,
                     &iAbort );

            if(iAbort != 0)
                KRATOS_ERROR << "Force calculation stop after calculating material stiffness, iAbort = " << iAbort;

            UDSM::Vector1DToMatrix(D, rOutput.AlgorithmicTangent, NonSym);
        }

        // compute for the undrained state
        if(IsUndr != 0)
        {
            // TODO add the bulk stiffness of water to the material stiffness matrix
        }
    }

//    if(PlasticState)
//    {
       // KRATOS_WATCH(StressVector)
       // KRATOS_WATCH(AlgorithmicTangent)
//        KRATOS_WATCH(CurrentStateVariables)
//    }
}

int UDSMImplicit::ComputeNumericalTangent( const UDSM::Input& rInput, const UDSM::Output& rRefOutput,
                                           Matrix& Tangent,
                                           const udsm_t UserFun, const UDSM::DebugInfo& pinfo,
                                           int ModelNumber, int IsUndr, double BulkW, const Vector& Props,
                                           double epsilon )
{
    int IDTask;
    Vector D(36);
    int iPl;
    int NonSym;
    int iStrsDep, iTang = 0, iTimeDep, iAbort = 0;
    int iPrjDir, iPrjLen;
    int nStat = rInput.LastStateVariables.size();
    double dummyX, dummyY, dummyZ;
    double time0 = pinfo.time, dtime = pinfo.delta_time;
    int iStep = pinfo.step, iTer = pinfo.iteration,
        Iel = pinfo.element_id, Int = pinfo.point_id;
    double* props = const_cast<double*>(&Props[0]);

    Vector Sig0(20), dEps(12), Sig(6);
    double Swp0, Swp;
    Vector StVar0(nStat), StVar(nStat);
    int iPl_diff_num = 0;

    for (unsigned int i = 0; i < 6; ++i)
    {
        // initialize the vector of previous effective stress component
        noalias(Sig0) = ZeroVector(20);
        noalias(subrange(Sig0, 0, 6)) = rInput.LastStress;

        // initialize the previous excess pore pressure
        Swp0 = rInput.LastExcessPorePressure;

        // initialize the previous value of state variables
        noalias(StVar0) = rInput.LastStateVariables;

        // initialize the vector of incremental strain
        noalias(subrange(dEps, 0, 6)) = rRefOutput.CurrentStrain - rInput.LastStrain;
        noalias(subrange(dEps, 6, 12)) = rInput.LastStrain;
        dEps(i) += epsilon;

        // initialize the stress and other state values
        noalias(Sig) = rRefOutput.CurrentStress;
        noalias(StVar) = rRefOutput.CurrentStateVariables;
        Swp = rRefOutput.CurrentExcessPorePressure;

        // initialize D with elastic matrix
        IDTask = 6;
        UserFun( &IDTask,
                 &ModelNumber,
                 &IsUndr,
                 &iStep,        // iStep
                 &iTer,         // iter
                 &Iel,          // Iel
                 &Int,          // Int
                 &dummyX, &dummyY, &dummyZ, // X, Y, Z
                 &time0,        // Time0
                 &dtime,        // dTime
                 props,
                 &Sig0[0],
                 &Swp0,
                 &StVar0[0],
                 &dEps[0],
                 &D[0],
                 &BulkW,
                 &Sig[0],
                 &Swp,
                 &StVar[0],
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
        UserFun( &IDTask,
                 &ModelNumber,
                 &IsUndr,
                 &iStep,        // iStep
                 &iTer,         // iter
                 &Iel,          // Iel
                 &Int,          // Int
                 &dummyX, &dummyY, &dummyZ, // X, Y, Z
                 &time0,        // Time0
                 &dtime,        // dTime
                 props,
                 &Sig0[0],
                 &Swp0,
                 &StVar0[0],
                 &dEps[0],
                 &D[0],
                 &BulkW,
                 &Sig[0],
                 &Swp,
                 &StVar[0],
                 &iPl,
                 &nStat,
                 &NonSym,
                 &iStrsDep,
                 &iTimeDep,
                 &iTang,
                 &iPrjDir,
                 &iPrjLen,
                 &iAbort );

        // estimate the tangent
        noalias( column(Tangent, i) ) = (Sig - rRefOutput.CurrentStress) / epsilon;

        // check the plastic state
        if (iPl != rRefOutput.PlasticState)
            ++iPl_diff_num;
    }

    return iPl_diff_num;
}

} // Namespace Kratos
