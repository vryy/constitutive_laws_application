#include "includes/define.h"
#include "constitutive_laws/umat3.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/ublas_interface.h"
#include "structural_application/structural_application_variables.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
#include "custom_utilities/shared_library_handle_shield.h"

namespace Kratos
{

// REF: https://code.google.com/p/parafem/wiki/UMATs

typedef Kratos::ConstitutiveLaw::GeometryType GeometryType;

const int Umat3::A2K[] = {0, 1, 2, 3, 5, 4};
const int Umat3::K2A[] = {0, 1, 2, 3, 5, 4};
const int Umat3::PS[]  = {0, 1, 3};

#ifdef KRATOS_UMAT_LIBRARY_IS_PROVIDED
extern "C" void umat_( double* STRESS, double* STATEV, double* DDSDDE, double* SSE, double* SPD, double* SCD,
                       double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                       double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                       char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                       double* COORDS, double* DROT, double* PNEWDT, double* CELENT, double* DFGRD0,
                       double* DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC );
#endif

Umat3::Umat3()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    mp_umat_handle = nullptr;
    #endif
    Umat = nullptr;
}

Umat3::~Umat3()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    mp_umat_handle = nullptr;
    #endif
    Umat = nullptr;
}

int Umat3::Check( const Kratos::Properties& props, const GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo ) const
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    if(props.Has( ABAQUS_LIBRARY_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define ABAQUS_LIBRARY_NAME";
    }
    if(props.Has( UMAT_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_NAME";
    }
    #endif
    if(props.Has( UMAT_NDI ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_NDI";
    }
    if(props.Has( UMAT_NSHR ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_NSHR";
    }
    if(props.Has( UMAT_NSTATV ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_NSTATV";
    }
    if(props.Has( UMAT_CMNAME ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_CMNAME";
    }
    if(props.Has( MATERIAL_PARAMETERS ) == false)
    {
        KRATOS_ERROR << "Properties must define MATERIAL_PARAMETERS";
    }
    return 0;
}

bool Umat3::Has( const Variable<int>& rThisVariable )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        return true;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        return true;
    return false;
}

bool Umat3::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Umat3::Has( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
        return true;
    if(rThisVariable == STRAIN)
        return true;
    return false;
}

bool Umat3::Has( const Variable<Matrix>& rThisVariable )
{
    if(rThisVariable == CAUCHY_STRESS_TENSOR)
        return true;
    return false;
}

int& Umat3::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& Umat3::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == PRESTRESS_FACTOR)
    {
        rValue = mPrestressFactor;
    }
    else if (rThisVariable == PRESSURE_P)
    {
        if( mOldStress.size() == 6 ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            rValue = -1.0/3 * (mOldStress[0] + mOldStress[1] + mOldStress[2]);
        }
        else if( mOldStress.size() == 3 ) // 2D case    [o_xx  o_yy  o_xy]
        {
            rValue = -1.0/3 * (mOldStress[0] + mOldStress[1]);
        }
        else if( mOldStress.size() == 4 ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            rValue = -1.0/3 * (mOldStress[0] + mOldStress[1] + mOldStressZZ);
        }
    }
    else if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double p, sxx, syy, szz, sxy, sxz, syz;
        if( mOldStress.size() == 6 ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            p = 1.0/3 * (mOldStress[0] + mOldStress[1] + mOldStress[2]);
            sxx = mOldStress[0] - p;
            syy = mOldStress[1] - p;
            szz = mOldStress[2] - p;
            sxy = mOldStress[3];
            syz = mOldStress[4];
            sxz = mOldStress[5];
        }
        else if( mOldStress.size() == 3 ) // 2D case    [o_xx  o_yy  o_xy]
        {
            p = 1.0/3 * (mOldStress[0] + mOldStress[1]);
            sxx = mOldStress[0] - p;
            syy = mOldStress[1] - p;
            szz = -p;
            sxy = mOldStress[2];
            syz = 0.0;
            sxz = 0.0;
        }
        else if( mOldStress.size() == 4 ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            p = 1.0/3 * (mOldStress[0] + mOldStress[1] + mOldStressZZ);
            sxx = mOldStress[0] - p;
            syy = mOldStress[1] - p;
            szz = mOldStressZZ - p;
            sxy = mOldStress[2];
            syz = 0.0;
            sxz = 0.0;
        }

        rValue = sqrt( 3.0 * ( 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + syz*syz + sxz*sxz ) );
    }
    return rValue;
}

Vector& Umat3::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable == INTERNAL_VARIABLES)
    {
        if(rValue.size() != mCurrentStateVariables.size())
            rValue.resize(mCurrentStateVariables.size(), false);
        noalias(rValue) = mCurrentStateVariables;
    }
    if(rThisVariable == STRESSES)
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);

        if( mOldStress.size() == 6 ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            noalias(rValue) = mOldStress;
        }
        else if( mOldStress.size() == 3 ) // 2D case    [o_xx  o_yy  o_xy]
        {
            rValue[0] = mOldStress[0];
            rValue[1] = mOldStress[1];
            rValue[2] = 0.0;
            rValue[3] = mOldStress[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
        else if( mOldStress.size() == 4 ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            rValue[0] = mOldStress[0];
            rValue[1] = mOldStress[1];
            rValue[2] = mOldStressZZ;
            rValue[3] = mOldStress[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
    }
    if(rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS)
    {
        if(rValue.size() != mPrestress.size())
            rValue.resize(mPrestress.size(), false);
        noalias(rValue) = mPrestress;
    }
    if(rThisVariable == STRAIN)
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);

        if( mOldStrain.size() == 6 ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            noalias(rValue) = mOldStrain;
        }
        else if( mOldStrain.size() == 3 ) // 2D case    [o_xx  o_yy  o_xy]
        {
            rValue[0] = mOldStrain[0];
            rValue[1] = mOldStrain[1];
            rValue[2] = 0.0; // TODO is it zero?
            rValue[3] = mOldStrain[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
        else if( mOldStrain.size() == 4 ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            rValue[0] = mOldStrain[0];
            rValue[1] = mOldStrain[1];
            rValue[2] = 0.0;
            rValue[3] = mOldStrain[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
    }

    return rValue;
}

Matrix& Umat3::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == CAUCHY_STRESS_TENSOR)
    {
        if (rValue.size1() != 3 || rValue.size2() != 3)
            rValue.resize(3, 3, false);
        SD_MathUtils<double>::StressVectorToTensor(mOldStress, rValue);
        return rValue;
    }
    return rValue;
}

void Umat3::SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PARENT_ELEMENT_ID)
        mElementId = rValue;
    if(rVariable == INTEGRATION_POINT_INDEX )
        mIntPointIndex = rValue;
}

void Umat3::SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
}

void Umat3::SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        if (mPrestress.size() != rValue.size())
            mPrestress.resize(rValue.size());
        noalias(mPrestress) = rValue;
    }
    if ( rVariable == UMAT_STATEV )
    {
        noalias(mCurrentStateVariables) = rValue;
    }
}

void Umat3::SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{
}

void Umat3::InitializeMaterial( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues )
{
    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];

    unsigned int strain_size;
    if (NDI == 3 && NSHR == 3)
    {
        strain_size = 6;
    }
    else if (NDI == 3 && NSHR == 1)
    {
        bool is_axisymmetric = false;
        if (props.Has(IS_AXISYMMETRIC))
            if (props[IS_AXISYMMETRIC])
                is_axisymmetric = true;

        if (is_axisymmetric)
            strain_size = 4;
        else
            strain_size = 3;
    }
    else if (NDI == 2 && NSHR == 1)
    {
        strain_size = 3;
    }

    mCurrentStress.resize(strain_size);
    mOldStress.resize(strain_size);
    mPrestress.resize(strain_size);
    mCurrentStrain.resize(strain_size);
    mOldStrain.resize(strain_size);
    mCurrentStateVariables.resize(NSTATV);
    mOldStateVariables.resize(NSTATV);
    mPrestressFactor = 1.0;
    mOldStressZZ = 0.0;
    mCurrentStressZZ = 0.0;

    noalias(mCurrentStress) = ZeroVector(strain_size);
    noalias(mOldStress) = ZeroVector(strain_size);
    noalias(mPrestress) = ZeroVector(strain_size);
    noalias(mCurrentStrain) = ZeroVector(strain_size);
    noalias(mOldStrain) = ZeroVector(strain_size);
    noalias(mCurrentStateVariables) = ZeroVector(NSTATV);
    noalias(mOldStateVariables) = ZeroVector(NSTATV);

    mStep = 0;
    mIncrement = 0;

    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    {
        // get the library name and load the udsm subroutine
        std::string lib_name = props[ABAQUS_LIBRARY_NAME];
        mp_umat_handle = DLL::GetSharedLibraryHandle(lib_name);
        if(mp_umat_handle == nullptr)
        {
            KRATOS_ERROR << "Error loading Abaqus material library " << lib_name;
        }

        std::string umat_name = props[UMAT_NAME];
        Umat = (umat_t) DLL::GetSymbol(mp_umat_handle, umat_name);
        auto error = DLL::GetError();
        if(error)
        {
            KRATOS_ERROR << "Error loading subroutine " << umat_name << " in the " << lib_name << " library"
                         << ", error message = " << DLL::GetErrorMessage(error);
        }
        else
        {
            std::cout << "Successfully load Umat from " << lib_name << std::endl;
        }
    }
    #endif
}

void Umat3::ResetMaterial( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    if (mPrestress.size() == 6)
    {
        if (mOldStress.size() == 6)
        {
            noalias(mOldStress) = -mPrestressFactor*mPrestress;
        }
        else if (mOldStress.size() == 3)
        {
            mOldStress[0] = -mPrestressFactor*mPrestress[0];
            mOldStress[1] = -mPrestressFactor*mPrestress[1];
            mOldStress[2] = -mPrestressFactor*mPrestress[3];
            mOldStressZZ = -mPrestressFactor*mPrestress[2];
        }
    }
    else if (mPrestress.size() == 3)
    {
        if (mOldStress.size() == 6)
        {
            mOldStress[0] = -mPrestressFactor*mPrestress[0];
            mOldStress[1] = -mPrestressFactor*mPrestress[1];
            mOldStress[2] = 0.0;
            mOldStress[3] = -mPrestressFactor*mPrestress[2];
            mOldStress[4] = 0.0;
            mOldStress[5] = 0.0;
            mOldStressZZ = 0.0;
        }
        else if (mOldStress.size() == 3)
        {
            noalias(mOldStress) = -mPrestressFactor*mPrestress;
            mOldStressZZ = 0.0;
        }
    }

    noalias(mCurrentStress) = mOldStress;
    mCurrentStressZZ = mOldStressZZ;
    mCurrentStrain.clear();
    mOldStrain.clear();
    mCurrentStateVariables.clear();
    mOldStateVariables.clear();
}

void Umat3::InitializeSolutionStep( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{
    ++mStep;
    mIncrement = 0;
}

void Umat3::InitializeNonLinearIteration ( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo )
{
    ++mIncrement;
}

void Umat3::CalculateMaterialResponseCauchy( Parameters& parameters )
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
        , parameters.IsSetStrainVector()
        , parameters.IsSetDeformationGradientF()
    );
}

void Umat3::CalculateMaterialResponse( const Vector& StrainVector,
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
    this->CalculateMaterialResponse( StrainVector, DeformationGradient,
            StressVector, AlgorithmicTangent, CurrentProcessInfo, props, geom,
            ShapeFunctionsValues, CalculateStresses, CalculateTangent,
            SaveInternalVariables, true, false );
}

void Umat3::CalculateMaterialResponse( const Vector& StrainVector,
                                       const Matrix& DeformationGradient,
                                       Vector& StressVector,
                                       Matrix& AlgorithmicTangent,
                                       const ProcessInfo& CurrentProcessInfo,
                                       const Properties& props,
                                       const GeometryType& geom,
                                       const Vector& ShapeFunctionsValues,
                                       bool CalculateStresses,
                                       int CalculateTangent,
                                       bool SaveInternalVariables,
                                       bool IsSetStrain,
                                       bool IsSetDeformationGradient )
{
    if (CurrentProcessInfo[SET_CALCULATE_REACTION])
    {
        if (CalculateStresses)
        {
            noalias(StressVector) = mCurrentStress;
            return;
        }
    }

    Vector PROPS = props[MATERIAL_PARAMETERS];
    int NPROPS = PROPS.size();

    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];
    char* cmname = (char*)(props[UMAT_CMNAME].c_str());

    // TODO fix the variable size array issue
    std::vector<double> STRAN(NTENS);
    std::vector<double> STRES(NTENS);
    std::vector<double> DSTRAN(NTENS);
    std::vector<double> DDSDDE(NTENS*NTENS);
    std::vector<double> STATEV(NSTATV);
    std::vector<double> DDSDDT(NTENS);
    std::vector<double> DRPLDE(NTENS);
    std::vector<double> DRPLDT(NTENS);
    double DFGRD0[9];
    double DFGRD1[9];
    double DROT[9];
    double TIM[2];
    double DTIM;
    double SSE, SPD, SCD, RPL;
    double TEMP, DTEMP;
    int KSTEP = (int) mStep, KINC = (int) mIncrement;

    bool is_axisymmetric = false;

    if (IsSetStrain)
    {
        if( (NDI == 3) && (NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                STRAN[i] = StrainVector[K2A[i]];
                DSTRAN[i] = StrainVector[K2A[i]] - mOldStrain[K2A[i]];
                STRES[i] = mOldStress[K2A[i]] - mPrestressFactor*mPrestress[K2A[i]];
                for(int j = 0; j < NTENS; ++j)
                    DDSDDE[i*NTENS + j] = 0.0;
            }
        }
        else if( (NDI == 2) && (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                STRAN[i] = StrainVector[i];
                DSTRAN[i] = StrainVector[i] - mOldStrain[i];
                STRES[i] = mOldStress[i];
                for(int j = 0; j < NTENS; ++j)
                    DDSDDE[i*NTENS + j] = 0.0;
            }
        }
        else if( (NDI == 3) && (NSHR == 1) )
        {
            if (props.Has(IS_AXISYMMETRIC))
                if (props[IS_AXISYMMETRIC])
                    is_axisymmetric = true;

            if (is_axisymmetric)  // axisymmetric  Abaqus [o_xx  o_yy  o_zz  o_xy] <=> Kratos [o_xx  o_yy  o_xy  o_zz]
            {
                for(int i = 0; i < 3; ++i)
                {
                    STRAN[PS[i]] = StrainVector[i];
                    DSTRAN[PS[i]] = StrainVector[i] - mOldStrain[i];
                    STRES[PS[i]] = mOldStress[i];
                }
                STRAN[2] = StrainVector[3];
                DSTRAN[2] = StrainVector[3] - mOldStrain[3];
                for(int i = 0; i < 4; ++i)
                    for(int j = 0; j < 4; ++j)
                        DDSDDE[i*NTENS + j] = 0.0;
                STRES[2] = mOldStress[3];
            }
            else // 2D case, plane strain   Abaqus [o_xx  o_yy  o_zz  o_xy] <=> Kratos [o_xx  o_yy  o_xy]
            {
                for(int i = 0; i < 3; ++i)
                {
                    STRAN[PS[i]] = StrainVector[i];
                    DSTRAN[PS[i]] = StrainVector[i] - mOldStrain[i];
                    STRES[PS[i]] = mOldStress[i];
                }
                STRAN[2] = 0.0;
                DSTRAN[2] = 0.0;
                for(int i = 0; i < 4; ++i)
                    for(int j = 0; j < 4; ++j)
                        DDSDDE[i*NTENS + j] = 0.0;
                STRES[2] = mOldStressZZ;
            }
        }
    }

    if (IsSetDeformationGradient)
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                DFGRD0[i + 3*j] = 0.0;
                DFGRD1[i + 3*j] = DeformationGradient(i, j);
            }
        }
    }

    for(int i = 0; i < NSTATV; ++i)
    {
        STATEV[i] = mOldStateVariables[i];
    }

    TIM[0] = CurrentProcessInfo[TIME];
    TIM[1] = CurrentProcessInfo[TIME];
    DTIM = CurrentProcessInfo[DELTA_TIME];

    TEMP = 0.0;
    DTEMP = 0.0; // TODO: take into account the temperature

    // identity rotation matrix
    DROT[0] = 1.0; DROT[3] = 0.0; DROT[6] = 0.0;
    DROT[1] = 0.0; DROT[4] = 1.0; DROT[7] = 0.0;
    DROT[2] = 0.0; DROT[5] = 0.0; DROT[8] = 1.0;

    #ifdef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    umat_( STRES,
          STATEV,
          DDSDDE.data(),
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT.data(), DRPLDE.data(), DRPLDT.data(),
          STRAN.data(), DSTRAN.data(),
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          cmname,
          &NDI, &NSHR, &NTENS, &NSTATV, &PROPS[0], &NPROPS,
          NULL,
          DROT,
          NULL, NULL,
          DFGRD0, DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );
    #else
    Umat( STRES.data(),
          STATEV.data(),
          DDSDDE.data(),
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT.data(), DRPLDE.data(), DRPLDT.data(),
          STRAN.data(), DSTRAN.data(),
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          cmname,
          &NDI, &NSHR, &NTENS, &NSTATV, &PROPS[0], &NPROPS,
          NULL,
          DROT,
          NULL, NULL,
          DFGRD0, DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );
    #endif

    // KRATOS_WATCH(StrainVector)
    noalias(mCurrentStrain) = StrainVector;

    if(NDI == 3 && NSHR == 3)
    {
        for(int i = 0; i < NTENS; ++i)
        {
            mCurrentStress[i] = STRES[A2K[i]];
        }
    }
    else if(NDI == 2 && NSHR == 1)
    {
        for(int i = 0; i < NTENS; ++i)
        {
            mCurrentStress[i] = STRES[i];
        }
    }
    else if(NDI == 3 && NSHR == 1)
    {
        if (is_axisymmetric)
        {
            for(int i = 0; i < 3; ++i)
            {
                mCurrentStress[i] = STRES[PS[i]];
            }
            mCurrentStress[3] = STRES[2];
            mCurrentStressZZ = STRES[2];
        }
        else
        {
            for(int i = 0; i < 3; ++i)
            {
                mCurrentStress[i] = STRES[PS[i]];
            }
            mCurrentStressZZ = STRES[2];
        }
    }

    if(CalculateStresses)
        noalias(StressVector) = mCurrentStress; // TODO: take into account the temperature

    for(int i = 0; i < NSTATV; ++i)
    {
        mCurrentStateVariables[i] = STATEV[i];
    }

    if (CalculateTangent)
    {
        if( (NDI == 3) && (NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                for(int j = 0; j < NTENS; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[A2K[j]*NTENS + A2K[i]];
                }
            }
        }
        else if( (NDI == 2) && (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                for(int j = 0; j < NTENS; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[j*NTENS + i];
                }
            }
        }
        else if( (NDI == 3) && (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_zz  o_xy]
        {
            for(int i = 0; i < 3; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[PS[j]*NTENS + PS[i]];
                }
            }
        }
    //    KRATOS_WATCH(AlgorithmicTangent)
        // TODO: export the variable SSE, SPD, SCD, RPL
    }

    // std::cout << "Stress and tangent at element " << mElementId << ", point " << mIntPointIndex << std::endl;
    // KRATOS_WATCH(StressVector)
    // KRATOS_WATCH(AlgorithmicTangent)
}

void Umat3::FinalizeNonLinearIteration( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues,
                                        const ProcessInfo& CurrentProcessInfo )
{
}

void Umat3::FinalizeSolutionStep( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ,
                                  const ProcessInfo& CurrentProcessInfo )
{
    noalias(mOldStrain) = mCurrentStrain;
    noalias(mOldStress) = mCurrentStress;
    mOldStressZZ = mCurrentStressZZ;
    noalias(mOldStateVariables) = mCurrentStateVariables;
}

} // Namespace Kratos
