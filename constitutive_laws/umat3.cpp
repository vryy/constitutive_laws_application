#include <iostream>
#include <dlfcn.h>

#include "includes/define.h"
#include "constitutive_laws/umat3.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "constitutive_laws_application.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

// REF: https://code.google.com/p/parafem/wiki/UMATs

typedef Kratos::ConstitutiveLaw::GeometryType GeometryType;

const int Umat3::A2K[] = {0, 1, 2, 3, 5, 4};
const int Umat3::K2A[] = {0, 1, 2, 3, 5, 4};
const int Umat3::PS[]  = {0, 1, 3};

#ifdef KRATOS_UMAT_LIBRARY_IS_PROVIDED
extern "C" void umat_( double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                       double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                       double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                       char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                       double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                       double** DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC );
#endif

Umat3::Umat3()
{
}

Umat3::~Umat3()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    if(mp_umat_handle != 0)
        dlclose(mp_umat_handle);
    #endif
}

int Umat3::Check( const Kratos::Properties& props, const GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo )
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    if(props.Has( ABAQUS_LIBRARY_NAME ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define ABAQUS_LIBRARY_NAME", "")
    }
    if(props.Has( UMAT_NAME ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define UMAT_NAME", "")
    }
    #endif
    if(props.Has( UMAT_NDI ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define UMAT_NDI", "")
    }
    if(props.Has( UMAT_NSHR ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define UMAT_NSHR", "")
    }
    if(props.Has( UMAT_NSTATV ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define UMAT_NSTATV", "")
    }
    if(props.Has( UMAT_CMNAME ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define UMAT_CMNAME", "")
    }
    if(props.Has( MATERIAL_PARAMETERS ) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties must define MATERIAL_PARAMETERS", "")
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

    mCurrentStress.resize(NTENS);
    mOldStress.resize(NTENS);
    mPrestress.resize(NTENS);
    mCurrentStrain.resize(NTENS);
    mOldStrain.resize(NTENS);
    mCurrentStateVariables.resize(NSTATV);
    mOldStateVariables.resize(NSTATV);
    mPrestressFactor = 1.0;
    mOldStressZZ = 0.0;
    mCurrentStressZZ = 0.0;

    noalias(mCurrentStress) = ZeroVector(NTENS);
    noalias(mOldStress) = ZeroVector(NTENS);
    noalias(mPrestress) = ZeroVector(NTENS);
    noalias(mCurrentStrain) = ZeroVector(NTENS);
    noalias(mOldStrain) = ZeroVector(NTENS);
    noalias(mCurrentStateVariables) = ZeroVector(NSTATV);
    noalias(mOldStateVariables) = ZeroVector(NSTATV);

    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    // get the library name and load the udsm subroutine
    std::string lib_name = props[ABAQUS_LIBRARY_NAME];
//    mp_umat_handle = dlopen(lib_name.c_str(), RTLD_LAZY);
    mp_umat_handle = dlopen(lib_name.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if(mp_umat_handle == 0)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "The Abaqus material library does not exist:", lib_name)
    }

    std::string umat_name = props[UMAT_NAME];
    char* error;
    Umat = (void (*)(double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                double** DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC))
           dlsym(mp_umat_handle, umat_name.c_str());
    error = dlerror();
    if(error != NULL)
    {
        std::stringstream ss;
        ss << "Error loading subroutine " << umat_name << " in the " << lib_name << " library, error message = " << error;
        KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
    }
    #endif
}

void Umat3::ResetMaterial( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];

    noalias(mCurrentStress) = ZeroVector(NTENS);
    noalias(mOldStress) = ZeroVector(NTENS);
    noalias(mPrestress) = ZeroVector(NTENS);
    noalias(mCurrentStrain) = ZeroVector(NTENS);
    noalias(mOldStrain) = ZeroVector(NTENS);
    noalias(mCurrentStateVariables) = ZeroVector(NSTATV);
    noalias(mOldStateVariables) = ZeroVector(NSTATV);
    mPrestressFactor = 1.0;
    mOldStressZZ = 0.0;
    mCurrentStressZZ = 0.0;
}

void Umat3::InitializeSolutionStep( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{
}

void Umat3::InitializeNonLinearIteration ( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo )
{
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
    Vector PROPS = props[MATERIAL_PARAMETERS];
    int NPROPS = PROPS.size();

    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];
    char* cmname = (char*)(props[UMAT_CMNAME].c_str());

    // TODO fix the variable size array issue
    double STRAN[NTENS];
    double STRES[NTENS];
    double DSTRAN[NTENS];
    double DDSDDE[NTENS][NTENS];
    double STATEV[NSTATV];
    double DDSDDT[NTENS];
    double DRPLDE[NTENS];
    double DRPLDT[NTENS];
    double DFGRD0[3][3];
    double DFGRD1[3][3];
    double TIM[2];
    double DTIM;
    double SSE, SPD, SCD, RPL;
    double TEMP, DTEMP;
    int KSTEP, KINC;

    if( (NDI == 3) && (NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            STRAN[i] = StrainVector[K2A[i]];
            DSTRAN[i] = StrainVector[K2A[i]] - mOldStrain[K2A[i]];
            STRES[i] = mOldStress[K2A[i]] - mPrestressFactor*mPrestress[K2A[i]];
            for(int j = 0; j < NTENS; ++j)
                DDSDDE[i][j] = 0.0;
        }
    }
    else if( (NDI == 2) & (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            STRAN[i] = StrainVector[i];
            DSTRAN[i] = StrainVector[i] - mOldStrain[i];
            STRES[i] = mOldStress[i] - mPrestressFactor*mPrestress[i];
            for(int j = 0; j < NTENS; ++j)
                DDSDDE[i][j] = 0.0;
        }
    }
    else if( (NDI == 3) & (NSHR == 1) ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
    {
        for(int i = 0; i < 3; ++i)
        {
            STRAN[PS[i]] = StrainVector[i];
            DSTRAN[PS[i]] = StrainVector[i] - mOldStrain[i];
            STRES[PS[i]] = mOldStress[i] - mPrestressFactor*mPrestress[i];
        }
        STRAN[2] = 0.0;
        DSTRAN[2] = 0.0;
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
                DDSDDE[i][j] = 0.0;
        STRES[2] = mOldStressZZ;
    }

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            DFGRD0[i][j] = 0.0;
            DFGRD1[i][j] = DeformationGradient(i, j);
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

    #ifdef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    umat_( STRES,
          STATEV,
          (double**)DDSDDE,
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT, DRPLDE, DRPLDT,
          STRAN, DSTRAN,
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          cmname,
          &NDI, &NSHR, &NTENS, &NSTATV, &PROPS[0], &NPROPS,
          NULL, NULL, NULL, NULL,
          (double**)DFGRD0, (double**)DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );
    #else
    Umat( STRES,
          STATEV,
          (double**)DDSDDE,
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT, DRPLDE, DRPLDT,
          STRAN, DSTRAN,
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          cmname,
          &NDI, &NSHR, &NTENS, &NSTATV, &PROPS[0], &NPROPS,
          NULL, NULL, NULL, NULL,
          (double**)DFGRD0, (double**)DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );
    #endif

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
        for(int i = 0; i < 3; ++i)
        {
            mCurrentStress[i] = STRES[PS[i]];
        }
        mCurrentStressZZ = STRES[2];
    }

    noalias(mCurrentStress) += mPrestressFactor*mPrestress;
    noalias(StressVector) = mCurrentStress; // TODO: take into account the temperature

    for(int i = 0; i < NSTATV; ++i)
    {
        mCurrentStateVariables[i] = STATEV[i];
    }

    if( (NDI == 3) && (NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            for(int j = 0; j < NTENS; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[A2K[j]][A2K[i]];
            }
        }
    }
    else if( (NDI == 2) && (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            for(int j = 0; j < NTENS; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[j][i];
            }
        }
    }
    else if( (NDI == 3) && (NSHR == 1) ) // 2D case    [o_xx  o_yy  o_zz  o_xy]
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[PS[j]][PS[i]];
            }
        }
    }
//    KRATOS_WATCH(AlgorithmicTangent)
    // TODO: export the variable SSE, SPD, SCD, RPL
}

void Umat3::FinalizeNonLinearIteration ( const Properties& props,
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
