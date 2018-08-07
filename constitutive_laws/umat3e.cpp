#include <iostream>
#include <dlfcn.h>

#include "includes/define.h"
#include "constitutive_laws/umat3e.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "constitutive_laws_application.h"
#include "includes/ublas_interface.h"

#define DEBUG_UMAT
#define DEBUG_ELEMENT_ID    1
#define DEBUG_POINT_ID      0

namespace Kratos
{

// REF: https://code.google.com/p/parafem/wiki/UMATs

typedef Kratos::ConstitutiveLaw::GeometryType GeometryType;

const int Umat3e::A2K[] = {0, 1, 2, 3, 5, 4};
const int Umat3e::K2A[] = {0, 1, 2, 3, 5, 4};
const int Umat3e::PS[]  = {0, 1, 3};

#ifdef KRATOS_UMAT_LIBRARY_IS_PROVIDED
extern "C" void umat_( double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                       double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                       double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                       char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                       double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                       double** DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC );
#endif

Umat3e::Umat3e()
{
}

Umat3e::~Umat3e()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    if(mp_umat_handle != 0)
        dlclose(mp_umat_handle);
    #endif
}

int Umat3e::Check( const Kratos::Properties& props, const GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo )
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
    if((props.Has( UMAT_NDI ) == false) && (Has(UMAT_NDI) == false))
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties or ConstitutiveLaw must define UMAT_NDI", "")
    }
    if((props.Has( UMAT_NSHR ) == false) && (Has(UMAT_NSHR) == false))
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties or ConstitutiveLaw must define UMAT_NSHR", "")
    }
    if((props.Has( UMAT_NSTATV ) == false) && (Has(UMAT_NSTATV) == false))
    {
        KRATOS_THROW_ERROR(std::logic_error, "Properties or ConstitutiveLaw must define UMAT_NSTATV", "")
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

bool Umat3e::Has( const Variable<int>& rThisVariable )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        return true;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        return true;
    if(rThisVariable == UMAT_NDI)
        return true;
    if(rThisVariable == UMAT_NSHR)
        return true;
    if(rThisVariable == UMAT_NSTATV)
        return true;
    return false;
}

bool Umat3e::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Umat3e::Has( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
        return true;
    if(rThisVariable == STRAIN)
        return true;
    if(rThisVariable == MATERIAL_PARAMETERS)
        return true;
    return false;
}

bool Umat3e::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

int& Umat3e::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& Umat3e::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == PRESTRESS_FACTOR)
    {
        rValue = mPrestressFactor;
    }
    return rValue;
}

Vector& Umat3e::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if(rThisVariable == INTERNAL_VARIABLES)
    {
        if(rValue.size() != mOldStateVariables.size())
            rValue.resize(mOldStateVariables.size(), false);
        noalias(rValue) = mOldStateVariables;
    }
    if(rThisVariable == STRESSES)
    {
        if(rValue.size() != 6)
            rValue.resize(6, false);

        if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            noalias(rValue) = mOldStress;
        }
        else if( (mNDI == 2) & (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
        {
            rValue[0] = mOldStress[0];
            rValue[1] = mOldStress[1];
            rValue[2] = 0.0;
            rValue[3] = mOldStress[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
        else if( (mNDI == 3) & (mNSHR == 1) ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
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

        if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            noalias(rValue) = mOldStrain;
        }
        else if( (mNDI == 2) & (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
        {
            rValue[0] = mOldStrain[0];
            rValue[1] = mOldStrain[1];
            rValue[2] = 0.0; // TODO is it zero?
            rValue[3] = mOldStrain[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
        else if( (mNDI == 3) & (mNSHR == 1) ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            rValue[0] = mOldStrain[0];
            rValue[1] = mOldStrain[1];
            rValue[2] = 0.0;
            rValue[3] = mOldStrain[2];
            rValue[4] = 0.0;
            rValue[5] = 0.0;
        }
    }
    if(rThisVariable == MATERIAL_PARAMETERS)
    {
        if(rValue.size() != mPROPS.size())
            rValue.resize(mPROPS.size(), false);
        noalias(rValue) = mPROPS;
    }

    return rValue;
}

Matrix& Umat3e::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

std::string& Umat3e::GetValue( const Variable<std::string>& rThisVariable, std::string& rValue )
{
    if(rThisVariable == UMAT_CMNAME)
        rValue = mCMNAME;
    return rValue;
}

void Umat3e::SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PARENT_ELEMENT_ID)
        mElementId = rValue;
    if(rVariable == INTEGRATION_POINT_INDEX )
        mIntPointIndex = rValue;
    if(rVariable == UMAT_NDI )
        mNDI = rValue;
    if(rVariable == UMAT_NSHR )
        mNSHR = rValue;
    if(rVariable == UMAT_NSTATV )
        mNSTATV = rValue;
}

void Umat3e::SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
}

void Umat3e::SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        if(mPrestress.size() != rValue.size())
            mPrestress.resize(rValue.size(), false);
        noalias(mPrestress) = rValue;

        ResetState(); // WARNING: we must set PRESTRESS_FACTOR before setting PRESTRESS
    }
    if ( rVariable == UMAT_STATEV )
    {
        if(mCurrentStateVariables.size() != rValue.size())
            mCurrentStateVariables.resize(rValue.size(), false);
        noalias(mCurrentStateVariables) = rValue;
    }
    if ( rVariable == MATERIAL_PARAMETERS )
    {
        if(mPROPS.size() != rValue.size())
            mPROPS.resize(rValue.size(), false);
//            KRATOS_THROW_ERROR(std::logic_error, "The size of the MATERIAL_PARAMETERS variable is incompatible", "")
        noalias(mPROPS) = rValue;
        #ifdef DEBUG_UMAT
        if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
            std::cout << "At element " << mElementId << ", point " << mIntPointIndex << ", MATERIAL_PARAMETERS is set to " << mPROPS << std::endl;
        #endif
    }
    if ( rVariable == INTERNAL_VARIABLES )
    {
        if(mOldStateVariables.size() != rValue.size())
            mOldStateVariables.resize(rValue.size());
//            KRATOS_THROW_ERROR(std::logic_error, "The size of the INTERNAL_VARIABLES variable is incompatible", "")
        noalias(mOldStateVariables) = rValue;
    }
}

void Umat3e::SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{
}

void Umat3e::SetValue( const Variable<std::string>& rVariable, const std::string& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == UMAT_CMNAME )
    {
        mCMNAME = rValue;
    }
}

void Umat3e::InitializeMaterial( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues )
{
    mNDI = props[UMAT_NDI];
    mNSHR = props[UMAT_NSHR];
    mNSTATV = props[UMAT_NSTATV];
    mPROPS = props[MATERIAL_PARAMETERS];
    mCMNAME = props[UMAT_CMNAME];
    int NTENS = mNDI + mNSHR;

    mCurrentStress.resize(NTENS);
    mOldStress.resize(NTENS);
    mPrestress.resize(NTENS);
    mCurrentStrain.resize(NTENS);
    mOldStrain.resize(NTENS);
    mCurrentStateVariables.resize(mNSTATV);
    mOldStateVariables.resize(mNSTATV);
    mPrestressFactor = 1.0;
    mOldStressZZ = 0.0;
    mCurrentStressZZ = 0.0;

    noalias(mCurrentStress) = ZeroVector(NTENS);
    noalias(mOldStress) = ZeroVector(NTENS);
    noalias(mPrestress) = ZeroVector(NTENS);
    noalias(mCurrentStrain) = ZeroVector(NTENS);
    noalias(mOldStrain) = ZeroVector(NTENS);
    noalias(mCurrentStateVariables) = ZeroVector(mNSTATV);
    noalias(mOldStateVariables) = ZeroVector(mNSTATV);

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

void Umat3e::ResetMaterial( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    int NTENS = mNDI + mNSHR;

    if(mCurrentStress.size() != NTENS)
        mCurrentStress.resize(NTENS, false);
    noalias(mCurrentStress) = ZeroVector(NTENS);

    if(mPrestress.size() != NTENS)
        mPrestress.resize(NTENS, false);
    noalias(mPrestress) = ZeroVector(NTENS);

    if(mCurrentStrain.size() != NTENS)
        mCurrentStrain.resize(NTENS, false);
    noalias(mCurrentStrain) = ZeroVector(NTENS);

    if(mCurrentStateVariables.size() != mNSTATV)
        mCurrentStateVariables.resize(mNSTATV, false);
    noalias(mCurrentStateVariables) = ZeroVector(mNSTATV);

    mPrestressFactor = 1.0;
    mCurrentStressZZ = 0.0;

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
        KRATOS_WATCH(mPrestress)
    #endif

    ResetState();
}

void Umat3e::ResetState()
{
    int NTENS = mNDI + mNSHR;

    if(mOldStress.size() != NTENS)
        mOldStress.resize(NTENS, false);
//    noalias(mOldStress) = ZeroVector(NTENS);
    noalias(mOldStress) = -mPrestressFactor*mPrestress;

    mOldStressZZ = 0.0;

    if(mOldStrain.size() != NTENS)
        mOldStrain.resize(NTENS, false);
    noalias(mOldStrain) = ZeroVector(NTENS);

    if(mOldStateVariables.size() != mNSTATV)
        mOldStateVariables.resize(mNSTATV, false);
    noalias(mOldStateVariables) = ZeroVector(mNSTATV);
}

void Umat3e::InitializeSolutionStep( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{
}

void Umat3e::InitializeNonLinearIteration ( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo )
{
}

void Umat3e::CalculateMaterialResponse( const Vector& StrainVector,
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
    int NPROPS = mPROPS.size();
    int NTENS = mNDI + mNSHR;

    // TODO fix the variable size array issue
    double STRAN[NTENS];
    double STRES[NTENS];
    double DSTRAN[NTENS];
    double DDSDDE[NTENS][NTENS];
    double STATEV[mNSTATV];
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

    if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            STRAN[i] = StrainVector[K2A[i]];
            DSTRAN[i] = StrainVector[K2A[i]] - mOldStrain[K2A[i]];
            STRES[i] = mOldStress[K2A[i]];
            for(int j = 0; j < NTENS; ++j)
                DDSDDE[i][j] = 0.0;
        }
    }
    else if( (mNDI == 2) & (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            STRAN[i] = StrainVector[i];
            DSTRAN[i] = StrainVector[i] - mOldStrain[i];
            STRES[i] = mOldStress[i];
            for(int j = 0; j < NTENS; ++j)
                DDSDDE[i][j] = 0.0;
        }
    }
    else if( (mNDI == 3) & (mNSHR == 1) ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
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

    for(int i = 0; i < mNSTATV; ++i)
    {
        STATEV[i] = mOldStateVariables[i];
    }

    TIM[0] = CurrentProcessInfo[TIME];
    TIM[1] = CurrentProcessInfo[TIME];
    DTIM = CurrentProcessInfo[DELTA_TIME];

    TEMP = 0.0;
    DTEMP = 0.0; // TODO: take into account the temperature
//if(mElementId == 1)
//{
//    KRATOS_WATCH(mCMNAME)
//    KRATOS_WATCH(mPROPS)
//    KRATOS_WATCH(StrainVector)
//    KRATOS_WATCH(mOldStrain)
//    KRATOS_WATCH(mOldStress)
//    KRATOS_WATCH(mPrestress)
//    KRATOS_WATCH(mPrestressFactor)
//    std::cout << "STRAN:";
//    for(int i = 0; i < 6; ++i)
//        std::cout << " " << STRAN[i];
//    std::cout << std::endl;
//    std::cout << "DSTRAN:";
//    for(int i = 0; i < 6; ++i)
//        std::cout << " " << DSTRAN[i];
//    std::cout << std::endl;
//    std::cout << "STRES:";
//    for(int i = 0; i < 6; ++i)
//        std::cout << " " << STRES[i];
//    std::cout << std::endl;
//}
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
          (char*) mCMNAME.c_str(),
          &mNDI, &mNSHR, &NTENS, &mNSTATV, &mPROPS[0], &NPROPS,
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
          (char*) mCMNAME.c_str(),
          &mNDI, &mNSHR, &NTENS, &mNSTATV, &mPROPS[0], &NPROPS,
          NULL, NULL, NULL, NULL,
          (double**)DFGRD0, (double**)DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );
    #endif

    noalias(mCurrentStrain) = StrainVector;

    if(mNDI == 3 && mNSHR == 3)
    {
        for(int i = 0; i < NTENS; ++i)
        {
            mCurrentStress[i] = STRES[A2K[i]];
        }
    }
    else if(mNDI == 2 && mNSHR == 1)
    {
        for(int i = 0; i < NTENS; ++i)
        {
            mCurrentStress[i] = STRES[i];
        }
    }
    else if(mNDI == 3 && mNSHR == 1)
    {
        for(int i = 0; i < 3; ++i)
        {
            mCurrentStress[i] = STRES[PS[i]];
        }
        mCurrentStressZZ = STRES[2];
    }

    noalias(StressVector) = mCurrentStress; // TODO: take into account the temperature
//if(mElementId == 1)
//{
//    KRATOS_WATCH(StressVector)
//    KRATOS_WATCH("-----------------")
//}
    for(int i = 0; i < mNSTATV; ++i)
    {
        mCurrentStateVariables[i] = STATEV[i];
    }
//TODO verify if DDSDDE need to be transpose because of Fortran (!!!!!!IMPORTANT!!!!)
    if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            for(int j = 0; j < NTENS; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[A2K[j]][A2K[i]];
            }
        }
    }
    else if( (mNDI == 2) && (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
    {
        for(int i = 0; i < NTENS; ++i)
        {
            for(int j = 0; j < NTENS; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[j][i];
            }
        }
    }
    else if( (mNDI == 3) && (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_zz  o_xy]
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

void Umat3e::FinalizeNonLinearIteration ( const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues,
                                         const ProcessInfo& CurrentProcessInfo )
{
}

void Umat3e::FinalizeSolutionStep( const Properties& props,
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
