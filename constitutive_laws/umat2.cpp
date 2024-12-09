#include <iostream>

#ifdef _MSC_VER
#include <Windows.h>
#else
#include <dlfcn.h>
#endif

#include "includes/define.h"
#include "constitutive_laws/umat2.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/ublas_interface.h"
#include "constitutive_laws_application_variables.h"
#include "structural_application/structural_application_variables.h"

namespace Kratos
{

// REF: https://code.google.com/p/parafem/wiki/UMATs

typedef Kratos::ConstitutiveLaw::GeometryType GeometryType;

const int Umat2::A2K[] = {0, 1, 2, 3, 5, 4};
const int Umat2::K2A[] = {0, 1, 2, 3, 5, 4};
const int Umat2::PS[] = {0, 1, 3};

Umat2::Umat2()
{
}

Umat2::~Umat2()
{
    if (mp_umat_handle != 0)
    {
#ifdef _MSC_VER
        // TODO
#else
        dlclose(mp_umat_handle);
#endif
    }
}

int Umat2::Check( const Kratos::Properties& props, const GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo ) const
{
    if(props.Has( ABAQUS_LIBRARY_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define ABAQUS_LIBRARY_NAME";
    }
    if(props.Has( UMAT_NAME ) == false)
    {
        KRATOS_ERROR << "Properties must define UMAT_NAME";
    }
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

bool Umat2::Has( const Variable<int>& rThisVariable )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
        return true;
    if(rThisVariable == INTEGRATION_POINT_INDEX)
        return true;
    return false;
}

bool Umat2::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Umat2::Has( const Variable<Vector>& rThisVariable )
{
    if(rThisVariable == INTERNAL_VARIABLES)
        return true;
    if(rThisVariable == STRESSES)
        return true;
    if(rThisVariable == STRAIN)
        return true;
    return false;
}

bool Umat2::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

int& Umat2::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& Umat2::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == PRESTRESS_FACTOR)
    {
        rValue = mPrestressFactor;
    }
    return rValue;
}

Vector& Umat2::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
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
    if(rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS)
    {
        if(rValue.size() != mPrestress.size())
            rValue.resize(mPrestress.size(), false);
        noalias(rValue) = mPrestress;
    }
    if(rThisVariable == STRAIN)
    {
        if(rValue.size() != mCurrentStrain.size())
            rValue.resize(mCurrentStrain.size(), false);
        noalias(rValue) = mCurrentStrain;
    }

    return rValue;
}

Matrix& Umat2::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    return rValue;
}

void Umat2::SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable == PARENT_ELEMENT_ID)
        mElementId = rValue;
    if(rVariable == INTEGRATION_POINT_INDEX )
        mIntPointIndex = rValue;
}

void Umat2::SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
}

void Umat2::SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
}

void Umat2::SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{
}

void Umat2::InitializeMaterial( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues )
{

    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];

    mCurrentStress.resize(NTENS);
    mPrestress.resize(NTENS);
    mCurrentStrain.resize(NTENS);
    mCurrentStateVariables.resize(NSTATV);
    mPrestressFactor = 1.0;
    mCurrentStressZZ = 0.0;

    noalias(mCurrentStress) = ZeroVector(NTENS);
    noalias(mPrestress) = ZeroVector(NTENS);
    noalias(mCurrentStrain) = ZeroVector(NTENS);
    noalias(mCurrentStateVariables) = ZeroVector(NSTATV);

    // get the library name and load the udsm subroutine
    std::string lib_name = props[ABAQUS_LIBRARY_NAME];
#ifdef _MSC_VER
    // TODO
#else
//    mp_umat_handle = dlopen(lib_name.c_str(), RTLD_LAZY);
    mp_umat_handle = dlopen(lib_name.c_str(), RTLD_NOW | RTLD_GLOBAL);
#endif
    if(mp_umat_handle == 0)
    {
        KRATOS_ERROR << "The Abaqus material library " << lib_name << " does not exist";
    }

    std::string umat_name = props[UMAT_NAME];
#ifdef _MSC_VER
    // TODO
#else
    char* error;
    Umat = (void (*)(double* STRESS, double* STATEV, double* DDSDDE, double* SSE, double* SPD, double* SCD,
                double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                double* COORDS, double* DROT, double* PNEWDT, double* CELENT, double* DFGRD0,
                double* DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC))
           dlsym(mp_umat_handle, umat_name.c_str());
    error = dlerror();
    if(error != NULL)
    {
        KRATOS_ERROR << "Error loading subroutine " << umat_name << " in the " << lib_name << " library, error message = " << error;
    }
#endif
}

void Umat2::ResetMaterial( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    int NDI = props[UMAT_NDI];
    int NSHR = props[UMAT_NSHR];
    int NTENS = NDI + NSHR;
    int NSTATV = props[UMAT_NSTATV];

    noalias(mCurrentStress) = ZeroVector(NTENS);
    noalias(mPrestress) = ZeroVector(NTENS);
    noalias(mCurrentStrain) = ZeroVector(NTENS);
    noalias(mCurrentStateVariables) = ZeroVector(NSTATV);
    mPrestressFactor = 1.0;
}

void Umat2::InitializeSolutionStep( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo )
{
}

void Umat2::InitializeNonLinearIteration ( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues,
                                           const ProcessInfo& CurrentProcessInfo )
{
}

void Umat2::CalculateMaterialResponse( const Vector& StrainVector,
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
//    KRATOS_WATCH(PROPS(1))
//    PROPS(1) = mPrestressFactor * (mPrestress(0) + mPrestress(1) + mPrestress(2)) / 3; // turn it on if debugging the Hypoplasticity model
//    if(mElementId == 1 && mIntPointIndex == 0)
//        KRATOS_WATCH(PROPS(1))
    int NPROPS = PROPS.size();

    int my_NDI = props[UMAT_NDI];
    int my_NSHR = props[UMAT_NSHR];
    int NDI = 3;
    int NSHR = 3;
    int NTENS = NDI + NSHR; // 6
    int NSTATV = props[UMAT_NSTATV];
    char* cmname = (char*)(props[UMAT_CMNAME].c_str());

    std::vector<double> STRAN(NTENS);
    std::vector<double> STRES(NTENS);
    std::vector<double> DSTRAN(NTENS);
    std::vector<double> DDSDDE(NTENS*NTENS);
    std::vector<double> DDSDDT(NTENS);
    std::vector<double> DRPLDE(NTENS);
    std::vector<double> DRPLDT(NTENS);
    double DFGRD0[9];
    double DFGRD1[9];
    double TIM[2];
    double DTIM;
    double SSE, SPD, SCD, RPL;
    double TEMP, DTEMP;
    int KSTEP, KINC;

    if( (my_NDI == 3) && (my_NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < 6; ++i)
        {
            STRAN[i] = StrainVector[K2A[i]];
            DSTRAN[i] = StrainVector[K2A[i]] - mCurrentStrain[K2A[i]];
            STRES[i] = mCurrentStress[K2A[i]];
            for(int j = 0; j < 6; ++j)
                DDSDDE[i*6 + j] = 0.0;
        }
    }
    else if( ( (my_NDI == 2) && (my_NSHR == 1) )     // 2D case, plane strain    [o_xx  o_yy  o_xy]
          || ( (my_NDI == 3) && (my_NSHR == 1) ) )   // 2D case, plane strain    [o_xx  o_yy  o_zz o_xy]
    {
        for(int i = 0; i < 6; ++i)
        {
            STRAN[i] = 0.0;
            DSTRAN[i] = 0.0;
            STRES[i] = 0.0;
            for(int j = 0; j < 6; ++j)
                DDSDDE[i*6 + j] = 0.0;
        }
        for(int i = 0; i < 3; ++i)
        {
            STRAN[PS[i]] = StrainVector[i];
            DSTRAN[PS[i]] = StrainVector[i] - mCurrentStrain[i];
            STRES[PS[i]] = mCurrentStress[i];
        }
        STRES[2] = mCurrentStressZZ;
    }

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            DFGRD0[i + 3*j] = 0.0;
            DFGRD1[i + 3*j] = DeformationGradient(i, j);
        }
    }

    TIM[0] = CurrentProcessInfo[TIME];
    TIM[1] = CurrentProcessInfo[TIME];
    DTIM = CurrentProcessInfo[DELTA_TIME];

    TEMP = 0.0;
    DTEMP = 0.0; // TODO: take into account the temperature

    Umat( STRES.data(),
          &mCurrentStateVariables[0],
          DDSDDE.data(),
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT.data(), DRPLDE.data(), DRPLDT.data(),
          STRAN.data(), DSTRAN.data(),
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          cmname,
          &NDI, &NSHR, &NTENS, &NSTATV, &PROPS[0], &NPROPS,
          NULL, NULL, NULL, NULL,
          DFGRD0, DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );

    noalias(mCurrentStrain) = StrainVector;

    if(my_NDI == 3 && my_NSHR == 3)
    {
        for(int i = 0; i < 6; ++i)
        {
            mCurrentStress[i] = STRES[A2K[i]];
        }
    }
    else if( ( (my_NDI == 2) && (my_NSHR == 1) )     // 2D case, plane strain    [o_xx  o_yy  o_xy]
          || ( (my_NDI == 3) && (my_NSHR == 1) ) )   // 2D case, plane strain    [o_xx  o_yy  o_zz o_xy]
    {
        for(int i = 0; i < 3; ++i)
        {
            mCurrentStress[i] = STRES[PS[i]];
        }
        mCurrentStressZZ = STRES[2];
    }

    noalias(StressVector) = mCurrentStress + mPrestressFactor * mPrestress; // TODO: take into account the temperature

    if( (my_NDI == 3) && (my_NSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
    {
        for(int i = 0; i < 6; ++i)
        {
            for(int j = 0; j < 6; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[A2K[i] + A2K[j]*6];
            }
        }
    }
    else if( ( (my_NDI == 2) && (my_NSHR == 1) )     // 2D case, plane strain    [o_xx  o_yy  o_xy]
          || ( (my_NDI == 3) && (my_NSHR == 1) ) )   // 2D case, plane strain    [o_xx  o_yy  o_zz o_xy]
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                AlgorithmicTangent(i, j) = DDSDDE[PS[i] + PS[j]*6];
            }
        }
    }
//    KRATOS_WATCH(AlgorithmicTangent)

    // TODO: export the variable SSE, SPD, SCD, RPL
}

void Umat2::FinalizeNonLinearIteration ( const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues,
                                         const ProcessInfo& CurrentProcessInfo )
{
}

void Umat2::FinalizeSolutionStep( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues ,
                                 const ProcessInfo& CurrentProcessInfo )
{
}

} // Namespace Kratos
