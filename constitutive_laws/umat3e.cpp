#include <iostream>
#include <dlfcn.h>

#include "includes/define.h"
#include "constitutive_laws/umat3e.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry.h"
#include "constitutive_laws_application.h"
#include "structural_application/structural_application_variables.h"

// #define DEBUG_UMAT
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

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN )
    {
        if (mOldStress.size() == 3)
            rValue = (mOldStrain[0] + mOldStrain[1]);
        else if (mOldStress.size() == 6)
            rValue = (mOldStrain[0] + mOldStrain[1] + mOldStrain[2]);
    }

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN )
    {
        // REF: https://dianafea.com/manuals/d944/Analys/node405.html
        double exx, eyy, ezz, exy, eyz, exz, ev;
        if (mOldStrain.size() == 3)
        {
            ev = (mOldStrain[0] + mOldStrain[1]);
            exx = mOldStrain[0] - ev/3;
            eyy = mOldStrain[1] - ev/3;
            ezz = -ev/3;
            exy = mOldStrain[2];
            eyz = 0.0;
            exz = 0.0;
        }
        else if (mOldStrain.size() == 6)
        {
            ev = (mOldStrain[0] + mOldStrain[1] + mOldStrain[2]);
            exx = mOldStrain[0] - ev/3;
            eyy = mOldStrain[1] - ev/3;
            ezz = mOldStrain[2] - ev/3;
            exy = mOldStrain[3];
            eyz = mOldStrain[4];
            exz = mOldStrain[5];
        }

        rValue = (2.0/3) * std::sqrt(1.5*(pow(exx, 2) + pow(eyy, 2) + pow(ezz, 2) + 3*(pow(exy, 2) + pow(eyz, 2) + pow(exz, 2))));
    }

    if (rThisVariable == PRESSURE_P) // this follows the soil mechanics's hydrostatic pressure convention
    {
        if (mOldStress.size() == 3)
            rValue = -(mOldStress[0] + mOldStress[1] + mOldStressZZ) / 3;
        else if (mOldStress.size() == 6)
            rValue = -(mOldStress[0] + mOldStress[1] + mOldStress[2]) / 3;
    }

    if (rThisVariable == PRESSURE_Q || rThisVariable == VON_MISES_STRESS)
    {
        double s00, s11, s22, s01, s12, s02, p;
        double j2;
        if (mOldStress.size() == 3)
        {
            p = (mOldStress[0] + mOldStress[1] + mOldStressZZ) / 3;
            s00 = mOldStress[0] - p;
            s11 = mOldStress[1] - p;
            s22 = mOldStressZZ - p;
            s01 = mOldStress[2];
            s12 = 0.0;
            s02 = 0.0;
        }
        else if (mOldStress.size() == 6)
        {
            p = (mOldStress[0] + mOldStress[1] + mOldStress[2]) / 3;
            s00 = mOldStress[0] - p;
            s11 = mOldStress[1] - p;
            s22 = mOldStress[2] - p;
            s01 = mOldStress[3];
            s12 = mOldStress[4];
            s02 = mOldStress[5];
        }
        j2 = 0.5*( pow(s00, 2) + pow(s11, 2) + pow(s22, 2) ) + pow(s01, 2) + pow(s12, 2) + pow(s02, 2);
        rValue = sqrt( 3.0 * j2 );
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
    {
        mElementId = rValue;
    }
    if(rVariable == INTEGRATION_POINT_INDEX)
    {
        mIntPointIndex = rValue;
    }
    if(rVariable == UMAT_NDI)
    {
        mNDI = rValue;
    }
    if(rVariable == UMAT_NSHR)
    {
        mNSHR = rValue;
    }
    if(rVariable == UMAT_NSTATV)
    {
        mNSTATV = rValue;
    }
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
        if (mPrestress.size() != rValue.size())
            mPrestress.resize(rValue.size());
        noalias(mPrestress) = rValue;

        #ifdef DEBUG_UMAT
        if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(__LINE__)
            KRATOS_WATCH(rValue)
            KRATOS_WATCH(mPrestress)
        }
        #endif

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
    if ( rVariable == ABAQUS_LIBRARY_NAME )
    {
        mLibName = rValue;
    }
    if ( rVariable == UMAT_NAME )
    {
        mName = rValue;
    }
}

void Umat3e::InitializeMaterial( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues )
{
    if (props.Has(UMAT_NDI))
        mNDI = props[UMAT_NDI];

    if (props.Has(UMAT_NSHR))
        mNSHR = props[UMAT_NSHR];

    if (props.Has(UMAT_NSTATV))
        mNSTATV = props[UMAT_NSTATV];

    if (props.Has(MATERIAL_PARAMETERS))
        mPROPS = props[MATERIAL_PARAMETERS];

    if (props.Has(UMAT_NAME))
        mName = props[UMAT_NAME];

    if (props.Has(UMAT_CMNAME))
        mCMNAME = props[UMAT_CMNAME];

    if (props.Has(ABAQUS_LIBRARY_NAME))
        mLibName = props[ABAQUS_LIBRARY_NAME];

    int NTENS = mNDI + mNSHR;

    unsigned int strain_size;
    if (mNDI == 3 && mNSHR == 3)
    {
        strain_size = 6;
    }
    else if (mNDI == 3 && mNSHR == 1)
    {
        strain_size = 3;
    }
    else if (mNDI == 2 && mNSHR == 1)
    {
        strain_size = 3;
    }

    mCurrentStress.resize(strain_size);
    mOldStress.resize(strain_size);
    mPrestress.resize(strain_size);
    mCurrentStrain.resize(strain_size);
    mOldStrain.resize(strain_size);
    mCurrentStateVariables.resize(mNSTATV);
    mOldStateVariables.resize(mNSTATV);
    mPrestressFactor = 1.0;
    mOldStressZZ = 0.0;
    mCurrentStressZZ = 0.0;

    noalias(mCurrentStress) = ZeroVector(strain_size);
    noalias(mOldStress) = ZeroVector(strain_size);
    noalias(mPrestress) = ZeroVector(strain_size);
    noalias(mCurrentStrain) = ZeroVector(strain_size);
    noalias(mOldStrain) = ZeroVector(strain_size);
    noalias(mCurrentStateVariables) = ZeroVector(mNSTATV);
    noalias(mOldStateVariables) = ZeroVector(mNSTATV);

    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    // get the library name and load the udsm subroutine
//    mp_umat_handle = dlopen(mLibName.c_str(), RTLD_LAZY);
    mp_umat_handle = dlopen(mLibName.c_str(), RTLD_NOW | RTLD_GLOBAL);
    if(mp_umat_handle == 0)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "The Abaqus material library does not exist:", mLibName)
    }

    char* error;
    Umat = (void (*)(double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                double** DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC))
           dlsym(mp_umat_handle, mName.c_str());
    error = dlerror();
    if(error != NULL)
    {
        std::stringstream ss;
        ss << "Error loading subroutine " << mName << " in the " << mLibName << " library, error message = " << error;
        KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
    }
    #endif
}

void Umat3e::ResetMaterial( const Properties& props,
                           const GeometryType& geom,
                           const Vector& ShapeFunctionsValues )
{
    unsigned int strain_size;
    if (mNDI == 3 && mNSHR == 3)
    {
        strain_size = 6;
    }
    else if (mNDI == 3 && mNSHR == 1)
    {
        strain_size = 3;
    }
    else if (mNDI == 2 && mNSHR == 1)
    {
        strain_size = 3;
    }

    if(mCurrentStrain.size() != strain_size)
        mCurrentStrain.resize(strain_size, false);
    noalias(mCurrentStrain) = ZeroVector(strain_size);

    if(mCurrentStateVariables.size() != mNSTATV)
        mCurrentStateVariables.resize(mNSTATV, false);
    noalias(mCurrentStateVariables) = ZeroVector(mNSTATV);

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "Element " << mElementId << ", point " << mIntPointIndex << " at ResetMaterial:" << std::endl;
        KRATOS_WATCH(mPrestress)
        KRATOS_WATCH(mNDI)
        KRATOS_WATCH(mNSHR)
        KRATOS_WATCH(mNSTATV)
    }
    #endif

    ResetState();
}

void Umat3e::ResetState()
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

    unsigned int strain_size;
    if (mNDI == 3 && mNSHR == 3)
    {
        strain_size = 6;
    }
    else if (mNDI == 3 && mNSHR == 1)
    {
        strain_size = 3;
    }
    else if (mNDI == 2 && mNSHR == 1)
    {
        strain_size = 3;
    }

    if(mOldStrain.size() != strain_size)
        mOldStrain.resize(strain_size, false);
    noalias(mOldStrain) = ZeroVector(strain_size);

    if(mOldStateVariables.size() != mNSTATV)
        mOldStateVariables.resize(mNSTATV, false);
    noalias(mOldStateVariables) = ZeroVector(mNSTATV);

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "Element " << mElementId << ", point " << mIntPointIndex << " at ResetState:" << std::endl;
        KRATOS_WATCH(mPrestress)
        KRATOS_WATCH(mOldStress)
        KRATOS_WATCH(mOldStressZZ)
    }
    #endif
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
    if (CalculateStresses && !CalculateTangent)
    {
        noalias(StressVector) = mCurrentStress;
        return;
    }

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

    noalias(StressVector) = mCurrentStress;

    // TODO: take into account the temperature

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
