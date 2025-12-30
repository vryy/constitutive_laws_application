#include "includes/define.h"
#include "constitutive_laws/umat3.h"
#include "constitutive_laws/umat3e.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry.h"
#include "structural_application/structural_application_variables.h"
#include "custom_utilities/shared_library_handle_shield.h"

// #define DEBUG_UMAT
#define DEBUG_ELEMENT_ID    1
#define DEBUG_POINT_ID      0

namespace Kratos
{

// REF: https://code.google.com/p/parafem/wiki/UMATs

typedef Kratos::ConstitutiveLaw::GeometryType GeometryType;

Umat3e::Umat3e()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    mp_umat_handle = nullptr;
    #endif
    Umat = nullptr;
}

Umat3e::~Umat3e()
{
    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    mp_umat_handle = nullptr;
    #endif
    Umat = nullptr;
}

int Umat3e::Check( const Kratos::Properties& props, const GeometryType& geom, const Kratos::ProcessInfo& CurrentProcessInfo ) const
{
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

// bool Umat3e::Has( const Variable<std::string>& rThisVariable )
// {
//     if (rThisVariable == ABAQUS_LIBRARY_NAME)
//         return true;
//     if (rThisVariable == UMAT_NAME)
//         return true;
//     if (rThisVariable == UMAT_CMNAME)
//         return true;
//     return false;
// }

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

        #ifdef DEBUG_UMAT
        if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
        {
            KRATOS_WATCH(__LINE__)
            KRATOS_WATCH(mOldStateVariables)
        }
        #endif
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
    if(rThisVariable == ABAQUS_LIBRARY_NAME)
        rValue = mLibName;
    if(rThisVariable == UMAT_NAME)
        rValue = mName;
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
        mPresetInternalVariables = true;
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

    mStep = 0;
    mIncrement = 0;
    mPresetInternalVariables = false;

    {
        mp_umat_handle = DLL::GetSharedLibraryHandle(mLibName);
        if(mp_umat_handle == nullptr)
        {
            KRATOS_ERROR << "Error loading Abaqus material library " << mLibName;
        }

        Umat = (umat_t) DLL::GetSymbol(mp_umat_handle, mName);
        const char* error = DLL::GetError();
        if(error != nullptr)
        {
            KRATOS_ERROR << "Error loading subroutine " << mName << " in the " << mLibName << " library"
                         << ", error message = " << DLL::GetErrorMessage(error);
        }
        else
        {
            std::cout << "Successfully load Umat from " << mLibName << std::endl;
        }
    }
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
    noalias(mCurrentStateVariables) = mOldStateVariables;

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

    this->ResetState();
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

    if (!mPresetInternalVariables)
    {
        if(mOldStateVariables.size() != mNSTATV)
            mOldStateVariables.resize(mNSTATV, false);
        noalias(mOldStateVariables) = ZeroVector(mNSTATV);
    }
    else
        mPresetInternalVariables = false;

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "Element " << mElementId << ", point " << mIntPointIndex << " at ResetState:" << std::endl;
        KRATOS_WATCH(mPrestress)
        KRATOS_WATCH(mOldStress)
        KRATOS_WATCH(mOldStressZZ)
        KRATOS_WATCH(mOldStateVariables)
    }
    #endif
}

void Umat3e::InitializeSolutionStep( const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues ,
                                     const ProcessInfo& CurrentProcessInfo )
{
    ++mStep;
    mIncrement = 0;
}

void Umat3e::InitializeNonLinearIteration ( const Properties& props,
                                            const GeometryType& geom,
                                            const Vector& ShapeFunctionsValues,
                                            const ProcessInfo& CurrentProcessInfo )
{
    ++mIncrement;
}

void Umat3e::CalculateMaterialResponseCauchy( Parameters& parameters )
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
    this->CalculateMaterialResponse( StrainVector, DeformationGradient,
            StressVector, AlgorithmicTangent, CurrentProcessInfo, props, geom,
            ShapeFunctionsValues, CalculateStresses, CalculateTangent,
            SaveInternalVariables, true, false );
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

    int NPROPS = mPROPS.size();
    int NTENS = mNDI + mNSHR;

    // TODO fix the variable size array issue
    std::vector<double> STRAN(NTENS);
    std::vector<double> STRES(NTENS);
    std::vector<double> DSTRAN(NTENS);
    std::vector<double> DDSDDE(NTENS*NTENS);
    std::vector<double> STATEV(mNSTATV);
    std::vector<double> DDSDDT(NTENS);
    std::vector<double> DRPLDE(NTENS);
    std::vector<double> DRPLDT(NTENS);
    double DFGRD0[9];
    double DFGRD1[9];
    double TIM[2];
    double DTIM;
    double SSE, SPD, SCD, RPL;
    double TEMP, DTEMP;
    int KSTEP = (int) mStep, KINC = (int) mIncrement;

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
    {
        std::cout << "Element " << mElementId << ", point " << mIntPointIndex << " at CalculateMaterialResponse:" << std::endl;
        KRATOS_WATCH(mOldStress)
        KRATOS_WATCH(StrainVector)
        KRATOS_WATCH(mOldStrain)
        KRATOS_WATCH(StrainVector - mOldStrain)
        KRATOS_WATCH(mOldStateVariables)
        KRATOS_WATCH(CurrentProcessInfo[DELTA_TIME])
        // KRATOS_THROW_ERROR(std::logic_error, "stop here", "")
    }
    #endif

    if (IsSetStrain)
    {
        if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                STRAN[i] = StrainVector[Umat3::K2A[i]];
                DSTRAN[i] = StrainVector[Umat3::K2A[i]] - mOldStrain[Umat3::K2A[i]];
                STRES[i] = mOldStress[Umat3::K2A[i]];
                for(int j = 0; j < NTENS; ++j)
                    DDSDDE[i*NTENS + j] = 0.0;
            }
        }
        else if( (mNDI == 2) && (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
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
        else if( (mNDI == 3) && (mNSHR == 1) ) // 2D case, plane strain    [o_xx  o_yy  o_zz  o_xy]
        {
            for(int i = 0; i < 3; ++i)
            {
                STRAN[Umat3::PS[i]] = StrainVector[i];
                DSTRAN[Umat3::PS[i]] = StrainVector[i] - mOldStrain[i];
                STRES[Umat3::PS[i]] = mOldStress[i];
            }
            STRAN[2] = 0.0;
            DSTRAN[2] = 0.0;
            for(int i = 0; i < 4; ++i)
                for(int j = 0; j < 4; ++j)
                    DDSDDE[i*NTENS + j] = 0.0;
            STRES[2] = mOldStressZZ;
        }
    }

    if (IsSetDeformationGradient)
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                DFGRD0[i + j*3] = 0.0;
                DFGRD1[i + j*3] = DeformationGradient(i, j);
            }
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

    Umat( STRES.data(),
          STATEV.data(),
          DDSDDE.data(),
          &SSE, &SPD, &SCD, &RPL,
          DDSDDT.data(), DRPLDE.data(), DRPLDT.data(),
          STRAN.data(), DSTRAN.data(),
          TIM, &DTIM,
          &TEMP, &DTEMP,
          NULL, NULL,
          (char*) mCMNAME.c_str(),
          &mNDI, &mNSHR, &NTENS, &mNSTATV, &mPROPS[0], &NPROPS,
          NULL, NULL, NULL, NULL,
          DFGRD0, DFGRD1,
          &mElementId, &mIntPointIndex,
          NULL, NULL,
          &KSTEP, &KINC );

    noalias(mCurrentStrain) = StrainVector;

    if(mNDI == 3 && mNSHR == 3)
    {
        for(int i = 0; i < NTENS; ++i)
        {
            mCurrentStress[i] = STRES[Umat3::A2K[i]];
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
            mCurrentStress[i] = STRES[Umat3::PS[i]];
        }
        mCurrentStressZZ = STRES[2];
    }

    if (CalculateStresses)
        noalias(StressVector) = mCurrentStress;

    // TODO: take into account the temperature

    for(int i = 0; i < mNSTATV; ++i)
    {
        mCurrentStateVariables[i] = STATEV[i];
    }

    if (CalculateTangent)
    {
        if( (mNDI == 3) && (mNSHR == 3) ) // 3D case    [o_xx  o_yy  o_zz  o_xy  o_yz  o_xz]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                for(int j = 0; j < NTENS; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[Umat3::A2K[j]*NTENS + Umat3::A2K[i]];
                }
            }
        }
        else if( (mNDI == 2) && (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_xy]
        {
            for(int i = 0; i < NTENS; ++i)
            {
                for(int j = 0; j < NTENS; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[j*NTENS + i];
                }
            }
        }
        else if( (mNDI == 3) && (mNSHR == 1) ) // 2D case    [o_xx  o_yy  o_zz  o_xy]
        {
            for(int i = 0; i < 3; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    AlgorithmicTangent(i, j) = DDSDDE[Umat3::PS[j]*NTENS + Umat3::PS[i]];
                }
            }
        }
    //    KRATOS_WATCH(AlgorithmicTangent)
        // TODO: export the variable SSE, SPD, SCD, RPL
    }

    #ifdef DEBUG_UMAT
    if(mElementId == DEBUG_ELEMENT_ID && mIntPointIndex == DEBUG_POINT_ID)
    {
        if (CalculateStresses)
            KRATOS_WATCH(mCurrentStress)
        KRATOS_WATCH(mCurrentStateVariables)
        if (CalculateTangent)
            KRATOS_WATCH(AlgorithmicTangent)
    }
    #endif
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

#undef DEBUG_UMAT
#undef DEBUG_ELEMENT_ID
#undef DEBUG_POINT_ID
