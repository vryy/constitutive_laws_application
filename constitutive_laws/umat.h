#if !defined(KRATOS_UMAT_H_INCLUDED )
#define  KRATOS_UMAT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

class Umat : public ConstitutiveLaw
{

public:

    // Type Definitions
    typedef ConstitutiveLaw BaseType;
    typedef array_1d<double, 81> MaterialTensorType;
    KRATOS_CLASS_POINTER_DEFINITION( Umat );
    typedef array_1d<double, 3 > PlaneArrayType;
    typedef array_1d<double, 6 > SpaceArrayType;

    //zero constructor
    Umat();

    //destructor
    virtual ~Umat();

    //clone
    BaseType::Pointer Clone() const final
    {
        BaseType::Pointer p_clone ( new Umat() );
        return p_clone;
    }

    size_t WorkingSpaceDimension() final
    {
        return 3;
    }

    size_t GetStrainSize() const final
    {
        return 6;
    }

    //operators
    bool Has ( const Variable<double>& rThisVariable ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::Has(double)", "" );
    }

    bool Has ( const Variable<Vector>& rThisVariable ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::Has(Vector)", "" );
    }

    bool Has ( const Variable<Matrix>& rThisVariable ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::Has(Matrix)", "" );
    }

    bool Has ( const Variable<PlaneArrayType>& rThisVariable ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::Has(2DVector)", "" );
    }

    bool Has ( const Variable<SpaceArrayType>& rThisVariable ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::Has(3DVector)", "" );
    }

    double& GetValue ( const Variable<double>& rThisVariable,
                               double& rValue ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::GetValue(double)", "" );
    }

    Vector& GetValue ( const Variable<Vector>& rThisVariable,
                               Vector& rValue ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::GetValue(Vector)", "" );
    }

    Matrix& GetValue ( const Variable<Matrix>& rThisVariable,
                               Matrix& rValue ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::GetValue(Matrix)", "" );
    }

    PlaneArrayType& GetValue ( const Variable<PlaneArrayType>& rVariable,
                                       PlaneArrayType& rValue ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::GetValue(2DVector)", "" );
    }

    SpaceArrayType & GetValue ( const Variable<SpaceArrayType>& rVariable,
                                        SpaceArrayType& rValue ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::GetValue(3DVector)", "" );
    }


    void SetValue ( const Variable<double>& rThisVariable,
                    const double& rValue,
                    const ProcessInfo& rCurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::SetValue(Double)", "" );
    }

    void SetValue ( const Variable<Vector>& rThisVariable,
                    const Vector& rValue,
                    const ProcessInfo& rCurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::SetValue(Vector)", "" );
    }

    void SetValue ( const Variable<Matrix>& rThisVariable,
                    const Matrix& rValue,
                    const ProcessInfo& rCurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::SetValue(Matrix)", "" );
    }

    void SetValue ( const Variable<PlaneArrayType>& rThisVariable,
                    const PlaneArrayType& rValue,
                    const ProcessInfo& rCurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::SetValue(2DVector)", "" );
    }

    void SetValue ( const Variable<SpaceArrayType>& rVariable,
                    const SpaceArrayType& Value,
                    const ProcessInfo& rCurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::SetValue(3DVector)", "" );
    }

    bool ValidateInput ( const Properties& props ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::ValidateInput called", "" );
    }

    StrainMeasure GetStrainMeasure() final
    {
        return StrainMeasure_Infinitesimal;
    }

    StressMeasure GetStressMeasure() final
    {
        return StressMeasure_PK1;
    }

    bool IsIncremental() final
    {
        return false;
    }



    //Material parameters inizialization
    void InitializeMaterial ( const Properties& props,
                              const GeometryType& geom,
                              const Vector& ShapeFunctionsValues ) final;


    void InitializeSolutionStep ( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ,
                                  const ProcessInfo& CurrentProcessInfo ) final;

    void FinalizeSolutionStep ( const Properties& props,
                                const GeometryType& geom,
                                const Vector& ShapeFunctionsValues ,
                                const ProcessInfo& CurrentProcessInfo ) final;

    void InitializeNonLinearIteration ( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues,
                                        const ProcessInfo& CurrentProcessInfo ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::InitializeNonLinearIteration called", "" );
    }

    void CalculateMaterialResponse ( const Vector& StrainVector,
            const Matrix& DeformationGradient,
            Vector& StressVector,
            Matrix& AlgorithmicTangent,
            const ProcessInfo& CurrentProcessInfo,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            bool CalculateStresses = true,
            int CalculateTangent = true,
            bool SaveInternalVariables = true ) final;


    void CalculateVolumetricResponse ( const double VolumetricStrain,
            const Matrix& DeformationGradient,
            double& VolumetricStress,
            double& AlgorithmicBulk,
            const ProcessInfo& CurrentProcessInfo,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            bool CalculateStresses = true,
            int CalculateTangent = true,
            bool SaveInternalVariables = true ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::CalculateVolumetricResponse called", "" );
    }

    void CalculateDeviatoricResponse ( const Vector& StrainVector,
            const Matrix& DeformationGradient,
            Vector& StressVector,
            Matrix& AlgorithmicTangent,
            const ProcessInfo& CurrentProcessInfo,
            const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues,
            bool CalculateStresses = true,
            int CalculateTangent = true,
            bool SaveInternalVariables = true ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::CalculateDeviatoricResponse called", "" );
    }

    void ResetMaterial ( const Properties& props,
                         const GeometryType& geom,
                         const Vector& ShapeFunctionsValues ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function Umat::ResetMaterial called", "" );
    }

    void CalculateCauchyStresses ( Vector& Cauchy_StressVector,
                                   const Matrix& F,
                                   const Vector& PK2_StressVector,
                                   const Vector& GreenLagrangeStrainVector ) final;

    int Check ( const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo ) const final;

private:

    //member variables for umat
    double* STRESS;
    double* STATEV;
    double* PROPS;
    int* NSTATV;
    int* NPROPS;

    //static member variables for umat
    double* STRAN;
    double* DSTRAN;
    int* NDI;
    int* NSHR;
    int* NTENS;


    //variable for wrapper
    int* MaterialNumber;

    //serialization

    friend class Serializer;

    void save ( Serializer& rSerializer ) const override
    {
        rSerializer.save ( "name", "Umat" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

    void load ( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

}; // Class Umat

}  // namespace Kratos.

#endif // KRATOS_UMAT_H_INCLUDED  defined

