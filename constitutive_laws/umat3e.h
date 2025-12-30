//
//   Project Name:        KratosConstitutiveLawsApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Nov 18, 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_UMAT3_E_H_INCLUDED )
#define  KRATOS_UMAT3_E_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "constitutive_laws_application_variables.h"


namespace Kratos
{

/**
Umat interface KRATOS-ABAQUS
THis interface will call respective Umat version depending on the parameters given. Therefore Umat must implement both 2D and 3D.
This version of Umat is the same as Umat3 but is customized for ekate, where the material parameters are read from matfile
*/
class Umat3e : public ConstitutiveLaw
{

public:

    // Counted pointer definition
    KRATOS_CLASS_POINTER_DEFINITION( Umat3e );

    // Type Definitions
    typedef ConstitutiveLaw BaseType;

    // Default constructor
    Umat3e();

    // Destructor
    ~Umat3e() override;

    //clone
    BaseType::Pointer Clone() const final
    {
        BaseType::Pointer p_clone ( new Umat3e() );
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

    StrainMeasure GetStrainMeasure() final
    {
        return StrainMeasure_Infinitesimal;
    }

    StressMeasure GetStressMeasure() final
    {
        return StressMeasure_Cauchy;
    }

    void GetLawFeatures(Features& rFeatures) final
    {
        rFeatures.SetStrainMeasure(this->GetStrainMeasure());
    }

    bool IsIncremental() final
    {
        return true;
    }

    bool Has ( const Variable<int>& rThisVariable ) final;

    bool Has ( const Variable<double>& rThisVariable ) final;

    bool Has ( const Variable<Vector>& rThisVariable ) final;

    bool Has ( const Variable<Matrix>& rThisVariable ) final;

    // bool Has ( const Variable<std::string>& rThisVariable ) final;

    int& GetValue ( const Variable<int>& rThisVariable, int& rValue ) final;

    double& GetValue ( const Variable<double>& rThisVariable, double& rValue ) final;

    Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue ) final;

    Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue ) final;

    std::string& GetValue ( const Variable<std::string>& rThisVariable, std::string& rValue ) final;

    void SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<std::string>& rVariable, const std::string& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    virtual bool ValidateInput ( const Properties& props )
    {
        KRATOS_ERROR << "virtual function Umat3e::ValidateInput called";
    }

    void InitializeMaterial ( const Properties& props,
                              const GeometryType& geom,
                              const Vector& ShapeFunctionsValues ) final;

    void ResetMaterial ( const Properties& props,
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
                                        const ProcessInfo& CurrentProcessInfo ) final;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy( Parameters& parameters ) final;

    /// DEPRECATED interface
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

    void FinalizeNonLinearIteration ( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo ) final;

    int Check ( const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo ) const final;

private:

    Vector mCurrentStrain;              // the current strain
    Vector mOldStrain;                  // the strain in previous "converged" load step
    Vector mCurrentStress;              // the current stress
    Vector mOldStress;                  // the stress in previous "converged" load step
    Vector mCurrentStateVariables;      // the current state variables
    Vector mOldStateVariables;          // the state variables in previous "converged" load step
    Vector mPrestress;
    double mPrestressFactor;
    double mCurrentStressZZ;            // the current stress ZZ
    double mOldStressZZ;                // the stress ZZ in previous "converged" load step

    // REMARKS
    // The UMAT convention for stress and strain is [sigma_11, sigma_22, sigma_33, sigma_12, {sigma_13}, {sigma_23}].
    // Meanwhile, KRATOS convention is [sigma_11, sigma_22, sigma_33, sigma_12, {sigma_23}, {sigma_13}].
    // Hence we have to swap the value when calling UMAT

    int mElementId;
    int mIntPointIndex;
    int mNDI;
    int mNSHR;
    int mNSTATV;
    Vector mPROPS;
    std::string mLibName;
    std::string mName;
    std::string mCMNAME;

    std::size_t mStep;  // step counter
    int mIncrement;     // increment counter
    bool mPresetInternalVariables; // flag to control the preset of the internal variables
        // In principle, the internal variables are always set to zero when ResetMaterial is called
        // However, one can override this behaviour by setting INTERNAL_VARIABLES before calling ResetMaterial
        // This flag will be reset whenever ResetMaterial is called

    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
    // handle to Umat library
    void* mp_umat_handle;
    #endif
    // wrapper function for calling the UMAT fortran subroutine
    umat_t Umat;

    /// Reset only the material state to initial state
    void ResetState();

    ///
    void CalculateMaterialResponse ( const Vector& StrainVector,
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
            bool IsSetDeformationGradient );

    // serialization
    friend class Serializer;

    void save ( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
        rSerializer.save("mCurrentStrain", mCurrentStrain);
        rSerializer.save("mOldStrain", mOldStrain);
        rSerializer.save("mCurrentStress", mCurrentStress);
        rSerializer.save("mOldStress", mOldStress);
        rSerializer.save("mCurrentStateVariables", mCurrentStateVariables);
        rSerializer.save("mOldStateVariables", mOldStateVariables);
        rSerializer.save("mPrestress", mPrestress);
        rSerializer.save("mPrestressFactor", mPrestressFactor);
        rSerializer.save("mCurrentStressZZ", mCurrentStressZZ);
        rSerializer.save("mOldStressZZ", mOldStressZZ);

        rSerializer.save("mElementId", mElementId);
        rSerializer.save("mIntPointIndex", mIntPointIndex);
        rSerializer.save("mNDI", mNDI);
        rSerializer.save("mNSHR", mNSHR);
        rSerializer.save("mNSTATV", mNSTATV);
        rSerializer.save("mPROPS", mPROPS);
        rSerializer.save("mCMNAME", mCMNAME);
    }

    void load ( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
        rSerializer.load("mCurrentStrain", mCurrentStrain);
        rSerializer.load("mOldStrain", mOldStrain);
        rSerializer.load("mCurrentStress", mCurrentStress);
        rSerializer.load("mOldStress", mOldStress);
        rSerializer.load("mCurrentStateVariables", mCurrentStateVariables);
        rSerializer.load("mOldStateVariables", mOldStateVariables);
        rSerializer.load("mPrestress", mPrestress);
        rSerializer.load("mPrestressFactor", mPrestressFactor);
        rSerializer.load("mCurrentStressZZ", mCurrentStressZZ);
        rSerializer.load("mOldStressZZ", mOldStressZZ);

        rSerializer.load("mElementId", mElementId);
        rSerializer.load("mIntPointIndex", mIntPointIndex);
        rSerializer.load("mNDI", mNDI);
        rSerializer.load("mNSHR", mNSHR);
        rSerializer.load("mNSTATV", mNSTATV);
        rSerializer.load("mPROPS", mPROPS);
        rSerializer.load("mCMNAME", mCMNAME);
    }

}; // Class Umat3e

}  // namespace Kratos.

#endif // KRATOS_Umat3e_H_INCLUDED  defined

