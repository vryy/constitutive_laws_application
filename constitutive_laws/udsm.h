//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Apr 27, 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_UDSM_H_INCLUDED )
#define  KRATOS_UDSM_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "constitutive_laws_application_variables.h"

/**
    Interface class to incorporate the user-defined soil models from Plaxis to Kratos
*/
namespace Kratos
{

class UDSM : public ConstitutiveLaw
{
    public:
        // Counted pointer for UDSM
        KRATOS_CLASS_POINTER_DEFINITION( UDSM );

        // Type Definitions
        typedef ConstitutiveLaw BaseType;

        // Default constructor
        UDSM();

        // Destructor
        virtual ~UDSM();

        // clone
        BaseType::Pointer Clone() const override
        {
            BaseType::Pointer p_clone ( new UDSM() );
            return p_clone;
        }

        ConstitutiveLaw::StrainMeasure GetStrainMeasure() final
        {
            return StrainMeasure_Infinitesimal;
        }

        ConstitutiveLaw::StressMeasure GetStressMeasure() final
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

        bool Has ( const Variable<double>& rThisVariable ) final;

        bool Has ( const Variable<Vector>& rThisVariable ) final;

        double& GetValue ( const Variable<double>& rThisVariable, double& rValue ) final;

        Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue ) final;

        void SetValue ( const Variable<int>& rThisVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void SetValue ( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void ResetMaterial ( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) final;

        int Check ( const Properties& props,
                    const GeometryType& geom,
                    const ProcessInfo& CurrentProcessInfo ) const final;

        void InitializeMaterial ( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ) final;

        void InitializeSolutionStep ( const Properties& props,
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
                                         bool SaveInternalVariables = true ) override;

        void FinalizeNonLinearIteration ( const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues,
                                          const ProcessInfo& CurrentProcessInfo ) final;

        void FinalizeSolutionStep ( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo ) final;

    protected:

        unsigned int mElemId;
        unsigned int mGaussId;
        unsigned int mStep;
        unsigned int mIter;

        Vector mCurrentStrain;              // to store the current strain
        Vector mLastStrain;                 // to store the converged strain in the last step
        Vector mCurrentStress;              // link with Sig variable   // this is constitutive stress
        Vector mLastStress;                 // link with Sig0 variable  // this is constitutive stress
        double mCurrentExcessPorePressure;  // link with Swp variable
        double mLastExcessPorePressure;     // link with Swp0 variable
        Vector mCurrentStateVariables;      // link with StVar variable
        Vector mLastStateVariables;         // link with StVar0 variable
        int mPlasticState;                  // link with iPl variable
        int mModelNumber;                   // to choose the soil model number, link with iMod variable

        #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
        static unsigned long long minstances; // values to store the instances of this constitutive law
        static void* mp_udsm_handle; // handle to udsm library
        static udsm_t UserMod;
        #else
        udsm_t UserMod;
        #endif

        //serialization
        friend class Serializer;

        void save ( Serializer& rSerializer ) const override
        {
            rSerializer.save ( "name", "UDSM" );
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
        }

        void load ( Serializer& rSerializer ) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
        }

        void VectorTo3DVector(const Vector& vector, Vector& vector_3d) const;
        void Vector3DToVector(const Vector& vector_3d, Vector& vector) const;
        void Vector1DToMatrix(const Vector& D, Matrix& A, const int& non_sym) const;
}; // Class UDSM

class UDSMImplex : public UDSM
{
    public:
        // Counted pointer for UDSMImplex
        KRATOS_CLASS_POINTER_DEFINITION( UDSMImplex );

        // Type Definitions
        typedef UDSM BaseType;

        // Default constructor
        UDSMImplex() : BaseType() {}

        // Destructor
        virtual ~UDSMImplex() {}

        // clone
        ConstitutiveLaw::Pointer Clone() const final
        {
            ConstitutiveLaw::Pointer p_clone ( new UDSMImplex() );
            return p_clone;
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
};

class UDSMImplicit : public UDSM
{
    public:
        // Counted pointer for UDSMImplicit
        KRATOS_CLASS_POINTER_DEFINITION( UDSMImplicit );

        // Type Definitions
        typedef UDSM BaseType;

        // Default constructor
        UDSMImplicit() : BaseType() {}

        // Destructor
        virtual ~UDSMImplicit() {}

        // clone
        ConstitutiveLaw::Pointer Clone() const final
        {
            ConstitutiveLaw::Pointer p_clone ( new UDSMImplicit() );
            return p_clone;
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
};

}  // namespace Kratos.

#endif // KRATOS_UDSM_H_INCLUDED  defined

