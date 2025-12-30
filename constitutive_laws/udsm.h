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
 * Interface class to incorporate the user-defined soil models from Plaxis to Kratos
 */
namespace Kratos
{

class UDSM : public ConstitutiveLaw
{
    public:
        // Counted pointer for UDSM
        KRATOS_CLASS_POINTER_DEFINITION( UDSM );

        // Debug Info for UDSM
        struct DebugInfo
        {
            double time, delta_time;
            int step, iteration, element_id, point_id;
        };

        struct Input
        {
            Input (Vector& v1, Vector& v2, Vector& v3, double& d1)
            : LastStrain(v1), LastStress(v2), LastStateVariables(v3),
              LastExcessPorePressure(d1)
            {}

            Vector &LastStrain, &LastStress, &LastStateVariables;
            double &LastExcessPorePressure;
        };

        struct Output
        {
            Output (Vector& v1, Vector& v2, Vector& v3, double& d1, int& i1, Matrix& m1)
            : CurrentStrain(v1), CurrentStress(v2), CurrentStateVariables(v3),
              CurrentExcessPorePressure(d1), PlasticState(i1),
              AlgorithmicTangent(m1)
            {}

            Vector &CurrentStrain, &CurrentStress, &CurrentStateVariables;
            double &CurrentExcessPorePressure;
            int &PlasticState;
            Matrix &AlgorithmicTangent;
        };

        // Type Definitions
        typedef ConstitutiveLaw BaseType;

        // Default constructor
        UDSM();

        // Destructor
        ~UDSM() override;

        // clone
        BaseType::Pointer Clone() const override
        {
            BaseType::Pointer p_clone ( new UDSM() );
            return p_clone;
        }

        ConstitutiveLaw::StrainMeasure GetStrainMeasure() override
        {
            return StrainMeasure_Infinitesimal;
        }

        ConstitutiveLaw::StressMeasure GetStressMeasure() override
        {
            return StressMeasure_Cauchy;
        }

        void GetLawFeatures(Features& rFeatures) override
        {
            rFeatures.SetStrainMeasure(this->GetStrainMeasure());
        }

        bool IsIncremental() override
        {
            return true;
        }

        bool Has ( const Variable<double>& rThisVariable ) override;

        bool Has ( const Variable<Vector>& rThisVariable ) override;

        double& GetValue ( const Variable<double>& rThisVariable, double& rValue ) override;

        Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue ) override;

        Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;

        void SetValue ( const Variable<int>& rThisVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

        void SetValue ( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

        void SetValue ( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

        void ResetMaterial ( const Properties& props,
                             const GeometryType& geom,
                             const Vector& ShapeFunctionsValues ) override;

        int Check ( const Properties& props,
                    const GeometryType& geom,
                    const ProcessInfo& CurrentProcessInfo ) const override;

        void InitializeMaterial ( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ) override;

        void InitializeSolutionStep ( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues ,
                                      const ProcessInfo& CurrentProcessInfo ) override;

        void InitializeNonLinearIteration ( const Properties& props,
                                            const GeometryType& geom,
                                            const Vector& ShapeFunctionsValues,
                                            const ProcessInfo& CurrentProcessInfo ) override;

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
                                          const ProcessInfo& CurrentProcessInfo ) override;

        void FinalizeSolutionStep ( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues ,
                                    const ProcessInfo& CurrentProcessInfo ) override;

    protected:

        unsigned int mElemId;
        unsigned int mGaussId;
        unsigned int mStep;
        unsigned int mIter;

        Vector mCurrentStrain;              // to store the current strain
        Vector mLastStrain;                 // to store the converged strain in the last step
        Vector mCurrentStress;              // link with Sig variable   // this is constitutive stress
        Vector mLastStress;                 // link with Sig0 variable  // this is constitutive stress
        Vector mPrestress;
        double mPrestressFactor;
        double mCurrentExcessPorePressure;  // link with Swp variable
        double mLastExcessPorePressure;     // link with Swp0 variable
        Vector mCurrentStateVariables;      // link with StVar variable
        Vector mLastStateVariables;         // link with StVar0 variable
        int mPlasticState;                  // link with iPl variable
        int mModelNumber;                   // to choose the soil model number, link with iMod variable
        int mIsUndr;                        // flag for undrained state
        double mBulkW;                      // bulk modulus of water
        Vector mProps;                      // vector of material properties

        #ifndef KRATOS_UDSM_LIBRARY_IS_PROVIDED
        void* mp_udsm_handle; // handle to udsm library
        #endif
        udsm_t UserMod; // pointer to UDSM function

        void InitializeMaterial( int ModelNumber, int IsUndrained, const Vector& Props, bool need_determine_internal_params = false );

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

        static void VectorTo3DVector(const Vector& vector, Vector& vector_3d);
        static void Vector3DToVector(const Vector& vector_3d, Vector& vector);
        static void Vector1DToMatrix(const Vector& D, Matrix& A, int non_sym);
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
        ~UDSMImplex() override {}

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

        static void StressIntegration ( const Vector& StrainVector, Vector& StressVector,
                                        const UDSM::Input& rInput, UDSM::Output& rOutput,
                                        const udsm_t UserFun, double omega, const UDSM::DebugInfo& pinfo,
                                        int ModelNumber, int IsUndr, double BulkW, const Vector& Props );
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
        ~UDSMImplicit() override {}

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

        static void StressIntegration ( const UDSM::Input& rInput, UDSM::Output& rOutput,
                                        const udsm_t UserFun, const UDSM::DebugInfo& pinfo,
                                        int ModelNumber, int IsUndr, double BulkW, const Vector& Props,
                                        int CalculateTangent = 1 );

        static int ComputeNumericalTangent( const UDSM::Input& rInput, const UDSM::Output& rRefOutput,
                                            Matrix& Tangent,
                                            const udsm_t UserFun, const UDSM::DebugInfo& pinfo,
                                            int ModelNumber, int IsUndr, double BulkW, const Vector& Props,
                                            double epsilon );
};

}  // namespace Kratos.

#endif // KRATOS_UDSM_H_INCLUDED  defined
