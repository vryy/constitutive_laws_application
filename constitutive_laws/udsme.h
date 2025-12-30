//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Apr 27, 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_UDSM_E_H_INCLUDED )
#define  KRATOS_UDSM_E_H_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_laws/udsm.h"

/**
    Interface class to incorporate the user-defined soil models from Plaxis to Kratos
    This version of Udsm is the same as UDSM but is customized for ekate, where the material parameters are read from matfile
*/
namespace Kratos
{

class UDSMe : public UDSM
{
    public:
        // Counted pointer for UDSMe
        KRATOS_CLASS_POINTER_DEFINITION( UDSMe );

        // Type Definitions
        typedef UDSM BaseType;

        // Default constructor
        UDSMe() : BaseType() {}

        // Destructor
        virtual ~UDSMe() {}

        // clone
        ConstitutiveLaw::Pointer Clone() const override
        {
            BaseType::Pointer p_clone ( new UDSMe() );
            return p_clone;
        }

        bool Has ( const Variable<int>& rThisVariable ) final;

        bool Has ( const Variable<bool>& rThisVariable ) final;

        bool Has ( const Variable<double>& rThisVariable ) final;

        bool Has ( const Variable<Vector>& rThisVariable ) final;

        // bool Has ( const Variable<std::string>& rThisVariable ) final;

        int& GetValue ( const Variable<int>& rThisVariable, int& rValue ) final;

        bool& GetValue ( const Variable<bool>& rThisVariable, bool& rValue ) final;

        double& GetValue ( const Variable<double>& rThisVariable, double& rValue ) final;

        Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue ) final;

        std::string& GetValue ( const Variable<std::string>& rThisVariable, std::string& rValue ) final;

        void SetValue ( const Variable<int>& rThisVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void SetValue ( const Variable<bool>& rThisVariable, const bool& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void SetValue ( const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void SetValue ( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo ) final;

        void SetValue( const Variable<std::string>& rVariable, const std::string& rValue, const ProcessInfo& rCurrentProcessInfo) final;

        void InitializeMaterial ( const Properties& props,
                                  const GeometryType& geom,
                                  const Vector& ShapeFunctionsValues ) final;

        int Check ( const Properties& props,
                    const GeometryType& geom,
                    const ProcessInfo& CurrentProcessInfo ) const override;

    protected:

        int mModelNumber;                   // to choose the soil model number, link with iMod variable

        std::string mLibName;
        std::string mName;
        Vector mProps;
        int mIsUndr;
        double mBulkW;

        //serialization
        friend class Serializer;

        void save ( Serializer& rSerializer ) const override
        {
            rSerializer.save ( "name", "UDSMe" );
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, BaseType )
        }

        void load ( Serializer& rSerializer ) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, BaseType )
        }
}; // Class UDSMe

class UDSMeImplex : public UDSMe
{
    public:
        // Counted pointer for UDSMeImplex
        KRATOS_CLASS_POINTER_DEFINITION( UDSMeImplex );

        // Type Definitions
        typedef UDSMe BaseType;

        // Default constructor
        UDSMeImplex() : BaseType() {}

        // Destructor
        virtual ~UDSMeImplex() {}

        // clone
        ConstitutiveLaw::Pointer Clone() const final
        {
            ConstitutiveLaw::Pointer p_clone ( new UDSMeImplex() );
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

class UDSMeImplicit : public UDSMe
{
    public:
        // Counted pointer for UDSMeImplicit
        KRATOS_CLASS_POINTER_DEFINITION( UDSMeImplicit );

        // Type Definitions
        typedef UDSMe BaseType;

        // Default constructor
        UDSMeImplicit() : BaseType() {}

        // Destructor
        virtual ~UDSMeImplicit() {}

        // clone
        ConstitutiveLaw::Pointer Clone() const final
        {
            ConstitutiveLaw::Pointer p_clone ( new UDSMeImplicit() );
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

#endif // KRATOS_UDSM_E_H_INCLUDED  defined
