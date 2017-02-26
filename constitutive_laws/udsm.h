/*
 ==============================================================================
 KratosStructuralApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 janosch.stascheit@rub.de
 nagel@sd.rub.de
 - CIMNE (International Center for Numerical Methods in Engineering),
 Gran Capita' s/n, 08034 Barcelona, Spain
 - Ruhr-University Bochum, Institute for Structural Mechanics, Germany


 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 ==============================================================================
 */
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
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"

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
            virtual BaseType::Pointer Clone() const
            {
                BaseType::Pointer p_clone ( new UDSM() );
                return p_clone;
            }

            virtual std::size_t WorkingSpaceDimension()
            {
                return 3;
            }

            virtual std::size_t GetStrainSize()
            {
                return 6;
            }

            virtual bool Has ( const Variable<double>& rThisVariable );

            virtual bool Has ( const Variable<Vector>& rThisVariable );

            virtual double& GetValue ( const Variable<double>& rThisVariable, double& rValue );

            virtual Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue );

            virtual void SetValue ( const Variable<int>& rThisVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo );

            virtual void SetValue ( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo );

            StrainMeasure GetStrainMeasure()
            {
                return StrainMeasure_Infinitesimal;
            }

            StressMeasure GetStressMeasure()
            {
                return StressMeasure_Cauchy;
            }

            bool IsIncremental()
            {
                return true;
            }

            virtual void ResetMaterial ( const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues );

            virtual int Check ( const Properties& props,
                                const GeometryType& geom,
                                const ProcessInfo& CurrentProcessInfo );

            virtual void InitializeMaterial ( const Properties& props,
                                              const GeometryType& geom,
                                              const Vector& ShapeFunctionsValues );

            virtual void InitializeSolutionStep ( const Properties& props,
                                                  const GeometryType& geom,
                                                  const Vector& ShapeFunctionsValues ,
                                                  const ProcessInfo& CurrentProcessInfo );

            virtual void InitializeNonLinearIteration ( const Properties& props,
                                                        const GeometryType& geom,
                                                        const Vector& ShapeFunctionsValues,
                                                        const ProcessInfo& CurrentProcessInfo );

            virtual void CalculateMaterialResponse ( const Vector& StrainVector,
                                                     const Matrix& DeformationGradient,
                                                     Vector& StressVector,
                                                     Matrix& AlgorithmicTangent,
                                                     const ProcessInfo& CurrentProcessInfo,
                                                     const Properties& props,
                                                     const GeometryType& geom,
                                                     const Vector& ShapeFunctionsValues,
                                                     bool CalculateStresses = true,
                                                     int CalculateTangent = true,
                                                     bool SaveInternalVariables = true );

            virtual void FinalizeNonLinearIteration ( const Properties& props,
					                                  const GeometryType& geom,
					                                  const Vector& ShapeFunctionsValues,
					                                  const ProcessInfo& CurrentProcessInfo );

            virtual void FinalizeSolutionStep ( const Properties& props,
                                                const GeometryType& geom,
                                                const Vector& ShapeFunctionsValues ,
                                                const ProcessInfo& CurrentProcessInfo );

            virtual void CalculateCauchyStresses( Vector& Cauchy_StressVector,
                                                  const Matrix& F,
                                                  const Vector& PK2_StressVector,
                                                  const Vector& GreenLagrangeStrainVector );

        private:

            Vector mCurrentStrain;              // to store the current strain (for computation of incremental strain)
            Vector mCurrentStress;              // link with Sig variable
            double mCurrentExcessPorePressure;  // link with Swp variable
            Vector mCurrentStateVariables;      // link with StVar variable
            int mPlasticState;                  // link with iPl variable
            int mModelNumber;                   // to choose the soil model number, link with iMod variable

            void* mp_udsm_handle; // handle to udsm library
            void (*UserMod)(int* IDTask, int* iMod, int* IsUndr, int* iStep, int* iTer, int* iEl,
                            int* Int, double* X, double* Y, double* Z, double* Time0, double* dTime,
                            double* Props, double* Sig0, double* Swp0, double* StVar0, double* dEps,
                            double* D, double* BulkW, double* Sig, double* Swp, double* StVar,
                            int* ipl, int* nStat, int* NonSym, int* iStrsDep, int* iTimeDep,
                            int* iTang, int* iAbort);

            //serialization
            friend class Serializer;

            virtual void save ( Serializer& rSerializer ) const
            {
                rSerializer.save ( "name", "UDSM" );
                KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
            }

            virtual void load ( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
            }

    }; // Class UDSM

}  // namespace Kratos.

#endif // KRATOS_UDSM_H_INCLUDED  defined 

