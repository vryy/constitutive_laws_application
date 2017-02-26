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
//   Date:                $Date: Oct 31, 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_Umat2_2_H_INCLUDED )
#define  KRATOS_Umat2_2_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

/**
Umat interface KRATOS-ABAQUS
THis interface will always call the 3D version of Umat, even in 2D
*/
class Umat2 : public ConstitutiveLaw
{

public:

    static const int A2K[];
    static const int K2A[];
    static const int PS[];

    // Counted pointer definition
    KRATOS_CLASS_POINTER_DEFINITION( Umat2 );

    // Type Definitions
    typedef ConstitutiveLaw BaseType;

    // Default constructor
    Umat2();

    // Destructor
    virtual ~Umat2();

    //clone
    virtual BaseType::Pointer Clone() const
    {
        BaseType::Pointer p_clone ( new Umat2() );
        return p_clone;
    }

    size_t WorkingSpaceDimension()
    {
        return 3;
    }

    size_t GetStrainSize()
    {
        return 6;
    }

    virtual StrainMeasure GetStrainMeasure()
    {
        return StrainMeasure_Infinitesimal;
    }

    virtual StressMeasure GetStressMeasure()
    {
        return StressMeasure_PK1;
    }

    virtual bool IsIncremental()
    {
        return true;
    }

    virtual bool Has ( const Variable<int>& rThisVariable );

    virtual bool Has ( const Variable<double>& rThisVariable );

    virtual bool Has ( const Variable<Vector>& rThisVariable );

    virtual bool Has ( const Variable<Matrix>& rThisVariable );

    virtual int& GetValue ( const Variable<int>& rThisVariable, int& rValue );

    virtual double& GetValue ( const Variable<double>& rThisVariable, double& rValue );

    virtual Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue );

    virtual Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue );

    virtual void SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual bool ValidateInput ( const Properties& props )
    {
        KRATOS_THROW_ERROR ( std::logic_error, "virtual function Umat2::ValidateInput called", "" );
    }

    virtual void InitializeMaterial ( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues );

    virtual void ResetMaterial ( const Properties& props,
                                 const GeometryType& geom,
                                 const Vector& ShapeFunctionsValues );

    virtual void InitializeSolutionStep ( const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues ,
                                          const ProcessInfo& CurrentProcessInfo );

    virtual void FinalizeSolutionStep ( const Properties& props,
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

    virtual int Check ( const Properties& props,
                        const GeometryType& geom,
                        const ProcessInfo& CurrentProcessInfo );

protected:

private:

    Vector mCurrentStrain;              // store the current strain
    Vector mCurrentStress;              // store the current stress
    Vector mCurrentStateVariables;      // store the current state variables
    Vector mPrestress;
    double mPrestressFactor;
    double mCurrentStressZZ;

    int mElementId;
    int mIntPointIndex;

    // handle to umat library
    void* mp_umat_handle;

    /**
     * wrapper function for calling the UMAT fortran subroutine
     * @param STRESS ......... the vector of stresses
     * @param STATEV ......... the vector of state variables
     * @param DDSDDE ......... the material tangent
     * @param SSE ............
     * @param SPD ............
     * @param SCD ............
     * @param RPL ............
     * @param DDSDDT .........
     * @param DRPLDE .........
     * @param DRPLDT .........
     * @param STRAN .......... the vector of total strains
     * @param DSTRAN ......... the vector of incremental strains
     * @param TIME ........... current time
     * @param DTIME .......... current time increment
     * @param TEMP ........... current temperature
     * @param DTEMP .......... current increment of temperature
     * @param PREDEF .........
     * @param DPRED ..........
     * @param CMNAME ......... material name
     * @param NDI ............ number of direct strain components (3 in 3D)
     * @param NSHR ........... number if shear strain components (3 in 3D)
     * @param NTENS .......... number of stress components (6 in 3D)
     * @param NSTATV ......... number of state variables (size of STATEV)
     * @param PROPS .......... material parameters
     * @param NPROPS ......... number of material paramters (size of PROPS)
     * @param COORDS ......... coordinates of the integration point
     * @param DROT ........... rotation increment (3, 3)
     * @param PNEWDT .........
     * @param CELENT ......... characteristic element length
     * @param DFGRD0 .........
     * @param DFGRD1 .........
     * @param NOEL ........... element number
     * @param NPT ............ integration point number
     * @param KSLAY ..........
     * @param KSPT ...........
     * @param KSTEP .......... step number
     * @param KINC ........... increment number
     */
    void (*Umat)(double* STRESS, double* STATEV, double** DDSDDE, double* SSE, double* SPD, double* SCD,
                 double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
                 double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                 char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
                 double* COORDS, double** DROT, double* PNEWDT, double* CELENT, double** DFGRD0,
                 double** DFGRD1, int* NOEL, int* NPT, int* KSLAY, int* KSPT, int* KSTEP, int* KINC);

    // serialization
    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const
    {
        rSerializer.save ( "name", "Umat2" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

    virtual void load ( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

}; // Class Umat2

}  // namespace Kratos.

#endif // KRATOS_Umat2_H_INCLUDED  defined 

