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


namespace Kratos
{

/**
Umat interface KRATOS-ABAQUS
THis interface will call respective Umat version depending on the parameters given. Therefore Umat must implement both 2D and 3D.
This version of Umat is the same as Umat3 but is customized for ekate
*/
class Umat3e : public ConstitutiveLaw
{

public:

    static const int A2K[];
    static const int K2A[];
    static const int PS[];

    // Counted pointer definition
    KRATOS_CLASS_POINTER_DEFINITION( Umat3e );

    // Type Definitions
    typedef ConstitutiveLaw BaseType;

    // Default constructor
    Umat3e();

    // Destructor
    virtual ~Umat3e();

    //clone
    virtual BaseType::Pointer Clone() const
    {
        BaseType::Pointer p_clone ( new Umat3e() );
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

    virtual std::string& GetValue ( const Variable<std::string>& rThisVariable, std::string& rValue );

    virtual void SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual void SetValue( const Variable<std::string>& rVariable, const std::string& rValue, const ProcessInfo& rCurrentProcessInfo);

    virtual bool ValidateInput ( const Properties& props )
    {
        KRATOS_THROW_ERROR ( std::logic_error, "virtual function Umat3e::ValidateInput called", "" );
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

    #ifndef KRATOS_UMAT_LIBRARY_IS_PROVIDED
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
    #endif

    /// Reset only the material state to initial state
    void ResetState();

    // serialization
    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const
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

    virtual void load ( Serializer& rSerializer )
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

