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
//   Date:                $Date: Nov 3, 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_OPENSEES_MAT_H_INCLUDED )
#define  KRATOS_OPENSEES_MAT_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

#include "SRC/OPS_Globals.h"
#include "SRC/handler/StandardStream.h"
#include "SRC/actor/channel/Channel.h"
#include "SRC/material/Material.h"
#include "SRC/material/nD/NDMaterial.h"
#include "SRC/matrix/Vector.h"
#include "SRC/matrix/Matrix.h"
#include "SRC/recorder/response/MaterialResponse.h"

class MyOpenSeesChannel : public Channel
{
public:
    MyOpenSeesChannel() : Channel()
    {}

    virtual ~MyOpenSeesChannel()
    {}

    virtual char *addToProgram(void) {return 0;}
    virtual int setUpConnection(void) {return 0;}
    virtual int setNextAddress(const ChannelAddress &theAddress) {return 0;}
    virtual ChannelAddress *getLastSendersAddress(void) {return 0;}

    virtual int sendObj(int commitTag, MovableObject &theObject, ChannelAddress *theAddress =0) {return 0;}
    virtual int recvObj(int commitTag, MovableObject &theObject, FEM_ObjectBroker &theBroker, ChannelAddress *theAddress =0) {return 0;}
    virtual int sendMsg(int dbTag, int commitTag, const Message &theMessage, ChannelAddress *theAddress =0) {return 0;}
    virtual int recvMsg(int dbTag, int commitTag, Message &theMessage, ChannelAddress *theAddress =0) {return 0;}
    virtual int recvMsgUnknownSize(int dbTag, int commitTag, Message &theMessage, ChannelAddress *theAddress =0) {return 0;}
    virtual int sendMatrix(int dbTag, int commitTag, const Matrix &theMatrix, ChannelAddress *theAddress =0) {return 0;}
    virtual int recvMatrix(int dbTag, int commitTag, Matrix &theMatrix, ChannelAddress *theAddress =0) {return 0;}

    virtual int sendVector(int dbTag, int commitTag, const Vector& theVector, ChannelAddress* theAddress = 0)
    {
        mTheVector = theVector;
        return 0;
    }

    virtual int recvVector(int dbTag, int commitTag, Vector& theVector, ChannelAddress* theAddress =0)
    {
        theVector = mTheVector;
        return 0;
    }

    virtual int sendID(int dbTag, int commitTag, const ID &theID, ChannelAddress *theAddress =0) {return 0;}
    virtual int recvID(int dbTag, int commitTag, ID &theID, ChannelAddress *theAddress =0) {return 0;}

    void copyVector(Kratos::Vector& rValue)
    {
        if(rValue.size() != mTheVector.Size())
            rValue.resize(mTheVector.Size());
        for(int i = 0; i < mTheVector.Size(); ++i)
            rValue(i) = mTheVector(i);
    }

private:
    Vector mTheVector;
};

void OpenSees_setTrialStrain(boost::shared_ptr<NDMaterial> pMat, const Kratos::Vector& rValue);

void OpenSees_getStress(boost::shared_ptr<NDMaterial> pMat, Kratos::Vector& rValue);

void OpenSees_getTangent(boost::shared_ptr<NDMaterial> pMat, Kratos::Matrix& rValue);

void OpenSees_getInitialTangent(boost::shared_ptr<NDMaterial> pMat, Kratos::Matrix& rValue);

namespace Kratos
{

class OpenSeesMat : public ConstitutiveLaw
{

public:

    // Counted pointer definition
    KRATOS_CLASS_POINTER_DEFINITION( OpenSeesMat );

    // Type Definitions
    typedef ConstitutiveLaw BaseType;

    // Default constructor
    OpenSeesMat();

    // Destructor
    virtual ~OpenSeesMat();

    //clone
    BaseType::Pointer Clone() const final
    {
        BaseType::Pointer p_clone ( new OpenSeesMat() );
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
        return StressMeasure_PK1;
    }

    bool IsIncremental() final
    {
        return false;
    }

    bool Has ( const Variable<int>& rThisVariable ) final;

    bool Has ( const Variable<double>& rThisVariable ) final;

    bool Has ( const Variable<Vector>& rThisVariable ) final;

    bool Has ( const Variable<Matrix>& rThisVariable ) final;

    int& GetValue ( const Variable<int>& rThisVariable, int& rValue ) final;

    double& GetValue ( const Variable<double>& rThisVariable, double& rValue ) final;

    Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue ) final;

    Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue ) final;

    void SetValue( const Variable<int>& rVariable, const int& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    void SetValue( const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo) final;

    bool ValidateInput ( const Properties& props ) final
    {
        KRATOS_THROW_ERROR ( std::logic_error, "function OpenSeesMat::ValidateInput called", "" );
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
            bool SaveInternalVariables = true );

    void FinalizeNonLinearIteration ( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo ) final;

    int Check ( const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo ) const final;

private:

    Vector mCurrentStrain;              // store the current strain
    Vector mCurrentStress;              // store the current stress
    Vector mCurrentStateVariables;      // store the current state variables

    boost::shared_ptr<NDMaterial> mpConstitutiveLaw;
    boost::shared_ptr<MyOpenSeesChannel> mpChannel;

    std::size_t mParentElementId;
    int mIntegrationPointId;

    // serialization
    friend class Serializer;

    void save ( Serializer& rSerializer ) const override
    {
        rSerializer.save ( "name", "OpenSeesMat" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

    void load ( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
    }

}; // Class OpenSeesMat

}  // namespace Kratos.

#endif // KRATOS_OpenSeesMat_H_INCLUDED  defined

