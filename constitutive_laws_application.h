//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:54:44 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "structural_application.h"
#include "constitutive_laws/umat3e.h"

namespace Kratos
{

///@name Kratos Globals
///@{

// Variables definition

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) KratosConstitutiveLawsApplication : public KratosApplication
{

public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosExternalSolversApplication
    KRATOS_CLASS_POINTER_DEFINITION ( KratosConstitutiveLawsApplication );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConstitutiveLawsApplication();

    /// Destructor.
    ~KratosConstitutiveLawsApplication() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosConstitutiveLawsApplication";
    }

    /// Print information about this object.
    void PrintInfo ( std::ostream& rOStream ) const override
    {
        rOStream << Info();
        PrintData ( rOStream );
    }

    ///// Print object's data.
    void PrintData ( std::ostream& rOStream ) const override
    {
        KRATOS_WATCH ( "in KratosConstitutiveLawsApplication" );
        KRATOS_WATCH ( KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData ( rOStream );
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData ( rOStream );
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData ( rOStream );
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    const Umat3e mUmat3e;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosConstitutiveLawsApplication& operator= ( KratosConstitutiveLawsApplication const& rOther );

    /// Copy constructor.
    KratosConstitutiveLawsApplication ( KratosConstitutiveLawsApplication const& rOther );


    ///@}

}; // Class KratosConstitutiveLawsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_EXTERNAL_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED  defined 


