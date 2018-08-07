//
//   Project Name:        Kratos
//   Last modified by:    $Author: nagel $
//   Date:                $Date: 2009-01-12 08:17:36 $
//   Revision:            $Revision: 1.10 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_constitutive_laws_to_python.h"

#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/condition.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"

//constitutive laws
#include "includes/constitutive_law.h"
#include "constitutive_laws/umat.h"
#include "constitutive_laws/udsm.h"
#include "constitutive_laws/umat2.h"
#include "constitutive_laws/umat3.h"
#include "constitutive_laws/umat3e.h"

#ifdef KRATOS_USE_OPENSEES
#include "constitutive_laws/opensees_mat.h"
#endif

namespace Kratos
{

namespace Python
{

using namespace boost::python;

typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;


void  AddConstitutiveLawsToPython()
{
    class_< Umat, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Umat", init<>() )
    ;

    class_< UDSMImplex, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "UDSMImplex", init<>() )
    ;

    class_< UDSMImplicit, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "UDSMImplicit", init<>() )
    ;

    class_< Umat2, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Umat2", init<>() )
    ;

    class_< Umat3, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Umat3", init<>() )
    ;

    class_< Umat3e, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "Umat3e", init<>() )
    ;

    #ifdef KRATOS_USE_OPENSEES
    class_< OpenSeesMat, bases< ConstitutiveLawBaseType >, boost::noncopyable >
    ( "OpenSeesMat", init<>() )
    ;
    #endif
}
}  // namespace Python.
} // Namespace Kratos
