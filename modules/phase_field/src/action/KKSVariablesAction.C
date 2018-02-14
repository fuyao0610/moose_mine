/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "KKSVariablesAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"
#include "libmesh/string_to_enum.h"


template <>
InputParameters
validParams<KKSVariablesAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up order parameter variables for a polycrystal sample");
  params.addParam<std::string>(
      "family", "LAGRANGE", "Specifies the family of FE shape functions to use for this variable");
  params.addParam<std::string>(
      "order", "FIRST", "Specifies the order of the FE shape function to use for this variable");
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");
  params.addRequiredParam<unsigned int>("op_num",
                                        "specifies the number of order parameters to create");
  return params;
}

KKSVariablesAction::KKSVariablesAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num"))
{
}

void
KKSVariablesAction::act()
{
#ifdef DEBUG
  Moose::err << "Inside the KKSVariablesAction Object\n"
             << "\torder: " << getParam<std::string>("order")
             << "\tfamily: " << getParam<std::string>("family") << std::endl;
#endif

  _problem->addVariable(
			"c",
			FEType(Utility::string_to_enum<Order>(getParam<std::string>("order")),
			       Utility::string_to_enum<FEFamily>(getParam<std::string>("family"))),
			getParam<Real>("scaling"));

  _problem->addVariable(
			"lambda",
			FEType(Utility::string_to_enum<Order>(getParam<std::string>("order")),
			       Utility::string_to_enum<FEFamily>(getParam<std::string>("family"))),
			getParam<Real>("scaling"));


  // Loop through the number of order parameters
  for (unsigned int op = 0; op < _op_num; op++)
  {
    // Create variable names
    std::string eta_name = "eta" + Moose::stringify(op + 1);;
    std::string c_name = "c" + Moose::stringify(op + 1);;
    //    std::stringstream out;
    //out << op + 1;
    //eta_name.append(out.str());
    //c_name.append(out.str());

    _problem->addVariable(
        eta_name,
        FEType(Utility::string_to_enum<Order>(getParam<std::string>("order")),
               Utility::string_to_enum<FEFamily>(getParam<std::string>("family"))),
        getParam<Real>("scaling"));

    _problem->addVariable(
        c_name,
        FEType(Utility::string_to_enum<Order>(getParam<std::string>("order")),
               Utility::string_to_enum<FEFamily>(getParam<std::string>("family"))),
        getParam<Real>("scaling"));

  }
}

