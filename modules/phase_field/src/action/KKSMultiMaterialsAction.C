/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "KKSMultiMaterialsAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"
#include "libmesh/string_to_enum.h"

template <>
InputParameters
validParams<KKSMultiMaterialsAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up KKS multiphase parameters");
  params.addRequiredParam<unsigned int>(
      "op_num", "Specifies the total number of phases");
  params.addRequiredParam<std::vector<Real>>(
      "eq_concentration", "Specifies the equilibrium phase concentration of all phases in the binary alloy");
  return params;
}

KKSMultiMaterialsAction::KKSMultiMaterialsAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _eq_c(getParam<std::vector<Real>>("eq_concentration"))
{
}

void
KKSMultiMaterialsAction::act()
{

  std::vector<VariableName> v_etanames;
  v_etanames.resize(_op_num); 

  for (unsigned int j = 0; j < _op_num; j++) {
      v_etanames[j] = "eta" + Moose::stringify(j + 1);
  }

  std::vector<std::string> var_name_mater;
  var_name_mater.resize(2); 
  var_name_mater[0] = "D";

  std::vector<VariableName> var_name_c;
  var_name_c.resize(1); 
  std::vector<VariableName> var_name_eta;
  var_name_eta.resize(1); 
  std::string func_expression;

  for (unsigned int j = 0; j < _op_num; ++j)
  {

    var_name_eta[0] = "eta" + Moose::stringify(j + 1);

    //
    // Set up switching function
    //
    {
      InputParameters params = _factory.getValidParams("SwitchingFunctionMultiPhaseMaterial");
      std::string h_name = "h" + Moose::stringify(j + 1);
      params.set<MaterialPropertyName>("h_name") = h_name;
      //params.set<std::string>("f_name") = h_name;
      params.set<std::vector<VariableName>>("all_etas") = v_etanames;
      params.set<std::vector<VariableName>>("phase_etas") = var_name_eta;
      std::string kernel_name = "h" + Moose::stringify(j + 1);
      _problem->addMaterial("SwitchingFunctionMultiPhaseMaterial", kernel_name, params);
    }

    //
    // Set up free energy material parameters
    //

    {
      InputParameters params = _factory.getValidParams("DerivativeParsedMaterial");
      //      std::string f_name = "F" + Moose::stringify(j + 1);
      params.set<std::string>("f_name") = "F" + Moose::stringify(j + 1);
      var_name_c[0] = "c" + Moose::stringify(j + 1);
      params.set<std::vector<VariableName>>("args") = var_name_c;
      if (j == 0)  func_expression = "20*(c" + Moose::stringify(j + 1) + "-" + Moose::stringify(_eq_c[0]) + ")^2";
      else func_expression = "20*(c" + Moose::stringify(j + 1) + "-" + Moose::stringify(_eq_c[1]) + ")^2";
      params.set<std::string>("function") = func_expression;
      std::string kernel_name = "f" + Moose::stringify(j + 1);
      _problem->addMaterial("DerivativeParsedMaterial", kernel_name, params);
    }

    //
    // Set up coefficients for diffusion equation
    //

    {
      InputParameters params = _factory.getValidParams("DerivativeParsedMaterial");
      std::string h_name = "h" + Moose::stringify(j + 1);
      var_name_mater[1] = h_name;
      params.set<std::vector<std::string>>("material_property_names") = var_name_mater;
      params.set<std::string>("function") = "D*" + h_name;
      //std::string f_name = "Dh" + Moose::stringify(j + 1);
      params.set<std::string>("f_name") = "Dh" + Moose::stringify(j + 1);
      std::string kernel_name = "Dh" + Moose::stringify(j + 1);
      _problem->addMaterial("DerivativeParsedMaterial", kernel_name, params);
    }


    //
    // Set up coefficients for Barrier functions
    //
    {
      InputParameters params = _factory.getValidParams("BarrierFunctionMaterial");
      params.set<MooseEnum>("g_order") = "SIMPLE";
      params.set<std::vector<VariableName>>("eta") = var_name_eta;
      params.set<std::string>("function_name") = "g" + Moose::stringify(j + 1);
      std::string kernel_name = "g" + Moose::stringify(j + 1);
      _problem->addMaterial("BarrierFunctionMaterial", kernel_name, params);
    }


  }

}
