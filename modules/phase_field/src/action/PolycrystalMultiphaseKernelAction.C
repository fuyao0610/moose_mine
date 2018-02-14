/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PolycrystalMultiphaseKernelAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

template <>
InputParameters
validParams<PolycrystalMultiphaseKernelAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription(
      "Set up Multiphase and TimeDerivative kernels");
  params.addRequiredParam<unsigned int>(
      "op_num", "specifies the total number of grains (deformed + recrystallized) to create");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addParam<VariableName>("c", "Name of coupled concentration variable");
  params.addParam<Real>("en_ratio", 1.0, "Ratio of surface to GB energy");
  params.addParam<unsigned int>("ndef", 0, "specifies the number of deformed grains to create");
  params.addParam<bool>("implicit", true, "Whether kernels are implicit or not");
  params.addParam<VariableName>("T", "Name of temperature variable");
  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  return params;
}

PolycrystalMultiphaseKernelAction::PolycrystalMultiphaseKernelAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _implicit(getParam<bool>("implicit"))
{
}

void
PolycrystalMultiphaseKernelAction::act()
{
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    //
    // Create variable names
    //

    std::string var_name = _var_name_base + Moose::stringify(op);
    std::vector<VariableName> v;
    v.resize(_op_num);

    unsigned int ind = 0;
    for (unsigned int j = 0; j < _op_num; ++j)
        v[ind++] = _var_name_base + Moose::stringify(j);

    //
    // Set up Multiphase kernels
    //

    {
      InputParameters params = _factory.getValidParams("Multiphase");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<std::vector<VariableName>>("v") = v;
      params.set<unsigned int>("op") = op;
      params.set<bool>("implicit") = _implicit;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");
      if (isParamValid("T"))
        params.set<std::vector<VariableName>>("T") = {getParam<VariableName>("T")};

      std::string kernel_name = "ACBulk_" + var_name;
      _problem->addKernel("Multiphase", kernel_name, params);
    }

    //
    // Set up ACInterface kernels
    //

    {
      InputParameters params = _factory.getValidParams("ACInterface");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<bool>("implicit") = getParam<bool>("implicit");
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "ACInt_" + var_name;
      _problem->addKernel("ACInterface", kernel_name, params);
    }

   
    //
    // Set up TimeDerivative kernels
    //

    {
      InputParameters params = _factory.getValidParams("TimeDerivative");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<bool>("implicit") = true;
      params.set<bool>("use_displaced_mesh") = getParam<bool>("use_displaced_mesh");

      std::string kernel_name = "IE_" + var_name;
      _problem->addKernel("TimeDerivative", kernel_name, params);
    }

  }
}
