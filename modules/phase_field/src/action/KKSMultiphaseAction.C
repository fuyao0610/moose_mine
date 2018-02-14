/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "KKSMultiphaseAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

template <>
InputParameters
validParams<KKSMultiphaseAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription(
      "Set up concentration diffusion, Allen-Cahn, Lagrange multiplier, phase constraint, and concentration constraints kernels");
  params.addRequiredParam<unsigned int>(
      "op_num", "specifies the total number of phases to create");
  //params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  
  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  return params;
}

KKSMultiphaseAction::KKSMultiphaseAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num"))
{
}

void
KKSMultiphaseAction::act()
{
  std::vector<MaterialPropertyName> v_Fjnames;
  v_Fjnames.resize(_op_num);
  std::vector<MaterialPropertyName> v_hjnames;
  v_hjnames.resize(_op_num); 
  std::vector<VariableName> v_cjnames;
  v_cjnames.resize(_op_num); 
  std::vector<VariableName> v_etanames;
  v_etanames.resize(_op_num); 

  unsigned int ind = 0;
  for (unsigned int j = 0; j < _op_num; ++j) {
      v_Fjnames[ind] = "F" + Moose::stringify(j);
      v_hjnames[ind] = "h" + Moose::stringify(j);
      v_cjnames[ind] = "c" + Moose::stringify(j);
      v_etanames[ind] = "eta" + Moose::stringify(j);
      ind ++;
  }

  std::vector<VariableName> var_name_c0;
  var_name_c0.resize(1); 
  var_name_c0[0] = "c";


  for (unsigned int j = 0; j < _op_num; ++j)
  {
    //
    // Create variable names
    //

    std::string var_name_g = "g" + Moose::stringify(j);
    std::string var_name_D = "Dh" + Moose::stringify(j);

    std::vector<VariableName> var_name_c;
    var_name_c.resize(1); 
    var_name_c[0] = "c" + Moose::stringify(j);


    std::vector<VariableName> var_name_eta;
    var_name_eta.resize(1); 
    var_name_eta[0] = "eta" + Moose::stringify(j);

    std::vector<VariableName> var_name_lambda;
    var_name_lambda.resize(1); 
    var_name_lambda[0] = "lambda";

    std::vector<VariableName> var_name_cb;
    var_name_cb.resize(1); 
    var_name_cb[0] = "c" + Moose::stringify(j + 1);


    std::vector<VariableName> v_bulkF;
    v_bulkF.resize(2 * _op_num - 1); //notice this is the _op_num - 1, so that the op_index itself isn't included
    std::vector<VariableName> v_bulkC;
    v_bulkC.resize(_op_num - 1);

    unsigned int ind_F = 0;
    unsigned int ind_C = 0;
    for (unsigned int k = 0; k < _op_num; ++k) {
      v_bulkF[ind_F++] = "c" + Moose::stringify(k);
      if (k != j)
        v_bulkC[ind_C++] = "eta" + Moose::stringify(k);
    }
    for (unsigned int k = 0; k < _op_num; ++k) 
      if (k != j)
	v_bulkF[ind_F++] = "eta" + Moose::stringify(k);


    //
    // Set up Matdiffusion kernels
    //
    if (j == 0) {
      InputParameters params = _factory.getValidParams("TimeDerivative");
      params.set<NonlinearVariableName>("variable") = "c";
      std::string kernel_name = "diff_time";
      _problem->addKernel("TimeDerivative", kernel_name, params);
    }

    {
      InputParameters params = _factory.getValidParams("MatDiffusion");
      params.set<NonlinearVariableName>("variable") = "c";
      params.set<MaterialPropertyName>("D_name") = var_name_D;
      //params.set<std::vector<VariableName>>("v") = v;
      params.set<std::vector<VariableName>>("conc") = var_name_c;
      std::string kernel_name = "diff_c" + Moose::stringify(j);
      _problem->addKernel("MatDiffusion", kernel_name, params);
    }

    //
    // Set up Allen-Cahn equation
    //
    if ( j < _op_num - 1)
    {
      {
	InputParameters params = _factory.getValidParams("TimeDerivative");
	params.set<NonlinearVariableName>("variable") = var_name_eta[0];
	std::string kernel_name = "delta" + Moose::stringify(j) + "dt";
	_problem->addKernel("TimeDerivative", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("KKSMultiACBulkF");
	params.set<NonlinearVariableName>("variable") = var_name_eta[0];
	params.set<std::vector<MaterialPropertyName>>("Fj_names") = v_Fjnames;
	params.set<std::vector<MaterialPropertyName>>("hj_names") = v_hjnames;
	params.set<MaterialPropertyName>("gi_name") = var_name_g;
	params.set<std::vector<VariableName>>("eta_i") = var_name_eta;
	params.set<Real>("wi") = 1.0;
	params.set<std::vector<VariableName>>("args") = v_bulkF;
	std::string kernel_name = "ACBulkF" + Moose::stringify(j);
	_problem->addKernel("KKSMultiACBulkF", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("KKSMultiACBulkC");
	params.set<NonlinearVariableName>("variable") = var_name_eta[0];
	params.set<std::vector<MaterialPropertyName>>("Fj_names") = v_Fjnames;
	params.set<std::vector<MaterialPropertyName>>("hj_names") = v_hjnames;
	params.set<std::vector<VariableName>>("cj_names") = v_cjnames;
	params.set<std::vector<VariableName>>("eta_i") = var_name_eta;
	params.set<std::vector<VariableName>>("args") = v_bulkC;
	std::string kernel_name = "ACBulkC" + Moose::stringify(j);
	_problem->addKernel("KKSMultiACBulkC", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("ACInterface");
	params.set<NonlinearVariableName>("variable") = var_name_eta[0];
	params.set<MaterialPropertyName>("kappa_name") = "kappa";
	std::string kernel_name = "ACInterface" + Moose::stringify(j);
	_problem->addKernel("ACInterface", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("MatReaction");
	params.set<NonlinearVariableName>("variable") = var_name_eta[0];
	params.set<std::vector<VariableName>>("v") = var_name_lambda;
	params.set<MaterialPropertyName>("mob_name") = "L";
	std::string kernel_name = "multipler" + Moose::stringify(j);
	_problem->addKernel("MatReaction", kernel_name, params);
      }
    }

    //
    // Set up Lagrange multiplier equation
    //
    if (j == 0) {
      InputParameters params = _factory.getValidParams("MatReaction");
      params.set<NonlinearVariableName>("variable") = "lambda";
      params.set<MaterialPropertyName>("mob_name") = "3";
      std::string kernel_name = "mult_lamda";
      _problem->addKernel("MatReaction", kernel_name, params);
    }

    {
      {
	InputParameters params = _factory.getValidParams("KKSMultiACBulkF");
	params.set<NonlinearVariableName>("variable") = "lambda";
	params.set<std::vector<MaterialPropertyName>>("Fj_names") = v_Fjnames;
	params.set<std::vector<MaterialPropertyName>>("hj_names") = v_hjnames;
	params.set<MaterialPropertyName>("gi_name") = var_name_g;
	params.set<std::vector<VariableName>>("eta_i") = var_name_eta;
	params.set<Real>("wi") = 1.0;
	params.set<MaterialPropertyName>("mob_name") = "1";
	params.set<std::vector<VariableName>>("args") = v_bulkF;
	std::string kernel_name = "mult_ACBulkF_" + Moose::stringify(j);
	_problem->addKernel("KKSMultiACBulkF", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("KKSMultiACBulkC");
	params.set<NonlinearVariableName>("variable") = "lambda";
	params.set<std::vector<MaterialPropertyName>>("Fj_names") = v_Fjnames;
	params.set<std::vector<MaterialPropertyName>>("hj_names") = v_hjnames;
	params.set<std::vector<VariableName>>("cj_names") = v_cjnames;
	params.set<std::vector<VariableName>>("eta_i") = var_name_eta;
	params.set<std::vector<VariableName>>("args") = v_bulkC;
	params.set<MaterialPropertyName>("mob_name") = "1";
	std::string kernel_name = "mult_ACBulkC" + Moose::stringify(j);
	_problem->addKernel("KKSMultiACBulkC", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("SimpleCoupledACInterface");
	params.set<NonlinearVariableName>("variable") = "lambda";
	params.set<std::vector<VariableName>>("v") = var_name_eta;
	params.set<MaterialPropertyName>("kappa_name") = "kappa";
	params.set<MaterialPropertyName>("mob_name") = "1";
	std::string kernel_name = "multi_CoupledACInt_" + Moose::stringify(j);
	_problem->addKernel("SimpleCoupledACInterface", kernel_name, params);
      }
    }

    //
    // Set up kernels for constraint equation \sum {eta_i} = 1
    //
    if ( j < _op_num - 1) {
      InputParameters params = _factory.getValidParams("MatReaction");
      params.set<NonlinearVariableName>("variable") = "eta" + Moose::stringify(_op_num);
      params.set<std::vector<VariableName>>("v") = var_name_eta;
      params.set<MaterialPropertyName>("mob_name") = "1";
      std::string kernel_name = "eta" + Moose::stringify(j) + "reaction";
      _problem->addKernel("MatReaction", kernel_name, params);
    }
    if ( j == _op_num - 1) {
      {
	InputParameters params = _factory.getValidParams("MatReaction");
	params.set<NonlinearVariableName>("variable") = "eta" + Moose::stringify(_op_num);
	params.set<MaterialPropertyName>("mob_name") = "1";
	std::string kernel_name = "eta" + Moose::stringify(_op_num) + "reaction";
	_problem->addKernel("MatReaction", kernel_name, params);
      }

      {
	InputParameters params = _factory.getValidParams("BodyForce");
	params.set<NonlinearVariableName>("variable") = "eta" + Moose::stringify(_op_num);
	params.set<Real>("value") = -1.0;
	std::string kernel_name = "one";
	_problem->addKernel("BodyForce", kernel_name, params);
      }
    }


    //
    // Set up phase concentration constraints
    //

    if ( j < _op_num - 1) {
      InputParameters params = _factory.getValidParams("KKSPhaseChemicalPotential");
      params.set<NonlinearVariableName>("variable") = "c" + Moose::stringify(j);
      params.set<std::vector<VariableName>>("cb") = var_name_cb;
      params.set<MaterialPropertyName>("fa_name") = "F" + Moose::stringify(j);
      params.set<MaterialPropertyName>("fb_name") = "F" + Moose::stringify(j + 1);
      std::string kernel_name = "chempot" + Moose::stringify(j) + Moose::stringify(j + 1);
      _problem->addKernel("KKSPhaseChemicalPotential", kernel_name, params);
    }
    if ( j == _op_num - 1) {
      InputParameters params = _factory.getValidParams("KKSMultiPhaseConcentration");
      params.set<NonlinearVariableName>("variable") = "c" + Moose::stringify(j);
      params.set<std::vector<VariableName>>("cj") = v_cjnames;
      params.set<std::vector<MaterialPropertyName>>("hj_names") = v_hjnames;
      params.set<std::vector<VariableName>>("etas") = v_etanames;
      params.set<std::vector<VariableName>>("c") = var_name_c0;
      std::string kernel_name = "phaseconcentration";
      _problem->addKernel("KKSMultiPhaseConcentration", kernel_name, params);
    }

  }
}
