/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "KKSMultiCircleICAction.h"
#include "Factory.h"
#include "MooseMesh.h"
#include "FEProblem.h"
#include "Conversion.h"
#include "MooseVariable.h"
#include "MooseRandom.h"

template <>
InputParameters
validParams<KKSMultiCircleICAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Random multicircle ics for KKS multiphase action");
  params.addRequiredParam<unsigned int>("op_num", "number of order parameters to create");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addRequiredParam<Real>("circlespac",
                                "minimum spacing of circles, measured from center to center");
  params.addRequiredParam<Real>("int_width",
                                "interfacial width of the created width");
  params.addRequiredParam<std::vector<Real>>("twophase_eq_concentration",
					     "equilibrium concentration of two phases");
  params.addRequiredParam<Real>("xmin", "mininum value in x axis");
  params.addRequiredParam<Real>("xmax", "maxinum value in x axis");
  params.addRequiredParam<Real>("ymin", "mininum value in y axis");
  params.addRequiredParam<Real>("ymax", "maxinum value in y axis");
  params.addRequiredParam<Real>("zmin", "mininum value in z axis");
  params.addRequiredParam<Real>("zmax", "maxinum value in z axis");
  params.addParam<unsigned int>("numtries", 1000, "The number of tries");
  params.addParam<unsigned int>("rand_seed", 12444, "The random seed");
  params.addRequiredParam<Real>("avg_radius", "Mean radius value for the circles");
  params.addParam<Real>("radius_variation",
                        0.0,
                        "Plus or minus fraction of random variation in "
                        "the grain radius for uniform, standard "
                        "deviation for normal");
  MooseEnum rand_options("uniform none", "none");
  params.addParam<MooseEnum>("radius_variation_type",
                             rand_options,
                             "Type of distribution that random circle radii will follow");
  params.addParam<bool>(
      "columnar_3D", false, "3D microstructure will be columnar in the z-direction?");
  return params;
}

KKSMultiCircleICAction::KKSMultiCircleICAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _circlespac(getParam<Real>("circlespac")),
    _int_width(getParam<Real>("int_width")),
    _eq_concentration(getParam<std::vector<Real>>("twophase_eq_concentration")),
    _xmin(getParam<Real>("xmin")),
    _xmax(getParam<Real>("xmax")),
    _ymin(getParam<Real>("ymin")),
    _ymax(getParam<Real>("ymax")),
    _zmin(getParam<Real>("zmin")),
    _zmax(getParam<Real>("zmax")),
    _max_num_tries(getParam<unsigned int>("numtries")),
    _rand_seed(getParam<unsigned int>("rand_seed")),
    _avg_radius(getParam<Real>("avg_radius")),
    _radius_variation(getParam<Real>("radius_variation"))
    //    _radius_variation_type(getParam<MooseEnum>("radius_variation_type")),
    //_columnar_3D(getParam<MooseEnum>("columnar_3D"))
{
}

void
KKSMultiCircleICAction::act()
{
#ifdef DEBUG
  Moose::err << "Inside the KKSMultiCircleICAction Object\n";
#endif

  //obtain the dimension of the simulation box
  _bottom_left(0) = _xmin;
  _bottom_left(1) = _ymin;
  _bottom_left(2) = _zmin;
  _range(0) = _xmax - _xmin;
  _range(1) = _ymax - _ymin;
  _range(2) = _zmax - _zmin;

  // Randomly generate the radius of the individual grains
  computeCircleRadii();
  // Randomly generate the centers of the individual grains
  computeCircleCenters();

  std::vector<Real> _x_pos, _y_pos, _z_pos, _radius;
  _x_pos.resize(_op_num - 1);
  _y_pos.resize(_op_num - 1);
  _z_pos.resize(_op_num - 1);
  _radius.resize(_op_num - 1);

  // Loop through the number of order parameters
  for (unsigned int op = 1; op < _op_num; op++) {
    _x_pos[op - 1] = _centers[op](0);
    _y_pos[op - 1] = _centers[op](1);
    _z_pos[op - 1] = _centers[op](2);
    _radius[op - 1] = _radii[op];
  }

  // Loop through the number of order parameters
  for (unsigned int op = 0; op < _op_num; op++) {
    if (op == 0) {
      {
	InputParameters poly_params = _factory.getValidParams("SpecifiedSmoothCircleIC");
	poly_params.set<VariableName>("variable") = _var_name_base + Moose::stringify(op + 1);
	poly_params.set<std::vector<Real>>("x_positions") = _x_pos;
	poly_params.set<std::vector<Real>>("y_positions") = _y_pos;
	poly_params.set<std::vector<Real>>("z_positions") = _z_pos;
	poly_params.set<std::vector<Real>>("radii") = _radius;
	poly_params.set<Real>("invalue") = 0.0;
	poly_params.set<Real>("outvalue") = 1.0;
	poly_params.set<Real>("int_width") = _int_width;
	// Add initial condition
	_problem->addInitialCondition("SpecifiedSmoothCircleIC", _var_name_base + Moose::stringify(op + 1), poly_params);
      }

      {
	InputParameters poly_params = _factory.getValidParams("SpecifiedSmoothCircleIC");
	poly_params.set<VariableName>("variable") = "c";
	poly_params.set<std::vector<Real>>("x_positions") = _x_pos;
	poly_params.set<std::vector<Real>>("y_positions") = _y_pos;
	poly_params.set<std::vector<Real>>("z_positions") = _z_pos;
	poly_params.set<std::vector<Real>>("radii") = _radius;
	poly_params.set<Real>("invalue") = _eq_concentration[1];
	poly_params.set<Real>("outvalue") = _eq_concentration[0];
	poly_params.set<Real>("int_width") = _int_width;
	// Add initial condition
	_problem->addInitialCondition("SpecifiedSmoothCircleIC", "c", poly_params);
      }
    }
    else {
      InputParameters poly_params = _factory.getValidParams("SmoothCircleIC");
      poly_params.set<VariableName>("variable") = _var_name_base + Moose::stringify(op + 1);
      poly_params.set<Real>("x1") = _centers[op](0);
      poly_params.set<Real>("y1") = _centers[op](1);
      poly_params.set<Real>("z1") = _centers[op](2);
      poly_params.set<Real>("radius") = _radii[op];
      poly_params.set<Real>("invalue") = 1.0;
      poly_params.set<Real>("outvalue") = 0.0;
      poly_params.set<Real>("int_width") = _int_width;
      // Add initial condition
      _problem->addInitialCondition("SmoothCircleIC", _var_name_base + Moose::stringify(op + 1), poly_params);
    }  //end else

  } //end for


}

void
KKSMultiCircleICAction::computeCircleRadii()
{
  MooseRandom::seed(_rand_seed);

  _radii.resize(_op_num);

  for (unsigned int i = 0; i < _op_num; i++) {
      // Vary bubble radius
    _radii[i] = _avg_radius * (1.0 + (1.0 - 2.0 * MooseRandom::rand()) * _radius_variation);
    _radii[i] = std::max(_radii[i], 0.0);
  }
  
}

void
KKSMultiCircleICAction::computeCircleCenters()
{

  MooseRandom::seed(_rand_seed);

  _centers.resize(_op_num);

  for (unsigned int phase = 0; phase < _op_num; phase++)
  {
    // Vary circle center positions
    unsigned int num_tries = 0;
    while (num_tries < _max_num_tries)
    {
      num_tries++;

      for (unsigned int i = 0; i < LIBMESH_DIM; i++)
	_centers[phase](i) = _bottom_left(i) + _range(i) * MooseRandom::rand();
      if (_columnar_3D)
	_centers[phase](2) = _bottom_left(2) + _range(2) * 0.5;

      for (unsigned int j = 0; j < phase; j++) {
	double _dx = _centers[j](0) - _centers[phase](0);
	double _dy = _centers[j](1) - _centers[phase](1);
	double _dz = _centers[j](2) - _centers[phase](2);
	double _dist = std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
        if (_dist < _circlespac)
          goto fail;
      }
      // accept the position of the new center
      goto accept;

    // retry a new position until tries are exhausted
    fail:
      continue;
    }

    if (num_tries == _max_num_tries)
      mooseError("Too many tries in KKSMultiCircleIC");

  accept:
    continue;
  }
}

