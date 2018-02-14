/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PolycrystalMultiCircleIC.h"
#include "IndirectSort.h"
#include "MooseRandom.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<PolycrystalMultiCircleIC>()
{
  InputParameters params = validParams<InitialCondition>();

  params.addClassDescription(
      "Random Voronoi tesselation polycrystal (used by PolycrystalVoronoiICAction)");
  params.addRequiredParam<unsigned int>("op_num", "Number of order parameters");
  params.addRequiredParam<unsigned int>(
      "grain_num", "Number of grains being represented by the order parameters");
  params.addRequiredParam<unsigned int>("op_index", "The index for the current order parameter");

  params.addRequiredParam<Real>("circlespac",
                                "minimum spacing of circles, measured from center to center");
  params.addParam<unsigned int>("numtries", 1000, "The number of tries");
  params.addParam<unsigned int>("rand_seed", 12444, "The random seed");
  params.addRequiredParam<Real>("radius", "Mean radius value for the circles");
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

PolycrystalMultiCircleIC::PolycrystalMultiCircleIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _mesh(_fe_problem.mesh()),
    _dim(_mesh.dimension()),
    _op_num(getParam<unsigned int>("op_num")),
    _grain_num(getParam<unsigned int>("grain_num")),
    _op_index(getParam<unsigned int>("op_index")),
    _rand_seed(getParam<unsigned int>("rand_seed")),
    _columnar_3D(getParam<bool>("columnar_3D")),
    _max_num_tries(getParam<unsigned int>("numtries")),
    _circlespac(getParam<Real>("circlespac")),
    _radius(getParam<Real>("radius")),
    _radius_variation(getParam<Real>("radius_variation")),
    _radius_variation_type(getParam<MooseEnum>("radius_variation_type"))
{
}

void
PolycrystalMultiCircleIC::initialSetup()
{
  // Set up domain bounds with mesh tools
  for (unsigned int i = 0; i < LIBMESH_DIM; i++)
  {
    _bottom_left(i) = _mesh.getMinInDimension(i);
    _top_right(i) = _mesh.getMaxInDimension(i);
  }
  _range = _top_right - _bottom_left;

  if (_radius_variation > 0.0 && _radius_variation_type == 1)
    mooseError("If radius_variation > 0.0, you must pass in a radius_variation_type in "
               "MultiSmoothCircleIC");

  if (_op_num > _grain_num)
    mooseError("ERROR in PolycrystalMultiCircleIC: Number of order parameters (op_num) can't be larger "
               "than the number of grains (grain_num)");

  // Randomly generate the radius of the individual grains
  computeCircleRadii();

  // Randomly generate the centers of the individual grains
  computeCircleCenters();

  _assigned_op = PolycrystalICTools::assignPointsToVariables_solution(_centers, _op_num, _mesh, _var);


}


void
PolycrystalMultiCircleIC::computeCircleRadii()
{
  MooseRandom::seed(_rand_seed);

  _radii.resize(_grain_num);

  for (unsigned int i = 0; i < _grain_num; i++)
  {
    // Vary bubble radius
    switch (_radius_variation_type)
    {
      case 0: // Uniform distribution
        _radii[i] = _radius * (1.0 + (1.0 - 2.0 * MooseRandom::rand()) * _radius_variation);
        break;
      case 1: // no variation
        _radii[i] = _radius;
    }

    _radii[i] = std::max(_radii[i], 0.0);
  }
}

void
PolycrystalMultiCircleIC::computeCircleCenters()
{

  MooseRandom::seed(_rand_seed);

  _centers.resize(_grain_num);
  for (unsigned int grain = 0; grain < _grain_num; ++grain)
  {
    // Vary circle center positions
    unsigned int num_tries = 0;
    while (num_tries < _max_num_tries)
    {
      num_tries++;

      RealTensorValue ran;
      for (unsigned int i = 0; i < LIBMESH_DIM; i++)
	_centers[grain](i) = _bottom_left(i) + _range(i) * MooseRandom::rand();
      if (_columnar_3D)
	_centers[grain](2) = _bottom_left(2) + _range(2) * 0.5;

      for (unsigned int j = 0; j < grain; ++j)
        if (_mesh.minPeriodicDistance(_var.number(), _centers[j], _centers[grain]) < _circlespac)
          goto fail;

      // accept the position of the new center
      goto accept;

    // retry a new position until tries are exhausted
    fail:
      continue;
    }

    if (num_tries == _max_num_tries)
      mooseError("Too many tries in PolycrystalMultiCircleIC");

  accept:
    continue;
  }
}


Real
PolycrystalMultiCircleIC::value(const Point & p)
{
  unsigned int min_index =
    PolycrystalICTools::assignPointToGrain_multicircle(p, _centers, _radii, _mesh, _var, _range.norm());

  // If the current order parameter index (_op_index) is equal to the min_index, set the value to
  // 1.0
  if (_assigned_op[min_index] == _op_index) // Make sure that the _op_index goes from 0 to _op_num-1
    return 0.99;
  else
    return 0.01 / (_op_num - 1);
}
