/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POLYCRYSTALMULTICIRCLEIC_H
#define POLYCRYSTALMULTICIRCLEIC_H

#include "InitialCondition.h"
#include "PolycrystalICTools.h"

// Forward Declarations
class PolycrystalMultiCircleIC;

template <>
InputParameters validParams<PolycrystalMultiCircleIC>();

/**
 * PolycrystalMultiCircleIC creates a polycrystal initial condition.
 */
class PolycrystalMultiCircleIC : public InitialCondition
{
public:
  PolycrystalMultiCircleIC(const InputParameters & parameters);

  virtual Real value(const Point & p);
  virtual void initialSetup();

protected:
  MooseMesh & _mesh;

  /// mesh dimension
  unsigned int _dim;

  unsigned int _op_num;
  unsigned int _grain_num;
  unsigned int _op_index;
  unsigned int _rand_seed;

  bool _cody_test;
  bool _columnar_3D;

  const unsigned int _max_num_tries;
  const Real _circlespac;

  Point _bottom_left;
  Point _top_right;
  Point _range;

  const Real _radius;
  const Real _radius_variation;
  const MooseEnum _radius_variation_type;

  std::vector<Point> _centers;
  std::vector<Real> _radii;

  //MooseRandom _random;

  std::vector<unsigned int> _assigned_op;

  virtual void computeCircleRadii();
  virtual void computeCircleCenters();

};

#endif // POLYCRYSTALMULTICIRCLEIC_H
