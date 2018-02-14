/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef KKSMULTICIRCLEICACTION_H
#define KKSMULTICIRCLEICACTION_H

#include "InputParameters.h"
#include "Action.h"

/**
 * Random multicircle phase action
 */
class KKSMultiCircleICAction : public Action
{
public:
  KKSMultiCircleICAction(const InputParameters & params);

  virtual void act();

private:
  MooseMesh & _mesh;
  const unsigned int _op_num;
  const std::string _var_name_base;
  const Real _circlespac;
  const Real _int_width;
  const std::vector<Real> _eq_concentration;

  const unsigned int _max_num_tries;
  const unsigned int _rand_seed;

  const Real _avg_radius;
  const Real _radius_variation;
  //const MooseEnum _radius_variation_type;

  bool _columnar_3D;

  Point _bottom_left;
  Point _top_right;
  Point _range;

  std::vector<Point> _centers;
  std::vector<Real> _radii;


};

template <>
InputParameters validParams<KKSMultiCircleICAction>();

#endif // KKSMULTICIRCLEICACTION_H
