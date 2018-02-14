/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POLYCRYSTAL_SOLUTIONKERNELACTION_H
#define POLYCRYSTAL_SOLUTIONKERNELACTION_H

#include "Action.h"

/**
 * Action that sets up ACGrGrPoly, ACInterface, TimeDerivative, and ACGBPoly
 * kernels.
 */
class Polycrystal_solutionKernelAction : public Action
{
public:
  Polycrystal_solutionKernelAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _var_name_base;

  /// kernels are implicit?
  const bool _implicit;
};

template <>
InputParameters validParams<Polycrystal_solutionKernelAction>();

#endif // POLYCRYSTAL_SOLUTIONKERNELACTION_H
