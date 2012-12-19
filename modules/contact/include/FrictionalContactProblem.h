/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef FRICTIONALCONTACTPROBLEM_H
#define FRICTIONALCONTACTPROBLEM_H

#include "ReferenceResidualProblem.h"

class FrictionalContactProblem;

struct InteractionParams
{
  Real _friction_coefficient;
  Real _slip_factor;
  Real _slip_too_far_factor;
};

enum ContactState
{
  STICKING,
  SLIPPING,
  SLIPPED_TOO_FAR
};

template<>
InputParameters validParams<FrictionalContactProblem>();

/**
 * FEProblem derived class for frictional contact-specific callbacks
 */
class FrictionalContactProblem : public ReferenceResidualProblem
{
public:
  FrictionalContactProblem(const std::string & name, InputParameters params);
  virtual ~FrictionalContactProblem();

  class SlipData
  {
    public:
      SlipData(const Node* node, unsigned int dof, Real slip);
      SlipData(const SlipData & sd);
      ~SlipData();

      const Node* _node;
      unsigned int _dof;
      Real _slip;
  };

  virtual void initialSetup();
  virtual void timestepSetup();
  virtual bool shouldUpdateSolution();
  virtual bool updateSolution(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution);
  virtual void predictorCleanup(NumericVector<Number>& ghosted_solution);
  bool enforceRateConstraint(NumericVector<Number>& vec_solution, NumericVector<Number>& ghosted_solution);
  bool calculateSlip(const NumericVector<Number>& ghosted_solution,
                     std::vector<SlipData> * iterative_slip);
  static ContactState calculateInteractionSlip(RealVectorValue &slip,
                                    Real &slip_residual,
                                    const RealVectorValue &normal,
                                    const RealVectorValue &residual,
                                    const RealVectorValue &incremental_slip,
                                    const RealVectorValue &stiffness,
                                    const Real friction_coefficient,
                                    const Real slip_factor,
                                    const Real slip_too_far_factor,
                                    const int dim);
  void applySlip(NumericVector<Number>& vec_solution,
                 NumericVector<Number>& ghosted_solution,
                 std::vector<SlipData> & iterative_slip);
  unsigned int numLocalFrictionalConstraints();
  void updateContactReferenceResidual();
  virtual MooseNonlinearConvergenceReason checkNonlinearConvergence(std::string &msg,
                                                                    const int it,
                                                                    const Real xnorm,
                                                                    const Real snorm,
                                                                    const Real fnorm,
                                                                    Real &ttol,
                                                                    const Real rtol,
                                                                    const Real stol,
                                                                    const Real abstol,
                                                                    const int nfuncs,
                                                                    const int max_funcs,
                                                                    const Real ref_resid,
                                                                    const Real div_threshold);
  void updateContactPoints(NumericVector<Number>& ghosted_solution,
                           bool update_incremental_slip);
  void updateIncrementalSlip();

protected:
  std::map<std::pair<int,int>,InteractionParams> _interaction_params;
  NonlinearVariableName _disp_x;
  NonlinearVariableName _disp_y;
  NonlinearVariableName _disp_z;
  AuxVariableName _residual_x;
  AuxVariableName _residual_y;
  AuxVariableName _residual_z;
  AuxVariableName _diag_stiff_x;
  AuxVariableName _diag_stiff_y;
  AuxVariableName _diag_stiff_z;
  AuxVariableName _inc_slip_x;
  AuxVariableName _inc_slip_y;
  AuxVariableName _inc_slip_z;

  std::vector<std::string> _contactRefResidVarNames;
  std::vector<unsigned int> _contactRefResidVarIndices;
  Real _refResidContact;

  Real _slip_residual;
  bool _do_slip_update;
  int _num_slip_iterations;
  int _min_slip_iters;
  int _max_slip_iters;
  int _slip_updates_per_iter;
  Real _target_contact_residual;
  Real _target_relative_contact_residual;
  Real _contact_slip_tol_factor;
  int _num_nl_its_since_contact_update;
  int _num_contact_nodes;
  int _num_slipping;
  int _num_slipped_too_far;
  Real _inc_slip_norm;
  Real _it_slip_norm;
};

#endif /* FRICTIONALCONTACTPROBLEM_H */
