
#include "FullJacobianPreconditioner.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"

template<>
InputParameters validParams<FullJacobianPreconditioner>()
{
  InputParameters params = validParams<MoosePreconditioner>();

  return params;
}


FullJacobianPreconditioner::FullJacobianPreconditioner(const InputParameters & params) :
		MoosePreconditioner(params)
{
	 NonlinearSystem & nl = _fe_problem.getNonlinearSystem();
	 unsigned int n_vars = nl.nVariables();

	  CouplingMatrix * cm = new CouplingMatrix(n_vars);


	   for (unsigned int i = 0; i < n_vars; i++)
		   for (unsigned int j = 0; j < n_vars; j++)
			   (*cm)(i,j) = 1;


	   std::cout << "33" <<std::endl;
	  _fe_problem.setCouplingMatrix(cm);
//	  _fe_problem.setCoupling(Moose::COUPLING_DIAG);
}
