#include "ParakeetApp.h"

#include "ProblemElementalL2Error.h"
#include "Moose.h"
#include "AppFactory.h"

#include "ActionFactory.h"
#include "Syntax.h"

/// Action
#include "AddMultiDGKernel.h"
#include "AddMultiKernel.h"
#include "AddMultiBC.h"
#include "AddMultiIC.h"
#include "AddMultiVariable.h"

#include "CLawAuxVariablesAction.h"
#include "CommonPostProcessorAction.h"
#include "AddMultiAuxVariableAction.h"

/// 单元积分
#include "CFDCellKernel.h"
#include "EmptyTimeDerivative.h"

/// 面积分
#include "CFDFaceKernel.h"

/// 初始条件
#include "CFDInitialCondition.h"

/// 边界条件
#include "CFDBC.h"
#include "SlipWall.h"
#include "AdiabaticWall.h"
#include "IsoThermalWall.h"
#include "Symmetric.h"
#include "FarFieldPressure.h"
#include "FarFieldRiemann.h"
#include "IsoVortexExact.h"
#include "CouetteFlowExact.h"

/// 辅助kernel
#include "NSAuxVariable.h"
#include "NearestWallDistance.h"

/// 时间积分
#include "MultiTimeDerivative.h"

/// 时间步长增加策略
#include "RatioTimeStepper.h"

#include "EulerProblem.h"
#include "SodProblem.h"
#include "Riemann1DProblem.h"
#include "ShockVortexProblem.h"
#include "DoubleMachProblem.h"

#include "CFDProblem.h"
#include "NavierStokesProblem.h"
#include "IsoVortexProblem.h"
#include "CouetteFlowProblem.h"

/// PostProcessor
#include "CFDResidual.h"
#include "ElementExtremeTimeDerivative.h"
#include "NumTimeStep.h"
#include "VariableResidual.h"
#include "CouetteFlowElementL2Error.h"

#include "VariableJumpIndicator.h"
#include "FluxJumpIndicator.h"

template<>
InputParameters validParams<ParakeetApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

ParakeetApp::ParakeetApp(const InputParameters &parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  //ModulesApp::registerObjects(_factory);
  ParakeetApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  //ModulesApp::associateSyntax(_syntax, _action_factory);
  ParakeetApp::associateSyntax(_syntax, _action_factory);
}

ParakeetApp::~ParakeetApp()
{
}

// External entry point for dynamic application loading
extern "C" void ParakeetApp__registerApps() { ParakeetApp::registerApps(); }
void
ParakeetApp::registerApps()
{
#undef  registerApp
#define registerApp(name) AppFactory::instance().reg<name>(#name)

  registerApp(ParakeetApp);


#undef  registerApp
#define registerApp(name) AppFactory::instance().regLegacy<name>(#name)
}

// External entry point for dynamic object registration
extern "C" void ParakeetApp__registerObjects(Factory & factory) { ParakeetApp::registerObjects(factory); }
void
ParakeetApp::registerObjects(Factory & factory)
{
#undef registerObject
#define registerObject(name) factory.reg<name>(stringifyName(name))
	/// 注册单元积分
		registerKernel(CFDCellKernel);
		registerKernel(EmptyTimeDerivative);

		/// 注册面积分
		registerDGKernel(CFDFaceKernel);

		/// 注册初始条件
		registerInitialCondition(CFDInitialCondition);

		/// 注册边界条件
		registerBoundaryCondition(CFDBC);
		registerBoundaryCondition(AdiabaticWall);
		registerBoundaryCondition(SlipWall);
		registerBoundaryCondition(FarFieldPressure);
		registerBoundaryCondition(FarFieldRiemann);
		registerBoundaryCondition(IsoThermalWall);
		registerBoundaryCondition(Symmetric);

//		registerBoundaryCondition(CouetteFlowBC);
//		registerBoundaryCondition(NavierStokesBC);
//		registerBoundaryCondition(KOmegaBC);

		/// 注册函数
		registerFunction(IsoVortexExact);
//		registerFunction(CouetteFlowExact);

		/// 注册辅助kernel
		registerAux(NSAuxVariable);
		registerAux(NearestWallDistance);

		registerExecutioner(RatioTimeStepper);

		registerProblem(CFDProblem);
		registerProblem(EulerProblem);
		registerProblem(SodProblem);
		registerProblem(Riemann1DProblem);
		registerProblem(ShockVortexProblem);
		registerProblem(DoubleMachProblem);
		registerProblem(IsoVortexProblem);
		registerProblem(NavierStokesProblem);
		registerProblem(CouetteFlowProblem);

		registerPostprocessor(CFDResidual);
		registerPostprocessor(ElementExtremeTimeDerivative);
		registerPostprocessor(NumTimeStep);
		registerPostprocessor(VariableResidual);
		registerPostprocessor(ProblemElementalL2Error);
		registerPostprocessor(CouetteFlowElementL2Error);

		registerIndicator(VariableJumpIndicator);
		registerIndicator(FluxJumpIndicator);

#undef registerObject
#define registerObject(name) factory.regLegacy<name>(stringifyName(name))
}

// External entry point for dynamic syntax association
extern "C" void ParakeetApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ParakeetApp::associateSyntax(syntax, action_factory); }
void
ParakeetApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
#undef registerAction
#define registerAction(tplt, action) action_factory.reg<tplt>(stringifyName(tplt), action)

//	syntax.registerActionSyntax("CLawAuxVariablesAction", "AuxVariables");
//	registerAction(CLawAuxVariablesAction, "add_aux_variable");
//	registerAction(CLawAuxVariablesAction, "add_aux_kernel");

	syntax.registerActionSyntax("CommonPostProcessorAction", "Postprocessors", "add_postprocessor");
	registerAction(CommonPostProcessorAction, "add_postprocessor");

	syntax.registerActionSyntax("AddMultiVariable", "Problem/Variables");
	registerAction(AddMultiVariable, "add_variable");

	syntax.registerActionSyntax("AddMultiKernel", "Problem/Kernels");
	registerAction(AddMultiKernel, "add_kernel");

	syntax.registerActionSyntax("AddMultiDGKernel", "Problem/DGKernels");
	registerAction(AddMultiDGKernel, "add_dg_kernel");

	syntax.registerActionSyntax("AddMultiBC", "Problem/BCs/*");
	registerAction(AddMultiBC, "add_bc");

	syntax.registerActionSyntax("AddMultiIC", "ICs");
	registerAction(AddMultiIC, "add_ic");

	syntax.registerActionSyntax("AddMultiAuxVariableAction", "Problem/AuxVariables/*");
	registerAction(AddMultiAuxVariableAction, "add_aux_variable");
	registerAction(AddMultiAuxVariableAction, "add_aux_kernel");

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)
}
