#include "ParakeetApp.h"
#include "Moose.h"
#include "AppFactory.h"
//#include "ModulesApp.h"

#include "ActionFactory.h" 	 // <- Actions are special (they have their own factory)
#include "Syntax.h"

/// 单元积分
#include "EulerCellKernel.h"
#include "EmptyTimeDerivative.h"
#include "NavierStokesCellKernel.h"
#include "KOmegaCellKernel.h"
#include "BurgersCellKernel.h"

/// 面积分
#include "EulerFaceKernel.h"
#include "NavierStokesFaceKernel.h"
#include "KOmegaFaceKernel.h"
#include "BurgersFaceKernel.h"

/// 初始条件
#include "IsoVortexIC.h"
#include "CouetteFlowIC.h"
#include "SodIC.h"
#include "CFDPassFlowIC.h"
#include "SinIC.h"
#include "KOmegaIC.h"

/// 边界条件
#include "ConservationLawBC.h"
#include "CouetteFlowBC.h"
#include "IsoVortexBC.h"
#include "CFDBC.h"
#include "KOmegaBC.h"

/// 函数
#include "IsoVortexExact.h"
#include "CouetteFlowExact.h"

/// 辅助kernel
#include "NSAuxVariable.h"
#include "NearestWallDistance.h"

/// Action
#include "CFDAction.h"

/// 材料属性
#include "FluidMaterial.h"

/// 时间积分
#include "MultiTimeDerivative.h"

/// 时间步长增加策略
#include "RatioTimeStepper.h"

template<>
InputParameters validParams<ParakeetApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
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
		registerKernel(EulerCellKernel);
		registerKernel(EmptyTimeDerivative);
//		registerKernel(NavierStokesCellKernel);
//		registerKernel(KOmegaCellKernel);
//		registerKernel(BurgersCellKernel);
	//	registerKernel(MultiTimeDerivative);

		/// 注册面积分
		registerDGKernel(EulerFaceKernel);
//		registerDGKernel(NavierStokesFaceKernel);
//		registerDGKernel(KOmegaFaceKernel);
//		registerDGKernel(BurgersFaceKernel);

		/// 注册初始条件
		registerInitialCondition(IsoVortexIC);
		registerInitialCondition(CouetteFlowIC);
		registerInitialCondition(SodIC);
		registerInitialCondition(CFDPassFlowIC);
		registerInitialCondition(SinIC);
		registerInitialCondition(KOmegaIC);

		/// 注册边界条件
		registerBoundaryCondition(BurgersBC);
		registerBoundaryCondition(EulerBC);
		registerBoundaryCondition(IsoVortexBC);
//		registerBoundaryCondition(CouetteFlowBC);
//		registerBoundaryCondition(NavierStokesBC);
//		registerBoundaryCondition(KOmegaBC);

		/// 注册函数
		registerFunction(IsoVortexExact);
		registerFunction(CouetteFlowExact);

		/// 注册辅助kernel
		registerAux(NSAuxVariable);
		registerAux(NearestWallDistance);

		/// 注册材料属性
		registerMaterial(FluidMaterial);

		registerExecutioner(RatioTimeStepper);

//		registerNamedPreconditioner(FullJacobianPreconditioner, "FJP");


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

	registerAction(CFDAction, "add_kernel");

	syntax.registerActionSyntax("CFDAction", "CFDAction");

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)
}
