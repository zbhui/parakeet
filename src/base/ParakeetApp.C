#include "ParakeetApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<ParakeetApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

ParakeetApp::ParakeetApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ParakeetApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
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
  registerApp(ParakeetApp);
}

// External entry point for dynamic object registration
extern "C" void ParakeetApp__registerObjects(Factory & factory) { ParakeetApp::registerObjects(factory); }
void
ParakeetApp::registerObjects(Factory & factory)
{
}

// External entry point for dynamic syntax association
extern "C" void ParakeetApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ParakeetApp::associateSyntax(syntax, action_factory); }
void
ParakeetApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
