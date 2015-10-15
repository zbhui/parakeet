#ifndef PARAKEETAPP_H
#define PARAKEETAPP_H

#include "MooseApp.h"

class ParakeetApp;

template<>
InputParameters validParams<ParakeetApp>();

class ParakeetApp : public MooseApp
{
public:
  ParakeetApp(const InputParameters &parameters);

  virtual ~ParakeetApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);

private:
  void registerMooseObjects(Factory & factory);
};

#endif /* PARAKEETAPP_H */
