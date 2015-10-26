
// moose object

#include "libmesh/petsc_macro.h"

#include "Moose.h"
#include "Factory.h"
#include "NonlinearSystem.h"
#include "PetscSupport.h"
#include "ActionWarehouse.h"
#include "ActionFactory.h"
#include "Syntax.h"

// objects that can be created by MOOSE
// Mesh
#include "FileMesh.h"
#include "GeneratedMesh.h"
#include "TiledMesh.h"
// MeshModifiers
#include "MeshExtruder.h"
#include "SideSetsFromPoints.h"
#include "SideSetsFromNormals.h"
#include "AddExtraNodeset.h"
#include "BoundingBoxNodeSet.h"
#include "Transform.h"
#include "SideSetsAroundSubdomain.h"
#include "SideSetsBetweenSubdomains.h"
#include "AddAllSideSetsByNormals.h"
#include "SubdomainBoundingBox.h"
#include "OrientedSubdomainBoundingBox.h"
#include "RenameBlock.h"

// problems
#include "FEProblem.h"
#include "DisplacedProblem.h"

// kernels
#include "TimeDerivative.h"
#include "CoupledTimeDerivative.h"
#include "MassLumpedTimeDerivative.h"
#include "Diffusion.h"
#include "AnisotropicDiffusion.h"
#include "CoupledForce.h"
#include "UserForcingFunction.h"
#include "BodyForce.h"
#include "Reaction.h"
#include "MassEigenKernel.h"

// bcs
#include "ConvectiveFluxBC.h"
#include "DirichletBC.h"
#include "PenaltyDirichletBC.h"
#include "PresetBC.h"
#include "NeumannBC.h"
#include "FunctionDirichletBC.h"
#include "FunctionPenaltyDirichletBC.h"
#include "FunctionPresetBC.h"
#include "FunctionNeumannBC.h"
#include "MatchedValueBC.h"
#include "VacuumBC.h"
#include "SinDirichletBC.h"
#include "SinNeumannBC.h"
#include "VectorNeumannBC.h"
#include "WeakGradientBC.h"
#include "DiffusionFluxBC.h"
#include "PostprocessorDirichletBC.h"
#include "OneDEqualValueConstraintBC.h"

// auxkernels
#include "ConstantAux.h"
#include "FunctionAux.h"
#include "NearestNodeDistanceAux.h"
#include "NearestNodeValueAux.h"
#include "PenetrationAux.h"
#include "ProcessorIDAux.h"
#include "SelfAux.h"
#include "GapValueAux.h"
#include "MaterialRealAux.h"
#include "MaterialRealVectorValueAux.h"
#include "MaterialRealTensorValueAux.h"
#include "MaterialStdVectorAux.h"
#include "MaterialRealDenseMatrixAux.h"
#include "MaterialStdVectorRealGradientAux.h"
#include "DebugResidualAux.h"
#include "BoundsAux.h"
#include "SpatialUserObjectAux.h"
#include "SolutionAux.h"
#include "VectorMagnitudeAux.h"
#include "ConstantScalarAux.h"
#include "QuotientAux.h"
#include "NormalizationAux.h"
#include "VariableGradientComponent.h"
#include "ParsedAux.h"

// dirac kernels
#include "ConstantPointSource.h"

// DG kernels
#include "DGDiffusion.h"

// ics
#include "ConstantIC.h"
#include "BoundingBoxIC.h"
#include "FunctionIC.h"
#include "RandomIC.h"
#include "ScalarConstantIC.h"
#include "ScalarComponentIC.h"

// executioners
#include "Steady.h"
#include "Transient.h"
#include "InversePowerMethod.h"
#include "NonlinearEigen.h"
#include "PetscTSExecutioner.h"

// functions
#include "Axisymmetric2D3DSolutionFunction.h"
#include "ConstantFunction.h"
#include "CompositeFunction.h"
#include "MooseParsedFunction.h"
#include "MooseParsedVectorFunction.h"
#include "MooseParsedGradFunction.h"
#include "PiecewiseConstant.h"
#include "PiecewiseLinear.h"
#include "SolutionFunction.h"
#include "PiecewiseBilinear.h"
#include "SplineFunction.h"
#include "PiecewiseMultilinear.h"
#include "LinearCombinationFunction.h"

// materials
#include "GenericConstantMaterial.h"
#include "GenericFunctionMaterial.h"

// PPS
#include "AverageElementSize.h"
#include "AverageNodalVariableValue.h"
#include "NodalSum.h"
#include "ElementAverageValue.h"
#include "ElementAverageTimeDerivative.h"
#include "ElementW1pError.h"
#include "ElementH1Error.h"
#include "ElementH1SemiError.h"
#include "ElementIntegralVariablePostprocessor.h"
#include "ElementIntegralMaterialProperty.h"
#include "ElementL2Error.h"
#include "ElementVectorL2Error.h"
#include "EmptyPostprocessor.h"
#include "FunctionValuePostprocessor.h"
#include "NodalVariableValue.h"
#include "NumDOFs.h"
#include "TimestepSize.h"
#include "RunTime.h"
#include "PerformanceData.h"
#include "NumElems.h"
#include "NumNodes.h"
#include "NumNonlinearIterations.h"
#include "NumLinearIterations.h"
#include "Residual.h"
#include "ScalarVariable.h"
#include "NumVars.h"
#include "NumResidualEvaluations.h"
#include "Receiver.h"
#include "SideAverageValue.h"
#include "SideFluxIntegral.h"
#include "SideFluxAverage.h"
#include "SideIntegralVariablePostprocessor.h"
#include "NodalMaxValue.h"
#include "NodalProxyMaxValue.h"
#include "PlotFunction.h"
#include "ScalarL2Error.h"
#include "ElementalVariableValue.h"
#include "ElementL2Norm.h"
#include "NodalL2Norm.h"
#include "NodalL2Error.h"
#include "TotalVariableValue.h"
#include "VolumePostprocessor.h"
#include "AreaPostprocessor.h"
#include "PointValue.h"
#include "NodalExtremeValue.h"
#include "ElementExtremeValue.h"
#include "DifferencePostprocessor.h"
#include "NumPicardIterations.h"
#include "FunctionSideIntegral.h"
#include "ExecutionerAttributeReporter.h"
#include "PercentChangePostprocessor.h"
//#include "RealParameterReporter.h"

// vector PPS
#include "ConstantVectorPostprocessor.h"
#include "NodalValueSampler.h"
#include "SideValueSampler.h"
#include "PointValueSampler.h"
#include "LineValueSampler.h"
#include "VectorOfPostprocessors.h"
#include "LeastSquaresFit.h"
#include "ElementsAlongLine.h"
#include "LineMaterialRealSampler.h"

// user objects
#include "LayeredIntegral.h"
#include "LayeredAverage.h"
#include "LayeredSideIntegral.h"
#include "LayeredSideAverage.h"
#include "LayeredSideFluxAverage.h"
#include "NearestPointLayeredAverage.h"
#include "ElementIntegralVariableUserObject.h"
#include "NodalNormalsEvaluator.h"
#include "NodalNormalsCorner.h"
#include "NodalNormalsPreprocessor.h"
#include "SolutionUserObject.h"
#ifdef LIBMESH_HAVE_FPARSER
#include "Terminator.h"
#endif

// preconditioners
#include "PhysicsBasedPreconditioner.h"
#include "FiniteDifferencePreconditioner.h"
#include "SingleMatrixPreconditioner.h"

#include "SplitBasedPreconditioner.h"
#include "Split.h"
#include "ContactSplit.h"
#include "AddSplitAction.h"

// dampers
#include "ConstantDamper.h"
#include "MaxIncrement.h"

// DG
#include "DGDiffusion.h"
#include "DGFunctionDiffusionDirichletBC.h"

// Constraints
#include "TiedValueConstraint.h"
#include "CoupledTiedValueConstraint.h"
#include "AddBoundsVectorsAction.h"
#include "EqualValueConstraint.h"
#include "EqualValueBoundaryConstraint.h"

// ScalarKernels
#include "ODETimeDerivative.h"
#include "FunctionScalarAux.h"
#include "NodalEqualValueConstraint.h"
#include "ParsedODEKernel.h"
#include "QuotientScalarAux.h"

// indicators
#include "AnalyticalIndicator.h"
#include "LaplacianJumpIndicator.h"
#include "GradientJumpIndicator.h"

// markers
#include "ErrorToleranceMarker.h"
#include "ErrorFractionMarker.h"
#include "UniformMarker.h"
#include "BoxMarker.h"
#include "ComboMarker.h"
#include "ValueThresholdMarker.h"
#include "ValueRangeMarker.h"
#include "OrientedBoxMarker.h"

// time steppers
#include "ConstantDT.h"
#include "FunctionDT.h"
#include "IterationAdaptiveDT.h"
#include "SolutionTimeAdaptiveDT.h"
#include "DT2.h"
#include "PostprocessorDT.h"
#include "AB2PredictorCorrector.h"
// time integrators
#include "SteadyState.h"
#include "ImplicitEuler.h"
#include "BDF2.h"
#include "CrankNicolson.h"
#include "ExplicitEuler.h"
//#include "ExplicitMidpoint.h"
//#include "Dirk.h"
#include "LStableDirk2.h"
#include "ImplicitMidpoint.h"
#include "Heun.h"
#include "Ralston.h"
//
#include "SimplePredictor.h"
#include "AdamsPredictor.h"

// MultiApps
#include "TransientMultiApp.h"
#include "FullSolveMultiApp.h"
#include "AutoPositionsMultiApp.h"

// Transfers
#ifdef LIBMESH_HAVE_DTK
  #include "MultiAppDTKUserObjectTransfer.h"
  #include "MultiAppDTKInterpolationTransfer.h"
  #include "MoabTransfer.h"
#endif
#include "MultiAppPostprocessorInterpolationTransfer.h"
#include "MultiAppVariableValueSampleTransfer.h"
#include "MultiAppVariableValueSamplePostprocessorTransfer.h"
#include "MultiAppMeshFunctionTransfer.h"
#include "MultiAppUserObjectTransfer.h"
#include "MultiAppNearestNodeTransfer.h"
#include "MultiAppCopyTransfer.h"
#include "MultiAppInterpolationTransfer.h"
#include "MultiAppPostprocessorTransfer.h"
#include "MultiAppProjectionTransfer.h"
#include "MultiAppPostprocessorToAuxScalarTransfer.h"


// Actions
#include "AddBCAction.h"
#include "AddDiracKernelAction.h"
#include "AddICAction.h"
#include "AddInitialConditionAction.h"
#include "AddKernelAction.h"
#include "AddScalarKernelAction.h"
#include "AddDGKernelAction.h"
#include "AddPeriodicBCAction.h"
#include "AddVariableAction.h"
#include "AddAuxVariableAction.h"
#include "AddPostprocessorAction.h"
#include "AddVectorPostprocessorAction.h"
#include "AddDamperAction.h"
#include "AddFunctionAction.h"
#include "CreateExecutionerAction.h"
#include "DetermineSystemType.h"
#include "SetupTimePeriodsAction.h"
#include "EmptyAction.h"
#include "InitProblemAction.h"
#include "CopyNodalVarsAction.h"
#include "SetupMeshAction.h"
#include "AddMeshModifierAction.h"
#include "SetupMeshCompleteAction.h"
#include "AddOutputAction.h"
#include "CommonOutputAction.h"
#include "AddMaterialAction.h"
#include "GlobalParamsAction.h"
#include "AdaptivityAction.h"
#include "SetupDampersAction.h"
#include "CheckIntegrityAction.h"
#include "SetupQuadratureAction.h"
#include "SetupPreconditionerAction.h"
#include "SetupDebugAction.h"
#include "SetupResidualDebugAction.h"
#include "DeprecatedBlockAction.h"
#include "AddConstraintAction.h"
#include "CreateDisplacedProblemAction.h"
#include "CreateProblemAction.h"
#include "DynamicObjectRegistrationAction.h"
#include "AddUserObjectAction.h"
#include "AddControlAction.h"
#include "AddElementalFieldAction.h"
#include "AddIndicatorAction.h"
#include "AddMarkerAction.h"
#include "SetAdaptivityOptionsAction.h"
#include "AddMultiAppAction.h"
#include "AddTransferAction.h"
#include "AddNodalNormalsAction.h"
#include "SetupTimeStepperAction.h"
#include "SetupTimeIntegratorAction.h"
#include "SetupPredictorAction.h"
#include "AddMortarInterfaceAction.h"
#include "SetupPostprocessorDataAction.h"
#include "MaterialOutputAction.h"
#include "CheckOutputAction.h"
#include "SetupRecoverFileBaseAction.h"

// Outputs
#ifdef LIBMESH_HAVE_EXODUS_API
#include "Exodus.h"
#endif
#include "Nemesis.h"
#include "Console.h"
#include "CSV.h"
#include "VTKOutput.h"
#include "Checkpoint.h"
#include "XDA.h"
#include "GMVOutput.h"
#include "Tecplot.h"
#include "Gnuplot.h"
#include "SolutionHistory.h"
#include "MaterialPropertyDebugOutput.h"
#include "VariableResidualNormsDebugOutput.h"
#include "TopResidualDebugOutput.h"
#include "DOFMapOutput.h"

// Controls
#include "RealFunctionControl.h"
//

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
#include "ForwardStepProblem.h"

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
#include "ErrorMaxFractionMarker.h"

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
//  Moose::registerObjects(_factory);
  registerMooseObjects(_factory);
  ParakeetApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ParakeetApp::associateSyntax(_syntax, _action_factory);
}
void ParakeetApp::registerMooseObjects(Factory & factory)
{
	 // mesh
	registerMesh(FileMesh);
	registerMesh(GeneratedMesh);
	registerMesh(TiledMesh);

	// mesh modifiers
	registerMeshModifier(MeshExtruder);
	registerMeshModifier(SideSetsFromPoints);
	registerMeshModifier(SideSetsFromNormals);
	registerMeshModifier(AddExtraNodeset);
	registerMeshModifier(BoundingBoxNodeSet);
	registerMeshModifier(Transform);
	registerMeshModifier(SideSetsAroundSubdomain);
	registerMeshModifier(SideSetsBetweenSubdomains);
	registerMeshModifier(AddAllSideSetsByNormals);
	registerMeshModifier(SubdomainBoundingBox);
	registerMeshModifier(OrientedSubdomainBoundingBox);
	registerMeshModifier(RenameBlock);

	// problems
	registerProblem(FEProblem);
	registerProblem(DisplacedProblem);

	// kernels
	registerKernel(TimeDerivative);
	registerKernel(CoupledTimeDerivative);
	registerKernel(MassLumpedTimeDerivative);
	registerKernel(Diffusion);
	registerKernel(AnisotropicDiffusion);
	registerKernel(CoupledForce);
	registerKernel(UserForcingFunction);
	registerKernel(BodyForce);
	registerKernel(Reaction);
	registerKernel(MassEigenKernel);

	// bcs
	//	  registerBoundaryCondition(ConvectiveFluxBC);
	//	  registerBoundaryCondition(DirichletBC);
	//	  registerBoundaryCondition(PenaltyDirichletBC);
	//	  registerBoundaryCondition(PresetBC);
	//	  registerBoundaryCondition(NeumannBC);
	//	  registerBoundaryCondition(FunctionDirichletBC);
	//	  registerBoundaryCondition(FunctionPenaltyDirichletBC);
	//	  registerBoundaryCondition(FunctionPresetBC);
	//	  registerBoundaryCondition(FunctionNeumannBC);
	//	  registerBoundaryCondition(MatchedValueBC);
	//	  registerBoundaryCondition(VacuumBC);

	//	  registerBoundaryCondition(SinDirichletBC);
	//	  registerBoundaryCondition(SinNeumannBC);
	//	  registerBoundaryCondition(VectorNeumannBC);
	//	  registerBoundaryCondition(WeakGradientBC);
	//	  registerBoundaryCondition(DiffusionFluxBC);
	//	  registerBoundaryCondition(PostprocessorDirichletBC);
	//	  registerBoundaryCondition(OneDEqualValueConstraintBC);

	// dirac kernels
	registerDiracKernel(ConstantPointSource);

	// aux kernels
	registerAux(ConstantAux);
	registerAux(FunctionAux);
	registerAux(NearestNodeDistanceAux);
	registerAux(NearestNodeValueAux);
	registerAux(PenetrationAux);
	registerAux(ProcessorIDAux);
	registerAux(SelfAux);
	registerAux(GapValueAux);
	registerAux(MaterialRealAux);
	registerAux(MaterialRealVectorValueAux);
	registerAux(MaterialRealTensorValueAux);
	registerAux(MaterialStdVectorAux);
	registerAux(MaterialRealDenseMatrixAux);
	registerAux(MaterialStdVectorRealGradientAux);
	registerAux(DebugResidualAux);
	registerAux(BoundsAux);
	registerAux(SpatialUserObjectAux);
	registerAux(SolutionAux);
	registerAux(VectorMagnitudeAux);
	registerAux(ConstantScalarAux);
	registerAux(QuotientAux);
	registerAux(NormalizationAux);
	registerAux(FunctionScalarAux);
	registerAux(VariableGradientComponent);
	registerAux(ParsedAux);

	// Initial Conditions
	registerInitialCondition(ConstantIC);
	registerInitialCondition(BoundingBoxIC);
	registerInitialCondition(FunctionIC);
	registerInitialCondition(RandomIC);
	registerInitialCondition(ScalarConstantIC);
	registerInitialCondition(ScalarComponentIC);

	// executioners
	registerExecutioner(Steady);
	registerExecutioner(Transient);
	registerExecutioner(InversePowerMethod);
	registerExecutioner(NonlinearEigen);
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,4,0)
#if 0 // This seems to be broken right now -- doesn't work wiith petsc >= 3.4 either
	registerExecutioner(PetscTSExecutioner);
#endif
#endif

	// functions
	registerFunction(Axisymmetric2D3DSolutionFunction);
	registerFunction(ConstantFunction);
	registerFunction(CompositeFunction);
	registerNamedFunction(MooseParsedFunction, "ParsedFunction");
	registerNamedFunction(MooseParsedGradFunction, "ParsedGradFunction");
	registerNamedFunction(MooseParsedVectorFunction, "ParsedVectorFunction");
	registerFunction(PiecewiseConstant);
	registerFunction(PiecewiseLinear);
	registerFunction(SolutionFunction);
	registerFunction(PiecewiseBilinear);
	registerFunction(SplineFunction);
	registerFunction(PiecewiseMultilinear);
	registerFunction(LinearCombinationFunction);

	// materials
	registerMaterial(GenericConstantMaterial);
	registerMaterial(GenericFunctionMaterial);

	// PPS
	registerPostprocessor(AverageElementSize);
	registerPostprocessor(AverageNodalVariableValue);
	registerPostprocessor(NodalSum);
	registerPostprocessor(ElementAverageValue);
	registerPostprocessor(ElementAverageTimeDerivative);
	registerPostprocessor(ElementW1pError);
	registerPostprocessor(ElementH1Error);
	registerPostprocessor(ElementH1SemiError);
	registerPostprocessor(ElementIntegralVariablePostprocessor);
	registerPostprocessor(ElementIntegralMaterialProperty);
	registerPostprocessor(ElementL2Error);
	registerPostprocessor(ElementVectorL2Error);
	registerPostprocessor(ScalarL2Error);
	registerPostprocessor(EmptyPostprocessor);
	registerPostprocessor(NodalVariableValue);
	registerPostprocessor(NumDOFs);
	registerPostprocessor(TimestepSize);
	registerPostprocessor(RunTime);
	registerPostprocessor(PerformanceData);
	registerPostprocessor(NumElems);
	registerPostprocessor(NumNodes);
	registerPostprocessor(NumNonlinearIterations);
	registerPostprocessor(NumLinearIterations);
	registerPostprocessor(Residual);
	registerPostprocessor(ScalarVariable);
	registerPostprocessor(NumVars);
	registerPostprocessor(NumResidualEvaluations);
	registerDeprecatedObjectName(FunctionValuePostprocessor, "PlotFunction", "09/18/2015 12:00");
	registerPostprocessor(Receiver);
	registerPostprocessor(SideAverageValue);
	registerPostprocessor(SideFluxIntegral);
	registerPostprocessor(SideFluxAverage);
	registerPostprocessor(SideIntegralVariablePostprocessor);
	registerPostprocessor(NodalMaxValue);
	registerPostprocessor(NodalProxyMaxValue);
	registerPostprocessor(ElementalVariableValue);
	registerPostprocessor(ElementL2Norm);
	registerPostprocessor(NodalL2Norm);
	registerPostprocessor(NodalL2Error);
	registerPostprocessor(TotalVariableValue);
	registerPostprocessor(VolumePostprocessor);
	registerPostprocessor(AreaPostprocessor);
	registerPostprocessor(PointValue);
	registerPostprocessor(NodalExtremeValue);
	registerPostprocessor(ElementExtremeValue);
	registerPostprocessor(DifferencePostprocessor);
	registerPostprocessor(FunctionValuePostprocessor);
	registerPostprocessor(NumPicardIterations);
	registerPostprocessor(FunctionSideIntegral);
	registerPostprocessor(ExecutionerAttributeReporter);
	registerPostprocessor(PercentChangePostprocessor);
//	registerPostprocessor(RealParameterReporter);

	// vector PPS
	registerVectorPostprocessor(ConstantVectorPostprocessor);
	registerVectorPostprocessor(NodalValueSampler);
	registerVectorPostprocessor(SideValueSampler);
	registerVectorPostprocessor(PointValueSampler);
	registerVectorPostprocessor(LineValueSampler);
	registerVectorPostprocessor(VectorOfPostprocessors);
	registerVectorPostprocessor(LeastSquaresFit);
	registerVectorPostprocessor(ElementsAlongLine);
	registerVectorPostprocessor(LineMaterialRealSampler);

	// user objects
	registerUserObject(LayeredIntegral);
	registerUserObject(LayeredAverage);
	registerUserObject(LayeredSideIntegral);
	registerUserObject(LayeredSideAverage);
	registerUserObject(LayeredSideFluxAverage);
	registerUserObject(NearestPointLayeredAverage);
	registerUserObject(ElementIntegralVariableUserObject);
	registerUserObject(NodalNormalsPreprocessor);
	registerUserObject(NodalNormalsCorner);
	registerUserObject(NodalNormalsEvaluator);
	registerUserObject(SolutionUserObject);
#ifdef LIBMESH_HAVE_FPARSER
	registerUserObject(Terminator);
#endif

	// preconditioners
	registerNamedPreconditioner(PhysicsBasedPreconditioner, "PBP");
	registerNamedPreconditioner(FiniteDifferencePreconditioner, "FDP");
	registerNamedPreconditioner(SingleMatrixPreconditioner, "SMP");
#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,3,0)
	registerNamedPreconditioner(SplitBasedPreconditioner, "SBP");
#endif
	// dampers
	registerDamper(ConstantDamper);
	registerDamper(MaxIncrement);
	// DG
	//	  registerDGKernel(DGDiffusion);
	//	  registerBoundaryCondition(DGFunctionDiffusionDirichletBC);

	// Constraints
	registerConstraint(TiedValueConstraint);
	registerConstraint(CoupledTiedValueConstraint);
	registerConstraint(EqualValueConstraint);
	registerConstraint(EqualValueBoundaryConstraint);

	// Scalar kernels
	registerScalarKernel(ODETimeDerivative);
	registerScalarKernel(NodalEqualValueConstraint);
	registerScalarKernel(ParsedODEKernel);
	registerScalarKernel(QuotientScalarAux);

	// indicators
	registerIndicator(AnalyticalIndicator);
	registerIndicator(LaplacianJumpIndicator);
	registerIndicator(GradientJumpIndicator);

	// markers
	registerMarker(ErrorToleranceMarker);
	registerMarker(ErrorFractionMarker);
	registerMarker(UniformMarker);
	registerMarker(BoxMarker);
	registerMarker(OrientedBoxMarker);
	registerMarker(ComboMarker);
	registerMarker(ValueThresholdMarker);
	registerMarker(ValueRangeMarker);

	// splits
	registerSplit(Split);
	registerSplit(ContactSplit);

	// MultiApps
	registerMultiApp(TransientMultiApp);
	registerMultiApp(FullSolveMultiApp);
	registerMultiApp(AutoPositionsMultiApp);

	// time steppers
	registerTimeStepper(ConstantDT);
	registerTimeStepper(FunctionDT);
	registerTimeStepper(IterationAdaptiveDT);
	registerTimeStepper(SolutionTimeAdaptiveDT);
	registerTimeStepper(DT2);
	registerTimeStepper(PostprocessorDT);
	registerTimeStepper(AB2PredictorCorrector);
	// time integrators
	registerTimeIntegrator(SteadyState);
	registerTimeIntegrator(ImplicitEuler);
	registerTimeIntegrator(BDF2);
	registerTimeIntegrator(CrankNicolson);
	registerTimeIntegrator(ExplicitEuler);
//	registerDeprecatedObjectName(ExplicitMidpoint, "RungeKutta2", "09/25/2015 12:00");
//	registerTimeIntegrator(ExplicitMidpoint);
//	registerDeprecatedObjectName(Dirk, "Dirk", "09/22/2015 12:00");
	registerTimeIntegrator(LStableDirk2);
	registerTimeIntegrator(ImplicitMidpoint);
	registerTimeIntegrator(Heun);
	registerTimeIntegrator(Ralston);
	// predictors
	registerPredictor(SimplePredictor);
	registerPredictor(AdamsPredictor);

	// Transfers
#ifdef LIBMESH_HAVE_DTK
	registerTransfer(MultiAppDTKUserObjectTransfer);
	registerTransfer(MultiAppDTKInterpolationTransfer);
	registerTransfer(MoabTransfer);
#endif
	registerTransfer(MultiAppPostprocessorInterpolationTransfer);
	registerTransfer(MultiAppVariableValueSampleTransfer);
	registerTransfer(MultiAppVariableValueSamplePostprocessorTransfer);
	registerTransfer(MultiAppMeshFunctionTransfer);
	registerTransfer(MultiAppUserObjectTransfer);
	registerTransfer(MultiAppNearestNodeTransfer);
	registerTransfer(MultiAppCopyTransfer);
	registerTransfer(MultiAppInterpolationTransfer);
	registerTransfer(MultiAppPostprocessorTransfer);
	registerTransfer(MultiAppProjectionTransfer);
	registerTransfer(MultiAppPostprocessorToAuxScalarTransfer);

	// Outputs
#ifdef LIBMESH_HAVE_EXODUS_API
	registerOutput(Exodus);
#endif
#ifdef LIBMESH_HAVE_NEMESIS_API
	registerOutput(Nemesis);
#endif
	registerOutput(Console);
	registerOutput(CSV);
#ifdef LIBMESH_HAVE_VTK
	registerNamedOutput(VTKOutput, "VTK");
#endif
	registerOutput(Checkpoint);
	registerNamedOutput(XDA, "XDR");
	registerOutput(XDA);
	registerNamedOutput(GMVOutput, "GMV");
	registerOutput(Tecplot);
	registerOutput(Gnuplot);
	registerOutput(SolutionHistory);
	registerOutput(MaterialPropertyDebugOutput);
	registerOutput(VariableResidualNormsDebugOutput);
	registerOutput(TopResidualDebugOutput);
	registerNamedOutput(DOFMapOutput, "DOFMap");

	// Controls
	registerControl(RealFunctionControl);

	//	  Moose::registered = true;
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
	registerKernel(CFDCellKernel);
	registerKernel(EmptyTimeDerivative);

	registerDGKernel(CFDFaceKernel);

	registerInitialCondition(CFDInitialCondition);

	registerBoundaryCondition(CFDBC);
	registerBoundaryCondition(AdiabaticWall);
	registerBoundaryCondition(SlipWall);
	registerBoundaryCondition(FarFieldPressure);
	registerBoundaryCondition(FarFieldRiemann);
	registerBoundaryCondition(IsoThermalWall);
	registerBoundaryCondition(Symmetric);

	registerFunction(IsoVortexExact);

	registerAux(NSAuxVariable);
	registerAux(NearestWallDistance);

	registerTimeStepper(RatioTimeStepper);

	registerProblem(CFDProblem);
	registerProblem(EulerProblem);
	registerProblem(SodProblem);
	registerProblem(Riemann1DProblem);
	registerProblem(ShockVortexProblem);
	registerProblem(DoubleMachProblem);
	registerProblem(ForwardStepProblem);

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

	registerMarker(ErrorMaxFractionMarker);
}

// External entry point for dynamic syntax association
extern "C" void ParakeetApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ParakeetApp::associateSyntax(syntax, action_factory); }
void
ParakeetApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
#undef registerAction
#define registerAction(tplt, action) action_factory.reg<tplt>(stringifyName(tplt), action)

	syntax.registerActionSyntax("CommonPostProcessorAction", "Postprocessors", "add_postprocessor");
	registerAction(CommonPostProcessorAction, "add_postprocessor");

	syntax.registerActionSyntax("AddMultiVariable", "Problem/Variables");
	syntax.registerActionSyntax("AddMultiVariable", "Variables");
	registerAction(AddMultiVariable, "add_variable");

	syntax.registerActionSyntax("AddMultiIC", "ICs");
	syntax.registerActionSyntax("AddMultiIC", "Variables");
	registerAction(AddMultiIC, "add_ic");

	syntax.registerActionSyntax("AddMultiKernel", "Problem/Kernels");
	registerAction(AddMultiKernel, "add_kernel");

	syntax.registerActionSyntax("AddMultiDGKernel", "Problem/DGKernels");
	registerAction(AddMultiDGKernel, "add_dg_kernel");

	syntax.registerActionSyntax("AddMultiBC", "Problem/BCs/*");
	syntax.registerActionSyntax("AddMultiBC", "BoundaryCondition/*");
	registerAction(AddMultiBC, "add_bc");



	syntax.registerActionSyntax("AddMultiAuxVariableAction", "Problem/AuxVariables");
	registerAction(AddMultiAuxVariableAction, "add_aux_variable");
	registerAction(AddMultiAuxVariableAction, "add_aux_kernel");

	syntax.registerActionSyntax("AddElementalFieldAction", "ArtificialViscosity/Indicators/*");
	syntax.registerActionSyntax("AddIndicatorAction", "ArtificialViscosity/Indicators/*");

	syntax.registerActionSyntax("AddElementalFieldAction", "ArtificialViscosity/Markers/*");
	syntax.registerActionSyntax("AddMarkerAction", "ArtificialViscosity/Markers/*");

#undef registerAction
#define registerAction(tplt, action) action_factory.regLegacy<tplt>(stringifyName(tplt), action)
}
