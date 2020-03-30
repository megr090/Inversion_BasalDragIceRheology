function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar)


%%  This input file is used if Ua is run directly from the source code folder.
%  
% 
%
%


%%
UserVar.MisExperiment='ice0';            % This I use in DefineMassBalance
UserVar.Outputsdirectory='ResultsFiles'; % This I use in UaOutputs
UserVar.RunType = 'icestream';
%%

CtrlVar.Experiment=['Bindschadler',UserVar.MisExperiment];   
%% Types of run
%
CtrlVar.doInverseStep=1 ;

CtrlVar.TimeDependentRun=0; 
CtrlVar.TotalNumberOfForwardRunSteps=3;
CtrlVar.TotalTime=100;
CtrlVar.Restart=0;  


CtrlVar.dt=0.01; 
CtrlVar.time=0; 

CtrlVar.UaOutputsDt=0; % interval between calling UaOutputs. 0 implies call it at each and every run step.
                       % setting CtrlVar.UaOutputsDt=1; causes UaOutputs to be called every 1 years.
                       % This is a more reasonable value once all looks OK.

CtrlVar.ATStimeStepTarget=1;
%% Restart
CtrlVar.Restart=0;  CtrlVar.WriteRestartFile=1;
CtrlVar.NameOfRestartFiletoRead='InverseRestart.mat';
CtrlVar.NameOfRestartFiletoWrite=CtrlVar.NameOfRestartFiletoRead;

%% Inverse   -inverse

CtrlVar.Inverse.MinimisationMethod='MatlabOptimization'; % {'MatlabOptimization','UaOptimization'}
CtrlVar.Inverse.Iterations=100;

CtrlVar.Inverse.InvertFor = 'logAGlenlogC';
CtrlVar.CisElementBased=0;
CtrlVar.AGlenisElementBased=0;
CtrlVar.Inverse.CalcGradI=true;  

% UaOptimization parameters, start :
CtrlVar.Inverse.GradientUpgradeMethod='confjgrad' ; %{'SteepestDecent','conjgrad'}
CtrlVar.Inverse.InitialLineSearchStepSize=[];
CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize=1e-20; % minimum step size in backtracking
CtrlVar.Inverse.MinimumRelativelLineSearchStepSize=1e-5; % minimum fractional step size relative to initial step size
CtrlVar.Inverse.MaximumNumberOfLineSeachSteps=50;
% end, UaOptimization parameters

% MatlabOptimisation parameters, start :
CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','trust-region',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'FunctionTolerance',1e-10,...
    'StepTolerance',1e-20,...
    'PlotFcn',{@optimplotfval,@optimplotstepsize},...
    'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI,...
    'HessianFcn','objective');

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'MaxIterations',CtrlVar.Inverse.Iterations,...
    'MaxFunctionEvaluations',1000,...
    'Display','iter-detailed',...
    'OutputFcn',@fminuncOutfun,...
    'Diagnostics','on',...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',1e-20,...
    'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI);


% end, MatlabOptimisation parameters.

CtrlVar.Inverse.InfoLevel=1;  % Set to 1 to get some basic information, >=2 for additional info on backtrackgin,
                                 % >=100 for further info and plots

CtrlVar.InfoLevelNonLinIt=0; CtrlVar.InfoLevel=0;

CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint' ; % {'Adjoint','FixPointC'}
CtrlVar.Inverse.AdjointGradientPreMultiplier='I'; % {'I','M'}

% Testing adjoint parameters, start:
CtrlVar.Inverse.TestAdjoint.isTrue=0; % If true then perform a brute force calculation 
                                      % of the directinal derivative of the objective function.  
CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType='central' ; % {'central','forward'}

CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize=1e-8 ;
CtrlVar.Inverse.TestAdjoint.iRange=[] ;  % range of nodes/elements over which brute force gradient is to be calculated.
                                         % if left empty, values are calulated for every node/element within the mesh. 
                                         % If set to for example [1,10,45] values are calculated for these three
                                         % nodes/elements.
% end, testing adjoint parameters. 
                                                    
UserVar.AddDataErrors=0;
CtrlVar.LinSolveTol=1e-6;

load('effsr.mat');
load('InverseRestart.mat','GF');
e = eff;
e(abs(e)<5e-4)=5e-4;
e(abs(e)>1e-2)=1e-2;
f = GF.node;
f(f<0.9)=0;
f(f>=0.9)=1;
f = abs(f-1);


CtrlVar.Inverse.Regularize.Field=CtrlVar.Inverse.InvertFor;
CtrlVar.Inverse.Regularize.C.gs=1;
CtrlVar.Inverse.Regularize.C.ga=1;
CtrlVar.Inverse.Regularize.logC.ga=10; 
CtrlVar.Inverse.Regularize.logC.gs=1e5;  

CtrlVar.Inverse.Regularize.AGlen.gs=1;
CtrlVar.Inverse.Regularize.AGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.ga =1e-2./((1-f).*e+(1e-2.*f));
CtrlVar.Inverse.Regularize.logAGlen.gs=1e2/((1-f).*e+(1e-2.*f));

CtrlVar.Inverse.DataMisfit.HessianEstimate='0'; % {'0','I','MassMatrix'}

CtrlVar.Inverse.DataMisfit.Multiplier=1;
CtrlVar.Inverse.Regularize.Multiplier=1;

CtrlVar.Inverse.DataMisfit.FunctionEvaluation='integral';
CtrlVar.MUA.MassMatrix=1;
CtrlVar.MUA.StiffnessMatrix=1;
%% Reading in mesh
CtrlVar.ReadInitialMesh=1;    % if true then read FE mesh (i.e the MUA variable) directly from a .mat file
                              % unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='AdaptMesh.mat';
CtrlVar.SaveInitialMeshFileName='NewMeshFile.mat';
%% Plotting options
CtrlVar.doplots=0;
CtrlVar.PlotMesh=0; 
CtrlVar.PlotBCs=0;
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
CtrlVar.doRemeshPlots=1;
CtrlVar.PlotXYscale=1000; 
%%

CtrlVar.TriNodes=3;

% very coarse mesh resolution
CtrlVar.MeshSizeMax=10e3;    % max element size
CtrlVar.MeshSizeMin=2e3;     % min element size
CtrlVar.MeshSize = 5000;

CtrlVar.MaxNumberOfElements=250e3;           % max number of elements. If #elements larger then CtrlMeshSize/min/max are changed

CtrlVar.AdaptMesh=1;           % 
CtrlVar.SaveAdaptMeshFileName='AdaptMesh.mat'; 



CtrlVar.AdaptMeshInitial=0 ;       % if true, then a remeshing will always be performed at the inital step
CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doRemeshPlots=1 to get plots)
CtrlVar.doAdaptMeshPlots=0;       % if true and if CtrlVar.doplots true also, then do some extra plotting related to adapt meshing

CtrlVar.RefineCriteriaWeights=1;                %  
CtrlVar.RefineCriteriaFlotationLimit=NaN;     

CtrlVar.RefineCriteria={'effective strain rates'};
  
CtrlVar.AdaptMeshInterval=1;  % number of run-steps between mesh adaptation
CtrlVar.AdaptMeshMaxIterations=1;


CtrlVar.MeshAdapt.GLrange=[10000 5000 ; 3000 2000];

%% Pos. thickness constraints
CtrlVar.ThickMin=1; % minimum allowed thickness without (potentially) doing something about it
CtrlVar.ResetThicknessToMinThickness=0;  % if true, thickness values less than ThickMin will be set to ThickMin
CtrlVar.ThicknessConstraints=1  ;        % if true, min thickness is enforced using active set method
CtrlVar.ThicknessConstraintsItMax=5  ; 

load('InverseRestart_Works.mat','CtrlVarInRestartFile');
MeshBoundaryCoordinates = CtrlVarInRestartFile.MeshBoundaryCoordinates;

 
end
