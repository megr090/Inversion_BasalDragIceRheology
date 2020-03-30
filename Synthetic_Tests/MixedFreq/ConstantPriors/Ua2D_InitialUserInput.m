function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput(UserVar,CtrlVar)



if ~isfield(UserVar,'RunType')
    UserVar.RunType='IceStream';
end


CtrlVar.doplots=0;
CtrlVar.InverseIceStream=0;

xd = 100e3; xu = -100e3; yl = 10e3; yr = -10e3;
MeshBoundaryCoordinates=[xu yr ; xd yr ; xd yl ; xu yl];
MeshBoundaryCoordinates=flipud(MeshBoundaryCoordinates);
%% Types of runs
CtrlVar.DevelopmentVersion=1;  %
CtrlVar.doInverseStep=1 ;
CtrlVar.TriNodes=3;

%% Restart
CtrlVar.Restart=0;  CtrlVar.WriteRestartFile=1;
CtrlVar.NameOfRestartFiletoRead=['InverseRestart.mat'];
CtrlVar.NameOfRestartFiletoWrite=CtrlVar.NameOfRestartFiletoRead;
%% Inverse   -inverse

CtrlVar.Inverse.MinimisationMethod='MatlabOptimization'; % {'MatlabOptimization','UaOptimization'}
CtrlVar.Inverse.Iterations=500;
CtrlVar.Inverse.InvertFor = 'logAGlenlogC';
CtrlVar.CisElementBased=0;
CtrlVar.AGlenisElementBased=0;
CtrlVar.Inverse.CalcGradI=true;  

% UaOptimization parameters, start :
CtrlVar.Inverse.GradientUpgradeMethod='conjgrad' ; %{'SteepestDecent','conjgrad'}
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


CtrlVar.Inverse.Regularize.Field=CtrlVar.Inverse.InvertFor;
CtrlVar.Inverse.Regularize.C.gs=1;
CtrlVar.Inverse.Regularize.C.ga=1;
CtrlVar.Inverse.Regularize.logC.ga=10;
CtrlVar.Inverse.Regularize.logC.gs=1e4; % 1e6  works well with I

CtrlVar.Inverse.Regularize.AGlen.gs=1;
CtrlVar.Inverse.Regularize.AGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.gs=1e5;

CtrlVar.Inverse.DataMisfit.HessianEstimate='0'; % {'0','I','MassMatrix'}

CtrlVar.Inverse.DataMisfit.Multiplier=1;
CtrlVar.Inverse.Regularize.Multiplier=1;

CtrlVar.Inverse.DataMisfit.FunctionEvaluation='integral';
CtrlVar.MUA.MassMatrix=1;
CtrlVar.MUA.StiffnessMatrix=1;

%% Mesh generation and remeshing parameters

CtrlVar.meshgeneration=1;
CtrlVar.GmshMeshingAlgorithm=8;    % see gmsh manual
% 1=MeshAdapt
% 2=Automatic
% 5=Delaunay
% 6=Frontal
% 7=bamg
% 8=DelQuad (experimental)

CtrlVar.ReadInitialMesh=0;    % if true then read FE mesh (i.e the MUA variable) directly from a .mat file
% unless the adaptive meshing option is used, no further meshing is done.
CtrlVar.ReadInitialMeshFileName='AdaptMeshFile10k';
%CtrlVar.ReadInitialMeshFileName='UniformMesh';
CtrlVar.SaveInitialMeshFileName='NewMeshFile.mat';



CtrlVar.MeshSize = 1000;
CtrlVar.MeshSizeMin=0.05*CtrlVar.MeshSize;    % min element size
CtrlVar.MeshSizeMax=CtrlVar.MeshSize;

CtrlVar.GmshGeoFileAdditionalInputLines{1}='Periodic Line {1,2} = {3,4};';
CtrlVar.AdaptMesh=0;
CtrlVar.SaveAdaptMeshFileName='AdaptMeshFile';

CtrlVar.RefineCriteria='effective strain rates';
CtrlVar.RefineCriteriaFlotationLimit=NaN ;
CtrlVar.RefineCriteriaWeights=1;
CtrlVar.RefineDiracDeltaInvWidth=1000;

%% plotting
CtrlVar.PlotXYscale=1000;
CtrlVar.PlotLabels=0 ; CtrlVar.PlotMesh=1; CtrlVar.PlotBCs=0;

CtrlVar.Inverse.NameOfRestartInputFile=CtrlVar.Inverse.NameOfRestartOutputFile;
%CtrlVar.Inverse.NameOfRestartInputFile=filename;
%%

CtrlVar.Experiment=['I-',CtrlVar.Inverse.AdjointGradientPreMultiplier,'-',num2str(CtrlVar.TriNodes)];

end
