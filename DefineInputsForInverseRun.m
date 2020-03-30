function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)



x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
Lx=max(x)-min(x); Ly=max(y)-min(y);

if CtrlVar.CisElementBased
    xC=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
    yC=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
else
    xC=x; yC=y;
end


if CtrlVar.AGlenisElementBased
    xA=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
    yA=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
else
    xA=x ; yA=y;
end
%

%% Define boundary conditions of adjoint problem
% Generally there is nothing that needs to be done here If BCsAdjoint is not
% modified, then Ua will define the BCs of the adjoint problem based on the BCs
% of the forward problem. 
%BCsAdjoint=BCs; % periodic BCs of forward model -> periodic BCs of adjoint model



%%  Covariance matrices of priors
% 
if CtrlVar.AGlenisElementBased
    CAGlen=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
else
    CAGlen=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
end

if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
    Err=1e-2 ; Sigma=1e3 ; DistanceCutoff=10*Sigma;
    
    if CtrlVar.CisElementBased
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    else
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    end
    
else
    if CtrlVar.CisElementBased
        CC=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
    else
        CC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
    end
end

Priors.CovAGlen=CAGlen;
Priors.CovC=CC;

Priors.s=F.s;
Priors.b=F.b;
Priors.S=F.S;
Priors.B=F.B;


[UserVar,Priors.C,Priors.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[UserVar,Priors.AGlen,Priors.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);

[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);
[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);

Priors.rho=F.rho;
Priors.rhow=F.rhow;

%% Define start values
            sx = 10e3/2;
            sy = 10e3/2;

InvStartValues.C = 1 + zeros(MUA.Nnodes,1);
InvStartValues.AGlen = AGlenVersusTemp(0)-0.5*mean(AGlenVersusTemp(0)) + zeros(MUA.Nnodes,1);
InvStartValues.m=Priors.m;
InvStartValues.n=Priors.n;
Priors.C = 1;
Priors.AGlen = AGlenVersusTemp(0)-0.5*mean(AGlenVersusTemp(0));

%% Define measurements and measurement errors

fprintf(' Creating synthetic data for iC \n')

CtrlVar.doDiagnostic=1;
[UserVar,F.C,F.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[UserVar,F.AGlen,F.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[F.AGlen,F.n]=TestAGlenInputValues(CtrlVar,MUA,F.AGlen,F.n);
[F.C,F.m]=TestSlipperinessInputValues(CtrlVar,MUA,F.C,F.m);


[UserVar,RunInfo,F,l,Kuv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);

Priors.TrueC=F.C;
Priors.TrueAGlen=F.AGlen;



Meas.us=F.ub ;
Meas.vs=F.vb;
Meas.ws=F.ub*0;

VelScale=mean(F.ub);

usError = 1;
vsError = 1;
wsError = 1;

Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.wsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,wsError.^2,MUA.Nnodes,MUA.Nnodes);
 
% if add errors

Meas.us=Meas.us+UserVar.AddDataErrors*usError.*randn(MUA.Nnodes,1);
Meas.vs=Meas.vs+UserVar.AddDataErrors*vsError.*randn(MUA.Nnodes,1);
Meas.ws=Meas.ws+UserVar.AddDataErrors*wsError.*randn(MUA.Nnodes,1);

end
