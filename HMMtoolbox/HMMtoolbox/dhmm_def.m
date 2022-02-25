% ------------------------------------------------------------	
% 
% 
%                  %%%%%%%%%%%%
%                  % DHMM_DEF % 
%                  %%%%%%%%%%%%
% 
% 
% 
% This function is described only as example to set up the HMM
% 
% This function defines the parameters of the discrete HMM. We can split the parameters in 4 groups:
% 
% Entry:	fhmm is the HMM´s file name.
% 
% 	vDB defines the database:
% vDB=[' nc ng agrup Np'];
% where nc is the number of class, ng is the number of group,
% agrup defines how to gather the parameters together
% NP is a vector with the number of each group of parameter.
% We have the following relations:
% [nc,ng]=size(vl);
% agrup{ig}=[1 .... size(vl{1,ig}{1},2)+1];
% Np(ig)=length(agrup{ig})-1;
% 
% 	vHMM defines the HMM :
% 	vHMM=[' BAKIS salto maxiter umbral maxitermi TOPN Ne Ns A B Pi'];
% if Bakis = 1 implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% Maxiter is the maximum iteration number in the Bakis HMM (condition to stop the algorithm)
% Umbral is the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion) 
% TOPN number of labels to take into account for the multi-labelling during the training phase
% Maxitermi is the maximum iteration for the initial model
% Ne states number 
% Ns symbols number for  every state
% A, B and Pi are the matrices of the HMM (it will be calculated thanks to the HMM) but must de defined here.
% Salhmm is a matrix used to store the probabilities values
% We have the following relations:
% TOPN=cell(ng,1);
% 
% 	 vQ is the matrix parameter for the library :
% vVQ=[' LBG dpztoLBG maxiterVQ umbralVQ men biblio'];
% where LBG is the choice of the algorithm to compile the library
% LBG = 1 the algorithm is LBG, in the other case we choose kMedias
% dptzoLBG percentage maximum for the distance of the code vector in the algorithm LBG (see gen_bib) 
% maxiterVQ is the maximum number of iteration to compute the library (condition to stop the algorithm)
% umbralVQ is the threshold condition to stop the library computation (condition to stop the algorithm)
% men=1 or 0 if men =1 the library generates message else no.
% biblio is a matrix parameter for the library computation (we define it here but calculate in the gen_bib function)
% 
% 	 vTEST is the matrix parameter for the test :
% vTEST=[' TOPNtest salhmm'];
% where TOPNtest is the number of labels to take into account for the multi-labeling during the test phase
% Salhmm is a matrix used to store the probabilities values
% We have the following relations:
% TOPNtest=cell(ng,1);
% salhmm=cell(nc,ng);
% ---------------------------------------------------------
function dhmm_def(fhmm)

% Name of the set-up file
%fhmm='hmm';


% vDB defines the database:
vDB=[' nc ng agrup Np'];
% where nc is the number of class, ng is the number of group,
% agrup defines how to gather the parameters together
% NP is a vector with the number of each group of parameter.
% with [nc,ng]=size(vl);
nc=10;
ng=2;                                                   %seb 3

% Definition of the gathering
% with agrup{ig}=[1 .... size(vl{1,ig}{1},2)+1];
agrup=cell(ng,1);
Np=zeros(ng,1);
for ig=1:ng
   agrup{ig}=[1 2 3];% 1 2 3
end
agrup{1}=[1 2 3 4 ];% 1 2
agrup{2}=[1 2 3 4];% 1 2
for ig=1:ng
  Np(ig)=length(agrup{ig})-1;
end

% VARIABLeS OF THE HMM
vHMM=[' BAKIS salto maxiter umbral maxitermi TOPN Ne Ns A B Pi'];
% if Bakis = 1 implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% Maxiter is the maximum iteration number in the Bakis HMM (condition to stop the algorithm)
% Umbral is the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion) 
% Maxitermi is the maximum iteration for the initial model
% TOPN number of labels to take into account for the multi-labelling during the training phase
% Ne states number 
% Ns symbols number for  every state
% A, B and Pi are the matrices of the HMM (it will be calculated thanks to the HMM) but must de defined here.
% Salhmm is a matrix used to store the probabilities values
% We have the following relations:
% TOPN=cell(ng,1);
TOPN=cell(ng,1);
for ig=1:ng   
   TOPN{ig}=1.*ones(Np(ig),1);
end


% Number of sates.
Ne=10.*ones(nc,ng);
% Number of symbols per state.
Ns=cell(ng,1);
for ig=1:ng   
   Ns{ig}=32.*ones(Np(ig),1);
end

BAKIS=1;

salto=1;

maxiter=30;

umbral=0.005;

maxitermi=10;
% Matrices of the HMM
A=cell(nc,ng);
B=cell(nc,ng);
Pi=cell(nc,ng);
%for ic=1:nc,
%   for ig=1:ng
%      A{ic,ig}=zeros(Ne(ic,ig),Ne(ic,ig));
%      B{ic,ig}=cell(Np(ig));
%      for ip=1:Np(ig)
%         B{ic,ig}{ip}=zeros(Ne(ic,ig),Ns{ig}(ip));
%      end
%      P{ic,ig}=zeros(Ne(ic,ig),1);
%   end
%end


% vQ is the matrix parameter for the library :
vVQ=[' LBG dpztoLBG maxiterVQ umbralVQ men biblio'];

% where LBG is the choice of the algorithm to compile the library
% LBG = 1 the algorithm is LBG, in the other case we choose kMedias
% dptzoLBG percentage maximum for the distance of the code vector in the algorithm LBG (see gen_bib) 
% maxiterVQ is the maximum number of iteration to compute the library (condition to stop the algorithm)
% umbralVQ is the threshold condition to stop the library computation (condition to stop the algorithm)
% men=1 or 0 if men =1 the library generates message else no.
% biblio is a matrix parameter for the library computation (we define it here but calculate in the gen_bib function)

LBG=1;

dpztoLBG=0.1;

maxiterVQ=40;

umbralVQ=0.001;

men=1;
% we definate the varaiables for the Libray
biblio=cell(ng,1);
%for ig=1:ng   
%   biblio{ig}=cell(Np(ig),1);
%   for ip=1:Np(ig)
%      dinic=agrup{ig}(ip);dfinal=agrup{ig}(ip+1)-1;
%      dim=dfinal-dinic+1;
%      biblio{ig}{ip}=zeros(Ns{ig}(ip),dim);
%   end
%end

% VARIABLES FOR THE TEST
vTEST=[' TOPNtest salhmm'];
% vTEST is the matrix parameter for the test :
% where TOPNtest is the number of labels to take into account for the multi-labeling during the test phase
% Salhmm is a matrix used to store the probabilities values
% We have the following relations:
% TOPNtest=cell(ng,1);
% salhmm=cell(nc,ng);

TOPNtest=cell(ng,1);
for ig=1:ng   
   TOPNtest{ig}=1.*ones(Np(ig),1);
end
% Salhmm is a matrix used to store the probabilities values for the output per test/model
salhmm=cell(nc,ng);

% we save the file and the parameters of the HMM
guardar=['save ',fhmm,vDB,vHMM,vVQ,vTEST,' vDB vHMM vVQ vTEST guardar'];
eval([guardar]);
