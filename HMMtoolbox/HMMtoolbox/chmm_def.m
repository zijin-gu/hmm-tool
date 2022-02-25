%---------------------------------------------------------------------------
%
%                           %%%%%%%%%%% 	                    
%                           %CHMM_DEF %
%                           %%%%%%%%%%%	
% 
%
% This function is described only as example to set up the HMM.
% 
% function chmm_def(fhmm)
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
% 	vHMM =[' BAKIS salto maxiter umbral maxitermi Ngauss Ne A B Med Var Pi'];
% if Bakis = 1 it implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% Maxiter is the maximum iteration number in the Bakis HMM (condition to stop the algorithm)
% Umbral is the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion) 
% Maxitermi is the iteration number to find the initial model.
% Ne states number 
% Ngauss is the number of Gaussians in the mixture of Gaussians used to represent the distribution of the observation by state.
% A, B and Pi are the matrices of the HMM (it will be calculated thanks to the CHMM) but must de defined here.
% Med and Var are the matrices of the mean and variances of the Gaussians.
% 
%  vTEST is the matrix parameter for the test :
% vTEST=[' salhmm'];
% Salhmm is a matrix used to store the probabilities values
% We have the following relations:
% salhmm=cell(nc,ng);
%
%
%
% -------------------------------------------------------------------------


function chmm_def(fhmm)

% name of the file of the CHMM
%fhmm='hmm';

iniciar=1;

% VARIABLES of the Database
vDB=[' nc ng agrup Np'];
% classes number and groups number with 
% the relation [nc,ng]=size(vl);
nc=10;
ng=7;

% Definition and numbers of the inputs.
% each realisation is a group of vectors and every vector contains
% a sequence of components that we gather in oreder to create the paramenters inputs.
% We define the gathering of the components and the number of groups.
% we have the relation agrup{ig}=[1 .... size(vl{1,ig}{1},2)+1];
grup=cell(ng,1);
Np=zeros(ng,1);
for ig=1:ng
   agrup{ig}=[1 2 3];
end
for ig=1:ng
   Np(ig)=length(agrup{ig})-1;
end



% VARIABLES of the HMM
vHMM=[' BAKIS salto maxiter umbral maxitermi Ngauss Ne A B Med Var Pi men'];
% Number of gaussians by groups and parameter.
% All the states have the same number of gaussians
Ngauss=cell(ng,1);
for ig=1:ng   
   Ngauss{ig}=6.*ones(Np(ig),1);
end
men=1;
% Variables for the HMM
% States number.
Ne=40.*ones(nc,ng);
% if Bakis = 1 it implements a Bakis HMM (also called left-right), else an ergodic HMM is defined.
BAKIS=1;
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
salto=1;
% Maximum number of iterations in the Bakis HMM (condition to stop the algorithm).
maxiter=10;
% the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion) 
umbral=0.05;
% the iteration number to find the initial model.
maxitermi=5;
% Matrices of the HMM
A=cell(nc,ng);
B=cell(nc,ng);
Med=cell(nc,ng);
Var=cell(nc,ng);
Pi=cell(nc,ng);
for ic=1:nc,
   for ig=1:ng
      A{ic,ig}=zeros(Ne(ic,ig),Ne(ic,ig));
      B{ic,ig}=cell(Np(ig),1);
      Med{ic,ig}=cell(Np(ig),1);
      Var{ic,ig}=cell(Np(ig),1);      
      for ip=1:Np(ig)
         B{ic,ig}{ip}=cell(Ne(ic,ig),1);
         Med{ic,ig}{ip}=cell(Ne(ic,ig),1);
         Var{ic,ig}{ip}=cell(Ne(ic,ig),1);
         for ie=1:Ne(ic,ig),
            B{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),1);      
            Med{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),agrup{ig}(ip+1)-agrup{ig}(ip));
            Var{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),agrup{ig}(ip+1)-agrup{ig}(ip));
         end
      end
      Pi{ic,ig}=zeros(Ne(ic,ig),1);
   end
end


% VARIABLES for the TEST
vTEST=[' salhmm'];

% outputs vectors for each sequence of the test/model
salhmm=cell(nc,ng);

% we save the HMM
guardar=['save ',fhmm,vDB,vHMM,vTEST,' iniciar vDB vHMM vTEST guardar'];
eval([guardar]);
