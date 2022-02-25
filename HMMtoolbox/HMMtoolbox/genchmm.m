%-------------------------------------------------------------------
%
%                              %%%%%%%%%%%
%                              % GENCHMM %
%                              %%%%%%%%%%%
% 
% 
% This function generates a CHMM with N states and Ngauss Gaussians.
% function [A,B,Med,Var,Pi]=genchmm(N,vcombl,vmean,vstd,agrup,Ngauss,Np,BAKIS,salto);		
% Entry:	N is the state number
% vcombl is a parameter of the mixture of Gaussians. See the function gmminit of the netlab software.
% vmean is the vector of the initial means
% vstd is the vector of the initial diversion
% agrup(Np+1,1) defines how to cluster the parameters together
% Ngauss(Np,1) number of Gaussians
% Np number of parameters
% if Bakis = 1 implements a Bakis HMM (also called left-right),else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% 
% Result:	A, B, Med, Var and Pi are the matrices of the HMM defined within this function.
% 
% A(N, N) is the transition probability matrix with the following conditions:
% 	 (sum(A(i,:)=1))
% 	Aij >= 0  (A(i,j) >= 0)
% 	For a Bakis HMM, we have aij=0 if j<i (A(i,j)=0 for j<i)
% 
% B{Np}{N}(Ngauss{ip},1) is a structure with the weight of each Gaussian for each state and to obtain each parameter.
% 
% Med: Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the mean for each Gaussian in each state and for each gathering of parameters
% 
% Var: Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the variance for each Gaussian in each state and for each gathering of parameters
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% sum Pi(i) =1
% For a Bakis HMM Pi(i)= delta(i)
% With delta (i) is the Kronecker function equal to 1 for i=0 and null for i<>0
% 
% We define here the matrix and update it within the training process of the HMM like the Baum function.
% 
%--------------------------------------------------------------------

function [A,B,Med,Var,Pi]=genchmm(N,vcombl,vmean,vstd,agrup,Ngauss,Np,BAKIS,salto);
randn('seed',0);

% Loop index.
i=0;
% shift of the mean and variance
fdpzto=0.2;
% size of each group
dim=diff(agrup);

% Matrix of the Transition probability.
A=rand(N,N);
if BAKIS,
	A=triu(A)-triu(A,salto+1);%+eps.*ones(size(A));
end
A= A./(sum(A')'*ones(1,N));

% Matrix of the probability of symbol.
% B=cell(Np,1);
% for ip=1:Np,
%    B{ip}=cell(N,1);
%    for ie=1:N
%       B{ip}{ie}=rand(Ngauss(ip),1);
%       B{ip}{ie}= B{ip}{ie}./sum(B{ip}{ie})+1e-6;   
%    end
% end
B=cell(Np,1);
for ip=1:Np,
   B{ip}=cell(N,1);
   for ie=1:N
      B{ip}{ie}=vcombl{ip}.*(1+fdpzto.*rand(Ngauss(ip),1));
      B{ip}{ie}= B{ip}{ie}./sum(B{ip}{ie})+1e-6;   
   end
end

% Matrix of the means and variances.
Med=cell(Np,1);
Var=cell(Np,1);
for ip=1:Np,
   Med{ip}=cell(N,1);
   Var{ip}=cell(N,1);
   for ie=1:N
         Med{ip}{ie}=vmean{ip}.*(1+fdpzto*randn(Ngauss(ip),dim(ip)));
         Var{ip}{ie}=vstd{ip}.*(1+fdpzto*randn(Ngauss(ip),dim(ip)));
   end
end

% Distribution of the initial states.
if BAKIS,
   Pi=zeros(N,1);
   Pi(1)=1;
else
   Pi=rand(N,1);
	Pi=Pi./sum(Pi);
end
return