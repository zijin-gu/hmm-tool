%-------------------------------------------------------------------
% 
% 
%                       %%%%%%%%%%
%                       % genHmm %
%                       %%%%%%%%%%
% 
% 
% This function generates a HMM with N states and M symbols by state.
% function [A,B,Pi]=genhmm(N,M,Np,BAKIS,salto);
% 		
% Entry:	N is the state number
% M is the symbols number by state
% Np is the number of parameters for each symbol
% if Bakis = 1 implements a Bakis HMM (also called left-right),else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% 
% Result:	A, B and Pi are the matrices of the HMM defined within this function.
% 
% A(N, N) is the transition probability matrix with the following conditions:
% 	 (sum(A(i,:)=1))
% 	Aij >= 0  (A(i,j) >= 0)
% 	For a Bakis HMM, we have aij=0 if j<i (A(i,j)=0 for j<i)
% 
% B(N,M) is the distribution probability symbol matrix 
% B(i,k) is the probability to obtain the kth  symbol when we are at the state i.
% Bik >= 0  (B(i,k) >= 0) 
%  (sum(B(i,:)=1))
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% sum Pi(i) =1
% For a Bakis HMM Pi(i)= delta(i)
% With delta(i) is the Kronecker function equal to 1 for i=0 and null for i<>0
% 
% We define here the matrix and update it within the training process of the HMM like the Baum function.
% 
% 
%--------------------------------------------------------------------

function [A,B,Pi]=genhmm(N,M,Np,BAKIS,salto);
randn('seed',0);

% Loop index.
i=0;

% Matrix of the transition probability.
A=rand(N,N);
if BAKIS,
	A=triu(A)-triu(A,salto+1);%+eps.*ones(size(A));
end
A= A./(sum(A')'*ones(1,N));

% Matrix of the symbol probability.
B=cell(Np,1);
for ip=1:Np,
   B{ip}=rand(N,M(ip));
   B{ip}= B{ip}./(sum(B{ip}')'*ones(1,M(ip)))+1e-6;
end
% Distribution of the initial state.
if BAKIS,
   Pi=zeros(N,1);
   Pi(1)=1;
else
   Pi=rand(N,1);
	Pi=Pi./sum(Pi);
end
return