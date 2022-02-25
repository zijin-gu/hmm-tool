%------------------------------------------------------------------------
% 
%               %%%%%%%%%%%%%%
%               %    Alfa    %
%               %%%%%%%%%%%%%%
% 
% This function calculates the alpha from the HMM defined with the matrix A, B and Pi. 
% The alfa is scaled to avoid the problem of the precision.
% 	function [alfa,c]=Alfa(A,B,Pi,O)
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B{Np}(N,M) is the distribution probability symbol matrix  for each parameter.
% B(i,k) is the probability to obtain the kth  symbol when we are at the state i.
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% 
% O is a matrix with the sequence monitored used to estimate all the parameters.
% 
% Result:	alfa is the forward probability matrix.
% 
% alfa(i,t) when alfa is not scaled is the probability to observe the sequence O(1,:),...,O(t,:) 
% for the state I at the instant t for all the parameters.
% We have in particular P(O/(A, B, Pi))= sum(alfa(:,T)). 
% 
% When alfa is scaled, we have sum(alfa(:,t))=1.
% 
% c(T,1) is the vector where we store the scale value for the instant t.
% We have c(T,1)= ones(T,1) if the alfa is not scaled. 

%--------------------------------------------------------------------------
function [alfa,c]=Alfa(A,B,Pi,O)

% Variables for the hmm:
N=size(A,1);
% Variables for the vectores sequences.
% T: Length of the sequence of states and observations.
T=size(O{1},1);
% Matrix of the forward probability
alfa=zeros(N,T);
% definition of the sacaled variables alfa and beta.
c=ones(T,1);
% Loop index.
t=0;

% We calculate those variables in order to compute alfa and beta
prob=zeros(N,T);
for t=1:T
   prob(:,t)=prodBO(B,O,t);
end

% Scaled matrix forward:
alfa(:,1)=Pi.*(prob(:,1));
c(1)=1/sum(alfa(:,1));		% Put c(1)=1 if we do not scale
alfa(:,1)=alfa(:,1).*c(1);	% Scaled.
for t=2:T,
	alfa(:,t)=(alfa(:,t-1)'*A)'.*(prob(:,t));
	c(t)=1/sum(alfa(:,t));	% Put c(1)=1 if we do not scale
	alfa(:,t)=alfa(:,t).*c(t)+realmin;
end;

return