%---------------------------------------------------------------------
%
%
%
%                           %%%%%%%%%%% 	                    
%                           %AlfaBetac%
%                           %%%%%%%%%%%
%                       
%
%
% This function calculates the alpha and beta from the HMM defined with the matrix A, B and Pi. The alfa and Beta are scaled to avoid the problem of the precision.
% 	function [alfa,beta,c]=alfabetac(A,B,Med,Var,Pi,vectores,agrup)
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B{Np}{N}(Ngauss{ip},1) is the weight for the Gaussians of the HMM.
% 
% Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the medians for the Gaussians of the HMM.
% 
% Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the variances for the Gaussians of the HMM.
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% 
% Vectores: is a matrix with the sequence monitored used to estimate all the parameters.
% 
% Result:	alfa is the forward probability matrix.
%           alfa(i,t) when alfa is not scaled is the probability to observe the sequence O(1,:),....,O(t,:) for the state i at the instant t for all the parameters.
%           We have in particular P(O/(A, B, Pi))= sum(alfa(:,T)). 
%           When alfa is scaled, we have sum(alfa(:,t))=1.
% 
%         	beta is the backward probability matrix.
%           beta(i,t) when beta is not scaled is the probability to observe the sequence O(t+1,:),....,O(T,:) 
%           for the state i at the instant t for all the parameters.
% 
%           c(T,1) is the vector where we store the scale value for the instant t.
%           We have c(T,1)= ones(T,1) if the alfa and beta are not scaled. 
% 
%           We have the following relations on the alfa and beta:
% 	        Relations between alfa scaled and not scaled
% 		    h=cumprod(c);
% 		    alfa_scaled(:,t)= h(t)*alfa_no_scaled;
% 	        Relations between beta scaled and not scaled
% 		    g=cumprod(c(T:-1:1);
% 		    beta_scaled(:,t)=g(t)*beta_no_scaled(:,t);
% 
% If alfa and beta are not scaled, the probability to observe O is (this product is independent of t):
% 	alfa( :,t)'*beta( :,t) = sum (alfa(T,:))
% 
% If alfa and beta are scaled (it depends on t):
% 	alfa(t,:)'*beta(t,:) = c(t)*sum (alfa(T,:))
%
%-----------------------------------------------------------------------------------

function [alfa,beta,c]=alfabetac(A,B,Med,Var,Pi,vectores,agrup)

% Variables gfor the hmm.
% N: States Number of the HMM.
N=size(B{1},1);
% Variables of the observation sequence.
% T: length of the sequence of states and observations.
T=size(vectores,1);
% Forward probability matrix.
alfa=zeros(N,T);
% ¡backward probability matrix .
beta=zeros(N,T);
% scale variable used for the alfa and beta.
c=ones(T,1);
% Loop.
t=0;

% We calculate the intermediate variable to compute the alfa and beta
prob=zeros(N,T);
for t=1:T
   prob(:,t)=probsimb(B,Med,Var,vectores(t,:),agrup);
end

% Scaled forward matrix..
alfa(:,1)=Pi.*(prob(:,1));
c(1)=1/sum(alfa(:,1));		% c(1)=1 to avoid the scaling.
alfa(:,1)=alfa(:,1).*c(1);	% alfa scaled.
for t=2:T,
	alfa(:,t)=(alfa(:,t-1)'*A)'.*(prob(:,t));
	c(t)=1/sum(alfa(:,t));	% c(t)=1 to avoid the scaling.
	alfa(:,t)=alfa(:,t).*c(t)+realmin;
end;

% computation of the backward probability matrix.
% c=ones(T,1) to avoid the scaling.
beta(:,T)=ones(N,1)*c(T);
for t=T-1:-1:1,
	beta(:,t)=A*(prob(:,t+1).*beta(:,t+1))*c(t)+realmin;
end;
return