%------------------------------------------------------------------------
%
%                           %%%%%%%%%%%
%                           % ProbSec %
%                           %%%%%%%%%%% 
% 
% This function estimates the log of the probability to monitor a sequence O given an HMM. We use the log to avoid numerical problem.
% 
% 	function [logPO,logalfaT]=probsec(A,B,Pi,O)
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B{Np}(N,M) is the distribution probability symbol matrix  for each parameter.
% B(i,k) is the probability to obtain the kth  symbol when we are at the state i.
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% 
% O is a matrix with the sequence to estimate.
% 
% Result:	logPO is the log of the probability to monitor the sequence O.
% 	logalfaT: is the log (alfa(:, T)) where alfa is the forward probability matrix and T is the running time for the sequence O. 
% 
% We use in this function the prodBO function.
% 
%--------------------------------------------------------------------------
function [logPO,logalfaT]=probsec(A,B,Pi,O)

% Variables of the hmm.
% N: Number of states if the hmm.
[N,M]=size(B{1});
% Variables of the sequences of vectors.
% T: Length of the sequence of states and observations.
[T,M]=size(O{1});
% Matrix of the Forward probability.
alfa=zeros(N,T);
% variable scaled of the alfa.
c=ones(T,1);
% Logarithm of the sequence O probability.
logPO=0.;
% Loop index.
t=0;

% Calcul of the matrices alfa y c.
% Variable for for the calcul fo the alfa
prob=zeros(N,T);
for t=1:T
   prob(:,t)=prodBO(B,O,t);
end

% Calcul of the matrix forward scaled.
alfa(:,1)=Pi.*(prob(:,1));
c(1)=1/sum(alfa(:,1));		% Put c(1)=1 if we do not scale.
alfa(:,1)=alfa(:,1).*c(1);	% Scaled.
for t=2:T,
	alfa(:,t)=(alfa(:,t-1)'*A)'.*(prob(:,t));
	c(t)=1/sum(alfa(:,t));	% Put c(1)=1 if we do not scale
	alfa(:,t)=alfa(:,t).*c(t)+realmin;
end;


% Probability of the sequence PO.
% By definition, for the alfa and beta not scaled:
%	PO=sum(alfa(:,T))=alfa(:,t)'*beta(:,t);
% For alfa and beta scaled, PO is:
%	h=cumprod(c);
%	sum(alfa(:,T))=1=h(T)*sum(alfa_no_escalada(:,T))=h(T)*PO
%	and: PO=1/h(T);
% For both case, we use to calculate the alfa and beta scaled or not:
%	h=cumprod(c);
%	PO=sum(alfa(:,T))/h(T);
% In the case not scaled, we have: PO=sum(alfa(:,T)). And h(T)=1.
% Scaled: PO=1/h(T). And sum(alfa(:,T))=1 is correct.
% NOTA: logPO=log(sum(exp(logalfaT)))

logPO=-sum(log(c))+log(sum(alfa(:,T)));
logalfaT=log(alfa(:,T)')-sum(log(c));

return