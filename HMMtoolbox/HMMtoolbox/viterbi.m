%------------------------------------------------------------------------
% 
%                      %%%%%%%%%%%
%                      % Viterbi %
%                      %%%%%%%%%%%
% 
% 
% This function calculates the sequence of the more probable states given the HMM and the sequence monitored O. We use the algorithm of Viterbi with the log for the numerical precision problems.
% 
% 	function qP=viterbi(A,B,Pi,O)
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
% Result:	qP(T,1) is the sequence of the more probable states.

%--------------------------------------------------------------------------
function qP=viterbi(A,B,Pi,O)

% Variables of the HMM:
N=size(A,1);
% Variables of the sequence of vectors.
% T: length of the sequence of states and observations.
% Dim: size of O.
T=size(O{1},1);

% Most probable state sequence.
qP=zeros(T,1);
% Variable delta: delta(i,t).
% Highest probability to be in the state i at the instant t
% given the sequence of observation O.
delta=zeros(N,T);
% Variable Psi. Psi(i,t).
% Most probable sequence of states to be in the state i at the instant t
% given the sequence of observation O.
Psi=zeros(N,T);
% Loop index.
t=0;j=0;

% Log of the hmm.
A=log(A+realmin);
Pi=log(Pi+realmin);
% Calcul of the delta matrix.
delta(:,1)=Pi+log(prodBO(B,O,1));
Psi(:,1)=zeros(N,1);
for t=2:T,
	for j=1:N,
		[delta(j,t),Psi(j,t)]=max(delta(:,t-1)+A(:,j));
	end;
	delta(:,t)=delta(:,t)+log(prodBO(B,O,t));
end;
% Calcul of the most probable state sequence.
[aux,qP(T)]=max(delta(:,T));
for t=T-1:-1:1,
	qP(t)=Psi(qP(t+1),t+1);
end;
return