%------------------------------------------------------------------------
% 
%                      %%%%%%%%%%%%
%                      % Viterbic %
%                      %%%%%%%%%%%%
%
% This function calculates the sequence of the more probable states given the HMM and the sequence monitored O. We use the algorithm of Viterbi with the log for the numerical precision problems.
% 
% 	function qP=viterbic(A,B,Med,Var,Pi,O,agrup)
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B{Np}{N}(Ngauss{ip},1) is a structure with the weight of each Gaussian for each state and to obtain each parameter.
% 
% Med: Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the mean for each Gaussian in each state and for each gathering of parameters
% 
% Var: Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the variance for each Gaussian in each state and for each gathering of parameters
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% 
% O is a matrix with the sequence to estimate.
% 
% agrup(Np+1,1) defines how to cluster the parameters together
% 
% 
% Result:	qP(T,1) is the sequence of the most probable states.
% 
%--------------------------------------------------------------------------
function qP=viterbic(A,B,Med,Var,Pi,O,agrup)


% Variables of the hmm.
% Ne: Number of states of the HMM
Ne=size(B{1},1);

% Variables of the sequence of vectors.
T: length of the sequence of states and observations.
% Dim: size of O.
T=size(O,1);

% Most probable state sequence.
qP=zeros(T,1);
% Variable delta: delta(i,t).
% Highest probability to be in the state i at the instant t
% given the sequence of observation O.
delta=zeros(Ne,T);
% Variable Psi. Psi(i,t).
% Most probable sequence of states to be in the state i at the instant t
% given the sequence of observation O.
Psi=zeros(Ne,T);
% Loop index.
t=0;j=0;

% % Log of the hmm.
A=log(A+realmin);
Pi=log(Pi+realmin);
% Calcul of the delta matrix.
delta(:,1)=Pi+log(probsimb(B,Med,Var,O(1,:),agrup));
Psi(:,1)=zeros(Ne,1);
for t=2:T,
	for j=1:Ne,
		[delta(j,t),Psi(j,t)]=max(delta(:,t-1)+A(:,j));
	end;
	delta(:,t)=delta(:,t)+log(probsimb(B,Med,Var,O(t,:),agrup));
end;
% Calcul of the most probable state sequence.
[aux,qP(T)]=max(delta(:,T));
for t=T-1:-1:1,
	qP(t)=Psi(qP(t+1),t+1);
end;
return