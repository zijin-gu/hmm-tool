%------------------------------------------------------------------------
% 
%                     %%%%%%%%   
%                     % Baum %
%                     %%%%%%%%
% 
%
%
% This function computes the Baum Welch algorithm in order to estimate the HMM parameters.
% function [A1,B1,Pi1,logPOs]=baum(A,B,Pi,lrep)
% 
% A, B and Pi are the matrices of the HMM defined with the function genHMM for example.
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B(N,M) is the distribution probability symbol matrix 
% B(i,k) is the probability to obtain the kth  symbol when we are at the state i.
% 
% Pi(N,1) is the distribution probability for the initial state i.
% Pi(i) is the probability to start in the state i.
% 
% Lrep is a matrix with the vectors number and symbols number.
% 
% A1(N, N) is the transition probability matrix updated with the Baum Welch algorithm.
% B1(N, N) is the distribution probability symbol matrix updated with the Baum Welch algorithm.
% Pi(N,1) is the distribution probability for the initial state updated.
% logPOs(nr,1): probability that the HMM with the parameters (A,B,Pi) generates each realisation of the labels entries of the subroutine. 
% 
% This function makes a call to the alfabeta function.
% 
%--------------------------------------------------------------------------

function [A1,B1,Pi1,logPOs]=baum(A,B,Pi,lrep)

%   Variables for the HMM
% Np: Parameters to calculate.
% N: States number.
% M(ip): Number of synbols for each parameter
Np=size(B,1);
M=zeros(Np,1);
for ip=1:Np
   [N,M(ip)]=size(B{ip});
end
% nr: Number of repetition (Number of sequence observed for the training)
nr=length(lrep);

% Matrices of the HMM reestimated.
A1=zeros(N,N);Pi1=zeros(N,1);
B1=cell(Np,1);
for ip=1:Np
   B1{ip}=zeros(N,M(ip));
end
% Matrices used to calculate A1.
numA1=zeros(N,N);denA1=zeros(N,1);Eps=zeros(N,N);
% Matrices used to calculate B1.
numB1=cell(Np,1);denB1=cell(Np,1);Psi=cell(Np,1);
for ip=1:Np,
   numB1{ip}=zeros(N,M(ip));denB1{ip}=zeros(N,1);Psi{ip}=zeros(N,M(ip));
end
% Variable para calcular the probability.
POl=0.0;
% Vector of probability
logPOs=zeros(nr,1);
% Loop index.
i=0;k=0;t=0;l=0;

% Calcul of the improvement of each sequence of Observation
% from the training set.

fidpl=fopen('c:\temphmm\vpl.tmp','r');
for ir=1:nr,
	 % We recalculate the matrices that depends on T
    T=lrep(ir);
	% Forward probability matrix.
	alfa=zeros(N,T);
	% Backward probability matrix.
	beta=zeros(N,T);
	
	% Variable gama.
    % Gama(i,t) is the probability to be in the state i at the instant t
    % given the observations sequence O.
    % sum(gama(:,t))=1 for any t...
    
	gama=zeros(N,T);
	% Matrices intermediate to compute B1.
	punt=zeros(T,1);unos=ones(T,1);
	% Variables scaled for alfa and beta.
	c=ones(T,1);

	% calcul of the alfa and beta scaled.
   % The following algorithm can be used as well with the scaled alfa and beta 
   % as with the no scaled.
   
   % we read the labels for the repetition ir
   labels=cell(Np,1);
   for ip=1:Np
       labels{ip}=fread(fidpl,lrep(ir)*M(ip),'double');
       labels{ip}=reshape(labels{ip},lrep(ir),M(ip));
   end
   
   % The probability of the sequence is calculated in each state
   prob=zeros(N,T);
   for t=1:T
      prob(:,t)=prodBO(B,labels,t);
   end
   
   % Calcul of the Forward probability Matrix.
   alfa(:,1)=Pi.*(prob(:,1));
   c(1)=1/sum(alfa(:,1));		% Put c(1)=1 to not scale.
   alfa(:,1)=alfa(:,1).*c(1);	% Scaled.
   for t=2:T,
      alfa(:,t)=(alfa(:,t-1)'*A)'.*(prob(:,t));
      c(t)=1/sum(alfa(:,t));	% Pur c(t)=1 to not scale.
      alfa(:,t)=alfa(:,t).*c(t)+realmin;
   end;
   
   % Calcul of the Matrix Backward probability.
   % Put c=ones(T,1) to not scale.
   beta(:,T)=ones(N,1)*c(T);
   for t=T-1:-1:1,
      beta(:,t)=A*(prob(:,t+1).*beta(:,t+1))*c(t)+realmin;
   end;
   POl=sum(alfa(:,T));
   logPOs(ir)=log(POl)-sum(log(c));
	
	% Matrix gama.
	for t=1:T,
		gama(:,t)=(alfa(:,t).*beta(:,t))./(alfa(:,t)'*beta(:,t))+realmin;
	end;

	 % we calculate Psit to reestimate B, Mean and Variance.
    % Psi is the time sum of the Psit
	Eps=zeros(N,N);
	for t=1:T-1,
		for i=1:N,
			Eps(i,:)=Eps(i,:)+alfa(i,t)*(A(i,:).*(prob(:,t+1))'.*beta(:,t+1)');
		end;
	end;
	numA1=numA1+Eps/POl;
	denA1=denA1+sum(gama(:,1:T-1)')'/POl;

	% Calcul of the part of the numerator and denominator for the distribution matrix
	% of the symbols corresponding to the sequence of observation of the training.
	% Psi is initialized with very small values because, if there is a symbol that never occurs, the probability
    % would be small enough but not a zero (error would occur)
    %   Psi=eps.*ones(N,M);
   for ip=1:Np,
      Psi=zeros(N,M(ip));
      for i=1:N
         for k=1:M(ip),
            % Reestestimete B with the heuristic probability.
            Psi(i,k)=Psi(i,k)+gama(i,:)*labels{ip}(:,k);
            % Reestimate B with contribution normalized.
%            for t=0:T-1,
%               Psi(i,k)=Psi(i,k)+gama(i,t+1)*labels{ir}{ip}(t+1,k)*B{ip}(i,k)/sum(B{ip}(i,:).*labels{ir}{ip}(t+1,:));
%            end;
         end;
      end;
      numB1{ip}=numB1{ip}+Psi/POl;
      denB1{ip}=denB1{ip}+sum(gama')'/POl;
   end;
   	
   % Calcul of the contribution to the distribution of the initials
	% state corresponding to the observation sequence 1 of the training	
	Pi1=Pi1+gama(:,1);
end;
fclose(fidpl);
% Reestimation of the states transition matrix, of the
% symbol probability distribution, and the distribution
% of the initials states.
A1=numA1./(denA1*ones(1,N));
for ip=1:Np,
   B1{ip}=numB1{ip}./(denB1{ip}*ones(1,M(ip)))+1e-6;
end
Pi1=Pi1/nr;
%fclose(fidpl);
return
