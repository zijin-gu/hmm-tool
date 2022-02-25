%
%                           %%%%%%%%%%% 	                    
%                           %  BAUMC  %
%                           %%%%%%%%%%%
%
%
% This function computes the Baum Welch algorithm in order to estimate the HMM parameters.
% function [A1,B1,Med1,Var1,Pi1,logPOs]=baumc(A,B,Med,Var,Pi,nftrain,lrep,agrup)
% 
% A, B and Pi are the matrices of the HMM defined with the function genHMM for example.
% 
% Entry:	A(N, N) is the transition probability matrix from a state I to a state k.
% 
% B{Np}{N}(Ngauss{ip},1) is the weight for the Gaussians of the HMM.
% 
% Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the medians for the Gaussians of the HMM.
% 
% Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the variances for the Gaussians of the HMM.
% 
% Pi(i) is the probability to start in the state i.
% 
% Lrep is a matrix with the vectors number and symbols number.
% 
% A1(N, N) is the transition probability matrix updated with the Baum Welch algorithm.
% B1{Np}{N}(Ngauss{ip},1) is the distribution probability symbol matrix updated with the Baum Welch algorithm.
% Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the medians for the Gaussians of the HMM updated.
% Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)): Matrix with the values of the variances for the Gaussians of the HMM updated.
% Pi1(N,1) is the distribution probability for the initial state updated.
% logPOs(nr,1): probability that the HMM with the parameters (A,B,Pi) generates each realisation of the labels entries of the subroutine. 
% realisation of the labels entries of the subroutine. 
% 
% This function makes a call to the probsimb function.
%
%
%--------------------------------------------------------



function [A1,B1,Med1,Var1,Pi1,logPOs]=baumc(A,B,Med,Var,Pi,nftrain,lrep,agrup)

% Variables for the HMM
% Np: Parameters to calculate.
% N: States number.
Np=size(B,1);
Ngauss=zeros(Np,1);
for ip=1:Np
    Ngauss(ip)=size(B{ip}{1},1);
end
N=size(B{1},1);

% nr: Number of repetition (Number of sequence observed for the training)
nr=size(lrep,1);

% Matrices for the HMM reestimate.
A1=zeros(N,N);Pi1=zeros(N,1);
B1=cell(Np,1);
for ip=1:Np
    B1{ip}=cell(N,1);
    for ie=1:N
        B1{ip}{ie}=zeros(Ngauss(ip),agrup(ip+1)-agrup(ip));
    end
end
% Matrices used to calculate A1.
numA1=zeros(N,N);denA1=zeros(N,1);
% Matrices used to calculate B1.
numB1=cell(Np,1);denB1=cell(Np,1);
numMed1=cell(Np,1);numVar1=cell(Np,1);denMedVar=cell(Np,1);
for ip=1:Np,
    numB1{ip}=cell(N,1);
    denB1{ip}=cell(N,1);
    numMed1{ip}=cell(N,1);
    numVar1{ip}=cell(N,1);
    denMedVar{ip}=cell(N,1);
    for ie=1:N
        numB1{ip}{ie}=zeros(Ngauss(ip),1);
        denB1{ip}{ie}=0;
        numMed1{ip}{ie}=zeros(Ngauss(ip),agrup(ip+1)-agrup(ip));
        numVar1{ip}{ie}=zeros(Ngauss(ip),agrup(ip+1)-agrup(ip));
        denMedVar{ip}{ie}=zeros(Ngauss(ip),1);
    end
end
% Variable intermediate to calculate the probability.
POl=0.0;
% Vector of probability
logPOs=zeros(nr,1);

% Calcul of the improvement of each sequence of Observation
% from the training set.

fidvtrain=fopen(nftrain,'r');
for ir=1:nr,
    % the vector is readen
    dimrep=agrup(end)-agrup(1);
    vlt=fread(fidvtrain,lrep(ir)*dimrep,'double');
    vlt=reshape(vlt,lrep(ir),dimrep);
    
    
    % We recalculate the matrices that depends on T
    
    T=lrep(ir);
    % Forward probability Matrix.
    alfa=zeros(N,T);
    % Backward probability Matrix.
    beta=zeros(N,T);
    % Variable gama.
    % Gama(i,t) is the probability to be in the state i at the instant t
    % given the observations sequence O.
    % sum(gama(:,t))=1 for any t..
    gama=zeros(N,T);
    % Variable scaled for the alfa and beta.
    c=ones(T,1); 
    
    % We calculate the probabilities for each symbols of the sequence 
    % in each state and for each gaussian.
    prob=zeros(N,T);
    Probg=cell(T,1);
    for t=1:T
        [prob(:,t),probg{t}]=probsimb(B,Med,Var,vlt(t,:),agrup);
    end
    
    % calcul of the alfa and beta scaled.
    % The algorithms implemented can be used either for the variables scaled
    % or for the non scaled.
    % Intermediate variables for the alfa and beta
    % Forward scaled matrix
    alfa(:,1)=Pi.*(prob(:,1));
    c(1)=1/sum(alfa(:,1));		% c(1)=1 to avoid the scaling.
    alfa(:,1)=alfa(:,1).*c(1);	% alfa scaled.
    for t=2:T,
        alfa(:,t)=(alfa(:,t-1)'*A)'.*(prob(:,t));
        c(t)=1/sum(alfa(:,t));	% c(t)=1 to avoid the scaling.
        alfa(:,t)=alfa(:,t).*c(t)+realmin;
    end;
    % Backward scaled matrix.
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
    
    for t=1:T-1,
        for i=1:N,
            numA1(i,:)=numA1(i,:)+alfa(i,t)*(A(i,:).*prob(:,t+1)'.*beta(:,t+1)');
        end;
    end;
    denA1=denA1+sum(gama(:,1:T-1)')';
    
    % we calculate Psit to reestimate B, Mean and Variance.
    % Psi is the time sum of the Psit
    for ip=1:Np
        Psi=zeros(Ngauss(ip),1);
        Psit=zeros(Ngauss(ip),T);
        dimp=[agrup(ip):agrup(ip+1)-1];
        onesg=ones(Ngauss(ip),1);onesp=ones(1,length(dimp));
        for i=1:N,
            Psit(:,1)=gama(i,1).*probg{1}{ip}{i};
            Psi=Psit(:,1);
            for t=2:T
                Psit(:,t)=gama(i,t).*probg{t}{ip}{i};
                Psi=Psi+Psit(:,t);
            end
            % we calculate numB1 and denB1
            numB1{ip}{i}=numB1{ip}{i}+Psi;
            denB1{ip}{i}=denB1{ip}{i}+sum(Psi);         
            % we calculate the numerator of the mean and the variance
            % and the denominator of the mean and the variance
            for t=1:T,
                numMed1{ip}{i}=numMed1{ip}{i}+Psit(:,t)*vlt(t,dimp);
                cteA=onesg*vlt(t,dimp)-Med{ip}{i};
                numVar1{ip}{i}=numVar1{ip}{i}+(Psit(:,t)*onesp).*(cteA.*cteA);
            end
            denMedVar{ip}{i}=denMedVar{ip}{i}+Psi;      
        end
    end
    
    % We calculate the part from the distribution of the initial states
    % due to the training sequence l
    Pi1=Pi1+gama(:,1);
end;

% Reestimation of the transition states matrix, of the
% distribution probability of symbol, and of the distribution
% of the initials states.
A1=A;B1=B;Med1=Med;Var1=Var;
A1=numA1./(denA1*ones(1,N));
for ip=1:Np,
    for i=1:N,
        B1{ip}{i}=numB1{ip}{i}./denB1{ip}{i}+1e-6;
        Med1{ip}{i}=numMed1{ip}{i}./(denMedVar{ip}{i}*ones(1,agrup(ip+1)-agrup(ip)));
        Var1{ip}{i}=sqrt(numVar1{ip}{i}./(denMedVar{ip}{i}*ones(1,agrup(ip+1)-agrup(ip))))+1e-6;
    end
end
Pi1=Pi1/nr;
fclose(fidvtrain);
return
