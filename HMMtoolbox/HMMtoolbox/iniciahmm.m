%---------------------------------------------------
% 
% 
%               %%%%%%%%%%%%%
%           	% iniciaHMM %
%               %%%%%%%%%%%%%    
%     
% 
% This function calculates the HMM optimal initial model for the Baum-Welch training.
% 
% 
% Entry: 	if Bakis = 1 implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
% Salto is the maximum path authorized in the Bakis HMM from a state to another.
% Maxitermi is the maximum iteration to calculate the initial HMM.
% Ne state number 
% Ns symbol number for every state
% Np number parameters by symbol
% Lrep is a Matrix with the size from the training set data for every classes groups and parameters (see also the function formato_lectura_secuencial)
% 
% Result: 	Aa Ba and Pia are the matrices for the parameters of the HMM.
% 
% In this function, we will generate a number of symbols in every state bigger than twice Ns. Thus, we make sure to initialize correctly the HMM.
% 		
% We use the viterbi and the genhmm functions.
% 
% 
%---------------------------------------------------
function [Aa,Ba,Pia]=iniciahmm(Ne,Ns,Np,BAKIS,salto,lrep,maxitermi)

% Ir reads the database for the training of this class and group 
% to initialize the initial model
nr=length(lrep);
Ltd=sum(lrep);
labels=cell(Np,1);
OqP=cell(nr,1);
for ir=1:nr
    OqP{ir}=zeros(lrep(ir),1);
end


% It initializes the Markov models.
itermi=1;
% Initialization of the matrices of the Markov model
[A,B,Pi]=genhmm(Ne,Ns,Np,BAKIS,salto);
Aa=A;Ba=B;Pia=Pi;
if gt(maxitermi,1)
   % To check that initialization is good, the critarion
   % is that all the states have more than twice symbols 
   % generated than Ns.
   % Calcul of the number of symbols generated by each state of the HMM.
   fidpl=fopen('c:\temphmm\vpl.tmp','r');
   for ir=1:nr,
       for ip=1:Np
           labels{ip}=fread(fidpl,lrep(ir)*Ns(ip),'double');
           labels{ip}=reshape(labels{ip},lrep(ir),Ns(ip));
       end
       OqP{ir}=viterbi(A,B,Pi,labels);
   end;
   fclose(fidpl);
   % And we calculate the largest distance to the equi-occupation
   dmim=max(abs(Ltd/Ne-hist(cat(1,OqP{:}),[1:Ne])));
    % If all the states do not have the number bigger than
   %  the minimum, we reionitiate the model.
   % Nota Bene: we can be blocked in an infinite loop
   while and(dmim>0.3*Ltd/Ne,itermi<maxitermi),
      itermi=itermi+1;
      [A,B,Pi]=genhmm(Ne,Ns,Np,BAKIS,salto);
      fidpl=fopen('c:\temphmm\vpl.tmp','r');
      for ir=1:nr,
       for ip=1:Np
           labels{ip}=fread(fidpl,lrep(ir)*Ns(ip),'double');
           labels{ip}=reshape(labels{ip},lrep(ir),Ns(ip));
       end
          OqP{ir}=viterbi(A,B,Pi,labels);
      end;
      fclose(fidpl);
      dmi=max(abs(Ltd/Ne-hist(cat(1,OqP{:}),[1:Ne])));
      if dmi<dmim,
         dmim=dmi;
         Aa=A;Ba=B;Pia=Pi;
      end;
   end;
end
return