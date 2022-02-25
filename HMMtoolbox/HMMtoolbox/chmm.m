% ------------------------------------------------------------------------------------------------------
% This program calls the diffrent functions to design the HMMs (one for each class and for each group)
% designed in function CHMM_DEF.
% 
% function chmm(fhmm,fptrain,fptest,fhmmout)
% 
% Entry:  fhmm is the name of the HMMs set up in the function CHMM_DEF
%         
%         fptrain is the name of the file containing the sequences of parameters to train each HMM 
%         with its own group of training set. In each repetition, we find the a sequence of parameters
%         for a class and a group.
% 
%         fptest is the name of the file containing the sequences of parameters to test each HMM 
%         with its own group of test set. In each repetition, we find the a sequence of parameters
%         for a class and a group.
% 
%         fsalhmm: Name of the file containing the outputs of each classifier for each sample of the database of the test 
%		  in a cell array salhmm{class, group}(repetition).
%        
%       Variables defined in the fhmm with the function CHMM_DEF
% 
%     nc: number of classes.
%     ng: number of groups.
%     agrup{ng}: the way to gather the parameters together.
%     Np(ng): number of parameters by group.
%		The variables of the HMM:
%				cell array A=cell(nc,ng);
%				cell array B=cell(nc,ng);
% 				cell array Med=cell(nc,ng);
% 				cell array Var=cell(nc,ng);
%				cell array Pi=cell(nc,ng);
%		Ne(nc,ng): Number of states of the HMM for each class and group.
%       Maxitermi is the maximum iteration for the initial model
%		umbral: Umbral is the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion)
%		maxiter: Maxiter is the maximum iteration number in the Bakis HMM (condition to stop the algorithm)
%       salto: Salto is the maximum path authorized in the Bakis HMM from a state to another.
%		BAKIS: if Bakis = 1 implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
%
%           Salhmm is a matrix used to store the probabilities values
%
% NOTA: use the scripts chmm_def, and chmm_men.
% NOTA: use the functions: iniciacHMM, genchmm, alfabetac,  probsecC, viterbic, baumc
% NOTA: if fptrain='' we only realize the test
% NOTA: if fptest='' we only make the training
%------------------------------------------------------------------------------------------------------
function chmm(fhmm,fptrain,fptest,fhmmout)

% Name of the files ot the training set
%fptrain='vtrain.mat';
% Name of the file of the test set
%fptest='vtest.mat';

% We read the file of set up for the HMMs defined in the DHMM_DEF

eval(['load ',fhmm]);
if eq(nargin,4),fhmm=fhmmout;end
guardar=['save ',fhmm,vDB,vHMM,vTEST,' iniciar vDB vHMM vTEST guardar'];

% random variables
randn('state',sum(100*clock));
rand('state',sum(100*clock));

% We check that we have to make or not the training (of the length is 0, we don't train).
if ne(length(fptrain),0)
    % Here we make the training:
    
    % We change the format of the file to make the training
    [lrep]=formato_lectura_secuencial(fptrain,'vtrain');
    

   
   % We calculate the classifiers
   for ig=1:ng  
         for ic=1:nc, 
            % Messages for the design of the classifier
            if men
              chmm_men
          end
            nftrain=['c:\temphmm\vtrain_c',num2str(ic),'_g',num2str(ig),'.tmp'];         
            % we calculate the model initial of the HMM for each class and group
            if iniciar
                hora=clock;
                fprintf('Calcul of the initial model. Hour: %g:%g\n',hora(4),hora(5))
                [A{ic,ig},B{ic,ig},Med{ic,ig},Var{ic,ig},Pi{ic,ig}]=iniciachmm(Ne(ic,ig),Np(ig),BAKIS,salto,nftrain,lrep{ic,ig},agrup{ig},Ngauss{ig},maxitermi);
            % We start to train the HMM in order to initialize the model
            % with the Baum-Welch (EM) method.
            % The training is terminated when the distance between two models calculated is under the threshold
            % or when we have made the maximum iteration authorized.
                hora=clock;
                fprintf('End of the calcul of the initial model.  Hour: %g:%g\n',hora(4),hora(5))
            else
                hora=clock;
                fprintf('Start of Baum Welch. Hour: %g:%g\n',hora(4),hora(5))
            end
            
            iter=0;
            cfe=umbral+1;logPOs=zeros(size(lrep{ic,ig},1),1);
            while and(abs(cfe)>umbral,iter<maxiter),            
               
               iter=iter+1;
               logPOa=logPOs;
               % Baum-Welch algorithm.
               [A{ic,ig},B{ic,ig},Med{ic,ig},Var{ic,ig},Pi{ic,ig},logPOs]=baumc(A{ic,ig},B{ic,ig},Med{ic,ig},Var{ic,ig},Pi{ic,ig},nftrain,lrep{ic,ig},agrup{ig});
               % Criterion to stop.
               if eq(iter,1),
                   logPObase=mean(logPOs);
                   hora=clock;
                   fprintf('iter. Baum: %g,\tProb.: %g\tStd: %g\tHour: %g:%g\n',iter,mean(logPOs),std(logPOs),hora(4),hora(5));
               else;
                   cfe=1-(mean(logPOa-logPObase)/mean(logPOs-logPObase));
                   hora=clock;
                   fprintf('iter. Baum: %g,\tProb.: %g\tStd: %g\tCfe: %g\tHour: %g:%g\n',iter,mean(logPOs),std(logPOs),cfe,hora(4),hora(5));
               end
           end;
           fprintf('End of the HMM, group %g, class %g.\n',ig,ic);
          % We save a file with all the variables of the HMM
            eval([guardar]);
         end;
      fprintf('End of the training.\n')   
   end   
end
% We save in a file all the varaiables
eval([guardar]);

% We check if we have to make or no the test
if ne(length(fptest),0)
    % If the length>0 we make it:
    % We convert the test set into the adequat format for the HMM
    [lrep]=formato_lectura_secuencial(fptest,'vtest');

    % Messages :
    if men,
        fprintf('Input parameters of the HMM: %s\n',fptest);
        fprintf('Number of the class: %g\n',nc);
        fprintf('Number of the group: %g\n',ng);
    end;
    for ig=1:ng,
        for ic=1:nc,
            hora=clock;
            fprintf('Calcul of the output of the group:: %g, class: %g. Number of the repetition: %g. Hour: %g:%g\n',ig,ic,size(lrep{ic,ig},1),hora(4),hora(5))
            % We create a file to store the output (results) of the HMM:
            fidsalhmm=fopen(['c:\temphmm\salhmm_c',num2str(ic),'_g',num2str(ig),'.tmp'],'w');
           % We open the file containing the sequence of parameters for the test
            fidvtest=fopen(['c:\temphmm\vtest_c',num2str(ic),'_g',num2str(ig),'.tmp'],'r');               
            % we make the test
            dimrep=agrup{ig}(end)-agrup{ig}(1);
            nrep=size(lrep{ic,ig},1);

            for ir=1:nrep,
                %hora=clock;
                %fprintf('\tCalcul of the output of the HMM: %g, class: %g. Repetition: %g. Hour: %g:%g\n',ig,ic,ir,hora(4),hora(5))
                % We read the repetition
                vlt{1}=fread(fidvtest,lrep{ic,ig}(ir)*dimrep,'double');
                vlt{1}=reshape(vlt{1},lrep{ic,ig}(ir),dimrep);
                % We evaluate the probability for each model
                salida=zeros(nc,1);
                for ihmm=1:nc;
                    salida(ihmm,1)=probsecc(A{ihmm,ig},B{ihmm,ig},Med{ihmm,ig},Var{ihmm,ig},Pi{ihmm,ig},vlt{1},agrup{ig});
                end
                fwrite(fidsalhmm,salida,'double');
            end
            fclose(fidsalhmm);
            fclose(fidvtest);
        end;
    end;
     
    % Join the outputs of the HMM salhmm to keep all together
    fprintf('Union of the outputs salhmmm.\n');
    for ig=1:ng
        for ic=1:nc
            fid=fopen(['c:\temphmm\salhmm_c',num2str(ic),'_g',num2str(ig),'.tmp'],'r');
            nrep=size(lrep{ic,ig},1);
            salhmm{ic,ig}=cell(nrep,1);
            for ir=1:nrep
                salhmm{ic,ig}{ir}=fread(fid,nc,'double');
            end
            fclose(fid);
        end
    end

      
   fprintf('FIN test.\n');
    % We save the outputs of the HMM.
    eval([guardar]);
end
return
