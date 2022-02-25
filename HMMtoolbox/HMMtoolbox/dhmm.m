% ------------------------------------------------------------------------------------------------------
% This program calls the diffrent functions to design the HMMs (one for each class and for each group)
% designed in function DHMM_DEF.
% 
% function dhmm(fhmm,fptrain,fptest,fhmmout)
% 
% Entry:  fhmm is the name of the HMMs seted up in the function DHMM_DEF
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
%       Variables defined in the fhmm with the function DHMM_DEF
% 
%     nc: number of classes.
%     ng: number of groups.
%     agrup{ng}: the way to gather the parameters together.
%     Np(ng): number of parameters by group.
%		The variables of the HMM:
%				cell array A=cell(nc,ng);
%				cell array B=cell(nc,ng);
%				cell array Pi=cell(nc,ng);
%		Ne(nc,ng): Number of states of the HMM for each class and group.
%		Ns{ng}(Np): The number of symbols for each parameter of each group.
%       TOPN{ng}(Np): Number of labels for the multilabelling for each parameter of each group
%       Maxitermi is the maximum iteration for the initial model
%		umbral: Umbral is the threshold condition to stop the HMM (the error is calculated with the maximum likelihood criterion)
%		maxiter: Maxiter is the maximum iteration number in the Bakis HMM (condition to stop the algorithm)
%       salto: Salto is the maximum path authorized in the Bakis HMM from a state to another.
%		BAKIS: if Bakis = 1 implements a Bakis HMM (also called left-right), else an ergodic HMM is defined 
%
%	     
%     The libraries are calculated thanks to the set up made in the DHMM_DEF
%            We calculate the library for each group and class
%            The number of centroids is equal to the number of symbols
%               
%           We store the library in the cell array biblio:
%               biblio{group}
%           maxiterVQ is the maximum number of iteration to compute the library (condition to stop the algorithm)
%           umbralVQ is the threshold condition to stop the library computation (condition to stop the algorithm)
%           men=1 or 0 if men =1 the library generates message else no.
%           LBG = 1 the algorithm is LBG, in the other case we choose kMedias
%           dptzoLBG percentage maximum for the distance of the code vector in the algorithm LBG (see gen_bib) 
%    
%           TOPNtest is the number of labels to take into account for the multi-labeling during the test phase
%           Salhmm is a matrix used to store the probabilities values
%
% NOTA: use the scripts dhmm_def, and dhmm_men.
% NOTA: use the functions: etiquetado, iniciaHMM, genhmm, alfabeta, alfa, probsec, viterbi, baum,
%                              gen_bib and kmedias 
% NOTA: if fptrain='' we only realize the test
% NOTA: if fptest='' we only make the training
%------------------------------------------------------------------------------------------------------
function dhmm(fhmm,fptrain,fptest,fhmmout)

% Name of the files ot the training set
%fptrain='vtrain.mat';
% Name of the file of the test set
%fptest='vtest.mat';

% We read the file of set up for the HMMs defined in the DHMM_DEF
eval(['load ',fhmm]);
if eq(nargin,4),fhmm=fhmmout;end
guardar=['save ',fhmm,vDB,vHMM,vVQ,vTEST,' vDB vHMM vVQ vTEST guardar'];

% random variables
randn('state',sum(100*clock));
rand('state',sum(100*clock));

% We check that we have to make or not the training (of the length is 0, we don't train).
if ne(length(fptrain),0)
    % Here we make the training:
    
    % We change the format of the file to make the training
    [lrep,lgrupo]=formato_lectura_secuencial(fptrain,'vtrain',agrup,'vtrainVQ');
    % We generate the library for each parameter of each group
    biblio=gen_bib('vtrainVQ',lgrupo,agrup,Ns,LBG,dpztoLBG,maxiterVQ,umbralVQ,men,biblio);
    eval([guardar])
    
    
    % Once we have the libraries, we are going to label the input parameters and then pass the sequence of labels to train the HMM
    for ig=1:ng,
        for ic=1:nc,
            % We call the function dhmm_men to print some messages if men is true (=1)
            if men,dhmm_men;end

            % We open the files containing the sequences of parameters for the training
            fidvtrain=fopen(['c:\temphmm\vtrain_c',num2str(ic),'_g',num2str(ig),'.tmp'],'r');
            % We keep the sequence of labels in those files:
            fidpl=fopen('c:\temphmm\vpl.tmp','w');
            % It reads the sequence of data for this class and this group
            hora=clock;
            fprintf('Start the function etiquetado to label the parameters. Hour: %g:%g\n',hora(4),hora(5))
            dimrep=agrup{ig}(end)-agrup{ig}(1);
            nrep=size(lrep{ic,ig},1);
            vlt=cell(1,1);
            for ir=1:nrep,
                vlt=fread(fidvtrain,lrep{ic,ig}(ir)*dimrep,'double');
                vlt=reshape(vlt,lrep{ic,ig}(ir),dimrep);               
                % We label each sequence and store it into the temp file corresponding
                pl=etiquetado(vlt,agrup{ig},Ns{ig},biblio{ig},TOPN{ig});
                for ip=1:Np(ig),
                    fwrite(fidpl,pl{ip},'double');
                end
            end
            fclose(fidvtrain);
            fclose(fidpl);
            hora=clock;
            fprintf('End of label of the parameters. Calcul of the initial model. Hour: %g:%g\n',hora(4),hora(5))
            % We calculate the initial model (see iniciaHMM function)
            [A{ic,ig},B{ic,ig},Pi{ic,ig}]=iniciahmm(Ne(ic,ig),Ns{ig},Np(ig),BAKIS,salto,lrep{ic,ig},maxitermi);
            hora=clock;
            fprintf('End of the calcul of the initial model. Start of Baum Welch. Hour: %g:%g\n',hora(4),hora(5))
            % We start to train the HMM in order to initialize the model
            % with the Baum-Welch (EM) method.
            % The training is terminated when the distance between two models calculated is under the threshold
            % or when we have made the maximum iteration authorized.
            iter=0;
            cfe=umbral+1;logPOs=zeros(nrep(1,1),1);
            while and(cfe>umbral,iter<maxiter),            
                iter=iter+1;
                logPOa=logPOs;
                % Baum-Welch algorithm.
                [A{ic,ig},B{ic,ig},Pi{ic,ig},logPOs]=baum(A{ic,ig},B{ic,ig},Pi{ic,ig},lrep{ic,ig});          
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
            if men,fprintf('End of the HMM, group %g, class %g.\n',ig,ic);end;
            % We save a file with all the variables of the HMM
            eval([guardar]);
        end;
        fprintf('End of the training.\n')      
    end
end


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
            fprintf('Calcul of the output of the group: %g, class: %g. Number of the repetition: %g. Hour: %g:%g\n',ig,ic,size(lrep{ic,ig},1),hora(4),hora(5))
            % We create a file to store the output (results) of the HMM:
            fidsalhmm=fopen(['c:\temphmm\salhmm_c',num2str(ic),'_g',num2str(ig),'.tmp'],'w');
            % We open the file containing the sequence of parameters for the test
            fidvtest=fopen(['c:\temphmm\vtest_c',num2str(ic),'_g',num2str(ig),'.tmp'],'r');
            % We make the test
            dimrep=agrup{ig}(end)-agrup{ig}(1);
            nrep=size(lrep{ic,ig},1);

            for ir=1:nrep,
                %hora=clock;
                %fprintf('\tCalcul of the output of the HMM: %g, class: %g. Repetition: %g. Hour: %g:%g\n',ig,ic,ir,hora(4),hora(5))
                % We read the repetition
                vlt=fread(fidvtest,lrep{ic,ig}(ir)*dimrep,'double');
                vlt=reshape(vlt,lrep{ic,ig}(ir),dimrep);
                % It labels this repetition
                pl=etiquetado(vlt,agrup{ig},Ns{ig},biblio{ig},TOPNtest{ig});
                % Evaluate the probability for each label.
                salida=zeros(nc,1);
                for ihmm=1:nc;
                    salida(ihmm,1)=probsec(A{ihmm,ig},B{ihmm,ig},Pi{ihmm,ig},pl);
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
            nrep=size(lrep{ic,ig},1);

            fid=fopen(['c:\temphmm\salhmm_c',num2str(ic),'_g',num2str(ig),'.tmp'],'r');
            salhmm{ic,ig}=cell(nrep,1);
            for ir=1:nrep
                salhmm{ic,ig}{ir}=fread(fid,nc,'double');
            end
            fclose(fid);
        end
    end
    % We save the outputs of the HMM.
    eval([guardar]);
    fprintf('End of test.\n');

end
return
