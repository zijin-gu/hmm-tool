%------------------------------------------------------------
%
%
%           %%%%%%%%%%%
%           %DHMM_MEN %
%           %%%%%%%%%%%
% 
% 
% 
% 
% 
% This script is a small program that prints some messages 
% to describe the computation of the HMM.
% 
%-------------------------------------------------------------

fprintf('\n_______________________________________________________\n\n');
fprintf('Generation of the HMM for the class %g, group (gathering) %g of the file %s\n',ic,ig,fptrain)
fprintf('Number of states of the model: %g\n',Ne(ic,ig));
fprintf('Number of parameters: %g\n',Np(ig));
fprintf('Size of the vectors for each parameter: %g\n',diff(agrup{ig}))
fprintf('Number of labels per symbol (Multilabeling)/parameter: %g\n',TOPN{ig});
fprintf('maximum iterations, threshold of stop: %g %g\n',maxiter,umbral);
fprintf('Length of the training set: %g\n',length(lrep{ic,ig}));

