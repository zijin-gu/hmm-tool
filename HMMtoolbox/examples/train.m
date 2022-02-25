
%------------------------------------------------------------------------------------
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          %   This script separates the database between training   %
%          % and test. It also calls the function DHMM_DEF and DHMM. %
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------------------

% WARNING: before runing create c:\temphmm and add HMMtoolbox to Matlab path
% For Continuous HMM thye NETLAB toolbox is required

%create example
generatedatabasetest

% train HMM
% load database
load parampoligonos vlcp

% we separate the database in training and test.
[nc,ng]=size(vlcp);
ng=size(vlcp,2);
nr=size(vlcp{1,1},1);
% The training rate
ptrain=50;
nrt=ceil(nr*ptrain/100);
nrtest=nr-nrt;
vtrain=cell(nc,ng);
vtest=cell(nc,ng);
for ic=1:nc
    for ig=1:ng
        vtrain{ic,ig}=cell(nrt,1);
        vtest{ic,ig}=cell(nr-nrt,1);
    end
end
ind=randperm(nr);
indtr=ind(1:nrt);
indts=ind(nrt+1:nr);
for ic=1:nc
    for ir=1:nrt
        for ig=1:ng
            vtrain{ic,ig}{ir}=vlcp{ic,ig}{indtr(ir)};
        end
    end
    for ir=1:nrtest
        for ig=1:ng
            vtest{ic,ig}{ir}=vlcp{ic,ig}{indts(ir)};
        end
    end
end
vl=vtrain;
save vtrain vl
vl=vtest;
save vtest vl
clear vtest vtrain vl vlcp;


% define dhmm morphology
%chmm_def('hmmpoligonos2.mat');
dhmm_def('hmmpoligonos.mat');

% train and test DHMM
dhmm('hmmpoligonos.mat','vtrain','vtest');

%Train and test CHMM. 
%dhmm2chmm('hmmpoligonos.mat','hmmpoligonos2.mat');
%chmm('hmmpoligonos2.mat','vtrain','vtest');

%result
Mc=resulhmm('hmmpoligonos');

