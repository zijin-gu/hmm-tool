% ----------------------------------------------------------------
% 
%                      %%%%%%%
%                      % ROC %
%                      %%%%%%%
%
% 
% This functions analyses the results from the HMM. 
% In particular it deals with the problem of 
% false acceptance rate (FAR called here FMR false match rate) and 
% false rejection rate (FRR called here FNMR false non match rate).
% 
% NOTA: Use probsec and alfabeta
%---------------------------------------------------------------------------------
function [Fa,Fr,Umbral]=ROC(fhmm)



% Read the HMM, database and libraries
eval(['load ',fhmm])
% Number of iteration
nr=zeros(nc,1);
for ic=1:nc
   nr(ic)=size(salhmm{ic,1},1);
end
% Vectors with the output probability
logPO=zeros(nc,1);
%  genuine distance and impostor distance of the vector:
DGEN=zeros(sum(nr),1);
DIMP=zeros(sum(nr)*(nc-1),1);
% Number of points to calculate
npe=100;
% Definition of the matrices of false acceptance rate (FAR called here FMR false match rate) 
% and false rejection rate (FRR called here FNMR false non match rate)for the groups 2 to 2, 3 to 3 ... ng to ng
% example: FMR{2}, is the probability to have the false acceptance of the classifier
% where the outputs are combined 2 to 2
% and FMR{2}{4} is the probability to have the false acceptance of the classifier
% with the output of the group G(4) with G=nchoosek([1:ng],2);
FMR=cell(ng,1);
FNMR=cell(ng,1);
for ig=1:ng
   ncomb=size(nchoosek([1:ng],ig),1);
   FMR{ig}=cell(ncomb,1);
   FNMR{ig}=cell(ncomb,1);
   for icomb=1:ncomb
      FMR{ig}{icomb}=zeros(npe,1);
      FNMR{ig}{icomb}=zeros(npe,1);
   end 
end

% Messages of thhe analysis.
fprintf('Number of classes: %g\n',nc);
fprintf('Number of groups: %g\n',ng);

% The output of the classifiers are normalized [0 1]
% in an independent way
for ic=1:nc
   for ig=1:ng,
      for ir=1:nr(ic),
         margen=diff(minmax(salhmm{ic,ig}{ir}')');
         minimo=min(salhmm{ic,ig}{ir});
         salhmm{ic,ig}{ir}=(salhmm{ic,ig}{ir}-ones(nc,1)*minimo)./(ones(nc,1)*margen);
      end
   end
end

% The ROC curve of each combination of groups is calculated
% icc: the classifiers of the groups are taken into account as 1, 2, icc, until taken all together.
for icsc=1:ng,
   Mcomb=nchoosek([1:ng],icsc);
   ncomb=size(Mcomb,1);
  
   for icomb=1:ncomb,    
      %scan all the classes and repetitions
      DGEN=[];DIMP=[];
      for ic=1:nc,
         for ir=1:nr(ic),
            % Read all the outputs of the classifier selected and make the mean
            logPO=zeros(nc,1);
            for isc=1:icsc
               logPO=logPO+salhmm{ic,Mcomb(icomb,isc)}{ir}/icsc;
            end
            logPO=1-logPO;
            % Geniune distance and impostor distance
            DGEN=[DGEN; logPO(ic)];
            DIMP=[DIMP; logPO(1:ic-1); logPO(ic+1:nc)];
         end;
      end
      % Distribution of DGEN and DIMP
      tao=linspace(0,1,100);
      FGEN=hist(DGEN,tao);
      FIMP=hist(DIMP,tao);
      % We normalize
      paso=tao(2)-tao(1);
      FGEN=FGEN./sum(FGEN.*paso);
      FIMP=FIMP./sum(FIMP.*paso);
      % and plot
      %plot(tao,FGEN,'r',tao,FIMP,'b');
      %ylabel('Funcion densidad de probabilidad');
      %xlabel('distancia');
      % Calcul of the probability of False rejectiom
 
      FNMR{ig}{icomb}=1-cumsum(FGEN.*paso);
      % and false acceptation
      FMR{ig}{icomb}=cumsum(FIMP.*paso);
      %plot(tao,FNMR{ig}{icomb},'r',tao,FMR{ig}{icomb},'b')
      plot(FMR{ig}{icomb},FNMR{ig}{icomb},'-*');
      title(['Receiver Operating Characteristic Curve groups ',num2str(Mcomb(icomb,:))])
      xlabel('FALSE MATCH RATE');ylabel('FALSE NON-MATCH RATE')
      axis([0 0.5 0 0.5]);grid
      pause
   end   
end
return

