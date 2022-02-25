% ------------------------------------------------------------------------------
% 
%                      %%%%%%%%%%%%
%                      % Resulhmm %
%                      %%%%%%%%%%%%
%
% 
% This function calculates the confusion matrix group by group. 
% 	function [Mcmed,Mcvot,recmed]=resulhmm(fhmm)
% 
% Entry:	fhmm is the name of the HMM. We will use it to store the results from the HMM.
% 
% 
% Result:	Mcmed is the matrix of confusion calculated with a media criterion. 
% Mcvot is not yet implemented.
% recmed is the average of the recognition.
% 
% Example: Mc = resulhmm(fhmm);
% Mc{agroup}{ng}: Mc is the matrix with groups clustered agroup by agroup and the Ngth combination.
% Mc{2}{4} is the 4th combination for the groups joined 2 by 2.
% 
%---------------------------------------------------------------------------------
function [Mcmed,Mcvot,recmed]=resulhmm(fhmm)

% It reads the database, the libraries and the HMMs
eval(['load ',fhmm])

% Definition of the matrix of confusion
% Example: Mc = resulhmm(fhmm);
% Mc{agroup}{ng}: Mc is the matrix with groups clustered agroup by agroup and the Ngth combination.
% Mc{2}{4} is the 4th combination for the groups joined 2 by 2.

recmed=cell(ng,1);
Mcmed=cell(ng,1);
for ig=1:ng
   ncomb=size(nchoosek([1:ng],ig),1);
   Mcmed{ig}=cell(ncomb,1);
   recmed{ig}=cell(ncomb,1);
   for icomb=1:ncomb
      Mcmed{ig}{icomb}=zeros(nc,nc);
      recmed{ig}{icomb}=0;
   end 
end

% matrix of confusion (the most voted case)
Mcvot=Mcmed;
% The outputs vectors:
logPO=zeros(nc,ng);
% Mensajes del programa de analisis de resultados.
fprintf('Number of class: %g\n',nc);
fprintf('Number of groups: %g\n',ng);

for ic=1:nc,
   % it makes the test for this realisation:
   for ir=1:size(salhmm{ic,ig},1)
      for ig=1:ng
         % each realisation is made of a group of vectors
         logPO(:,ig)=salhmm{ic,ig}{ir};
      end
      % The outputs from each group is normalized between [0 1] in an independent way.
      logPO=(logPO-ones(nc,1)*min(logPO))./(ones(nc,1)*diff(minmax(logPO')'));
      % "The decision is taken" for different possible combination
      for ig=1:ng
         Mcomb=nchoosek([1:ng],ig);
         ncomb=size(Mcomb,1);
         for icomb=1:ncomb
            % decision by mean
            [aux,Modrec]=max(mean([logPO(:,Mcomb(icomb,:)) zeros(nc,1)]'));
            % matrix of confusion is fulfiled.
            Mcmed{ig}{icomb}(ic,Modrec)=Mcmed{ig}{icomb}(ic,Modrec)+1;
            % decision: the most voted case
            [Y,I]=sort(logPO(:,Mcomb(icomb,:)));
            [aux,Modrec]=max(hist(I(nc,:),1:nc));
            Mcvot{ig}{icomb}(ic,Modrec)=Mcvot{ig}{icomb}(ic,Modrec)+1;
         end
      end
      % Mean of each decision
      %fprintf('repetition: %g\tclass real: %g\tclass estimada: %g\n',ir,ic,Modrec);
   end;
end;

%  Matriz of confusion is printed combining all the classes together
if 0
   Mp=[Mc1 diag(Mc0)./sum(Mc0')'.*100];
   fprintf('CLASS');for ic=1:nc,fprintf('\t %2.0f',ic-1);end;fprintf('\t OK\n')
   for i=1:nc,
      fprintf('%2.0f',i-1);
      for j=1:nc,
         fprintf('\t%3.0f',Mp(i,j))
      end;
      fprintf('\t%4.2f %%\n',Mp(i,nc+1));
   end;
end
% Mean of the total recognition
for ig=1:ng
   Mcomb=nchoosek([1:ng],ig);
   ncomb=size(Mcomb,1);
   for icomb=1:ncomb
      recmed{ig}{icomb}=sum(diag(Mcmed{ig}{icomb}))./sum(sum(Mcmed{ig}{icomb})).*100;
      fprintf(['RATE OF RECOGNITION BY GROUP ',num2str(Mcomb(icomb,:)),': %g\n'],recmed{ig}{icomb})
      recvot=sum(diag(Mcvot{ig}{icomb}))./sum(sum(Mcvot{ig}{icomb})).*100;
%      fprintf(['RATE OF RECOGNITION BY VOTED GROUP ',num2str(Mcomb(icomb,:)),': %g\n'],recvot)
   end
end
return
