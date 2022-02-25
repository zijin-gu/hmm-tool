%-----------------------------------------------------------------
% 	
% 
%                 %%%%%%%%%%%%%%
%                 % Etiquetado %
%                 %%%%%%%%%%%%%%
%             
% 
% This function implements the multi-labelling for the discrete HMM. To do that, we use a clustering algorithm like the k-mean to quantify the library VQ (for further information on the k-mean algorithm see kmedia function). The multi-labelling enables to give various labels for a given parameter. The possible labels are directly linked to the number of symbol by state.
% 
% This function, in particular, labels all the data set given a class and a group.
% 
% 
% Entry:	v1 is the database.
% 	agrup defines how to cluster the parameters together
% 	Ns number of symbols by state
% 	Biblio is the matrix of the library given this class and this group.
% TOPN number of labels to take into account for the multi-labelling during the training phase
% 
% Result:	p1 is a matrix where we store the labels for all the vectors.
% 
% 
%------------------------------------------------------
function pl=etiquetado(vl,agrup,Ns,biblio,TOPN)
% To design a discrete HMM, we need to label every parameter 
% of the vetor input from the library calculated
% This program enables the Multilabelling.

% number of parameters
Np=length(agrup)-1;
% The label for each parameter
pl=cell(Np,1);
% we get the labels
for ip=1:Np
    Nslib=Ns(ip);TOPNl=TOPN(ip);
    dinic=agrup(ip);dfinal=agrup(ip+1)-1;
    pl{ip}=zeros(size(vl(:,dinic:dfinal),1),Nslib);
    for iv=1:size(pl{ip},1),
        vaux1=ones(Nslib,1)*vl(iv,dinic:dfinal)-biblio{ip};
        vaux2=1./sum(([realmin*Nslib*ones(Nslib,1) vaux1.*vaux1])');
        [Y,I]=sort(vaux2);
        pl{ip}(iv,I(Nslib-TOPNl+1:Nslib))=Y(Nslib-TOPNl+1:Nslib).*ones(1,TOPNl);
    end;
    % We nornalize.
    pl{ip}=pl{ip}./(ones(Nslib,1)*sum(pl{ip}'))';
end
return
