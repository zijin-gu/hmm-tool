%------------------------------------------------------------------------
% 
%                   %%%%%%%%%%%
%                   % verfdp1 %
%                   %%%%%%%%%%%
% 
% This function plots the pdf and the distribution for each state. We stand that the frp is Gaussian.
% 
% 	function [Ptotal,Ototal]=verfdp1(B,Med,Var,agrup)
%     
% B{Np}{N}(Ngauss{ip},1) is a structure with the weight of each Gaussian for each state and to obtain each parameter.
% 
% Med: Med{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the mean for each Gaussian in each state and for each gathering of parameters
% 
% Var: Var{Np}{N}(Ngauss(ip),agrup(ip+1)-agrup(ip)) is the variance for each Gaussian in each state and for each gathering of parameters
% 
% agrup(Np+1,1) defines how to cluster the parameters together
% 
%--------------------------------------------------------------------------
function [Ptotal,Ototal]=verfdp1(B,Med,Var,agrup)

% Variables of the hmm.
% Np: Number of parameters per symbol.
% N: Number of states of the hmm.

Np=size(B,1);
Ngauss=zeros(Np,1);
for ip=1:Np
   Ngauss(ip)=size(B{ip}{1},1);
end
Ne=size(B{1},1);
Ptotal=cell(Np,1);
Ototal=cell(Np,1);
for ip=1:Np,
    Ptotal{ip}=[];
    Ototal{ip}=[];
end
np=100;
for ie=1:Ne
   igraf=1;
   for ip=1:Np,
      dimp=agrup(ip+1)-agrup(ip);
      [desde,I]=min(Med{ip}{ie}-2*Var{ip}{ie});ind=[0:dimp-1]*Ngauss(ip)+I;desde=desde-1*Var{ip}{ie}(ind);
      [hasta,I]=max(Med{ip}{ie}+2*Var{ip}{ie});ind=[0:dimp-1]*Ngauss(ip)+I;hasta=hasta+1*Var{ip}{ie}(ind);
      ind=find(Var{ip}{ie}<1e-3);
      Var{ip}{ie}(ind)=sqrt(2*pi)\ones(size(ind));   
      for id=1:dimp
         O=linspace(desde(id),hasta(id),np);
         P=zeros(np,1);
         for ix=1:np;
             cte=normpdf(O(ix),Med{ip}{ie}(:,id),Var{ip}{ie}(:,id));
            P(ix)=sum(B{ip}{ie}.*cte)+realmin;
         end
         P=P./sum(P);
         Ptotal{ip}=[Ptotal{ip}, P];
         Ototal{ip}=[Ototal{ip}, O'];
         
%         subplot(agrup(Np+1)-1,1,igraf);plot(O,P);
%         title(['state ',num2str(ie),', parameter ',num2str(ip),', size ',num2str(id)])
      end
   end
%   pause
end
return
