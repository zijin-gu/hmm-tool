%---------------------------------------------------------------------------- 
% 
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     	    % formato_lectura_secuencial %
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         
% This function splits the database to smaller files in order to decrease the HMM computational time. 
% 
% Entry:	fptrain is a file containing the full database.
% 	fplecsec: the database split and stored at the direction c:\temphmm\fplecsec_cclase_ggrupo.mat 
% 	agrup: if we have agrup, we must prepare the training data set for VQ with a clustering algorithm.
% 	fpVQ: is a file containing the training database for VQ.
% 
% Result:	lgrupo: is a vector containing the number vectors for every group.
% 	Lrep: is a matrix with the length for each class, each group and repetition. lrep{ic,ig}(ir) contains the length of the icth class, igth group and irth vector.
%----------------------------------------------------------------------------
function [lrep,lgrupo]=formato_lectura_secuencial(fptrain,fplecsec,agrup,fpVQ)
% The DB is readen
eval(['load ',fptrain,' vl']);
% se definen variable
[nc,ng]=size(vl);
nrep=zeros(nc,ng);
dim=zeros(ng,1);
lrep=cell(nc,ng);
lgrupo=0;
% We split the DB by class en groups (gathering)
for ig=1:ng
    dim(ig)=size(vl{1,ig}{1},2);
    for ic=1:nc
        nrep(ic,ig)=size(vl{ic,ig},1);
        lrep{ic,ig}=zeros(nrep(ic,ig),1);
        fid=fopen(['c:\temphmm\',fplecsec,'_c',num2str(ic),'_g',num2str(ig),'.tmp'],'w');
        for ir=1:nrep(ic,ig)
            vlt=vl{ic,ig}{ir};
            lrep{ic,ig}(ir)=size(vlt,1);
            fwrite(fid,vlt,'double');
        end
        fclose(fid);
    end
end
% We store the DB for the training
if gt(nargin,2),
    lgrupo=zeros(ng,1);
    for ig=1:ng,
        Np=length(agrup{ig})-1;
        for ip=1:Np,
            lgrupo(ig)=0;%all the parameters from the same group have the same length
            fid=fopen(['c:\temphmm\',fpVQ,'_g',num2str(ig),'_p',num2str(ip),'.tmp'],'w');
            for ic=1:nc,
                for ir=1:nrep(ic,ig),
                    vlt=vl{ic,ig}{ir}(:,agrup{ig}(ip):agrup{ig}(ip+1)-1);
                    lgrupo(ig)=lgrupo(ig)+size(vlt,1); %number of vectors for each group
                    fwrite(fid,vlt','double');
                end
            end
            fclose(fid);
        end
    end
end
return
