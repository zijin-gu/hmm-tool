%-------------------------------------------------------------------
% 	
%                         %%%%%%%%%%%
%                         %dhmm2chmm%
%                         %%%%%%%%%%%
% 
% 
%
% This function is used to migrate from a Discrete HMM to a Continuous HMM. The idea is to use the DHMM trained to calculate the initials parameters of the CHMM. This function uses the Netlab utility. 
% function dhmm2chmm(filedhmm,filechmm)
% 
% 	filedhmm is the file with the parameters of the discrete HMM.
% 	filechmm is a file containing the parameters of the continuous HMM wanted.
% 
% This function eliminates the variable Ns (number of symbols) and initializes the variables Med Var maxitermi and Ngauss (see the CHMM_DEF function for further information).
% We generate 1000 samples to simulate each Gaussian with the variance and mean calculated.
% We then save the continuous HMM variables in the file 'filechmm'.
% % 
%--------------------------------------------------------------------

function dhmm2chmm(filedhmm,filechmm)

% it reads the file for the DHMM
eval(['load ',filedhmm]);
% it generates the variable for the CHMM
% it eliminates the Ns number of symbols.
vHMM=[' BAKIS salto maxiter umbral maxitermi Ngauss Ne A B Med Var Pi men'];
% it generathes the Gaussians
Ngauss=cell(ng,1);
for ig=1:ng   
   Ngauss{ig}=6.*ones(Np(ig),1);
end
% The mean, varaiances and coefficienct of the mixture of Gaussians
% for the CHMM
Bcont=cell(nc,ng);
Med=cell(nc,ng);
Var=cell(nc,ng);
for ic=1:nc,
   for ig=1:ng
      Bcont{ic,ig}=cell(Np(ig),1);
      Med{ic,ig}=cell(Np(ig),1);
      Var{ic,ig}=cell(Np(ig),1);      
      for ip=1:Np(ig)
         Bcont{ic,ig}{ip}=cell(Ne(ic,ig),1);
         Med{ic,ig}{ip}=cell(Ne(ic,ig),1);
         Var{ic,ig}{ip}=cell(Ne(ic,ig),1);
         for ie=1:Ne(ic,ig),
            Bcont{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),1);      
            Med{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),agrup{ig}(ip+1)-agrup{ig}(ip));
            Var{ic,ig}{ip}{ie}=zeros(Ngauss{ig}(ip),agrup{ig}(ip+1)-agrup{ig}(ip));
         end
      end
   end
end

for ig=1:ng  
    for ic=1:nc, 
        hora=clock;
        fprintf('Transformation Group %g class %g. Hour: %g:%g\n',ig,ic,hora(4),hora(5))
        Med{ic,ig}=cell(Np(ig),1);
        Var{ic,ig}=cell(Np(ig),1);
        Bcont{ic,ig}=cell(Np(ig),1);
        for ip=1:Np(ig)
            Med{ic,ig}{ip}=cell(Ne(ic,ig),1);
            Var{ic,ig}{ip}=cell(Ne(ic,ig),1);
            Bcont{ic,ig}{ip}=cell(Ne(ic,ig),1);
            for ie=1:Ne(ic,ig)
                fpdie=cumsum(B{ic,ig}{ip}(ie,:));
                % 1000 samples are used to simulate each Gaussian
                muestras=zeros(1000,agrup{ig}(ip+1)-agrup{ig}(ip));
                for im=1:1000
                    ind=find(le(rand-fpdie,0));
                    muestras(im,:)=biblio{ig}{ip}(ind(1),:);%.*(1+randn*0.2)
                end
                
                % it calculates the centers
                mix = gmm(agrup{ig}(ip+1)-agrup{ig}(ip),Ngauss{ig}(ip),'spherical');
                options=foptions;
                options(1)=0;		% Prints out error values.
                options(5)=1;        %check variance
                options(14) = 10;
                mix = gmminit(mix,muestras, options);
                mix = gmmem(mix, muestras, options);
                Med{ic,ig}{ip}{ie}=mix.centres;
%                Var{ic,ig}{ip}{ie}=mix.covars+1e-1;
                Var{ic,ig}{ip}{ie}=mix.covars'*ones(1,agrup{ig}(ip+1)-agrup{ig}(ip))+1e-1;
                Bcont{ic,ig}{ip}{ie}=mix.priors';                      
            end
        end
    end   
end
% We save the variables into the filechmm
B=Bcont;
iniciar=0;
umbral=0.0051;
%maxiter=1;
guardar=['save ',filechmm,vDB,vHMM,vTEST,' iniciar vDB vHMM vTEST guardar'];
eval([guardar]);
