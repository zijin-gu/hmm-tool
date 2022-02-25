%----------------------------------------------------------------------------- 
%
% 
%                         %%%%%%%%%%%
%                         % gen_bib %
%                         %%%%%%%%%%%
% 	
%     
%     
% This function computes a library for all the different models’ parameters (a library for) based on the training vectors.
% According the set up made in the DHMM_DEF, the library is generated with the LBG algorithm or with the k-means algorithm (see kmedias function for further information about the k-means).
% 
% 
% Entry: 	fpVQ is the matrix with the vectors to train the library (see DHMM and formato_lectura_secuencial functions)
% Nv number of vectors for the training
% agrup defines how to join the parameters together
% LBG is the choice of the algorithm to compile the library
% LBG = 1 the algorithm is LBM, in the other case we choose kMedias
% dptzoLBG percentage maximum for the distance of the code vector in the algorithm LBG 
% maxiterVQ is the maximum number of iteration to compute the library (condition to stop the algorithm)
% umbralVQ is the threshold condition to stop the library computation (condition to 	stop the algorithm)
% men=1 or 0 if men =1 the library generates message else no.
% biblio is a matrix parameter for the library computation (we define it in DHMM_DEF but calculate it here)
% 
% Result: 	biblio is now the library.
% 
% 
% 
% VQ codebook design:
% We want to design the codebook for each parameter in each class.
% We try to find the codebook size and vectors in order to have the overall distortion minimized.
% With fpVQ the training data set and s a vector from fpVQ
% R the centroids of the library with r a vector of R (in the function is an element from biblio).
% d(s, r)= ||s-r||^2
% 
% LBG (Linde, Buzo and Gray) algorithm:
% 1.	initialization of the library with the centroid r calculated with the vectors fpVQ from this class
% 2.	Define the vectors r1 = r + ?, r2 = r- ?
% 3.	The closest vectors from r1 (r2) are s1 (s2)
% 4.	Search for the centroids from s1 and s2
% 5.	Make the steps 3-4 several times. UmbralVQ and maxiterVQ are the conditions to stop.
% 6.	Make 1 to 5 up to obtain the desired numbers or this class
%
%-----------------------------------------------------------------------------

function biblio=gen_bib(fpVQ,Nv,agrup,Ns,LBG,dpztoLBG,maxiterVQ,umbralVQ,men,biblio)

% number of class and groups
ng=size(Ns,1);
% Number of parameters bu group
Np=zeros(ng,1);
for ig=1:ng
    Np(ig)=size(Ns{ig},1);
end

% Messages Generals
if men,fprintf('CALCUL OF THE VQ\n')
    fprintf('Number of groups to calculte the Library: %g\n',ng)
    if (LBG==1),
        fprintf('Algorithm LBG. Percentage of the shift used: %g\n',dpztoLBG);
    else
        fprintf('Algorithm Kmedias or K-mean\n');
    end;
    fprintf('maximum of iterations, threshold of the end of iteration: %g %g\n',maxiterVQ,umbralVQ);end

% Design of the library of the Quantifier
tiempo=0;
% We calculate the libraries of the parameters for each group
for ig=1:ng,  
    % We print the messages of the particular design
    if men,fprintf('CALCUL of the VQ for the GROUP %g de %g\n',ig,ng)
        fprintf('Number of Libraries (parameters) to generate: %g\n',Np(ig))
        fprintf('Number of centroids for each Library: %g\n',Ns{ig})
        fprintf('Size of the centroids for each Library: %g\n',diff(agrup{ig}));
        fprintf('Number of vectors for the training: %g\n',Nv(ig));end
    for ip=1:Np(ig),
        if men,fprintf('CALCUL OF THE VQ for the GRUP %g PARAMETER %g\n',ig,ip);end
        Nslib=Ns{ig}(ip);
        dim=agrup{ig}(ip+1)-agrup{ig}(ip);
        if (LBG==1),
            % we calculate the first library.
            bib=zeros(1,dim);
            fidVQ=fopen(['c:\temphmm\',fpVQ,'_g',num2str(ig),'_p',num2str(ip),'.tmp'],'r');
            for iv=1:Nv(ig);
                bib=bib+fread(fidVQ,dim,'double')';
            end
            bib=bib/Nv(ig);
            fclose(fidVQ);
            while (2*size(bib,1)<=Nslib),
                % we double the library
                bib=[bib; bib.*(1.+randn(size(bib)).*dpztoLBG)];
                % we increase the efficiency of the library with the kmedia (k-mean algorithm).
                bib=kmedias([fpVQ,'_g',num2str(ig),'_p',num2str(ip)],Nv(ig),bib,maxiterVQ,umbralVQ, men);
                if (men), fprintf('\tEND library of %g centroids.\n',size(bib,1)); end
            end;
            if (size(bib,1)<Nslib)
                addbib = bib(1:Nslib-size(bib,1),:).*(1.+randn(Nslib-size(bib,1),dim).*dpztoLBG);
                bib = [bib; addbib];
                bib=kmedias([fpVQ,'_g',num2str(ig),'_p',num2str(ip)],Nv(ig),bib,maxiterVQ,umbralVQ, men);
                if (men), fprintf('\tEND Library of %g centroids.\n',Nslib); end
            end
            biblio{ig}{ip}=bib;
            if (men), fprintf('End of the calcul of the library for the group %g, parameter %g.\n',ig,ip); end
        else
            % Algorithm of the kmedias (k-mean).
            % We ininciate the library tomando periodically Lb vectors from the 
            % training sequence
            tic
            bib=zeros(Nslib,dim);
            periodo=Nv(ig)/Nslib;
            fidVQ=fopen(['c:\temphmm\',fpVQ,'_g',num2str(ig),'_p',num2str(ip),'.tmp'],'r');
            for ib=1:Nslib;
                bib(ib,:)=fread(fidVQ,dim,'double')';
                fseek(fidVQ,periodo*dim*8,0);
            end
            fclose(fidVQ);
            % And we design the final Library
            bib=kmedias([fpVQ,'_g',num2str(ig),'_p',num2str(ip)],Nv(ig),bib,maxiterVQ,umbralVQ, men);
            biblio{ig}{ip}=bib;
            if men,fprintf('End of the calcul of the Library for the group %g, parameter %g. Time %g\n',ig,ip,toc);end
        end;
    end
end   
return
