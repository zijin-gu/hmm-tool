%-------------------------------------------------------
% 
%               %%%%%%%%%%%
%               % kMedias %
%               %%%%%%%%%%%
% 
% 
% This function computes a k-mean algorithm for the library VQ.
% function	[biblio]=kmedias(fpVQ,Lt,biblio,maxiter,umbral,men);
% 	
% Entry:	fpVQ: vectors for the training if the library
% 	Lt is the vectors numbers for the training
% 	Biblio: library initial
% 	Maxiter: number maximum for the iteration
% 	Umbral is the threshold condition to stop the library computation
% 	Men =1 this function generates messages else no.
% 
% Result:	Biblio: library after the training
% 
% The algorithm k-mean is used to compute the N centroids in each library associated to each parameter in each class. The k-mean problem is to minimize the mean square error (MSE).
% The idea, in this algorithm, is to build the N centroids adding N new training vectors by step. We update the centroids step by step.
% 
% 1.	 The N first training vectors are the centoids(0).
% 2.	 Each vector added is assigned to the closest centroid(t) and the centroids(t+1) are recalculated with the new vectors added
% 3.	The algorithm is terminated when centroid (t)= centroid(t+1)
% 		 (in our case it will be terminated when 
% 		? ||centroid(t)-centroid(t+1)|| < threshold 
% 		or when the iteration number > maxiter)
% 
%--------------------------------------------------------
function	[biblio]=kmedias(fpVQ,Lt,biblio,maxiter,umbral,men);

% Variables to increase the speed.
[Lb,dim]=size(biblio);
punt=0;
vaux1=zeros(size(biblio));
bibnueva=zeros(size(biblio));
unos=ones(Lb,1);
vectxcentro=zeros(Lb,1);;
vector=zeros(dim,1);
difac=0;


% It starts the loop of the intialization
% of the vectorial quantifier 
if men,fprintf('Library, size of: %g\n',Lb);end
cfe=umbral+1;

iter=1;
fidVQ=fopen(['c:\temphmm\',fpVQ,'.tmp'],'r');
for i=1:Lt,
    tabla=fread(fidVQ,dim,'double')';
    vaux1=unos*tabla-biblio;
    [dif,puntero]=min(sum([zeros(1,Lb); (vaux1.*vaux1)']));
    difac=difac+dif;
    bibnueva(puntero,:)=bibnueva(puntero,:)+tabla;
    vectxcentro(puntero)=vectxcentro(puntero)+1;
end;
diftotal=difac/Lt;
fclose(fidVQ);

ind=find(vectxcentro==0);
for ivp=1:length(ind),
    vinsertado=fix(rand*(Lt-1))+1;
    fidVQ=fopen(['c:\temphmm\',fpVQ,'.tmp'],'r');
    fseek(fidVQ,dim*vinsertado*8,0)';    
    bibnueva(ind(ivp),:)=fread(fidVQ,dim,'double')';    
    vectxcentro(ind(ivp))=1;
    fclose(fidVQ);   
end
for ib=1:Lb
    biblio(ib,:)=bibnueva(ib,:)./vectxcentro(ib);
end
hora=clock;
if men,fprintf('\titer.: %g\tdis.: %g\thora: %g:%g\n',iter,diftotal,hora(4),hora(5));end

while (iter<maxiter)&(cfe>umbral),
    
    iter=iter+1;
    bibnueva=zeros(size(biblio));
    vectxcentro=zeros(Lb,1);
    difac=0;
    difant=diftotal;
    
    fidVQ=fopen(['c:\temphmm\',fpVQ,'.tmp'],'r');
    for i=1:Lt,
        tabla=fread(fidVQ,dim,'double')';
        vaux1=unos*tabla-biblio;
        [dif,puntero]=min(sum([zeros(1,Lb); (vaux1.*vaux1)']));
        difac=difac+dif;
        bibnueva(puntero,:)=bibnueva(puntero,:)+tabla;
        vectxcentro(puntero)=vectxcentro(puntero)+1;
    end;
    diftotal=difac/Lt;
    fclose(fidVQ);
    ind=find(vectxcentro==0);
    for ivp=1:length(ind),
        vinsertado=fix(rand*(Lt-1))+1;
        fidVQ=fopen(['c:\temphmm\',fpVQ,'.tmp'],'r');
        fseek(fidVQ,dim*vinsertado*8,0)';    
        bibnueva(ind(ivp),:)=fread(fidVQ,dim,'double')';    
        vectxcentro(ind(ivp))=1;
        fclose(fidVQ);   
    end
    for ib=1:Lb
        biblio(ib,:)=bibnueva(ib,:)./vectxcentro(ib);
    end
    
    % Criterion of the end of iteration.
    cfe=difant/diftotal-1;
    hora=clock;
    if men,fprintf('\titer.: %g,\tdis.: %g\tCfe: %g\thour: %g:%g\n',iter,diftotal,cfe,hora(4),hora(5));end
end;