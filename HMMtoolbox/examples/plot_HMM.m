        load hmmpoligonos
        load parampoligonos vlcp
        for ic=1:nc
            for ig=1:ng,
                % we etiquate all the repetitions for each class
              for ir=1:size(vlcp{ic,ig},1)
               pl=etiquetado(vlcp{ic,ig}{ir},agrup{ig},Ns{ig},biblio{ig},TOPNtest{ig});
         
              
                % we evaluate the probability for the 4 HMM
                for ihmm=1:nc;
                    salida(ihmm)=probsec(A{ihmm,ig},B{ihmm,ig},Pi{ihmm,ig},pl);
                    [alfa,beta,c]=alfabeta(A{ihmm,ig},B{ihmm,ig},Pi{ihmm,ig},pl);
                    T=size(pl{1},1);
                    gama=zeros(Ne(ihmm),T);
                    for t=1:T,
                        gama(:,t)=(alfa(:,t).*beta(:,t))./(alfa(:,t)'*beta(:,t))+realmin;
                    end;
                    prob=zeros(Ne(ihmm),T);
                    for t=1:T
                        prob(:,t)=prodBO(B{ihmm,ig},pl,t);
                    end
                    Eps=zeros(Ne(ihmm),Ne(ihmm));
                    for t=1:T-1,
                        for i=1:Ne(ihmm),
                            Eps(i,:)=Eps(i,:)+alfa(i,t)*(A{ihmm,ig}(i,:).*(prob(:,t+1))'.*beta(:,t+1)');
                        end;
                    end;   
                    qP=viterbi(A{ihmm,ig},B{ihmm,ig},Pi{ihmm,ig},pl);
                    h3=figure(3)
                    subplot(211)
                    %plot(eval(['angulo',num2str(ic),'/pi']),eval(['radio',num2str(ic)]))
                    % we plot the gamma values for each repetition and classes
                    %plot(eval(['radio',num2str(ic)]))
                    title(['class real: ',num2str(ic),'HMM model of the class: ',num2str(ihmm),' group and repetition: ', num2str(ig),' ',num2str(ir)])
                    subplot(212)
                    imagesc(gama);
                    title(['graph 1: radius, graph 2: Gamma']);
                    h4=figure(4) 
                    subplot(211)
                    plot(qP);
                    title([' class real: ',num2str(ic),' HMM model of the class: ',num2str(ihmm),' group and repetition: ', num2str(ig),' ',num2str(ir),])
                    subplot(212)
                    imagesc(alfa);
                    title(['graph 1: most probable state sequence with viterbi, graph 2: Alpha']);
                    pause
%                    saveas(h3,'fig3.jpg','jpg');
%              saveas(h4,'fig4.jpg','jpg');
                end
            end;  
        end
    end
    
