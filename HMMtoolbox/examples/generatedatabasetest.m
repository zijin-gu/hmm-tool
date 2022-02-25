%------------------------------------------- 
% 
% This script is made as an example.
% We create 4 different ploygons and classify its.
% The second part plots the ploigons and the third one,
% plots some caracteristics of the HMM:
%  Alpha and gamma.
%
%-------------------------------------------

% generation for the polygons:
nvert=3;
%radius for the circle in order to generate the polygons
rcirc=15;
% radius for for a small ball in order to create the vertexes
rbola=rcirc*1.5/9;
% clase 1
% it is a triangular
% clase 2
%it is a 6 side forms
% clase 3
% it is the clase 2 rotated
angrot=pi/nvert;
% clase 4
% is a form with 23 sides.


% number of classes
nc=4;
% number of repetitions
nr=20;

% number of groups
ng=2;
% Matriz for the parameters
vlcp=cell(nc,ng);
for ic=1:nc,
   vlc{ic,1}=cell(nr,1);
   vlp{ic,1}=cell(nr,1);
   for ig=1:ng
      vlcp{ic,ig}=cell(nr,1);
   end
end


% generation nr repetitions
for ir=1:nr
    % generation classe 1
    [radio1,angulo1,np1,d1,xlin1,ylin1]=poligono(nvert,rcirc,0,rbola*2);
    %vlcp{1,1}{ir}=[diff(radio1(1:end-1)') radio1(1:end-2)' angulo1(1:end-2)'];
    vlcp{1,2}{ir}=[diff(radio1(1:end-1)') angulo1(1:end-2)'];
    vlcp{1,1}{ir}=[ radio1(1:end-2)' angulo1(1:end-2)'];
    I1=pol2contourcont(radio1,angulo1);

  
    %generation classe 2
    [radio2,angulo2,np2,d3,xlin1,ylin1]=poligono(nvert+3,rcirc,0,rbola);
    %vlcp{2,1}{ir}=[diff(radio2(1:end-1)') radio2(1:end-2)' angulo2(1:end-2)'];
    vlcp{2,2}{ir}=[diff(radio2(1:end-1)')  angulo2(1:end-2)'];
    vlcp{2,1}{ir}=[ radio2(1:end-2)' angulo2(1:end-2)'];
    I2=pol2contourcont(radio2,angulo2);
  
  
    
    %generation classe 3: 
    [radio3,angulo3]=poligono(nvert,rcirc,0,rbola*2,[0:nvert-1]*2*pi/(nvert)+pi/5);
    longcont=length(radio3);
   %vlcp{3,1}{ir}=[diff(radio3(1:end-1)') radio3(1:end-2)' angulo3(1:end-2)'];
    vlcp{3,2}{ir}=[diff(radio3(1:end-1)') angulo3(1:end-2)'];
    vlcp{3,1}{ir}=[ radio3(1:end-2)' angulo3(1:end-2)'];
    I3=pol2contourcont(radio3,angulo3);
    
    
     %generation classe 4: 
    [radio4,angulo4]=poligono(nvert+20,rcirc,0,rbola);
    %vlcp{4,1}{ir}=[diff(radio4(1:end-1)') radio4(1:end-2)' angulo4(1:end-2)'];
    vlcp{4,2}{ir}=[diff(radio4(1:end-1)')  angulo4(1:end-2)'];
    vlcp{4,1}{ir}=[ radio4(1:end-2)' angulo4(1:end-2)'];
    I4=pol2contourcont(radio4,angulo4);
    
    
    % representations
    if 1
        % This shows the value of the radius and angles for the ng classes abd repetitions
        for ic=1:nc
            h3=figure(1)
            subplot(2,2,ic)
            %[x,y]=pol2cart(eval(['angulo',num2str(ic)]),eval(['radio',num2str(ic)]));
            %plot(x,y);axis square; axis off
            eval(['imshow(1-I',num2str(ic),')'])
            axis square
            title(['class number ',num2str(ic)])
            h4=figure(2)
            subplot(2,2,ic)
            eval(['ang=angulo',num2str(ic),';']);
            ang=ang+(sign(ang)==-1)*2*pi;
            [ang,ind]=sort(ang);eval(['radio=radio',num2str(ic),'(ind);']);
            plot(ang/pi,radio);axis([-0.1 2.1 min(radio) max(radio)])
            eval(['plot(ang(2:end-1)/pi,radio',num2str(ic),'(2:end-1))'])
            title(['class ',num2str(ic)])
            xlabel('angle');ylabel('radius')        
        end
        drawnow
    end
end  
save parampoligonos vlcp
