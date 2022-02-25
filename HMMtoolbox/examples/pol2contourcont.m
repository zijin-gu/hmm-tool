% This function builds a image of the polygins thanks to the 
% cartesian coordinate of the polygon.
function [I]=pol2contourcont(radio,angulo)
% se pasa las coordenadas a cartesianas
[xlin,ylin]=pol2cart(angulo,radio);
xlin=round(xlin+15);ylin=round(ylin);
% We check that the coordinates are continuous
xlinc=[];ylinc=[];
for i=1:length(xlin)-1;
    np(i)=max(abs(xlin(i)-xlin(i+1)),abs(ylin(i)-ylin(i+1)))+1;
    xlinc=[xlinc round(linspace(xlin(i),xlin(i+1),np(i)))];
    ylinc=[ylinc round(linspace(ylin(i),ylin(i+1),np(i)))];
end
% we moove to the image
altura=diff(minmax(ylinc))+1;
anchura=diff(minmax(xlinc))+1;
minx=abs(min(xlinc))+1;
miny=abs(min(ylinc))+1;
I=zeros(altura,anchura);
for i=1:length(xlinc);
    I(xlinc(i)+minx,ylinc(i)+miny)=1;
end
return