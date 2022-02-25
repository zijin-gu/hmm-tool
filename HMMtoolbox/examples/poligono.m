% It creates a polygon with nvert vertex
% wiithin a circle of the radius rcirc
% and eccentricity ex
% the vertex has a variability defined with rbola
% and with the angle angvert
% The polygone has in center the origin
% 
% The outputs are the polar coordinate
% nomber of points between two vertexes are np
% and the vertexes in the cartesians coordinates d

function [radio,angulo,np,d,xlin,ylin]=poligono(nvert,rcirc,ex,rbola,avpol)

% We build the ellipse
[xelipe,yelipe]=ellipse1(0,0,[rcirc ex]);
[aelipe,relipe]=cart2pol(xelipe,yelipe);
% We construct the nvert vertexes
if lt(nargin,5)
    avpol=[0:nvert-1]*2*pi/(nvert);
end
coord=zeros(nvert,2);
for icoord=1:nvert;
    if gt(avpol(icoord),pi),avpol(icoord)=avpol(icoord)-2*pi; end;
    [aux,irpol]=min(abs(avpol(icoord)-aelipe));
    coord(icoord,:)=[relipe(irpol)*cos(avpol(icoord)) relipe(irpol)*sin(avpol(icoord))];
end
coord=coord+(rand(size(coord))*2*rbola-rbola);
d=round([coord;coord(1,:)]);
% We calculate the cartesians coordinates of the circle
xlin=[];ylin=[];
for icoord=1:nvert;
    np(icoord)=max(abs(d(icoord,1)-d(icoord+1,1)),abs(d(icoord,2)-d(icoord+1,2)))+1;
    xlin=[xlin round(linspace(d(icoord,1),d(icoord+1,1),np(icoord)))];
    ylin=[ylin round(linspace(d(icoord,2),d(icoord+1,2),np(icoord)))];
end
% we migrate to the polares
% the center is 0 0

[angulo,radio]=cart2pol(xlin-15,ylin-15);
return