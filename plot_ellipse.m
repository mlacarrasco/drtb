% PLOT ELLIPSE (auxiliary function)
%
% (c) Miguel Carrasco mlacarrasco@gmail.com

function [x, y]= plot_ellipse(x1,y1,x2,y2)
a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
b = a*sqrt(1-e^2);
t = linspace(0,2*pi);
X = a*cos(t);
Y = b*sin(t);
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);

end
