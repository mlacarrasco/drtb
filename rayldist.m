% RAYLDIST Rayleight distribution
%
% (c) Miguel Carrasco mlacarrasco@gmail.com

function data=rayldist(x, sigma, pant)
data=(x./sigma^2).* exp(-x.^2./(2*sigma^2));
if pant
    figure, plot(data), axis tight; drawnow;
end
end
