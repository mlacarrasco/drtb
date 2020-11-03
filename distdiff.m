% DISTDIFF difference between two distributions
%
% (c) Miguel Carrasco mlacarrasco@gmail.com
function d=distdiff(dA, dB)


%--> Normalization
ndA=(dA(:)-min(dA(:)))./(max(dA(:))-min(dA(:)));
ndB=(dB(:)-min(dB(:)))./(max(dB(:))-min(dB(:)));
%--> RMSD Distancenormalization
d=sum(sqrt((ndA-ndB).^2))/size(ndA,1);

end