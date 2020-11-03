function [sts, idx, area]= segmentation(ima, level, pant)
% SEGMENTATION is in charge of automatic segmentation  level
% (c) Miguel Carrasco mlacarrasco@gmail.com

idx=[];
bw= im2bw(ima, level);
sts = regionprops(bw, 'all');
area=zeros(1,length(sts));
psx = cat(1,sts.Centroid);

for i=1:length(sts)
    area(i)=sts(i).Area;
    orientation(i)= sts(i).Orientation;
end

if pant
    figure, imshow(ima), title ('Segmented algorithm');axis on;
    hold on;
end


if not(isempty(psx))
    %Based on paper Leys C, Ley C, Klein O, Bernard P, Licata L (2013)
    %Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. 
    %Journal of Experimental So- cial Psychology 49(4):764?766, DOI 10.1016/j.jesp. 2013.03.013
    MAD = 1.4826*median(abs(area-median(area)));
    idx_leys= abs((area-median(area))/MAD) < 3;
    idx= find(idx_leys==1);
end
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);


for i=1:length(idx)
    k= idx(i);
    
    xbar = sts(k).Centroid(1);
    ybar = sts(k).Centroid(2);
    a = sts(k).MajorAxisLength/2;
    b = sts(k).MinorAxisLength/2;
    theta = pi*sts(k).Orientation/180;
    R = [ cos(theta)   sin(theta)
        -sin(theta)   cos(theta)];
    xy = [a*cosphi; b*sinphi];
    xy = R*xy;
    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    if pant
        %str= sprintf('%2.2f',ley_mean_dst(k))
        plot(x,y,'r','LineWidth',2);
        plot(xbar, ybar, 'b*');
        %text(xbar,ybar,str, 'Color','red');
        
        axis equal;
    end
    
end
end