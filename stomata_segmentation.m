function stomata_segmentation()
%STOMATA_SEGMENTATION: This is an exmple code of the paper 
% "Automatic stomatal segmentation based on Delaunay-Rayleigh frequency distance".

% Authors: Miguel Carrasco,  Patricio Toledo, Ramiro Velazquez, Odemir Bruno.
% Corresponding Author: Miguel Carrasco, mlacarrasco@gmail.com

% version 1.0. 04/10/2020
% runs on Matlab R2015b

clear all;
clc
close all;

%>> PARAMETERS
delta_t      = 1/7;
kappa        = 20;
option       = 2;
distance_tol = 10; %distance pixels

%>> RAILEY distribution
dB=rayldist(0:1:99, 10, 0);

%>> IMAGE DATABASE
bd_images= dir('database/*.JPG');
bd_coordinates = dir('database/*.mat');
STATS=[];


for id_ima=1:size(bd_images,1) %imagen de prueba #6
    
    %>> DATABASE READING
    namefile  = bd_images(id_ima).name;
    namepos   = bd_coordinates(id_ima).name;
    str       = sprintf('jatoba2/%s',namefile);
    str_pos   = sprintf('jatoba2/%s',namepos);
    raw_image = imread(str);
    
    load(str_pos);
    image = imresize(raw_image, 0.2);
    fprintf('\nImage File:\t%s\n',namefile);
    
    
    %PHASE 1. PREPROCESSING
    HSV=rgb2hsv(image);
    S=HSV(:,:,2);
   
    %>> PERONA-MALIK
    fprintf('PHASE 1:\n');
    fprintf('Perona-Malik Iteration');
    anim= anisodiff2D(image,3, delta_t,kappa,option);
    dff = anim(:,:,1);   
    [m,n]=size(dff);
    ima=(dff(:)-min(dff(:)))./(max(dff(:))-min(dff(:)));
    ima=reshape(uint8(ima.*255),m,n);
 
    %>>  MEANSHIFT
    fprintf('\nMeanshift Iteration...\n');
    dff_ms = meanShiftPixCluster(ima,5,10, 0.25, 0);
    dff_ms= S./dff_ms;
    dff_ms(isinf(dff_ms))=1;
    
    %-->initialization
    steps_level=linspace(0,1,100);aux=1;
    diff_AB=zeros(size(steps_level,1),2);
     
    
    %PHASE 2. DRTB ALGORITHM
    fprintf('\n PHASE2:\n');
    fprintf('DBRT Algorithm\n');
    h=waitbar(0,'Please Wait..');
    
    for level=steps_level
               
        %>>SEGMENTATION
        [sts, idx]= segmentation(dff_ms, level,0);
        waitbar(level);
        if not( isempty(idx))
            try
                pos_all = cat(1,sts.Centroid);
                pos_sel= pos_all(idx,:);
                x= pos_sel(:,1); y = pos_sel(:,2);
                
                 %>>DELAUNAY TESELATION
                TRI= delaunay(x,y);
                TRIE= [TRI, TRI(:,1)];
                DST=[];
                
                for d=1:3
                    col1=TRIE(:,d);
                    col2=TRIE(:,d+1);
                    DST(:,d) = sqrt((x(col1,:)-x(col2,:)).^2+(y(col1,:)-y(col2,:)).^2);
                    
                end
                
                %>>DISTACE TESELATION
                data= hist(DST(:),100);
                diff_AB(aux,:)=[distdiff(data,dB),  size(sts,1)];
                aux=aux+1;
            catch
                fprintf('Error\n');
            end
        end
    end
    close(h);
    [~, or]=sort(diff_AB(:,1));
    level_best=mean(steps_level(or(1:5)));
    
    %> BEST LEVEL
    fprintf('Selected Level: %g\n', level_best)
    figure(1), plot(diff_AB(:,1)), title('Delauney center-to-center distance frequency');
    axis tight; drawnow;
    
    [sts, idx, area]=segmentation(dff_ms, level_best,1);
    figure, imshow(image); title('Original Image'); axis on;
    hold on;
    
    %-> Outlier Detection
    MAD = median(area)-0.95*median(area);  %# REVISAR!!!
    idx_leys= area>MAD;   
    
    pos_all = cat(1,sts.Centroid);
    pos_sel= pos_all(idx_leys,:);  %<<- idx no idx_leys
    x= pos_sel(:,1); y = pos_sel(:,2);
    
    %>> DELAUNAY TESELATION
    TRI= delaunay(x,y);
    trimesh(TRI,x,y);
    
    %best delauney triangle.
    TRIE= [TRI, TRI(:,1)];
    DST=[];
    
    for d=1:3
        col1=TRIE(:,d);
        col2=TRIE(:,d+1);
        DST(:,d) = sqrt((x(col1,:)-x(col2,:)).^2+(y(col1,:)-y(col2,:)).^2);
        
    end
   
      
    %>> PERFORMANCE
    % TP: True Positive Number
    % FN: False Negative
    % FP: False Positive
    
    TP=0;
    FN=0;
    for k=1:size(PTS,1)
        pK= PTS(k,:);
        inflat_K= repmat(pK,  size(pos_sel,1),1);
        dst_K=sum(sqrt((inflat_K-pos_sel).^2),2);
        [val, id]=sort(dst_K);
        
        if (val(1)<distance_tol)
            TP=TP+1;
            plot(pK(1),pK(2), 'bo','MarkerSize',10);drawnow;
        else
            FN=FN+1;
            
            plot(pK(1),pK(2), 'rx','MarkerSize',10);drawnow;
        end
        
        FP=size(pos_sel,1)-TP;
        
    end
    
    STATS(id_ima,:)=[TP FP FN];
    
    fprintf('id_image: %i\t TP: %i\t FP:%i\t FN:%i\n', id_ima, TP,FP,FN);
    
  
    
end


save('stats.mat', 'STATS');
TPR_stats= STATS(:,1)./(STATS(:,1)+STATS(:,3));
FPR_stats= STATS(:,1)./(STATS(:,1)+STATS(:,2));


end %EOF
