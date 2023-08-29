function J=regiongrowing1(I,x,y,reg_maxdist)
 
 
sigma =0;
if(exist('reg_maxdist','var')==0), reg_maxdist=100; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end

J = zeros(size(I)); % Output 
 
Isizes = size(I); % Dimensions of input image

reg_mean = I(x,y); % The mean of the segmented region
reg_size = 1; % Number of pixels in region

 neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 

pixdist=0; % Distance of the region newest pixel to the regio mean

 neigb=[-1 0; 1 0; 0 -1;0 1;-1 1;-1 -1;1 -1;1 1];

  while(pixdist<reg_maxdist&&reg_size<numel(I))

     for j=1:8,
         xn = x +neigb(j,1); yn = y +neigb(j,2);
        
         ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        
    
        if(ins&&(J(xn,yn)==0))%&&(I(xn,yn)> reg_mean-3*sigma)&&(I(xn,yn)< reg_mean+3*sigma)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

     if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
     dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
     reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    sigma = sqrt(((reg_size-1)*(sigma).^2+(((reg_size+1)/reg_size)*( neg_list(index,3)-reg_mean).^2))/reg_size);   
    
     x = neg_list(index,1); y = neg_list(index,2);
    
     neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end

[x1,y1] = find((I> reg_mean-3*sigma)&(I<reg_mean+3*sigma)&(J<1));
if isempty(x1)
    return
else
    while(~isempty(x1))
        x = x1(1);
        y = y1(1);
%--------------------     
 

pixdist=0; 
    while(pixdist<reg_maxdist)

     for j=1:8,
         xn = x +neigb(j,1); yn = y +neigb(j,2);
        
         ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
         
 
        if(ins&&(J(xn,yn)==0))%&&(I(xn,yn)> reg_mean-3*sigma)&&(I(xn,yn)< reg_mean+3*sigma)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

     if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
     dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
     reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    sigma = sqrt(((reg_size-1)*(sigma).^2+(((reg_size+1)/reg_size)*(neg_list(index,3)-reg_mean).^2))/reg_size);   
    
     x = neg_list(index,1); y = neg_list(index,2);
    
     neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end
clear x1 y1
[x1,y1] = find((I> reg_mean-2*sigma)&(I<reg_mean+2*sigma)&(J == 0));    
    %-------------------
    end
end
%
 J=J>1;

