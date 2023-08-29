function d = NeighborPts(x1,y1,I,eps)
 
X = (x1-eps):(x1+eps);
Y = (y1-eps):(y1+eps);

Z1 = repmat(X,1,size(X,2));
Z2 = repmat(Y,size(Y,2),1);
Z2 = Z2(:);
d = [Z1;Z2'];

 