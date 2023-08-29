function I = FCMCluster(IM,k)

IM = double(IM); 
I = zeros(size(IM));
l1 = size(IM,1);
l2 = size(IM,2);
data=reshape(IM,l1*l2,1);

[center, U, obj_fcn] = fcm(data, k); 
maxU = max(U); 
for i = 1:k
index = find(U(i, :) == maxU);
data(index,1) = i;
end

I = reshape(data, l1,l2);
% figure, imshow(I,[])