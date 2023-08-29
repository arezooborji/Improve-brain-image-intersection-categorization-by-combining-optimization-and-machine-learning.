function Cost = MyCostFunction(Pop,PopSize,NPar,I,GM,WM,CSF,Imn2,Imn3,Imn4)
Cost = zeros(PopSize,1);
for ii = 1:PopSize
   P = round(Pop(ii,:)) ;
x = [GM(P(1,1),1),WM(P(1,2),1),CSF(P(1,3),1)];
y = [GM(P(1,1),2),WM(P(1,2),2),CSF(P(1,3),2)];

 
	x1 = x(1);
	y1 = y(1);
    J1=regiongrowing1(I,x1,y1,15); 
	x1 = x(2);
	y1 = y(2);
    J2=regiongrowing1(I,x1,y1,15); 
    x1 = x(3);
	y1 = y(3);
    J3=regiongrowing1(I,x1,y1,15); 
    
 Imn2 = imresize(Imn2,size(I));
 Imn3 = imresize(Imn3,size(I));
 Imn4 = imresize(Imn4,size(I));
 
 Cost(ii) =(sum(sum((Imn2-J1).^2))+sum(sum((Imn3-J2).^2))+sum(sum((Imn4-J3).^2)))/(3*size(I,1)*size(I,2));
end



