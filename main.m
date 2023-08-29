clc
clear
close all
%--------------------------------------------- Load data
 
I = dicomread('subject04_t10090.dcm');
 I = rgb2gray(I);
 
% I = imread('fig1.tif');
I = double(I);
figure, imshow(I,[])
% I = rgb2gray(I);
%% -------------- GM & WM & CSF & BG separation
nbins = 50;
[nelements,centers] = hist(I(:),nbins);
figure, plot(centers,nelements)
%------------------------------------------------------
I1 = I<40;
se = strel('disk',10);
I1 = imclose(I1,se);
figure, imshow(I1,[])
I = I.*(~I1);
Imm = im2bw(I);
[x,y] = find(Imm == 1);
X_max= max (x);
X_min= min (x);
Y_max= max (y);
Y_min= min (y);
I = I(X_min:X_max,Y_min:Y_max);
figure, imshow(I,[])
 
figure, imshow(I,[])
s = 3;
[x,y] = ginput(s);
y = round(y);
x = round(x);
hold on
plot(x,y,'r*')
% -------------------------------------------- Deep Based Region Growing 

	x1 = x(1);
	y1 = y(1);
    J_rg1=regiongrowing(I,y1,x1,15);
	x1 = x(2);
	y1 = y(2);
    J_rg2=regiongrowing(I,y1,x1,15);
    x1 = x(3);
	y1 = y(3);
    J_rg3=regiongrowing(I,y1,x1,15);
    
%% Ground truth
Immm = dicomread('subject04_t10090.dcm');

I2= dicomread('subject04_gm_v0180.dcm');
I2 = rgb2gray(I2);
I2 = I2(1:400,:);
I3= dicomread('subject04_wm_v0180.dcm');
I3 = rgb2gray(I3);
I3 = I3(1:400,:);
I4= dicomread('subject04_csf_v0180.dcm');
I4 = rgb2gray(I4);
I4 = I4(1:400,:);
% figure, imshow(I2,[]), title('CSF')
% figure, imshow(I3,[]), title('GM')
% figure, imshow(I4,[]), title('WM')
% % 
Imm = im2bw(I4);
[x,y] = find(Imm == 1);
X_max= max (x);
X_min= min (x);
Y_max= max (y);
Y_min= min (y);
Imn2 = I2(X_min:X_max,Y_min:Y_max);
Imn2 = imresize(Imn2,size(I));
Imn3 = I3(X_min:X_max,Y_min:Y_max);
Imn3 = imresize(Imn3,size(I));
Imn4 = I4(X_min:X_max,Y_min:Y_max);
Imn4 = imresize(Imn4,size(I));
% 
% 
figure, imshow((Immm))
figure, imshow((Imn2))
figure, imshow((Imn3))
figure, imshow((Imn4))
Imn2 = im2bw(Imn2);
Imn3 = im2bw(Imn3);
Imn4 = im2bw(Imn4);

%% ---------------   Deep Based Region Growing 
%       Deep algorithm 
 bw = FCMCluster(I,4);
figure, imshow(bw,[])
 

for i =1:4
    m(i) = mean(I(bw ==i));
    s(i) = std(I(bw ==i));
    if (m(i) <130)&&(m(i) >100)
        [x,y] = find(bw == i);
        GM = [x,y];
        gm_fcm = (bw == i);
    elseif (m(i) <150)&&(m(i) >130)
        [x,y] = find(bw == i);
        WM = [x,y];
        wm_fcm = (bw == i);
    elseif (m(i) <100)&&(m(i) >50)
        [x,y] = find(bw == i);
        CSF = [x,y];
        csf_fcm = (bw == i);
    end
end
%------------ Deep parameters
FunName = 'MyCostFunction';
NPar = 3;
VarHigh1 = size(GM,1);
VarHigh2 =  size(WM,1);
VarHigh3 =  size(CSF,1);

PopSize = 5;
MaxGenerations = 20;
RecomPercent = 10/100;
CrossPercent = 50/100;
MutatPercent = 1 - RecomPercent - CrossPercent;

RecomNum = round(PopSize*RecomPercent);
CrossNum = round(PopSize*CrossPercent);
if mod(CrossNum,2)~=0;
    CrossNum = CrossNum-1;
end
MutatNum = PopSize - RecomNum - CrossNum;
%--------------------
Pop =[round( rand(PopSize,1)*(VarHigh1)),...
       round( rand(PopSize,1)*(VarHigh2)),...
       round( rand(PopSize,1)*(VarHigh3))];

Cost = feval(FunName,Pop,PopSize,NPar,I,GM,WM,CSF,Imn2,Imn3,Imn4);

[Cost, inx] = sort(Cost);
Pop = Pop(inx,:);
%--------------------
%% Main Loop
MinCostMat = [];
figure, 
for Iter = 1:MaxGenerations
    % Recombination
    RecomPop = Pop(1:RecomNum,:);
    % CrossOver
    SelectedParents = MySelectionFcn(Cost,CrossNum,2);
    CrossPop = [];
    
    for ii = 1:2:CrossNum
        Par1Index = SelectedParents(ii);
        Par2Index = SelectedParents(ii+1);
        
        Parent1 = Pop(Par1Index,:);
        Parent2 = Pop(Par2Index,:);
        
        [OffSpring1 OffSpring2] = MyCrossOverFcn(Parent1,Parent2);
        
        CrossPop = [CrossPop ; OffSpring1 ; OffSpring2];
    end
    % Mutation
    MutatPop =  [round( rand(MutatNum,1)*(VarHigh1)),...
        round( rand(MutatNum,1)*(VarHigh2)),...
        round( rand(MutatNum,1)*(VarHigh3))];
    
    % New Population
    Pop = [RecomPop ; CrossPop ; MutatPop];
    Cost =  feval(FunName,Pop,PopSize,NPar,I,GM,WM,CSF,Imn2,Imn3,Imn4);
    
    [Cost, inx] = sort(Cost);
    Pop = Pop(inx,:);
    
    % Result Display
    BestCost = Cost(1);
    BestSolution = Pop(1,:);
    [Iter BestCost]
    MinCostMat = [MinCostMat min(Cost)];
    
    plot(MinCostMat,'-r','Linewidth',3);
    xlim([1 MaxGenerations])
%     pause(.01)
end
%%
BestSolution = round(Pop(1,:))
BestCost = Cost(1)
xlabel('Generations')
ylabel('Cost')

	x1 = GM(BestSolution(1,1),1);
	y1 = GM(BestSolution(1,1),2);
    J1=regiongrowing1(I,x1,y1,15); 
	x1 = WM(BestSolution(1,2),1);
	y1 = WM(BestSolution(1,2),2);
    J2=regiongrowing1(I,x1,y1,15); 
	x1 = CSF(BestSolution(1,3),1);
	y1 = CSF(BestSolution(1,3),2);
    J3=regiongrowing1(I,x1,y1,15); 
figure, 
subplot(2,2,1),imshow(I,[])
subplot(2,2,2),imshow(J1)
subplot(2,2,3),imshow(J2)
subplot(2,2,4),imshow(J3)
figure, 
subplot(2,2,1),imshow(I,[])
subplot(2,2,2),imshow(Imn2)
subplot(2,2,3),imshow(Imn3)
subplot(2,2,4),imshow(Imn4)

MSE_proposed =sqrt((sum(sum((Imn2-J1).^2))+sum(sum((Imn3-J2).^2))+sum(sum((Imn4-J3).^2)))/(size(I,1)*size(I,2)))
MSE_RG = sqrt((sum(sum((Imn2-J_rg1).^2))+sum(sum((Imn3-J_rg2).^2))+sum(sum((Imn4-J_rg3).^2)))/(size(I,1)*size(I,2)))

