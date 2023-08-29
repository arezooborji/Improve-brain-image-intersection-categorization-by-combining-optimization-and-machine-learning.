function [Off1 Off2] = MyCrossOverFcn(Parent1,Parent2)
Beta1 = rand;
Beta2 = rand;
Off1 = Beta1*Parent1 + (1-Beta1)*Parent2;
Off2 = Beta2*Parent1 + (1-Beta2)*Parent2;

 