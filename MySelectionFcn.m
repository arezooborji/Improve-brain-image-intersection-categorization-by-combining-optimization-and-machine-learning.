function SelectionIndexes=MySelectionFcn(Cost,NumOfSelection,Mode)
if Mode == 1 % 
    NumOfAll = numel(Cost);
    
    for ii = 1:NumOfSelection;
        R = randperm(NumOfAll);
        R = R(1:3);
        CR = Cost(R);
        Inx = find(CR == min(CR)); Inx = Inx(1);
        SelectionIndexes(ii) = R(Inx);
    end
else
    if Mode == 2 % Cost Weighting
        MaxCost = max(Cost);
        Cost2 = MaxCost-Cost;
        PDF = Cost2/sum(Cost2);
        
    elseif Mode == 3 % Rank Weighting
        MaxCost = max(Cost);
        Cost2 = MaxCost-Cost;
        PDF = Cost2/sum(Cost2);
    end
    
    CDF = [0];
    
    for i=1:numel(Cost);
        CDF(i+1)=CDF(i)+PDF(i);
    end
    
    for ii = 1:NumOfSelection
        R = rand;
        D = CDF-R;
        D2 = D(D<0);
        Selection=numel(D2);
        SelectionIndexes(ii)=Selection;
    end
end
end