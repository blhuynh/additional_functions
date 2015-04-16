function [allfactors]=listfactors(x)
%find all factors of x

[b,m,n]=unique(factor(x),'last');

%b is all prime factors

occurences=[m(1) diff(m)];

current_factors=[b(1).^(0:occurences(1))]';

for index=2:length(b)
    
    currentprime=b(index).^(0:occurences(index));
    
    current_factors=current_factors*currentprime;
    
    current_factors=current_factors(:);
end

allfactors=sort(current_factors);



