function pcomb = stouffer(p)
% Stouffer et al's (1949) unweighted method for combination of 
% independent p-values via z's 
    if isempty(p)
        pcomb = 1;
    elseif(numel(p) ==1)
        pcomb = p;
    else                
        pcomb = 1-normpdf(sum(norminv(1-p)) / sqrt(numel(p)));
%         pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
    end