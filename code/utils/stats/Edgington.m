function [ combined_p ] = Edgington( p )
% Edgington method of combing p-values when the sum of p-values is less than or equal to 1: Edgington (1972) An additive method for combining probability values from independent experiments
    S = sum(p);
    K = numel(p);
    if(S == 0)
        combined_p = 1;
    elseif(S == 1)
        combined_p = S^K / factorial(K);
    else        
        combined_p = sum(arrayfun(@(j) (-1)^j * 10 .^ (K*log10(max(0, S-j)) -logfactorial(K-j) - logfactorial(j)), 0:floor(S)));                
    end    
    combined_p = min(1, combined_p);
    combined_p = max(0, combined_p);    
end

