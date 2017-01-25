function [Pval] = ER_DensityPval(n, p, rho, subgraphSize)
    function [kappa] = Kappa(p, rho)
        if(rho == 1)
            kappa = log(1/p);
        else            
            kappa = rho*log(rho/p)  + (1-rho)*log( (1-rho)/(1-p) );
        end
    end
    
    kappa = Kappa(p, rho);

    r0 = ( log(n) - log(log(n)) + log(kappa) )   /  kappa;        
    if(r0 < subgraphSize)            
        epsilon = ( subgraphSize - r0 )   /   r0;
        Pval = (  (1+epsilon)*log(n)  )  / (  n^(epsilon*(1+epsilon)*r0)  );
    else
        Pval   = 1; 
    end
    
    if(Pval > 1 || Pval < 0) 
        Pval = nan;
    end
end
    