function [p_kost,p_fisher,Ckost,df_kost] = KostsMethod_mvhg_Cov(p_values, covar_matrix)
    %% Kost combination
    m = size(covar_matrix,1);
    df_fisher = 2.0*m;
    Expected = 2.0*m;
    cov_sum = sum(sum(covar_matrix))-sum(diag(covar_matrix));
    Var = 4.0*m+cov_sum;
    Ckost = Var/(2*Expected);
    df_kost = 2*(Expected^2)/Var;
    if df_kost > df_fisher;
        df_kost = df_fisher;
        Ckost = 1;
    end

    x = 2*sum(-log(p_values));
    p_kost = chi2cdf(1.0*x/Ckost,df_kost,'upper');
    p_fisher = chi2cdf(1.0*x,df_fisher,'upper');
end

