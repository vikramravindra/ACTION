function [ Z ] = my_zscore( X )
%     col_mean = nanmean(X);
%     col_std = nanstd(X);

    n = size(X, 1);
    col_mean = sum(X) / n;
    centered_X = bsxfun(@minus, X, col_mean);
    col_std = sqrt(sum(centered_X.^2) / (n-1));
    %std(X);
    
    Z = bsxfun(@rdivide, centered_X, col_std);
end

