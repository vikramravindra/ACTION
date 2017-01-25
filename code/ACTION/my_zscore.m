function [ Z ] = my_zscore( X )
%     col_mean = nanmean(X);
%     col_std = nanstd(X);

    col_mean = mean(X);
    col_std = std(X);
    
    centered_X = bsxfun(@minus, X, col_mean);
    Z = bsxfun(@rdivide, centered_X, col_std);
end

