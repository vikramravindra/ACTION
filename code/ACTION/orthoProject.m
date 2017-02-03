function [ X_R ] = orthoProject( X, Z )   
    if(isvector(Z))
        X_R = X - Z*(Z'*X) ./ (Z'*Z);
    else
        [~, pivcol] = rref(Z);    
        A = Z(:, pivcol);    
        clear Z;

    %     P = A/(A'*A)*A';
    %     I = eye(size(P));    
    %     P_ortho = I - P;    
    %     X_R = P_ortho*X;
        X_R = X - A/(A'*A)*(A'*X);    
    end        
    
end

