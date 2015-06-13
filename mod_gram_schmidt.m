function[result] = mod_gram_schmidt(varargin)
    % MATLAB implementation of Modified Gram-Schmidt Algorithm.

    % Convert input arguments into A matrix.
    A_raw = cellfun(@transpose, varargin, 'UniformOutput', false);  
    A = cell2mat(A_raw);
    
    % Initialise the first 'page' of V matrices.
    for i=1:nargin,
        V(:,i,1) = A(:,i);
    end
    
    % Initialise the first column of Q.
    Q(:,1) = V(:,1,1) / norm(V(:,1,1));
    
    % Initialise the rest of the column vectors of Q.
    for i=2:nargin,
        for j=i:nargin,
            V(:,j,i) = V(:,j,i-1) - dot(V(:,j,i-1), Q(:,i-1)) * Q(:,i-1);
        end
        Q(:,i) = V(:,i,i) / norm(V(:,i,i));
    end    
    result = Q;
end
