%Function of Gram-Schmidt process
%Author: Xiawei Huai
%Date: 2015/6/3

function [Q,znorm]=GS(A)

%   Construct a new orthogonal basis
[nr,nc] = size(A);
Q = zeros(size(A));
znorm = zeros(nc,1);

%   Normalize first vector
znorm(1) = norm(A(:,1),2);
Q(:,1) = A(:,1)/znorm(1); 

for j=2:nc
    %   Project a(j) on Span(q1,...,q(j-1))
    p = zeros(nr,1);
    for i=1:j-1
        p = p + A(:,j)'*Q(:,i)*Q(:,i);
    end
    znorm(j) = norm(A(:,j) - p);
    Q(:,j) = (A(:,j) - p)/znorm(j);
end

