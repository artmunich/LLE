%对m*3矩阵的Gram-Schmidt正交化
function A = ThreeGS(V)
v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3);
a1 = v1;
a2 = v2-((a1'*v2)/(a1'*a1))*a1;
a3 = v3-((a1'*v3)/(a1'*a1))*a1-((a2'*v3)/(a2'*a2))*a2;
A = [a1,a2,a3];