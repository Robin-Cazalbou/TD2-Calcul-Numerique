
function [r] = Richardson(alpha, A, b, k)
    
    n = size(A,1);
    G = eye(n,n)-alpha*A;
    
    res = zeros(n,1); // initialisation de x^0
    r = b-A*res;
    
    for i=1:k
        res = res + alpha*r;
        r = b-A*res;
    end
    
    r = norm(r)/norm(b);
    
endfunction
