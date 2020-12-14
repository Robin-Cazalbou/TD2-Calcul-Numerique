
function [r] = Jacobi(A, b, k)
    
    n = size(A,1);
    
    res = zeros(n,1); // initialisation de x^0
    r = b-A*res; //résidu initial
    
    for i=1:k
        res = res + diag(1./diag(A))*r; //création de x^k+1
        r = b-A*res; //mise à jour du résidu
    end
    
    r = norm(r)/norm(b);
    
endfunction
