
function [r] = GaussSiedel(A, b, k)
    
    n = size(A,1);
    
    res = zeros(n,1); // initialisation de x^0
    r = b-A*res; //résidu initial
    
    for i=1:k
        res = res + inv(tril(A))*r; //création de x^k+1
        r = b-A*res; //mise à jour du résidu
    end
    
    r = norm(r)/norm(b);
    
endfunction



function [res] = sol_exacte(T0, T1, n)
    for i=1:n
        res(i) = T0 + (T1-T0)/(n+1);
    end
endfunction
