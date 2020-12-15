
function [xk, r] = Jacobi(A, b, k)
    
    xk = zeros(size(A,1),1); // initialisation de x^0
    r(1) = norm(b-A*xk); //résidu initial
    
    for i=1:k
        xk = xk + inv(diag(diag(A)))*(b-A*xk); //création de x^k+1
        r(i+1) = norm(b-A*xk); //mise à jour du résidu
    end
    
endfunction


function [r] = JacobiWhile(A, b, epsilon)
    
    xk = zeros(size(A,1),1);
    r(1) = norm(b-A*xk);
    i=1;
    
    while (r(i)>epsilon)
        i = i+1;
        xk = xk + inv(diag(diag(A)))*(b-A*xk);
        r(i) = norm(b-A*xk);
    end
    
endfunction

