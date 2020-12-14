clear n A ex_sol b alpha errJ errGS errR;

n = 100;
T0 = -5;
T1 = 5;

A = diag(2*ones(1,n))+diag(-1*ones(1,n-1),-1)+diag(-1*ones(1,n-1),1);
ex_sol = sol_exacte(T0, T1, n);
b(1,1) = T0;
for i=2:n-1
    b(i,1)=0;
end
b(n,1) = T1;




alpha = 0.45;

for k=1:150
    errJ(k) = Jacobi(A, b, k);
    errGS(k) = GaussSiedel(A, b, k);
    errR(k) = Richardson(alpha, A, b, k);
end

subplot(1,2,1)
plot([1:150], [log10(errJ) log10(errGS) log10(errR)]) //bleu : Jacobi, vert : Gauss-Siedel, rouge : Richardson





for j=1:5
    alpha = j/10;
    clear errR;
    for k=1:150
        errR(k) = Richardson(alpha, A, b, k);
    end
    
    subplot(1,2,2)
    plot([1:150], log10(errR))
end



