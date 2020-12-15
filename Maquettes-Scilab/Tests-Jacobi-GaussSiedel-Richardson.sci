clear n A ex_sol b alpha errJ errGS errR;

n = 100;
T0 = 5;
T1 = 20;


A = diag(2*ones(1,n))+diag(-1*ones(1,n-1),-1)+diag(-1*ones(1,n-1),1);
b=zeros(n,1);
b(1,1)=T0;
b(n,1)=T1;




alpha = 0.4;
epsilon = 10^(-12);



errJ = log10(JacobiWhile(A, b, epsilon));
errGS = log10(GaussSiedelWhile(A, b, epsilon));
errR = log10(RichardsonWhile(alpha, A, b, epsilon));


abscissesJ = [1:size(errJ,1)];
abscissesGS = [1:size(errGS,1)];
abscissesR = [1:size(errR,1)];

plot(abscissesJ, errJ, abscissesGS, errGS, abscissesR, errR)
title('Historique de convergence')
legend(['Jacobi'; 'Gauss-Siedel'; 'Richardson'])
xlabel("Nombre itérations")
ylabel("Log10 de l''erreur relative résiduelle")





/*
k = 1000;
[xJ, errJ] = Jacobi(A, b, k); //contient les erreurs relatives résiduelles de 0 à k+1
[xGS, errGS] = GaussSiedel(A, b, k);
[xR, errR] = Richardson(alpha, A, b, k);


plot([1:k+1], [log10(errJ) log10(errGS) log10(errR)])
legend(['Jacobi'; 'Gauss-Siedel'; 'Richardson'])
xlabel("Nombre itérations")
ylabel("Log10 de l''erreur relative résiduelle")
*/

