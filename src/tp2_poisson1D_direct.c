/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la; // nbpoints = n+2, la = n
  int ku, kl, kv, lab; //ku : nb diagonale supérieures, kl : nb diago inférieures
                      //lab : nb de lignes de AB, kv : ??
  int *ipiv;
  int info;
  int NRHS; //taille de RHS
  double T0, T1; //conditions aux bords en 0 et 1
  double *RHS, *EX_SOL, *X; //tableaux pour : RHS = b (second membre), EX_SOL = solution exacte, X = solution approchée
  double *AB; //tableau pour la matrice A en stockage bandes

  double temp, relres; //pour le calcul de l'erreur relative résiduelle

  NRHS=1; //nb colonnes du second membre
  nbpoints=102; //x_0, x_1, ..., x_101, i.e. x_0 à x_(n+1) avec n = 100
  la=nbpoints-2; //on enlève x_0 et x_(n+1) : la = n
  T0=-5.0; //valeur au bord gauche
  T1=5.0; //valeur au bord droit

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la); //allocation tableau de doubles, de "la" cases
  EX_SOL=(double *) malloc(sizeof(double)*la); //idem avec "la" cases
  X=(double *) malloc(sizeof(double)*la); //idem avec "la" cases

  set_grid_points_1D(X, &la); //remplit le tableau X avec les valeurs de x_1 à x_n
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1); //initiliase le second membre avec g=0, et T0 et T1
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  write_vec(RHS, &la, "RHS.dat"); //prend les valeurs dans RHS pour remplir le fichier RHS.dat
  write_vec(EX_SOL, &la, "EX_SOL.dat"); //idem mais avec EX_SOL
  write_vec(X, &la, "X_grid.dat"); //idem mais avec X_grid

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la); //tableau de doubles, de taille lab * la (matrice A en stockage bandes)

  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int)); //tableau de int initialisé à 0, de taille la

  int row = 1; //

  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");

    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  }


  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");


  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1); //produit scalaire du second membre avec lui-même
  temp = sqrt(temp); //calcul norme du second membre f
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1); //EX_SOL = EX_SOL - RHS
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp; //erreur relative résiduelle : valeur finale

  printf("\nThe relative residual error is relres = %e\n",relres);




  //=============================================================================
  // Test de dgbmv :
    double *CB;
    CB = (double *) malloc(sizeof(double)*(lab-1)*la);
    for (int i=0; i<la; i++){
      CB[i*(lab-1)]=-1.0;
      CB[i*(lab-1) + 1]=2.0;
      CB[i*(lab-1) + 2]=-1.0;
    }
    CB[0]=0.0;
    CB[(lab-1)*la - 1] = 0.0;

    double *RHS_2;
    RHS_2=(double *) malloc(sizeof(double)*la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    //dgbmv en col major :
    cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, CB, lab-1, EX_SOL, 1, 0, RHS_2, 1);

    // Erreur relative résiduelle
    temp = cblas_ddot(la, RHS, 1, RHS,1);
    temp = sqrt(temp);
    cblas_daxpy(la, -1.0, RHS, 1, RHS_2, 1);
    relres = cblas_ddot(la, RHS_2, 1, RHS_2,1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nL'erreur relative résiduelle sur le second membre (colmajor) est = %e\n",relres);

  //============================================================================
  //============================================================================
  // Test de dgbmv :
    double *DB;
    DB=(double *) malloc(sizeof(double)*(lab-1)*la);
    for (int i=0; i<la; i++){
      DB[i]=-1.0;
      DB[i + la]=2.0;
      DB[i + 2*la]=-1.0;
    }
    DB[0]=0.0;
    DB[(lab-1)*la -1]=0.0;

    double *RHS_3;
    RHS_3=(double *) malloc(sizeof(double)*la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    //dgbmv en row major :
    cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, 1.0, DB, la, EX_SOL, 1, 0.0, RHS_3, 1);
    write_vec(RHS_3, &la, "my_dgbmv.dat");

    // Erreur relative résiduelle
    temp = cblas_ddot(la, RHS, 1, RHS,1);
    temp = sqrt(temp);
    cblas_daxpy(la, -1.0, RHS, 1, RHS_3, 1);
    relres = cblas_ddot(la, RHS_3, 1, RHS_3,1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nL'erreur relative résiduelle sur le second membre (rowmajor) est = %e\n",relres);

  //============================================================================



  //------------ Exercice 5 ------------------------
  //============================================================================
  // Matrice du système (matrice de Poisson)
  double* EB;
  lab=3;
  EB = (double *) malloc(sizeof(double)*lab*la);
  EB[0]=0.0;
  for (int i=1; i<la; i++){
    EB[i]=-1.0;
  }
  for (int i=0; i<la; i++){
    EB[i+la]=2.0;
  }
  for (int i=0; i<la-1; i++){
    EB[i+2*la]=-1.0;
  }
  EB[lab*la-1]=0.0;

  // Second membre du système (conditions aux bords T0 et T1, le reste = 0)
  double* b;
  b = (double *) malloc(sizeof(double)*la);
  set_dense_RHS_DBC_1D(b,&la,&T0,&T1);

  // Factorisation LU
  mylu_tridiag_rowmajor(EB, la);
  myresol_LU_tridiag_rowmajor(EB, la, b);
  write_vec(b, &la, "exo5_b.dat");

  // Validation de la méthode : calcul de l'erreur relative
  temp = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, EX_SOL, 1, b, 1);
  relres = cblas_ddot(la, b, 1, b, 1);
  relres = sqrt(relres);
  relres = relres / temp;

  printf("\nExercice 5 :\nL'erreur relative résiduelle sur la solution avec mylu est = %e\n",relres);


  //============================================================================




  free(RHS); //libération premier malloc
  free(EX_SOL); //libération deuxième malloc
  free(X); //libération troisième malloc
  free(RHS_2);
  free(RHS_3);
  free(AB);
  free(CB);
  free(DB);
  free(EB);
  free(b);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
}
