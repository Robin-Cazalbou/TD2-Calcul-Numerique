/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"


void set_GB_operator_rowMajor_poisson1D(double* AB, int *lab, int *la){
  for (int i=0; i<(*la); i++){
    AB[i]=0.0;
    AB[i + (*la)]=-1.0;
    AB[i + (*la)*2]=2.0;
    AB[i + (*la)*3]=-1.0;
  }
  AB[(*la)]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
        AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}

  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
        AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

// Initialise le second membre, pour g=0, avec T0 et T1
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}

// Crée grâce à x_1 ,..., x_n la solution exacte aux points x_1, ..., x_n : T(x_1), ..., T(x_n)
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1)); //h = 1/(n+1)
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	      fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
        fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

// Prend les valeurs contenues dans le vecteur donné et remplit le fichier donné avec ces valeurs
void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

//=========================================================================
// Factorisation LU avec stockage en priorité ligne :
void mylu_tridiag_rowmajor(double* AB, int n){
  AB[2*n]=AB[2*n]/AB[n]; //calcul de m_1 (à part car n'a pas de m_i*c_i au dénominateur)
  for (int i=1; i<n-1; i++){ //calcul des m_i sauf m_1
    AB[2*n+i]=AB[2*n+i]/(AB[n+i]-AB[2*n+i-1]*AB[i]);
  }
  for (int i=1; i<n; i++){ //calcul des d_i (diagonale de U)
    AB[n+i]= AB[n+i]-AB[2*n+i-1]*AB[i];
  }
}

// Résolution d'un système linéaire Ax=b, avec A sous forme LU compact (A tridiagonale)
// /!\ le second membre b est modifié à la fin de l'appel, et contient la solution
// du système Ax=b
void myresol_LU_tridiag_rowmajor(double* AB, int n, double* b){
  //Descente avec L : stockage du résultat de Ly=b dans b
  for (int i=1; i<n; i++){
    b[i]=b[i]-AB[2*n-1 + i]*b[i-1];
  }
  //Remontée avec U : stockage du résultat de Ux=y (ou Ux=b avec le stockage précédent dans b) dans b
  b[n-1]=b[n-1]/AB[2*n-1];
  for (int i=n-2; i>-1; i--){
    b[i]=(b[i]-AB[i+1]*b[i+1])/AB[n+i];
  }
}

//=========================================================================

void eig_poisson1D(double* eigval, int *la){
  int ii;
  double scal;
  for (ii=0; ii< *la; ii++){
    scal=(1.0*ii+1.0)*M_PI_2*(1.0/(*la+1));
    eigval[ii]=sin(scal);
    eigval[ii]=4*eigval[ii]*eigval[ii];
  }
}

double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax=sin(*la *M_PI_2*(1.0/(*la+1)));
  eigmax=4*eigmax*eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin=sin(M_PI_2*(1.0/(*la+1)));
  eigmin=4*eigmin*eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la){
  //TODO
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit){
  //TODO
}
