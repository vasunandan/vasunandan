#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coeffs.h"
typedef double** Matrix;
Matrix createMat(int m,int n);
Matrix inverse(Matrix a);
double **matmul(double **a, double **b, int m, int n, int p);

int main(){
	
	Matrix a=createMat(3,3);
	a[0][0]=-1 ;a[0][1]=1 ;a[0][2]=-1 ;
	a[1][0]=8 ;a[1][1]=-6 ;a[1][2]= 2;
	a[2][0]=-5 ;a[2][1]=3 ;a[2][2]=-1 ;
	Matrix g=createMat(3,1);
	g[0][0]=1;g[1][0]=2;g[2][0]=3;
	
	Matrix I=createMat(3,3);
	I[0][0]=-2 ;I[0][1]=0 ;I[0][2]=0 ;
	I[1][0]=0 ;I[1][1]=-2 ;I[1][2]= 0;
	I[2][0]=0 ;I[2][1]=0 ;I[2][2]=-2 ;
	Matrix d=createMat(1,3);
	d[0][0]=-1/2;d[0][1]=1/2;d[0][2]=-1/2;
	Matrix f=matmul(d,a,1,3,3);
	 Matrix e=matmul(f,g,1,3,1);
	
	double determinant = a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1] * (a[1][0]* a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0]
	 * a[2][1] - a[2][0] * a[1][1]);
   double determinant1=-sqrt(determinant);
  
  printf("DETERMINANT OF A is %f \n ", determinant1 );
  
 Matrix k =inverse(a) ;

 //print(k,3,3 );
Matrix M =matmul(I,k , 3, 3, 3);
 print(M,3,3 );
 printf("a=%f \n",M[0][2]);
  printf("b=%f \n",M[2][1]);
  
  print(f,1,3);
 
  print(e,3,1);
  return 0;
}
Matrix inverse(Matrix a){
	int i;
	int j;
	Matrix inv=createMat(3,3);
	double determinant = a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1] * (a[1][0]* a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0]
	 * a[2][1] - a[2][0] * a[1][1]);
 for(i=0;i<3;i++){
      for(j=0;j<3;j++){
           inv[j][i]=(((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant);
          
	   }
	}
   
    return inv;
}
