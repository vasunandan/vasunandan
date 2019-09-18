#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coeffs.h"

typedef double** Matrix; typedef double** Vector;
Matrix meshgrid(int len, int start, int stepX, int stepY);
Matrix *createPlane(Vector n, double c, int meshLen);
Vector createVec(double x,double y,double z);
Vector cross_pdt(Vector a, Vector b);
void set(Matrix mat, int row, double x,double y,double z);
double** line_gen(double**A,double**m);

int main(){
	// p = (1 -2 4)^T, x = (1 2 2)^T (x is point on plane)
	Vector p = createVec(-2/5,1,0);
	
	
	// Assigning normal vectors of given planes
	Vector n1 = createVec(3,-1,1);
	Vector n2 = createVec(1,4,-2);
	double c1 = 1; double c2 = 2;
	
	Vector n = cross_pdt(n1, n2);
	
	double**line=line_gen(p,n);
	

	
	int meshLen = 10;
	
	Matrix *plane1 = createPlane(n1, c1, meshLen);
	savetxt(plane1[0], "meshX1.dat", meshLen, meshLen);
	savetxt(plane1[1], "meshY1.dat", meshLen, meshLen);
	savetxt(plane1[2], "meshZ1.dat", meshLen, meshLen);
	
	Matrix *plane2 = createPlane(n2, c2, meshLen);
	savetxt(plane2[0], "meshX2.dat", meshLen, meshLen);
	savetxt(plane2[1], "meshY2.dat", meshLen, meshLen);
	savetxt(plane2[2], "meshZ2.dat", meshLen, meshLen);
	
	
	
	
	savetxt(p,"P.dat",3,1);
	savetxt(line,"line.dat",3,100);
	return 0;
}

Matrix meshgrid(int len, int start, int stepX, int stepY){
	Matrix ret = createMat(len, len);
	for (int i=0; i<len; i++)
		for (int j=0; j<len; j++){
			ret[i][j] = start + i*stepY + j*stepX;
		}
	return ret;
}
Matrix *createPlane(Vector n, double c, int meshLen){
	Matrix* ret = (Matrix*)malloc(3*sizeof(Matrix));
	ret[0] = meshgrid(meshLen, -meshLen/2, 2, 0);
	ret[1] = meshgrid(meshLen, -meshLen/2, 0, 2);
	ret[2] = meshgrid(meshLen, -meshLen/2, 0, 0);
	
	for (int i=0; i<meshLen; i++)
		for (int j=0; j<meshLen; j++)
			ret[2][i][j] = ((c-*n[0]*ret[0][i][j]-*n[1]*ret[1][i][j])*1.0) / (*n[2]); //z = ((c-n[0]*x-n[1]*y)*1.0)/(n[2])
	return ret;
}
Vector createVec(double x,double y,double z){
	Vector vec = createMat(3,1);
	*vec[0]=x; *vec[1]=y; *vec[2]=z;
	return vec;
}
Vector cross_pdt(Vector a, Vector b){
	Matrix n1o = createMat(3,3);
	set(n1o, 0,      0, -*a[2],  *a[1]);
	set(n1o, 1,  *a[2],      0, -*a[0]);
	set(n1o, 2, -*a[1],  *a[0],      0);
	return matmul(n1o, b, 3,3,1);
}
void set(Matrix mat, int row, double x,double y,double z){
	mat[row][0]=x; mat[row][1]=y; mat[row][2]=z;
}
double** line_gen(double**A,double**m){
int len=100;
double**x_AB=createMat(3,len);
double**lam=linspace(0,50,len);
for (int i=0; i<len; i++){
x_AB[0][i]=*A[0]+*lam[i]*(*m[0]);
x_AB[1][i]=*A[1]+*lam[i]*(*m[1]);
x_AB[1][i]=*A[2]+*lam[i]*(*m[2]);
}
return x_AB;
}

