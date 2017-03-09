/*
Solving the 2D steady state heat equation using the Successive Over-Relaxation method
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define L 0.3
#define W 0.4
#define T1 40.0
#define T3 10.0
#define T2 0.0
#define T4 0.0
#define mmax 80
#define kmax 10000
#define xmax 100
#define ymax 100
#define pi 3.14159265359
#define conv 0.000005

void f_an(double**,double**,double***);
int cnorm(double,double);
double f_error(double***,double***);
double abso(double);
void tri(double**,double**,double**,double**,double**,double**);

//Analytical function
void f_an(double** x, double** y, double*** T_an){
	double Ta, Tb;
	int m,i,j;
	double test1,test2,test3;
	
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			Ta = 0.0;
			Tb = 0.0;
			for(m=1;m<mmax;m++){
				Ta = Ta + ((1-cos((m*pi)))/(m*pi))*\
				((sinh((m*pi*W*(1-(*y)[j]))/(L)))/(sinh((m*pi*W)/(L))))*(sin(m*pi*W*((*x)[i])/(L)));
				Tb = Tb + ((1-cos((m*pi)))/(m*pi))*\
				((sinh((m*pi*W*(*y)[j])/L))/(sinh((m*pi*W)/L)))*(sin(m*pi*W*((*x)[i])/L));
			}
			(*T_an)[i][j] = (2*Ta) + (2*(T3/T1)*Tb);
		}
	}
}

//error function
double f_error(double*** T, double*** T_an){
	double err=0.0,sum=0.0;
	int i,j;
	for(i=1;i<xmax;i++){
		for(j=1;j<ymax;j++){
			sum = sum + (((*T)[i][j] - (*T_an)[i][j]) * ((*T)[i][j] - (*T_an)[i][j]));			
		}
	}
	err = (sqrt(sum))/(((double)xmax)*((double)ymax));
	return err;
}

double abso(double n){
	if(n<0){
		return (-1*n);
	}else{
		return n;
	}
}

//solve tridiagonal system
void tri(double** a, double** b, double** c, double** d, double** r, double** s){
	int i;
	for(i=1;i<xmax;i++){
		if(i==1){
			(*r)[1]=(*c)[1]/(*b)[1];
			(*s)[1]=(*d)[1]/(*b)[1];
		}else{
			(*r)[i]=(*c)[i]/((*b)[i]-((*a)[i]*(*r)[i-1]));
			(*r)[xmax-1]=0;
			(*s)[i]=((*d)[i]-((*a)[i]*(*s)[i-1]))/((*b)[i]-((*a)[i]*(*r)[i-1]));
		}
	}		
}

//main function
int main(){
	int i,j,k,q;
	double x_nd, y_nd, T1_nd, T3_nd;
	double dx, dy;
	double B;
	double *x, *y;
	double **T_an, **T, **T_nxt;
	int cflag, cv;
	double Told;
	double error;
	double wf;
	double diff,summ;
	double *a,*b,*c,*d,*r,*s;
	
	x_nd = L/W;
	y_nd = W/W;
	T1_nd = T1/T1;
	T3_nd = T3/T1;
	
	//define increments
	dx = x_nd/(xmax+1);
	dy = y_nd/(ymax+1);
	
	//define \beta
	B = dx/dy;	
	
	//memory allocation
	x = (double*)malloc((xmax+1)*sizeof(double));
	y = (double*)malloc((ymax+1)*sizeof(double));
	
	T_an = (double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
		T_an[i] = (double*)malloc((ymax+1)*sizeof(double));
	}
	
	T = (double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
		T[i] = (double*)malloc((ymax+1)*sizeof(double));
	}
	
	T_nxt = (double**)malloc((xmax+1)*sizeof(double*));
	for(i=0;i<=xmax;i++){
		T_nxt[i] = (double*)malloc((ymax+1)*sizeof(double));
	}
	
	a = (double*)malloc((xmax+1)*sizeof(double));
	b = (double*)malloc((xmax+1)*sizeof(double));
	c = (double*)malloc((xmax+1)*sizeof(double));
	d = (double*)malloc((xmax+1)*sizeof(double));
	r = (double*)malloc((xmax+1)*sizeof(double));
	s = (double*)malloc((xmax+1)*sizeof(double));

	
	//define grid
	for(i=0;i<=xmax;i++){
		x[i] = i*dx;
	}
	for(j=0;j<=ymax;j++){
		y[j] = j*dy;
	}
	
	//calculate analytical solution
	f_an(&x,&y,&T_an);
	
	//relaxation factor
	wf = 1.2;
	
	//Gauss-Siedel
	//initialize for k=0;
	for(i=0;i<=xmax;i++){
		for(j=0;j<=ymax;j++){
			if(i==0 || j==0 || i==xmax || j==ymax){
				T[0][j] = T2;
				T[xmax][j] = T4;
				T[i][0] = T1_nd;
				T[i][ymax] = T3_nd;
				T[0][0] = (T1_nd+T2)/2;
				T[0][ymax] = (T2+T3_nd)/2;
				T[xmax][0] = (T1_nd+T4)/2;
				T[xmax][ymax] = (T3_nd+T4)/2;	
			}else{
				T[i][j] = 0.0;
			}
			T_nxt[i][j] = T[i][j];
		}
	}
	
		a[0]=0, a[xmax]=0;
		b[0]=0, b[xmax]=0;
		c[0]=0, c[xmax]=0;
		d[0]=0, d[xmax]=0;
		r[0]=0, r[xmax]=0;
		s[0]=0, s[xmax]=0;
	
	//iteration loop
	for(k=1;k<=kmax;k++){
		diff = 0.0;
		summ = 0.0;
		for(j=1;j<ymax;j++){
			for(i=1;i<xmax;i++){
				a[i]=wf;
				b[i]=-2*(1+B*B);
				c[i]=wf;
				d[i]=(-1*wf*B*B)*T[i][j-1] + (-2*(1+B*B)*(1-wf))*T[i][j] + (-1*wf*B*B)*T[i][j+1];
			}
			a[1] = 0.0;
			c[xmax-1] = 0.0;
			d[1] = (-1*wf*B*B)*T[1][j-1] + (-2*(1+B*B)*(1-wf))*T[1][j] + \
			(-1*wf*B*B)*T[1][j+1] - (wf)*T[0][j];
			d[xmax-1] = (-1*wf*B*B)*T[xmax-1][j-1] + (-2*(1+B*B)*(1-wf))*T[xmax-1][j] + \
			(-1*wf*B*B)*T[xmax-1][j+1] - (wf)*T[xmax][j];
			
			//solve the tridiagonal system of equations
			tri(&a, &b, &c, &d, &r, &s);
				
			//update T
			T[xmax-1][j]=s[xmax-1];
			for(i=xmax-2;i>=1;i--){
				Told = T[i][j];

				T[i][j]=s[i]-r[i]*T[i+1][j];

				diff = diff + abso(T[i][j]-Told);
				summ = summ + abso(Told);
			}
		}

		if((diff/summ)<conv){
			printf("converged\tk=%d\n",k);
			break;
		}
}
	//calculate error
	error = f_error(&T,&T_an);

	//free memory
	for(i=0;i<=xmax;i++){
		free(T[i]);
	}
	free(T);

	for(i=0;i<=xmax;i++){
		free(T_an[i]);
	}
	free(T_an);
	
	for(i=0;i<=xmax;i++){
		free(T_nxt[i]);
	}
	free(T_nxt);

	free(x);
	free(y);
	
	free(a);
	free(b);
	free(c);
	free(d);
	free(r);
	free(s);

	return 0;
}
