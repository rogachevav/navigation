#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "kalman.h"`
#define nn 2
void transp(double *S,double *St);
void multMV(double *a,double *b,double *c,int n);
void multMM(double *a,double *b,double *c,int n);
//void data(double *x_k0,double *x_k1);
void measurement(double *x,double dt,double *zet,double *H);
double white_noise(double mean, double variance);
int main(void)
{
    FILE *fp;
    fp = fopen("Rez.txt","w");
    int i,j;
    double dt;//шаг по времени
    dt = 1e0;
    //--------Задание системы---------
    double *A,*q,*Q;
    A = (double*)malloc((nn*nn)*sizeof(double));
    q = (double*)malloc(nn*sizeof(double));
    Q = (double*)malloc((nn*nn)*sizeof(double));
    //Матрица А
    for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            if(i==j) A[i*nn+j]=0;
            else A[i*nn+j]=1;
            A[2]=-1;

        }
    }
     //Матрица Q
       for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            if(i==j) Q[i*nn+j]=1e-4;
            else Q[i*nn+j]=0;
        }
    }
     //----------------
    double *H,r,R[1],z[1];
    H = (double*)malloc(nn*sizeof(double));
    H[0]=1; H[1]=1;
    R[0] = 1e-4;
    double P_0[nn*nn],P[nn*nn];
    double x_0[nn],x[nn];
    double S_0[nn*nn],S[nn*nn],St[nn*nn];
    //Матрица P(0)
    for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            if(i==j) P_0[i*nn+j] = 1;
            else P_0[i*nn+j] = 0;
        }
    }

    for(i=0;i<nn;i++) x_0[i] = 0;
    rsb(P_0,S_0,nn);
    double t;
    t=0; //время работа фильтра
    double x_real[nn];
    x_real[0]=1;
    x_real[1]=2;
    while (t < 50)
    {
        fprintf(fp,"%lf ",t);
        measurement(x_real,dt,z,H);
//        for(i=0;i<N;i++) printf("%lf ",x_real[i]);
//        printf("   |   ");

        Kalman_Filter(A,H,S_0,S,x_0,x,nn,1,R,Q,dt,z);

        for(i=0;i<nn;i++) printf("%lf ",x[i]-x_real[i]);
        printf("--------------   \n   ");

        for(i=0;i<nn;i++)
        {
            fprintf(fp,"%lf ",x[i]);
        }
        for(i=0;i<nn;i++)
        {
            fprintf(fp,"%lf ",x_real[i]);
        }
        transp(S,St);
        multMM(S,St,P,nn);
        fprintf(fp,"%lf %lf %lf",P[0],P[1*nn+1],P[1]);
        fprintf(fp,"\n");
        t = t+dt;
    }

    return 0;
}
double white_noise( double mean, double variance)
{
    double noise;
    noise = 2*((rand()/((double)RAND_MAX)) - 0.5)*sqrt(3);
    return mean+variance*noise;
}
void measurement(double *x,double dt,double *z,double *H)
{
    double x_next[nn];
    int i;
    double mes;
    mes=0;
    x_next[0] = x[0] +  (x[1] + white_noise(0,1e-2))*dt;
    x_next[1] = x[1] +  (-x[0] + white_noise(0,1e-2))*dt;
    for(i=0;i<nn;i++) x[i] = x_next[i];
    for(i=0;i<nn;i++) mes+=H[i]*x[i];
    z[0] = mes + white_noise(0,1e-2);
}
void transp(double *S,double *St)
{
    int i,j;
    for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            St[i*nn+j] = S[j*nn+i];
        }
    }
}

void multMM(double *a,double *b,double *c,int n)
{
    int i,j,k;
   for( i = 0; i < n; i++)
        for( j = 0; j < n; j++)
        {
            c[i*nn+j] = 0;
            for( k = 0; k < n; k++)
                c[i*n+j] += (a[i*n+k] * b[k*n+j]);
        }
}
void multMV(double *a,double *b,double *c,int n)
{
    c[0]=a[0]*b[0]+a[1]*b[1];
    c[1]=a[2]*b[0]+a[3]*b[1];
}
