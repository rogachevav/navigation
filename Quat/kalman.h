//#define CHECK_OUTPUT
//#define ALL_STEPS
void mult_MM(double *a,double *b,double *c,int m,int n,int v);
void makeVV_M(double *f,double *m,int N);
void Kalman_Filter(double *A,double *H,double *S_0,double *S,double *x_,double *xp, int N, int mes,double *R,double *Q,double dt,double *z);
void Kalman_Filter(double *A,double *H,double *S_0,double *S,double *x_,double *xp, int N, int mes,double *R,double *Q,double dt,double *z)
{
    int i,j;
    double E[N*N],Fm[N*N],St[N*N];
    double W[N*N];
    double HNN[N*N];
    double f1[N];
    double alpha1;
    double R1,z1,H_sk;
    double k1[N];
    double H_mat[N*N],H_mat2[N*N];
    int iteration;
    double h1[N];
    //формирование переходной матрицы F = E + A*dt
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(i==j) { Fm[i*N+j] = A[i*N+j]*dt + 1.; E[i*N+j]=1.;}
            else {Fm[i*N+j] = A[i*N+j]*dt; E[i*N+j]=0;}
        }
    }
    for(i=0;i<N*N;i++) S[i]=S_0[i];
    for(i=0;i<N;i++) xp[i]=x_[i];

     //----------------------------------
     /*-----Этап коррекции-----*/
     for(iteration = 0; iteration < mes; iteration++ )
     {
            R1 = R[iteration*mes + iteration];
            z1=z[iteration];
            for(i=0;i<N;i++) h1[i] = H[iteration*N + i];
            for(i=0;i<N;i++) x_[i]=xp[i];
            for(i=0;i<N*N;i++) S_0[i]=S[i];
     #ifdef ALL_STEPS
     printf("ITERATION %d ===============================================\n",iteration);
     #endif // ALL_STEPS

     #ifdef ALL_STEPS
     printf("h1--->");
     for(i=0;i<N;i++) printf(" %0.1e",h1[i]);
     printf("\n");
     #endif // ALL_STEPS
     //Шаг 1 Формирование  f1
     transp(S_0,St,N);
     #ifdef ALL_STEPS
     printf("Matrix St\n");
     for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%0.2e ",St[i*N+j]);
            }
        printf("\n");
    }
    printf("------------------------------------------------------------------------------\n");
     #endif // ALL_STEPS
     multMV(St,h1,f1,N);
     #ifdef ALL_STEPS
     printf("f1--->");
     for(i=0;i<N;i++) printf(" %0.1e",f1[i]);
     printf("\n");
     #endif // ALL_STEPS
     //----------------------------------
     //Шаг 2 Формирование aplha
     alpha1 = 0;
     for(i=0;i<N;i++) alpha1+=f1[i]*f1[i];
     alpha1 += R1;
     //alpha1=1./alpha1;
     #ifdef ALL_STEPS
     printf("alpha---> %0.1e\n",alpha1);
     #endif // ALL_STEPS
     //----------------------------------
     //Шаг 3 Формирование  K
     multMV(S_0,f1,k1,N);
     for(i=0;i<N;i++) k1[i] = k1[i]/alpha1;
     #ifdef ALL_STEPS
     printf("k1--->");
     for(i=0;i<N;i++) printf(" %0.1e",k1[i]);
     printf("\n");
     #endif // ALL_STEPS

     //----------------------------------
     //Шаг 4 Получение S+
     makeVV_M(f1,H_mat,N);
     for(i=0;i<N*N;i++) H_mat[i] = E[i] - H_mat[i]/alpha1;
          #ifdef ALL_STEPS
     printf("Matrix sqrt(E - f*f^T/alpha)\n");
     for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%0.2e ",H_mat[i*N+j]);
            }
        printf("\n");
    }
    printf("------------------------------------------------------------------------------\n");
     #endif // ALL_STEPS
     rsb(H_mat,H_mat2,N);
     multMM(S_0,H_mat2,S,N);
     #ifdef ALL_STEPS
     printf("Matrix S+\n");
     for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%0.2e ",S[i*N+j]);
            }
        printf("\n");
    }
    printf("------------------------------------------------------------------------------\n");
     #endif // ALL_STEPS
     //----------------------------------
     //Шаг 5 Получение x+
     H_sk = 0;
     for(i=0;i<N;i++) H_sk += h1[i]*x_[i];
     z1=z1-H_sk;
     for(i=0;i<N;i++) xp[i] = x_[i] + k1[i]*z1;
     //----------------------------------
     }

     /*-----Этап прогноза-----*/
     //Шаг 1 Получение x_
     multMV(Fm,xp,x_,N);
     //----------------------------------
     //Шаг 2 Получение S_
     multMM(Fm,S,W,N);
     transp(W,St,N);
     multMM(W,St,HNN,N);
     for(i=0;i<N*N;i++) HNN[i] = HNN[i] + Q[i];
     rsb(HNN,S_0,N);
     #ifdef CHECK_OUTPUT
     printf("Matrix S\n");
     for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
            {
                printf("%0.2e ",S_0[i*N+j]);
            }
        printf("\n");
    }
    printf("------------------------------------------------------------------------------\n");
    printf("Vector xp\n");
    for(j=0;j<N;j++)
    {
        printf("%0.2e ",xp[j]);
    }
    printf("\n");
    printf("------------------------------------------------------------------------------\n");
     #endif // CHECK_OUTPUT

     //----------------------------------
     //Ну вроде все
}

void rsb(double *p,double *s,int m);
void rsb(double *p,double *s,int m)
{
   int i,j,k,j1,jj,ij,jm,im;
   double ds,y,c;
   j1=m*m-1;
   *(s+j1)=sqrt(*(p+j1));
   ds=1.0/(*(s+j1));
   for(i=0;i < m-1;i++)
   {
       j1=(i+1)*m-1;
       *(s+j1)=*(p+j1)*ds;
   }
   if(m>2)
   {
       for(j=m-2;j>0;j--)
       {
           y=0.0;
           j1=j+1;
           jm=j*m;
           for(k=j1;k < m;k++)
           {
                c=*(s+jm+k);
                y+=c*c;
		   }
           jj=jm+j;
           *(s+jj)=sqrt(*(p+jj)-y);
           ds=1.0/(*(s+jj));
           for(i=j-1;i >= 0;i--)
           {
              im=i*m;
              y=0.0;
              for(k=j1;k < m;k++) y+=*(s+im+k)**(s+jm+k);
              ij=im+j;
              *(s+ij)=(*(p+ij)-y)*ds;
           }
	   }
   }
   y=0.0;
   for(i=1;i < m;i++)
       {
            c=*(s+i);
            y+=c*c;
       }
       *s=sqrt(*p-y);
}

void mult_MM(double *a,double *b,double *c,int m,int n,int v)
{
    int i,j,k;
    for( i = 0; i < m; i++)
        {
            for( j = 0; j < v; j++)
            {
                c[i*v+j] = 0;
                for( k = 0; k < n; k++) c[i*v+j] +=  (a[i*n+k] * b[k*v+j]);
            }
        }
}
void Inverse(double *a,int n);
void Inverse(double *a,int n)
{
    int i,j;
    double det, mat[9],mata[9];
    det =((a[0]*(a[8]*a[4]-a[7]*a[5]))-(a[1]*(a[8]*a[3]-a[6]*a[5]))+(a[2]*(a[7]*a[3]-a[6]*a[4])));
    mat[0] = pow(-1,2)*(a[8]*a[4]-a[7]*a[5]);
    mat[3] = pow(-1,3)*(a[8]*a[3]-a[6]*a[5]);
    mat[6] = pow(-1,4)*(a[7]*a[3]-a[6]*a[4]);

    mat[1] = pow(-1,3)*(a[8]*a[1]-a[7]*a[2]);
    mat[4] = pow(-1,4)*(a[8]*a[0]-a[6]*a[2]);
    mat[7] = pow(-1,5)*(a[7]*a[0]-a[6]*a[1]);

    mat[2] = pow(-1,4)*(a[5]*a[1]-a[4]*a[2]);
    mat[5] = pow(-1,5)*(a[5]*a[0]-a[3]*a[2]);
    mat[8] = pow(-1,6)*(a[4]*a[0]-a[3]*a[1]);
    for(i=0;i<9;i++) a[i] = mat[i]/det;
}
void make_B_m(double *b,double l,double phi);
void make_B_m(double *b,double l,double phi)
{
    b[0*3+0] = -sin(l);
    b[0*3+1] = cos(l);
    b[0*3+2] = 0;

    b[1*3+0] = -cos(l)*sin(phi);
    b[1*3+1] = -sin(l)*sin(phi);
    b[1*3+2] = cos(phi);

    b[2*3+0] = cos(l)*cos(phi);
    b[2*3+1] = sin(l)*cos(phi);
    b[2*3+2] = sin(phi);

}
void makeVV_M(double *f,double *m,int n)
{
    int i,j;
    for( i = 0; i < n; i++)
    {
        for( j = 0; j < n; j++) m[i*n+j] = f[i]*f[j];
    }
}
