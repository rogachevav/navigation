void multMV(double *a,double *b,double *c, int N);
void multMV(double *a,double *b,double *c,int N)
{
//    int i,j;
    c[0]=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    c[1]=a[3]*b[0]+a[4]*b[1]+a[5]*b[2];
    c[2]=a[6]*b[0]+a[7]*b[1]+a[8]*b[2];
//    for(i=0;i<3;i++)
//    {
//        c[i] = 0;
//        for(j=0;j<3;j++)
//        {
//
//            c[i] = c[i] + a[i*3+j]*b[j];
//        }
//    }
}
void multMM(double *a,double *b,double *c, int N);
void multMM(double *a,double *b,double *c,int N)
{
    int i,j,k;
   for( i = 0; i < N; i++)
        for( j = 0; j < N; j++)
        {
            c[i*N+j] = 0;
            for( k = 0; k < N; k++)
                c[i*N+j] += (a[i*N+k] * b[k*N+j]);
        }
}
void transp(double *S,double *St);
void transp(double *S,double *St)
{
    int i,j;
    int N=3;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            St[i*N+j] = 0;
            St[i*N+j] = S[j*N+i];
        }
    }
}
void makeMatrix(double *v,double *M_v);
void makeMatrix(double *v,double *M_v)
{
    M_v[0]=0;     M_v[1]=v[2]; M_v[2]=-v[1];
    M_v[3]=-v[2];  M_v[4]=0;    M_v[5]=v[0];
    M_v[6]=v[1];   M_v[7]=-v[0]; M_v[8]=0;
}
void matrixB(double l,double phi,double *b,int N);
void matrixB(double l,double phi,double *b,int N)
{
    b[0*N+0] = -cos(l)*sin(phi);//cos(l)*cos(phi);//b[1*N+1]*b[2*N+2] - b[1*N+2]*b[2*N+1];
    b[0*N+1] = -sin(l)*sin(phi);//-sin(l);//cos(l);
    b[0*N+2] = cos(phi);//cos(l);//0;

    b[1*N+0] = cos(l)*cos(phi);//sin(l)*cos(phi);//b[2*N+1]*b[0*N+2] - b[2*N+2]*b[0*N+1];
    b[1*N+1] = sin(l)*cos(phi);//0;//-sin(l)*sin(phi);
    b[1*N+2] = sin(phi);//-cos(l)*sin(phi);//cos(phi);

    b[2*N+0] = cos(phi);//sin(phi);//b[0*N+1]*b[1*N+2] - b[0*N+2]*b[1*N+1];
    b[2*N+1] = sin(phi);//-sin(l)*sin(phi);//sin(l)*cos(phi);
    b[2*N+2] = 0;//cos(phi);//sin(phi);

}
void matrixD(double *ang,double *d,int N);
void matrixD(double *ang,double *d,int N)
{
    double p,ga,t;
    p=ang[0];ga=ang[1];t=ang[2];
    d[0] = cos(p)*cos(ga) - sin(p)*sin(t)*sin(ga);//-sin(p)*cos(t);
    d[1] = sin(p)*cos(ga) + cos(p)*sin(t)*sin(ga);//cos(p)*sin(ga) + sin(p)*sin(t)*cos(ga);
    d[2] = -cos(t)*sin(ga);//cos(p)*cos(ga) - sin(p)*sin(t)*sin(ga);

    d[3] = -sin(p)*cos(t);//cos(p)*cos(t);
    d[4] = cos(p)*cos(t);//sin(p)*sin(ga) - cos(p)*sin(t)*cos(ga);
    d[5] = sin(t);//sin(p)*cos(ga) + cos(p)*sin(t)*sin(ga);

    d[6] = cos(p)*sin(ga) + sin(p)*sin(t)*cos(ga);//sin(t);
    d[7] = sin(p)*sin(ga) - cos(p)*sin(t)*cos(ga);//cos(t)*cos(ga);
    d[8] = cos(t)*cos(ga);//-cos(t)*sin(ga);
}
