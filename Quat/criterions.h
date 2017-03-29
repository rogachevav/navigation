#define sqr(x) ((x)*(x))
double norm(double *m,int k);
void krit_function(double *f,double *x,double *krit_av,double *sigma_sum,double *ix1,double *ix2,double *ix3,int N,int w,int *check,double *counst_for_krit);
void krit_function(double *f,double *x,double *krit_av,double *sigma_sum,double *ix1,double *ix2,double *ix3,int N,int w,int *check,double *counst_for_krit)
{
    int i;
    double krit[4],sigma[3],krit3_slag[3];
             krit[0] = sqrt(sqr(f[0])+sqr(f[1])+sqr(f[2]));//n1,n1av,g1s,g1
             krit[1]  = sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2]));
             krit[2]  = sqrt(fabs(x[0]*f[0])+fabs(x[1]*f[1])+fabs(x[2]*f[2]));
             for(i=0;i<3;i++) krit_av[i] = krit_av[i]*w/(w+1.) + krit[i]/(w+1.);
             for(i=0;i<3;i++) sigma_sum[i] = sigma_sum[i]*w/(w+1.) + sqr(krit[i]-krit_av[i])/(w+1.);
             for(i=0;i<3;i++) sigma[i]=sqrt(sigma_sum[i]);

             mass(ix1,x[0],N);
             mass(ix2,x[1],N);
             mass(ix3,x[2],N);

             krit3_slag[0] = norm(ix1,N);
             krit3_slag[1] = norm(ix2,N);
             krit3_slag[2] = norm(ix3,N);
             krit[3] = sqrt(sqr(krit3_slag[0])+sqr(krit3_slag[1])+sqr(krit3_slag[2]));
             for(i=0;i<3;i++)
            {
                if(fabs(krit[i]-krit_av[i])>counst_for_krit[i]*sigma[i]) {check[i]=1; }
                else {check[i]=-5;}
            }
            if(krit[3]>counst_for_krit[3]) check[3]=1;
            else  check[3]=-5;
}
void krit4_function(int len,double *st1,double *st2,double *st3,double *A,double *f,int *check,double *counst_for_krit,double gg);
void krit4_function(int len,double *st1,double *st2,double *st3,double *A,double *f,int *check,double *counst_for_krit,double gg)
{
    double krit_int;
    int i;
      for(i=0;i<len-1;i++)
{
    st1[i] = st1[i+1];
    st2[i] = st2[i+1];
    st3[i] = st3[i+1];
}
st1[len-1]= A[0]*f[0]+A[3]*f[1]+A[6]*f[2];
st2[len-1]= A[1]*f[0]+A[4]*f[1]+A[7]*f[2]-gg;
st3[len-1]= A[2]*f[0]+A[5]*f[1]+A[8]*f[2];
krit_int = sqrt(sqr(norm(st1,len))+sqr(norm(st2,len))+sqr(norm(st3,len)));
if(krit_int > counst_for_krit[4]) check[4]=1;
else check[4]=-5;
}

void mass(double *ix1, double a,int N);
void mass(double *ix1, double a,int N)
{
    int q;
    for(q=1;q<N;q++)
    {
        ix1[q-1]=ix1[q];
    }
    ix1[N-1] = a;
}

double norm(double *m,int k)
{
    double S=0,dt=0.004;
    int i;
    for(i=0;i<k-1;i++)
    {
        S = S + 0.5*(m[i] + m[i+1])*dt;
    }
    return S;
}
void allig(double a,double *b,int n);
void allig(double a,double *b,int n)
{
    *b = (*b)*n/(n+1.) + a/(n+1);
}
void makaAmatrix(double* A,double *quat_int);
void makaAmatrix(double* A,double *quat_int) //Из связаной в географическую
{
    int i;
    for(i=0;i<3;i++) A[i*3+i] = 2*quat_int[0]*quat_int[0] + 2*quat_int[i+1]*quat_int[i+1] - 1;
    A[1] = 2*quat_int[1]*quat_int[2]-2*quat_int[0]*quat_int[3];
    A[2] = 2*quat_int[1]*quat_int[3]+2*quat_int[1]*quat_int[2];
    A[5] = 2*quat_int[3]*quat_int[2]-2*quat_int[0]*quat_int[1];
    A[3] = 2*quat_int[1]*quat_int[2]+2*quat_int[0]*quat_int[3];
    A[6] = 2*quat_int[1]*quat_int[3]-2*quat_int[1]*quat_int[2];
    A[7] = 2*quat_int[3]*quat_int[2]+2*quat_int[0]*quat_int[1];

}
