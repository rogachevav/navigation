void integrate_V(double *v,double *OmB,double *u,double *f,double *g_m,double dt,int N);
void integrate_V(double *v,double *OmB,double *u,double *f,double *g_m,double dt,int N)
{
    int i;
    double right[3],help1[9];
    double MatOmB_m[9],Matu[9];
    makeMatrix(OmB,MatOmB_m);
    makeMatrix(u,Matu);
    for(i=0;i<9;i++) help1[i] = MatOmB_m[i] + 2.*Matu[i];
    multMV(help1,v,right,N);
    //for(i=0;i<3;i++) right[i] = right[i] + f[i] + g_m[i];
    for(i=0;i<3;i++) v[i] = v[i] + dt*(right[i]+f[i]+g_m[i]);

    //printf("%lf  %lf  %lf\n",f[0]+g_m[0],f[1]+g_m[1],f[2]+g_m[2]);
}
