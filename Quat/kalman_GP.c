#define sqr(x) ((x)*(x))
void Kalman_GP(double *Am,double *Ht,double *S_0,double *S,double *x_,double *xp, int N, int mes,double *R,double *Q,double dt,double *zt);
void Kalman_GP(double *Am,double *Ht,double *S_0,double *S,double *x_,double *xp, int N, int mes,double *R,double *Q,double dt,double *zt)
{
    int i,k;
    int iter;
    int kk;
    double h[N];
    double e0[N];
    double z;
    double d[N+1];
    double b[N],c[N];
    double S_[N*N];
    double f[N];
    double hv[N];
    for(i = 0;i <N*N; i++) S[i] = S_[i];
    for(i = 0;i < N; i++) xp[i] = x_[i];
    //Коррекция
    for(iter = 0; iter < mes; iter ++)
    {
        for(i = 0;i <N*N; i++) S_[i] = S[i];
        for(i = 0;i < N; i++) x_[i] = xp[i];
        for(i = 0;i < N; i++) h[i] = Ht[iter*mes + i];
        for(i = 0;i < N; i++) {e0[i] = 0; d[i+1] = 0;}
        d[0] = R[iter*mes + iter];
        z = zt[iter];
        multMV(S_,h,f,N);
        for(k=0; k < N; k++)
        {
            kk=k+1;
            d[kk] = d[kk-1] + sqr(f[kk]);
            b[k] = sqrt(d[kk-1]/d[kk]);
            c[k] = f[k]/sqrt(d[kk-1]*d[kk]);
            for(i = 0;i < N; i++) S[i*mes + k] = S_[i*mes+k] * b[k] - e0[i]*c[k];//переставить
            for(i = 0;i < N; i++) e0[i] += S_[i*mes+k] * f[k];
        }
        for(i = 0;i < N; i++) hv[i] = h[i]*x_[i];
        for(i = 0;i < N; i++) xp[i] = x_[i] + e0[i]/d[N]*(z - hv[i]);
    }
}
