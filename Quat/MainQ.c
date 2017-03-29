#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#define sqr(x) ((x)*(x))
#define pi            3.14159265358979323846
#define g             9.780327
#define time_skip 0
#define time_aligment 60
#include "QuaternionFunctions.h"
#define Q_d 1e-1
#include "kalman.h"
//#include "kalman_GP.c"
//#define CHECK_INPUT
#define COUNTER
//#define PRINT_FILTER
void transpose_any(double *a,double *a_t,int u,int v);
void Alignment(double *f,double *fav,double n);
void Ali_ang(double *fav,double *angles);
void make_L(double* a, double* L);
void Integrate_Matrix(double* A,double* v,double dt);
void Intergate_V(double* vec,double* V_m,double* g_m,double* f,double* L,double dt);
void makeMatrix(double *v,double *M_v);
void transp(double *S,double *St,int N);
void multMM(double *a,double *b,double *c, int N);
void multMV(double *a,double *b,double *c, int N);
void Nav(double *Qo,double *p,double *x,double dt,double *Rezn);
void makeQuatMat(double *qw,double *a);
void Make_quat_allig(double *q,double *a);
void Kalman_GP(double *Am,double *Ht,double *S_0,double *S,double *x_,double *xp, int N, int mes,double *R,double *Q,double dt,double *zt);
int main(void)
{
    FILE *fp;
    FILE *fpp;
    fp=fopen("newData.txt","r");
    fpp=fopen("RR1.txt","w");
    double b,dt;
    int i,j;
    dt=0.004;
    b=0;
    double x[3],f[3],temperature;
    double x_av[3],f_av[3],koef=pi/180;
    double Mwx,Mwy,Mwz,Mfx,Mfy,Mfz,Afx,Afy,Afz;
    Mwx=1.991; Mfx=2.002; Afx=-254.6;
    Mwy=1.993; Mfy=2.003; Afy=261.75;
    Mwz=2.0;   Mfz=2.002; Afz=-261.25;
    double angles[3];//0 - тангаж, 1 - крен, 2 - курс.
    double al_count; al_count = 0;
    double lyam_m,phi_m,h_m;
    h_m = 144; phi_m = 55.425395 *pi/180; lyam_m = 37.314170*pi/180;
    double u_m[3],g_m[3],u,g_0,e,a_sa;
    g_0 = 9.78030; e = 6.6943799901413e-3; a_sa = 6378137; u=7.2921158553e-5;
    double RE_m,RN_m,OmB_m[3];
    double L[9],Az[9],Ax[9],Ax_T[9],V_m[3];
    for(i=0;i<3;i++) {OmB_m[i] = 0; x_av[i] = 0; f_av[i] = 0;V_m[i] = 0;}
    double vec[3];
    for(i=0;i<9;i++) Ax[i] = 0;
    Ax[0]=1;Ax[4]=1;Ax[8]=1;
    int check;
    check = 0;
    //Кватернионы
    double p_ax[4],p_az[4];
    double q_ax_0[4],q_az_0[4];
    p_ax[0]=1; p_az[0]=1; q_ax_0[0]=1;
    for(i=1;i<4;i++) {p_ax[i]=0;p_az[i]=0;q_ax_0[i]=0;}
    double r_ax[4],r_az[4];
    //Для коррекции
    FILE *gps; gps=fopen("DataGPS.txt","r");
    FILE *filter; filter=fopen("Filter.txt","w");
    int N_m,mes;
    N_m=15;
    mes=3;
    double A_c[N_m*N_m],Om_Sh;
    Om_Sh = 1.54e-6;
    double H[N_m*mes];
    double P_0[N_m*N_m],P[N_m*N_m],S_0[N_m*N_m],S[N_m*N_m],St[N_m*N_m],x_0[N_m],xp[N_m],R[mes*mes],Q[N_m*N_m];
    for(i=0;i<N_m*N_m;i++) P_0[i]=0;
    P_0[0]=25; P_0[N_m+1]=25; P_0[N_m*2+2]=25;
    P_0[N_m*3+3]=1e-4;    P_0[N_m*4+4]=1e-4;    P_0[N_m*5+5]=1e-4;
    P_0[N_m*6+6]=25e-6;   P_0[N_m*7+7]=25e-6;   P_0[N_m*8+8]=1e-2;
    P_0[N_m*9+9]=4e-8;    P_0[N_m*10+10]=4e-8;  P_0[N_m*11+11]=4e-8;
    P_0[N_m*12+12]=25e-4; P_0[N_m*13+13]=25e-4; P_0[N_m*14+14]=25e-4;
    double P_proc[N_m*N_m];
    for(i=0;i<N_m*N_m;i++) P_proc[i]=0;
    for(i=0;i<mes*mes;i++) R[i]=0;
    for(i=0;i<mes;i++)
    {
        for(j=0;j<mes;j++)
            if(i==j) R[i*mes+j]=25/a_sa;
    }
    for(i=0;i<N_m*N_m;i++) Q[i]=0;
    for(i=0;i<N_m;i++) Q[i*N_m+i]=Q_d;
    P_0[0]=25; P_0[N_m+1]=25; P_0[N_m*2+2]=25;
    P_0[N_m*3+3]=1e-4;    P_0[N_m*4+4]=1e-4;    P_0[N_m*5+5]=1e-4;
    P_0[N_m*6+6]=25e-6;   P_0[N_m*7+7]=25e-6;   P_0[N_m*8+8]=1e-2;
    P_0[N_m*9+9]=4e-8;    P_0[N_m*10+10]=4e-8;  P_0[N_m*11+11]=4e-8;
    P_0[N_m*12+12]=25e-4; P_0[N_m*13+13]=25e-4; P_0[N_m*14+14]=25e-4;



    rsb(P_0,S_0,N_m);

    double S_cons_mat[N_m*N_m];
    for(i=0;i<N_m*N_m;i++) S_cons_mat[i]=S_0[i];
    //что тут
    for(i=0;i<N_m;i++) x_0[i] = 0;
    double z[mes];
    int gps_count;
    gps_count=0;
    double GPS_D_Mass[14];
    double gps_lyam,gps_phi,gps_h,dc1,cc1;
    double Cor_ang[9];
    double B_m[9];
    int count_filter;
    count_filter = 0;
    double dV[3],d_angles[3];
    double d_phi,d_lyam,d_h;

    while(b<747)
    {

        fscanf(fp,"%lf",&temperature);
        fscanf(fp,"%lf",&x[0]);
        fscanf(fp,"%lf",&x[1]);
        fscanf(fp,"%lf",&x[2]);
        fscanf(fp,"%lf",&f[0]);
        fscanf(fp,"%lf",&f[1]);
        fscanf(fp,"%lf",&f[2]);

        x[0] = x[0]/(3600.*Mwx);//201
        x[1] = x[1]/(3600.*Mwy);
        x[2] = x[2]/(3600.*Mwz);
        f[0] = (f[0] - Afx)/Mfx;
        f[1] = (f[1] - Afy)/Mfy;
        f[2] = (f[2] - Afz)/Mfz;

        x[0] = koef*(x[0]);
        x[1] = koef*(x[1]);
        x[2] = koef*(x[2]);

        if(gps_count%250==0)
             {
                for(i=0;i<14;i++)  fscanf(gps,"%lf",&GPS_D_Mass[i]);
             }
        #ifdef COUNTER
        if(gps_count%2500==0) printf("%lf\n",b);
        #endif // COUNTER
        dc1=modf(GPS_D_Mass[1]*0.01,&cc1);
        gps_phi = cc1 + dc1*100/60;
        dc1=modf(GPS_D_Mass[3]*0.01,&cc1);
        gps_lyam = cc1 + dc1*100/60;
        gps_h=GPS_D_Mass[8]+GPS_D_Mass[10];
        //printf("%lf %lf %lf\n",gps_phi,gps_lyam,gps_h);
        gps_phi = gps_phi*koef;
        gps_lyam = gps_lyam*koef;
        gps_count++;

        //ВЫСТАВКА
        if(b>=time_skip && b<=time_skip+time_aligment)
        {
            Alignment(f,f_av,al_count);
            Alignment(x,x_av,al_count);
            al_count++;
            Ali_ang(f_av,angles);
            b = b+dt;
           // fprintf(fpp,"%lf  %lf  %lf  %lf 0 0 0  %lf  %lf\n",b,angles[0],angles[1],angles[2],lyam_m,phi_m);
            make_L(angles,Az);

//---------------------------------------------------------
            //makeMatQuat(Az,q_az_0);
            Make_quat_allig(q_az_0,angles);
//---------------------------------------------------------
            for(i=0;i<9;i++) L[i]=Az[i];
            continue;
        }
        if(b > time_skip+time_aligment)
        {
            x[0] = (x[0] - x_av[0]);
            x[1] = (x[1] - x_av[1]);
            x[2] = (x[2] - x_av[2]);
            u_m[0] = 0; u_m[1] = u*cos(phi_m); u_m[2] = u*sin(phi_m);//Ux
            g_m[0] = 0; g_m[1] = 0; g_m[2] = -g_0*(1.+0.75*e*sqr(sin(phi_m)));
            RE_m = a_sa/sqrt(1-e*sqr(sin(phi_m)))+h_m; RN_m = a_sa*(1-e)/((1-e*sqr(sin(phi_m)))*sqrt(1-e*sqr(sin(phi_m))))+h_m;
            OmB_m[0] = -V_m[1]/RN_m; OmB_m[1] = V_m[0]/RE_m; OmB_m[2] = V_m[0]*tan(phi_m)/RE_m;

            phi_m = phi_m + dt*(V_m[1]/RN_m); lyam_m = lyam_m + dt*(V_m[0]/(RE_m*cos(phi_m))); h_m = h_m + dt*V_m[2];
            for(i=0;i<3;i++) vec[i] = OmB_m[i] + 2.*u_m[i];
            intV_G(vec,V_m,g_m,f,L,dt);
            //интегрирование Ax Az
//            Integrate_Matrix(Az,x,dt);
//--------------------------------------------------------------------------
            Nav(q_az_0,p_az,x,dt,r_az);
            makeQuatMat(r_az,Az);
            for(i=0;i<3;i++) vec[i] = OmB_m[i] + u_m[i];
            Nav(q_ax_0,p_ax,vec,dt,r_ax);
            makeQuatMat(r_ax,Ax);
//            Integrate_Matrix(Ax,vec,dt);
//---------------------------------------------------------------------------
            transp(Ax,Ax_T,3);
            multMM(Az,Ax_T,L,3);
            if(check < 1000) check++;
//            if(check == 1000 || check==1)
//            {
//                printf("AA %lf  %lf  %lf \n",x_av[0],x_av[1],x_av[2]);
//                for(i=0;i<3;i++)
//                {
//                    printf("| ");
//                    for(j=0;j<3;j++)
//                    {
//                        printf("%lf ",Ax[i*3+j]);
//                    }
//                    printf(" |\n");
//                }
//                check++;
//            }

            angles[0] = atan2(L[5],sqrt(sqr(L[2])+sqr(L[8])));
            angles[1] = -atan2(L[2],L[8]);
            angles[2] = -atan2(L[3],L[4]);

            //КОРРЕКЦИЯ
            count_filter++;
           // if(count_filter==1)
          //if(b > 60 && b < 70)//&& gps_count%250==0)
            {
//make_B_m(B_m,lyam_m,phi_m);
z[0] = (lyam_m - gps_lyam)*cos(gps_phi);
z[1] = phi_m - gps_phi;
z[2] = h_m - gps_h;
//fprintf(filter,"%lf ",b);
//fprintf(filter,"%lf %lf %lf",gps_lyam,gps_phi,gps_h);
//fprintf(filter,"\n");
for(i=0;i<15*15;i++) A_c[i] = 0;
//Координаты
A_c[0]=V_m[2]/RE_m; A_c[1]=OmB_m[2]; A_c[2]=-OmB_m[1]; A_c[3]=1.; A_c[8]=V_m[1]; A_c[7]=-V_m[2];

A_c[N_m+0]=-OmB_m[2];  A_c[N_m+1]=V_m[2]/RN_m; A_c[N_m+2]=OmB_m[0]; A_c[N_m+4]=1.; A_c[N_m+6]=V_m[2]; A_c[N_m+8]=-V_m[0];

A_c[N_m*2+0]=OmB_m[1]-V_m[0]/RE_m;  A_c[N_m*2+1]=-OmB_m[0]-V_m[1]/RN_m; A_c[N_m*2+5]=1.; A_c[N_m*2+6]=-V_m[1]; A_c[N_m*2+7]=V_m[0];
//Скорости
A_c[15*3+0]=u_m[0]*V_m[1]/RE_m; A_c[15*3+1]=V_m[1]*u_m[1]/RN_m+V_m[2]*u_m[2]/RN_m;
A_c[15*3+4]=OmB_m[2]+2.*u_m[2]; A_c[15*3+5]=-OmB_m[1]-2.*u_m[1];
A_c[15*3+6]=V_m[1]*u_m[1]+V_m[2]*u_m[2]; A_c[15*3+7]=-V_m[1]*u_m[0]-g_0; A_c[15*3+8]=-V_m[2]*u_m[0];
A_c[15*3+9]=V_m[2]*L[1] - V_m[1]*L[2]; A_c[15*3+10]=V_m[2]*L[4] - V_m[1]*L[5];  A_c[15*3+11]=V_m[2]*L[7] - V_m[1]*L[8];
A_c[15*3+12]=L[0]; A_c[15*3+13]=L[3];  A_c[15*3+14]=L[6];

A_c[15*4+0]=-V_m[0]*u_m[0]/RE_m-V_m[2]*u_m[2]/RE_m; A_c[15*4+1]=-V_m[0]*u_m[1]/RN_m;
A_c[15*4+3]=-OmB_m[2]-2.*u_m[2]; A_c[15*4+5]=OmB_m[0]+2.*u_m[0];
A_c[15*4+6]=-V_m[0]*u_m[1]+g_0; A_c[15*4+7]=V_m[0]*u_m[0]+V_m[2]*u_m[2]; A_c[15*4+8]=-V_m[3]*u_m[1];
A_c[15*4+9]=-V_m[2]*L[0] + V_m[1]*L[2];  A_c[15*4+10]=-V_m[2]*L[3] + V_m[1]*L[5];  A_c[15*4+11]=-V_m[2]*L[6] + V_m[1]*L[8];
A_c[15*4+12]=L[1];  A_c[15*4+13]=L[4];  A_c[15*4+14]=L[7];

A_c[15*5+0]=V_m[1]*u_m[2]/RE_m; A_c[15*5+1]=-V_m[0]*u_m[2]/RN_m; A_c[15*5+2]=2*Om_Sh;
A_c[15*5+3]=OmB_m[1]+2.*u_m[1]; A_c[15*5+4]=-OmB_m[0]-2.*u_m[0];
A_c[15*5+6]=-u_m[2]*V_m[0]; A_c[15*5+7]=-u_m[2]*V_m[1]; A_c[15*5+8]=u_m[0]*V_m[0]+u_m[1]*V_m[1];
A_c[15*5+9]=V_m[1]*L[0] - V_m[0]*L[1];  A_c[15*5+10]=V_m[1]*L[3] - V_m[0]*L[4];   A_c[15*5+11]=V_m[1]*L[6] - V_m[0]*L[7];
A_c[15*5+12]=L[2];  A_c[15*5+13]=L[5];  A_c[15*5+14]=L[8];

//Углы
A_c[15*6+0]=-u_m[2]/RE_m; A_c[15*6+4]=-1./RN_m; A_c[15*6+7]=u_m[2]; A_c[15*6+8]=-u_m[1];
A_c[15*6+9]=L[0]; A_c[15*6+10]=L[3]; A_c[15*6+11]=L[6];

A_c[15*7+1]=-u_m[2]/RN_m; A_c[15*7+3]=1./RE_m; A_c[15*7+6]=-u_m[2]; A_c[15*7+8]=u_m[0];
A_c[15*7+9]=L[1]; A_c[15*7+10]=L[4]; A_c[15*7+11]=L[7];

A_c[15*8+0]=(OmB_m[0]+u_m[0])/RE_m; A_c[15*8+1]=(OmB_m[1]+u_m[1])/RN_m; A_c[15*8+6]=OmB_m[1]+u_m[1];A_c[15*8+7]=-OmB_m[0]-u_m[0];
A_c[15*8+9]=L[2]; A_c[15*8+10]=L[5]; A_c[15*8+11]=L[8];
for(i=0;i<N_m*mes;i++) H[i] = 0;
H[0]=1./RN_m; H[15+1]=1./RE_m; H[15*2+2]=1;

//for(i=0;i<N_m;i++) x_0[i] = 0;
//for(i=0;i<N_m*N_m;i++) S_0[i]=S_cons_mat[i];
#ifdef CHECK_INPUT
printf("Matrix A\n");
for(i=0;i<N_m;i++)
{
    for(j=0;j<N_m;j++)
    {
        printf("%0.2e ",A_c[i*N_m + j]);
    }
    printf("\n");
}
printf("------------------------------------------------------------------------------\n");
printf("Matrix Q\n");
for(i=0;i<N_m;i++)
{
    for(j=0;j<N_m;j++)
    {
        printf("%0.2e ",Q[i*N_m + j]);
    }
    printf("\n");
}
printf("------------------------------------------------------------------------------\n");
printf("Matrix S_0\n");
for(i=0;i<N_m;i++)
{
    for(j=0;j<N_m;j++)
    {
        printf("%0.2e ",S_0[i*N_m + j]);
    }
    printf("\n");
}
printf("------------------------------------------------------------------------------\n");
printf("Matrix R\n");
for(i=0;i<mes;i++)
{
    for(j=0;j<mes;j++)
    {
        printf("%0.2e ",R[i*mes + j]);
    }
    printf("\n");
}
printf("------------------------------------------------------------------------------\n");
printf("Vector x_\n");
    for(j=0;j<N_m;j++)
    {
        printf("%0.2e ",x_0[j]);
    }
    printf("\n");
printf("------------------------------------------------------------------------------\n");
printf("Vector z\n");
    for(j=0;j<mes;j++)
    {
        printf("%0.2e ",z[j]);
    }
    printf("\n");
printf("------------------------------------------------------------------------------\n");
#endif // CHECK_INPUT
Kalman_Filter(A_c,H,S_0,S,x_0,xp,N_m,mes,R,Q,dt,z);
//for(i=0;i<N_m;i++) printf("%0.2e ",xp[i]);
//printf("\n");
transp(S,St,N_m);
multMM(S,St,P,N_m);
for(i=0;i<N_m;i++) xp[i] = x_0[i];
//Оценки вектора состояния МОДЕЛЬНЫХ уравнений
dV[0] = xp[3] + xp[8]*V_m[1] - (xp[7] - xp[0]/RE_m)*V_m[2];
dV[1] = xp[4] - xp[8]*V_m[0] + (xp[6] + xp[1]/RN_m)*V_m[2];
dV[2] = xp[5] + (xp[7] - xp[0]/RE_m)*V_m[0] - (xp[6] + xp[1]/RN_m)*V_m[1];
d_phi = xp[1]/RN_m;
d_lyam = xp[0]/(RE_m * cos(phi_m));
d_h = xp[2];
//0 - тангаж, 1 - крен, 2 - курс.
d_angles[0] = -xp[6]*cos(angles[3]) - xp[7]*sin(angles[3]);
d_angles[1] = (xp[6]*sin(angles[3]) - xp[7]*cos(angles[3]))/cos(angles[0]);
d_angles[2] = -xp[8] - d_lyam*sin(phi_m) - d_angles[1]*sin(angles[0]);

fprintf(filter,"%lf ",b);
for(i=0;i<3;i++) fprintf(filter,"%0.2e ",d_angles[i]);
for(i=0;i<3;i++) fprintf(filter,"%0.2e ",dV[i]);
fprintf(filter,"%0.2e %0.2e %0.2e\n",d_phi,d_lyam,d_h);
#ifdef PRINT_FILTER
fprintf(filter,"%lf ",b);
for(i=0;i<N_m;i++) fprintf(filter,"%lf ",xp[i]);
for(i=0;i<N_m;i++)
{
    for(j=0;j<N_m;j++)
    {
        if(i==j) fprintf(filter,"%lf ",P[i]);
    }
}
 fprintf(filter,"\n ");
#endif // PRINT_FILTER
 count_filter++;
            }//КОНЕЦ КОРРЕКЦИИ ВЫХОД P и xp
//    for(i=0;i<3;i++)
//    {
//        V_m[i] = V_m[i] - dV[i];
//        angles[i] = angles[i] - d_angles[i];
//    }
    fprintf(fpp,"%lf  %lf  %lf  %lf %lf  %lf  %lf %lf  %lf\n",b,angles[0],angles[1],angles[2],V_m[0],V_m[1],V_m[2],lyam_m,phi_m);
            b = b+dt;
    }
}
    return 0;
}
void Alignment(double *f,double *fav,double n)
{
    fav[0] = (fav[0])*n/(n+1.) + f[0]/(n+1.);
    fav[1] = (fav[1])*n/(n+1.) + f[1]/(n+1.);
    fav[2] = (fav[2])*n/(n+1.) + f[2]/(n+1.);
}
void Ali_ang(double *fav,double *angles)
{
    angles[0] = atan2(fav[1],sqrt(sqr(fav[0])+sqr(fav[2])));
    angles[1] = -atan2(fav[0],fav[2]);
    angles[2] = 0;

}
void make_L(double* a, double* L)
{
   //0 - тангаж, 1 - крен, 2 - курс.
   L[0]=cos(a[2])*cos(a[1])-sin(a[0])*sin(a[1])*sin(a[2]);   L[1]=sin(a[2])*cos(a[1])+cos(a[2])*sin(a[0])*sin(a[1]); L[2]=-cos(a[0])*sin(a[1]);
   L[3]=-sin(a[2])*cos(a[0]);                                L[4]=cos(a[2])*cos(a[0]);                               L[5]=sin(a[0]);
   L[6]=cos(a[2])*sin(a[1])+sin(a[0])*cos(a[1])*sin(a[2]);   L[7]=sin(a[2])*sin(a[1])-cos(a[2])*sin(a[0])*cos(a[1]); L[8]=cos(a[0])*cos(a[1]);

}
void Integrate_Matrix(double* A,double* v,double dt)
{
   int i,j;
   double Left[9],norm,v_M[9],E[9],H1[9],H2[9],H[9],R[9];
   for(i=0;i<9;i++)  E[i] = 0;
   E[0]=1.;E[4]=1.;E[8]=1.;
   makeMatrix(v,v_M);
   norm = 0;
   for(i=0;i<3;i++) norm = norm + sqr(v[i]);
   norm = sqrt(norm);
       //   printf("norm %lf\n",norm);
   for(i=0;i<9;i++) H1[i] = v_M[i] * sin(norm*dt) / norm;
   multMM(v_M,v_M,H2,3);
   for(i=0;i<9;i++) H2[i] = H2[i] * (1. - cos(norm*dt)) / sqr(norm);
   for(i=0;i<9;i++) H[i] = E[i] + H1[i] + H2[i];
   multMM(H,A,R,3);
   for(i=0;i<9;i++) A[i] = R[i];
}
void Intergate_V(double* vec,double* V_m,double* g_m,double* f,double* L,double dt)
{
    int i;
    double vec_m[9],H1[3],H2[3],L_T[9];
//    for(i=0;i<3;i++) printf("%lf ",vec[i]);
//    printf("\n");
    makeMatrix(vec,vec_m);
    multMV(vec_m,V_m,H1,3);
    transp(L,L_T,3);
    multMV(L_T,f,H2,3);
    for(i=0;i<3;i++) V_m[i] = V_m[i] + dt*(H1[i] + H2[i] + g_m[i]);
}
void makeMatrix(double *v,double *M_v)
{
    M_v[0]=0;     M_v[1]=v[2]; M_v[2]=-v[1];
    M_v[3]=-v[2];  M_v[4]=0;    M_v[5]=v[0];
    M_v[6]=v[1];   M_v[7]=-v[0]; M_v[8]=0;
}
void transp(double *S,double *St,int N)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            St[i*N+j] = 0;
            St[i*N+j] = S[j*N+i];
        }
    }
}
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
void multMV(double *a,double *b,double *c,int N)
{
   int i,k;
   for(i=0;i<N;i++)
   {
           c[i]=0;
           for(k=0;k<N;k++) c[i]+=a[i*N+k]*b[k];
   }
}
void Nav(double *Qo,double *p,double *x,double dt,double *Rezn)
{
    int i;
    double q[4],R[4],qn[4],Rn[4],Rez[4];
    QuaternionOfRotation(q,x,dt);
    QuaternionNormalization(q, qn);//2
    QuaternionMultiplication(p,qn,R);
    QuaternionNormalization(R, Rn);
    for(i=0;i<4;i++) p[i]=Rn[i];
    QuaternionMultiplication(Qo,Rn,Rez);
    QuaternionNormalization(Rez, Rezn);
}
void makeQuatMat(double *qw,double *a)
{
   a[0] = (2.*qw[0]*qw[0]+2.*qw[1]*qw[1]-1.);
   a[1] = (2.*qw[1]*qw[2]+2.*qw[0]*qw[3]);
   a[2] = (2.*qw[1]*qw[3]-2.*qw[0]*qw[2]);

   a[3] = (2.*qw[1]*qw[2]-2.*qw[0]*qw[3]);
   a[4] = (2.*qw[0]*qw[0]+2.*qw[2]*qw[2]-1.);
   a[5] = (2.*qw[2]*qw[3]+2.*qw[0]*qw[1]);

   a[6] = (2.*qw[1]*qw[3]+2.*qw[0]*qw[2]);
   a[7] = (2.*qw[2]*qw[3]-2.*qw[0]*qw[1]);
   a[8] = (2.*qw[0]*qw[0]+2.*qw[3]*qw[3]-1.);
}
void Make_quat_allig(double *q,double *a)//0 - тангаж, 1 - крен, 2 - курс.//ПОЧЕМУ ЗДЕСТЬ ТРАНСПОНИРОВАННЫЙ КВАТЕРНИОН А? А? А?
{
    int i;
    for(i=0;i<3;i++) a[i]=a[i]*0.5;
    q[0]=cos(a[2])*cos(a[0])*cos(a[1]) + sin(a[2])*sin(a[0])*sin(a[1]);
    q[1]=cos(a[2])*sin(a[0])*cos(a[1]) + sin(a[2])*cos(a[0])*sin(a[1]);
    q[2]=cos(a[2])*cos(a[0])*sin(a[1]) - sin(a[2])*sin(a[0])*cos(a[1]);
    q[3]=+sin(a[2])*cos(a[0])*cos(a[1]) + cos(a[2])*sin(a[0])*sin(a[1]);
}
void transpose_any(double *a,double *a_t,int u,int v)
{
    int i,j;
    for(i=0;i<v;i++)
    {
        for(j=0;j<u;j++)
        {
            a_t[i*u+j]=0;
            a_t[i*u+j]=a[j*v+i];
        }
    }
}
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
    double S_[N*N],S_T[N*N];
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
        transp(S_,S_T,N);
        multMV(S_T,h,f,N);
        for(k=0; k < N; k++)
        {
            kk=k+1;
            d[kk] = d[kk-1] + sqr(f[k]);
            b[k] = sqrt(d[kk-1]/d[kk]);
            c[k] = f[k]/sqrt(d[kk-1]*d[kk]);
            for(i = 0;i < N; i++) S[i*mes + k] = S_[i*mes+k] * b[k] - e0[i]*c[k];//переставить
            for(i = 0;i < N; i++) e0[i] += S_[i*mes+k] * f[k];
        }
        for(i = 0;i < N; i++) hv[i] = h[i]*x_[i];
        for(i = 0;i < N; i++) xp[i] = x_[i] + e0[i]/d[N]*(z - hv[i]);
    }
}
void intV_G(double* om,double* v,double* g_mod,double* f,double* L,double dt)
{
    int i = 0;
    for(i = 0;i < 3; i++) om[i] = om[i]*dt;
    double G_const = 0;
    for(i = 0;i < 3; i++) G_const+=sqr(om[i]);
    G_const = sqrt(G_const);

    double G_m[9];
    makeMatrix(om,G_m);
    double G_m2[9];
    multMM(G_m,G_m,G_m2,3);
    double E[9]; for(i = 0;i < 9; i++) E[i] = 0;
    E[0] = 1; E[4] = 1; E[8] = 1;

//------------------------------
//    double h_m[9];
//    double rez1[3];
//    double rez2[3];
//    for(i = 0;i < 9; i++) h_m[i] = E[i] + G_m[i] + 0.5*G_m2[i];
//    multMV(h_m,v,rez1,3);
//
//    double fx[3],L_T[9];
//    transp(L,L_T,3);
//    multMV(L_T,f,fx,3);
//    double hv[3];
//    for(i = 0;i < 3; i++) hv[i] = (fx[i] + g_mod[i])*dt;
//
//    for(i = 0;i < 9; i++) h_m[i] = E[i] + 0.5*G_m[i] + G_m2[i]*1./6.;
//    multMV(h_m,hv,rez2,3);
//    for(i = 0;i < 3; i++) v[i] = rez1[i] + rez2[i];



    //--------------
    double h_m[9];
    for(i = 0;i < 9; i++) h_m[i] = E[i] + sin(G_const)*G_m[i]/G_const + (1 - cos(G_const))*G_m2[i]/sqr(G_const);
    double rez1[3];
    multMV(h_m,v,rez1,3);


    for(i = 0;i < 9; i++) h_m[i] = E[i] + (1 - cos(G_const))*G_m[i]/sqr(G_const) + (G_const - sin(G_const))*G_m2[i]/(sqr(G_const)*G_const);
    double fx[3],L_T[9];
    transp(L,L_T,3);
    multMV(L_T,f,fx,3);
    double hv[3];
    for(i = 0;i < 3; i++) hv[i] = (fx[i] + g_mod[i])*dt;
    double rez2[3];
    multMV(h_m,hv,rez2,3);

    for(i = 0;i < 3; i++) v[i] = rez1[i] + rez2[i];
}
