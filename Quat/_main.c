#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include "QuaternionFunctions.h"
#include "Navigation.h"
#include "criterions.h"
#include "readcom.h"
#define sqr(x) ((x)*(x))
#define pi            3.14159265358979323846
#define g             9.780327
#define time_skip 0
#define time_aligment 60
#define krit1 10
#define krit2 20
#define krit3 10
#define krit4 0.01
#define krit5 0.02
#include "help_func.h"
#include "integrate.h"
  typedef struct
{
unsigned char num;
unsigned char invpac;
unsigned short serial_number;
unsigned short serial_number2;
unsigned short flags;
unsigned short tempK;
unsigned short gyro_X1;
unsigned short gyro_X2;
unsigned short gyro_Y1;
unsigned short gyro_Y2;
unsigned short gyro_Z1;
unsigned short gyro_Z2;
unsigned short accel_X1;
unsigned short accel_X2;
unsigned short accel_Y1;
unsigned short accel_Y2;
unsigned short accel_Z1;
unsigned short accel_Z2;
unsigned short CRC;
unsigned short CRC2;

} AIST_DATA;

int main(void)
{
    int i,j;
    double Mat_help[9],count;
    //переменные для 1ого уравнения
    double V_m[3],OmB_m[3],u_m[3]; //относительная угловая скорость трехгранника, скорость вращения Земли
    //double MatOmB_m[N*N],Matu_m[N*N];
    double g_m[3];
    double u;//модуль u_m
    //для модели g
    double g_0,betta,lyam_m,phi_m,a_sa,e,h_m;
    //для 2ого уравнения
    double R2_m,R1_m;
    //ориентация
    double D_xz[9],fx[3],B_m[9];//D - xz for f || B - xetta for OmB_m
    double Om_m[3];//MatOm_m[N*N] абсолютная угловая скорость трехгранника
    h_m = 144; phi_m = 55.425395 *pi/180; lyam_m = 37.314170*pi/180;
    a_sa = 6378137;
    e = 6.59e-3;
    u=2*pi/(24*3600);
    R1_m = a_sa/sqrt(1-e*sqr(sin(phi_m)))+h_m; R2_m = a_sa*(1-e)/((1-e*sqr(sin(phi_m)))*sqrt(1-e*sqr(sin(phi_m))))+h_m;
    for(i=0;i<3;i++)
        {
            V_m[i]=0;
            OmB_m[i]=0;
        }
    betta = 0;
    g_0 = 9.78030;
    double B_0[4],p_B[4],Rez_Om[4];
    p_B[0]=1;
    for(i=1;i<4;i++) {p_B[i] = 0;}
    matrixB(lyam_m,phi_m,B_m,3);
    quat_for_B0(lyam_m,phi_m,B_0);
double D_t[9];
//Кватернионы для интегрирования ориентации
double qAx_0[4],qAz_0[4],p_Ax[4],p_Az[4],qAx[4],qAz[4];
double Ax[9],Az[9],D[9],Help_mat[9],help_vec[3];
qAx_0[0] = 1; qAx_0[1] = 0; qAx_0[2] = 0; qAx_0[3] = 0;
p_Ax[0] = 1; p_Ax[1] = 0; p_Ax[2] = 0; p_Ax[3] = 0;
p_Az[0] = 1; p_Az[1] = 0; p_Az[2] = 0; p_Az[3] = 0;
double quat_D[4];

double L_for_all[9],u_almost[3],dr[3],u_et[3];
u_et[0]=0;u_et[1]=0;u_et[2]=u;
multMV(B_m,u_et,u_m,3);
help_vec[2] = 0; help_vec[0] = u*cos(phi_m); help_vec[1] = u*sin(phi_m);
for(i=0;i<3;i++) printf("sravn %lf %lf\n",u_m[i],help_vec[i]);
g_m[2] = g_m[0] = 0; g_m[1] = -g_0*(1+0.75*e*sqr(sin(phi_m)));







       float gyro_X,gyro_Y,gyro_Z,accel_X,accel_Y,accel_Z;
       //Для порта
       unsigned char buff[1];
       AIST_DATA dst;
       float Mwx,Mwy,Mwz,Mfx,Mfy,Mfz,Afx,Afy,Afz;
       Mwx=1.991;
       Mwy=1.993;
       Mwz=2.0;

       Mfx=2.002;
       Mfy=2.003;
       Mfz=2.002;

       Afx=-254.6;
       Afy=261.75;
       Afz=-261.25;
       DWORD  sizerd;
       int dstsize = 38;
       DCB dcbCommPort;
       HANDLE port;
       port = CreateFile("COM3", GENERIC_READ, 0, NULL,
		OPEN_EXISTING, 0, NULL);
//       if (port == INVALID_HANDLE_VALUE)
//        {
//            MessageBox(NULL, "Невозможно открыть последовательный порт", "Error", MB_OK);
//            ExitProcess(1);
//        }

dcbCommPort.DCBlength = sizeof(DCB);
GetCommState(port, &dcbCommPort);

BuildCommDCB("baud=115200 parity=N data=8 stop=1 xon=off octs=off rts=off", &dcbCommPort);//115200,8,N,1
SetCommState(port, &dcbCommPort);
PurgeComm(port, PURGE_RXABORT|PURGE_RXCLEAR);

    int k=0,n=0,w,first_in_nav=0;
    char buf[100];
    double tang ,kre ,kur,tangal,kreal,kural;
    double a,b,x[3],xo[3],xav[3],f[3],fav[3],dt,koef=pi/180,*ix1,*ix2,*ix3,temperature;
    double Qo[4],Rezn[4];
    for(i=0;i<3;i++) fav[i]=0;
    double p[4];
    double Nq=0;
    double krit_av[3],sigma_sum[3];
    for(i=0;i<3;i++) {krit_av[i]=0;sigma_sum[i]=0;}

    double counst_for_krit[5];
    counst_for_krit[0]=krit1;
    counst_for_krit[1]=krit2;
    counst_for_krit[2]=krit3;
    counst_for_krit[3]=krit4;
    counst_for_krit[4]=krit5;
    int check[5],check_print[5],MAIN_CHECK;
    for(i=0;i<5;i++) {check[i]=0;check_print[i]=0;}
    int N;
    MAIN_CHECK = 0;
    char ch='0';
    //переменные для ИМ
    int len = 25,cnt;
    double st1[len],st2[len],st3[len],quat_int[4],quat_int2[4],one[4],A[9],fav_forint[3];

    double krit_int,kre_int,tang_int;
    cnt=0;
    for(i=0;i<4;i++) one[i] = 0;
    one[0]=1;
    //конец переменных для ИМ
    w=0;
    FILE *fp;
    FILE *fpp;
    FILE *fppp;
    FILE *ind;
    fp=fopen("newData.txt","r");//AIST-29B88555_20131211_1154_250hz.cbr","r");
    fpp=fopen("RR1.txt","w");
    fppp=fopen("Index.txt","w");
    ind=fopen("Ind.txt","w");
//    while(ch != '\n') {fscanf(fp,"%c",&ch);}
//    ch='0';
//    fscanf(fp,"%lf",&a);
//    while(ch != '\n') {fscanf(fp,"%c",&ch);}
//    ch='0';
//    for(i=0;i<4;i++) printf("%lf    ",Qo[i]);
//    printf("\n");
    //QuaternionNormalization(p, pn);

    dt=0.004;
    b=dt;
    N = 25;
    quat_int2[0] = 1;
    for(i=1;i<4;i++) quat_int2[i]=0;
    p[0] = 1;
    for(i=1;i<4;i++) p[i]=0;
        ix1 = (double*)malloc(N*sizeof(double));
    for(i=0;i<N;i++) ix1[i] = 0;
        ix2 = (double*)malloc(N*sizeof(double));
    for(i=0;i<N;i++) ix2[i] = 0;
        ix3 = (double*)malloc(N*sizeof(double));
    for(i=0;i<N;i++) ix3[i] = 0;
    for(i=0;i<3;i++) xav[i]=0;
    for(i=0;i<3;i++) fav[i]=0;
    for(i=0;i<3;i++) fav_forint[i]=0;
    //для ИМ
    for(i=0;i<len;i++)
    {
        st1[i]=0;
        st2[i]=0;
        st3[i]=0;
    }
    //конец для ИМ

    while(b<757)
    {
        if(b < time_skip)
        {
            b = b+dt;
            fgets(buf,100,fp);
            continue;
        }
//      //Чтение из порта
//      {
//           memset((char*)&dst,0,dstsize);
//        while(1)
//        {
//            sizerd = 0;
//
//            ReadFile(port,buff,1,&sizerd,NULL);
//
//            if(buff[0]==0x47)
//            {
//
//
//                ReadFile(port,buff,1,&sizerd,NULL);
//
//                if(buff[0]==0x59)
//                {
//                    break;
//                }
//
//            }
//        }
//    	ReadFile(port,(unsigned char *)&dst,dstsize,&sizerd,NULL);
//
//    	*(unsigned long*)(&gyro_X) = (((unsigned long)(dst.gyro_X2))*65536+dst.gyro_X1);
//    	*(unsigned long*)(&gyro_Y) = (((unsigned long)(dst.gyro_Y2))*65536+dst.gyro_Y1);
//    	*(unsigned long*)(&gyro_Z) = (((unsigned long)(dst.gyro_Z2))*65536+dst.gyro_Z1);
//    	*(unsigned long*)(&accel_X) = (((unsigned long)(dst.accel_X2))*65536+dst.accel_X1);
//    	*(unsigned long*)(&accel_Y) = (((unsigned long)(dst.accel_Y2))*65536+dst.accel_Y1);
//    	*(unsigned long*)(&accel_Z) = (((unsigned long)(dst.accel_Z2))*65536+dst.accel_Z1);
//    	xo[2] = gyro_X/(3600*Mwx);
//    	xo[0] = gyro_Y/(3600*Mwy);
//    	xo[1] = gyro_Z/(3600*Mwz);
//    	f[2] = (accel_X - Afx)/Mfx;
//    	f[0] = (accel_Y - Afy)/Mfy;
//    	f[1] = (accel_Z - Afz)/Mfz;
//    	//printf(" %f  %f  %f | %f  %f  %f\n",gyro_X,gyro_Y,gyro_Z,accel_X,accel_Y,accel_Z);
//       // swapbytes(&dst.serial_number,sizeof dst.serial_number);
//}


             //Чтение из файла
             {
             fscanf(fp,"%lf",&temperature);

             fscanf(fp,"%lf",&xo[2]);//201
             fscanf(fp,"%lf",&xo[0]);
             fscanf(fp,"%lf",&xo[1]);

             fscanf(fp,"%lf",&f[2]);
             fscanf(fp,"%lf",&f[0]);
             fscanf(fp,"%lf",&f[1]);
             //if (b < 0.012) printf("-> %lf %lf %lf %lf %lf %lf %lf\n",temperature,xo[2],xo[0],xo[1],f[2],f[0],f[1]);
             }
            xo[2] = xo[2]/(3600*Mwx);//201
            xo[0] = xo[0]/(3600*Mwy);
            xo[1] = xo[1]/(3600*Mwz);
            f[2] = (f[2] - Afx)/Mfx;
            f[0] = (f[0] - Afy)/Mfy;
            f[1] = (f[1] - Afz)/Mfz;

             x[1] = koef*(xo[1] - xav[1]);
             x[0] = koef*(xo[0] - xav[0]);
             x[2] = koef*(xo[2] - xav[2]);
             //if (b < 11) printf("%lf  %lf  %lf\n",f[2],f[0],f[1]);
//КРИТЕРИИ ЗДЕСЬ
            krit_function(f,x,krit_av,sigma_sum,ix1,ix2,ix3,N,w,check,counst_for_krit);
            w++;
kural=0;//315*pi/180;
{
  if (b<time_skip+2)
  {
            Alignment(f,fav_forint,cnt);
            if(b>1.995) {
//                    Ali_ang(fav_forint,&tang_int,&kre_int);
            q(kre_int,kural,tang_int,quat_int);
            makaAmatrix(A,quat_int);
                    //for(i=0;i<3;i++){for(j=0;j<3;j++) {printf("%lf ",A[i*3+j]);}printf("\n");}
            }
  }
  else{ krit4_function(len,st1,st2,st3,A,f,check,counst_for_krit,g);}
}


cnt++;
//Конец КРИТЕРИЕВ
//fprintf(ind,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",b,fabs(n1-n1av),g1,fabs(n2-n2av),g2,fabs(n3-n3av),g3,n4,krit4,krit_int,krit5);


 fprintf(fppp,"%lf %d %d %d %d %d\n",b,check[0],check[1],check[2],check[3],check[4]);
 //ВЫСТАВКА
        if(b>=time_skip && b<=time_skip+time_aligment)
        {
            Alignment(f,fav,n);
            Alignment(xo,xav,n);
            n++;
            Ali_ang(fav,&tangal,&kreal,&kural,xav);
//            q(kreal,kural,tangal,Qo);
            fprintf(fpp,"%lf  %lf  %lf  %lf 0 0 0  %lf  %lf\n",b,tangal,kreal,kural,lyam_m,phi_m);
            b = b+dt;
//            multMV(B_m,u_et,u_m,3);//u_m[2] = 0; u_m[0] = u*cos(phi_m); u_m[1] = u*sin(phi_m);
//            def_kur(tangal,kreal,&kural,u_m,xav,phi_m,u,L_for_all);
//            multMV(L_for_all,u_m,u_almost,3);
//            transp(L_for_all,D);
//            multMV(D,u_almost,dr,3);
//            kural=kural*180/pi;
//            multMV(L_for_all,fav,help_vec,3);
//
//           for(i=0;i<3;i++) dr[i] = help_vec[i] + g_m[i];
//           for(i=0;i<3;i++) printf("%lf ",dr[i]);
//           printf("\n");



            continue;

        }
//НАВИГАЦИЯ
        if(b > time_skip+time_aligment)//&& MAIN_CHECK==1)//&& b <=443)
        {

            if(first_in_nav==0)
                {
                    printf("angles alligment %lf %lf %lf\n",tangal*180/pi,kreal*180/pi,kural*180/pi);
                    quat_for_Az0(kreal,kural,tangal,qAz_0);
                    makeQuatMat(qAz_0,D);
//                    printf("-----1----->%lf %lf %lf %lf\n",qAz_0[0],qAz_0[1],qAx_0[2],qAz_0[3]);
//                    makeQuatMat(qAz_0,Help_mat);
//                                kur = -atan2(Help_mat[2],Help_mat[0]);//20
//                                kre = -atan2(Help_mat[7],Help_mat[4]);//74
//                                tang = atan2(Help_mat[1],sqrt(sqr(Help_mat[4])+sqr(Help_mat[7])));//102
//
////                    makeMatQuat(Help_mat,qAz_0);
////                    QuaternionNormalization(qAz_0,qAx_0);
//                    printf("-----2----->%lf %lf %lf\n",tangal*180/pi,kreal*180/pi,kural*180/pi);

//                    for(i=0;i<3;i++)
//                        {
//                            printf("%lf ",fx[i]+g_m[i]);
//                            for(j=0;j<3;j++)
//                            {
//                                printf("%lf ",Help_mat[i*3+j]);
//                            }
//                            printf("\n");
//                        }

                    first_in_nav++;
                }

            u_m[2] = 0; u_m[0] = u*cos(phi_m); u_m[1] = u*sin(phi_m);
            g_m[2] = 0; g_m[0] = 0; g_m[1] = -g_0*(1.+0.75*e*sqr(sin(phi_m)));
            R1_m = a_sa/sqrt(1-e*sqr(sin(phi_m)))+h_m; R2_m = a_sa*(1-e)/((1-e*sqr(sin(phi_m)))*sqrt(1-e*sqr(sin(phi_m))))+h_m;
            OmB_m[2] = -V_m[0]/R2_m; OmB_m[0] = V_m[2]/R1_m; OmB_m[1] = V_m[2]*tan(phi_m)/R1_m;
            //интегрирование Ax Az
//            for(i=0;i<3;i++) help_vec[i] = OmB_m[i]+u_m[i];
            //if(first_in_nav<2)for(i=0;i<3;i++) printf("hv %lf \n",help_vec[i]);
            transp(D,Help_mat);
            multMV(Help_mat,f,fx,3);
            integrate_V(V_m,OmB_m,u_m,fx,g_m,dt,3);
            V_m[1]=0;

            //интергрирование координат
            Nav(B_0,p_B,OmB_m,dt,Rez_Om);
            makeQuatMat(Rez_Om,B_m);
            lyam_m = -atan2(B_m[2],B_m[5]);//25
            phi_m = atan2(B_m[7],B_m[6]);//76

            for(i=0;i<3;i++) help_vec[i] = OmB_m[i]+u_m[i];

            Nav(qAx_0,p_Ax,help_vec,dt,qAx);
//            for(i=0;i<4;i++) printf(" %lf ",p_Az[i]);
//            printf("\n");
            Nav(qAz_0,p_Az,x,dt,qAz);
//            if(first_in_nav<15)
//            {
//                for(i=0;i<3;i++)
//                {
//                    printf("%lf ",u_m[i]);
//                }
//                printf("\n");
//                first_in_nav++;
//            }

//              makaAmatrix(Ax,qAx);
//              makaAmatrix(Az,qAz);
            makeQuatMat(qAx,Ax);
            makeQuatMat(qAz,Az);
            transp(Az,Help_mat);
            multMM(Help_mat,Ax,D,3);
//            transp(D,Help_mat);
//            for(i=0;i<9;i++) D[i]=Help_mat[i];

            kur = -atan2(D[2],D[0]);//20
            kre = -atan2(D[7],D[4]);//74
            tang = atan2(D[1],sqrt(sqr(D[4])+sqr(D[7])));//102

//            transp(D,Help_mat);
//            multMV(Help_mat,f,fx,3);
//            if(b<time_skip+time_aligment+10)
//            {
//                Alignment(fx,fav,first_in_nav);
//                for(i=0;i<3;i++) dr[i] = fav[i]+g_m[i];
//                first_in_nav++;
//            }
            //for(i=0;i<3;i++) g_m[i] = g_m[i]-dr[i];

            if(first_in_nav<2)
            {
                //printf("%lf \n",sqrt(sqr(f[0])+sqr(f[1])+sqr(f[2])));
                for(i=0;i<3;i++) printf("sum %lf\n",g_m[i]+fx[i]);
                printf("\n");
                printf("angles %lf %lf %lf\n",tang*180/pi,kre*180/pi,kur*180/pi);
                for(i=0;i<3;i++)
                {
//                    printf("%lf ",fx[i]+g_m[i]);
                    for(j=0;j<3;j++)
                    {
                        printf("%lf ",D[i*3+j]);
                        if(j==2) printf("| %lf = %lf ",f[i],fx[i]);
                    }
                    printf("\n");

                }
                printf("----------------\n");
                printf("%lf  %lf  %lf\n",fx[0]+g_m[0],fx[1]+g_m[1],fx[2]+g_m[2]);
                first_in_nav++;
            }
//            integrate_V(V_m,OmB_m,u_m,fx,g_m,dt,3);
//            V_m[1]=0;
//
//            //интергрирование координат
//            Nav(B_0,p_B,OmB_m,dt,Rez_Om);
//            makeQuatMat(Rez_Om,B_m);
//            lyam_m = -atan2(B_m[2],B_m[5]);//25
//            phi_m = atan2(B_m[7],B_m[6]);//76


//            if(first_in_nav==0) {q(kreal,kural,tangal,Qo);first_in_nav++;}
//            Nav(Qo,p,x,dt,Rezn);
//            ang(Rezn,&tang,&kre,&kur,b);
            fprintf(fpp,"%lf  %lf  %lf  %lf %lf  %lf  %lf %lf  %lf\n",b,tang,kre,kur,V_m[0],V_m[1],V_m[2],lyam_m,phi_m);
        }
       b = b+dt;
for(i=0;i<5;i++)
    {
        if(check[i]==1)
            {
                check_print[i]=1;
                MAIN_CHECK = 1;
            }
        else
        {
                check_print[i]=0;
                MAIN_CHECK = 0;
        }
    }

//if(t==0) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\
//%7.3lf   %+6.2lf %+7.2lf %+7.2lf   %1d %1d %1d %1d %1d",b,tangal*180/pi,kreal*180/pi,kural*180/pi,check_print[0],check_print[1],check_print[2],check_print[3],check_print[4]);
//else
//printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\
//%7.3lf   %+6.2lf %+7.2lf %+7.2lf   %1d %1d %1d %1d %1d",b,tang*180/pi,kre*180/pi,kur*180/pi,check_print[0],check_print[1],check_print[2],check_print[3],check_print[4]);
////if(nm1==1 || nm2==1 || nm3==1 || nm4==1) fprintf(fppp,"%lf %d %d %d %d %lf\n",b,nm1,nm2,nm3,nm4,krit_int);
    }
    CloseHandle(port);
    return 0;
}
