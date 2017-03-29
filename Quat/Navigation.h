#define pi            3.14159265358979323846
#define sqr(x) ((x)*(x))
void Alignment(double *f,double *fav,int n);
void Alignment(double *f,double *fav,int n)
{
    fav[0] = (fav[0])*n/(n+1.) + f[0]/(n+1);
    fav[1] = (fav[1])*n/(n+1.) + f[1]/(n+1);
    fav[2] = (fav[2])*n/(n+1.) + f[2]/(n+1);
}
void Nav(double *Qo,double *p,double *x,double dt,double *Rezn);
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
void ang(double *qw,double *tang, double *kre,double *kur,double b);
void ang(double *qw,double *tang, double *kre,double *kur,double b)
{
    double a11,a13,a12,a22,a23;
   a12 = 2*qw[1]*qw[2]+2*qw[0]*qw[3];
   a11 = 2*qw[0]*qw[0]+2*qw[1]*qw[1]-1;
   a13 = 2*qw[1]*qw[3]-2*qw[0]*qw[2];
   a22 = 2*qw[0]*qw[0]+2*qw[2]*qw[2]-1;
   a23 = 2*qw[2]*qw[3]-2*qw[0]*qw[1];
   *tang = atan2(a12,sqrt(a11*a11+a13*a13));
   a22 = a22/cos(*tang);
   a23= -a23/cos(*tang);
   a11 = a11/cos(*tang);
   a13 = -a13/cos(*tang);
   *kre = atan2(a23,a22);
   *kur = -atan2(a13,a11);
      if(a12>0.996)
   {

//       QuaternionMultiplication(qw,quat,rez);
//       QuaternionNormalization(rez,rezn);
//       for(i=0;i<4;i++) qw[i]=rezn[i];
       *tang=pi/2;
       *kur=0;
       *kre=atan2(-2*qw[2]*qw[3]-2*qw[0]*qw[1],(2*qw[1]*qw[2]-2*qw[0]*qw[3]));//*kre=atan2(-(2*qw[2]*qw[3]+2*qw[0]*qw[1]),2*qw[1]*qw[2]-2*qw[0]*qw[3]);
   }
      if(a12<-0.996)
   {
       *tang=-pi/2;
       *kur=0;
       *kre=atan2(-2*qw[2]*qw[3]-2*qw[0]*qw[1],(2*qw[1]*qw[2]-2*qw[0]*qw[3]));//*kre=atan2(-(2*qw[0]*qw[0]+2*qw[3]*qw[3]-1),2*qw[3]*qw[3]+2*qw[0]*qw[1]);//*kre=atan2(-(2*qw[1]*qw[2]-2*qw[0]*qw[3]),2*qw[2]*qw[3]+2*qw[0]*qw[1]);
   }
   return;
}
//Выставка

void Ali_ang(double *fav,double *tangal,double *kreal,double *kural,double *xav);
void Ali_ang(double *fav,double *tangal,double *kreal,double *kural,double *xav)
{
    double xavp[3];
    *tangal = atan2(fav[0],sqrt(fav[1]*fav[1]+fav[2]*fav[2]));//012
    *kreal = -atan2(fav[2],fav[1]);//21

//    xavp[0]=xav[0]*cos(*tangal) +xav[1]*sin(*tangal) +xav[2]*0;
//    xavp[1]=-sin(*tangal)*cos(*kreal)*xav[0] +cos(*tangal)*cos(*kreal)*xav[1]+xav[2]*sin(*kreal);
//    xavp[2]=sin(*tangal)*sin(*kreal)*xav[0]* -cos(*tangal)*sin(*kreal)*xav[1]+xav[2]*cos(*kreal);
//    *kural = -atan2(xavp[2],xavp[0]);
}
void def_kur(double tangal,double kreal,double *kural,double *u_m,double *xav,double phi_m,double u,double *L_for_all);
void def_kur(double tangal,double kreal,double *kural,double *u_m,double *xav,double phi_m,double u,double *L_for_all)
{
    double l1[3],l2[3],l3[3],norma;
    norma=0;
    int i,j;
    l2[0]=sin(tangal);
    l2[1]=cos(tangal)*cos(kreal);
    l2[2]=-cos(tangal)*sin(kreal);
    for(i=0;i<3;i++)
    {
        l1[i] = (xav[i]-u*l2[i]*sin(phi_m))/cos(phi_m);
    }
    for(i=0;i<3;i++) norma = norma + sqr(l1[i]);
    for(i=0;i<3;i++) l1[i] = l1[i]/sqrt(norma);

    l3[0]=l1[1]*l2[2]-l1[2]*l2[1];
    l3[1]=-l1[0]*l2[2]+l1[2]*l2[0];
    l3[2]=l1[0]*l2[1]-l1[1]*l2[0];
    for(i=0;i<3;i++) norma = norma + sqr(l3[i]);
    for(i=0;i<3;i++) l3[i] = l3[i]/sqrt(norma);

    for(i=0;i<3;i++)
    {
        L_for_all[i*3+0] = l1[i];
        L_for_all[i*3+1] = l2[i];
        L_for_all[i*3+2] = l3[i];
    }
   // for(i=0;i<3;i++) printf("allig %lf %lf %lf \n",l1[i],l2[i],l3[i]);
    *kural = - atan2(l3[0],l1[0]);


}

void q(double kre,double kur,double tang,double *Qo);
void q(double kre,double kur,double tang,double *Qo)
{
    Qo[0]=cos(0.5*tang)*cos(0.5*kre)*cos(0.5*kur)-sin(0.5*tang)*sin(0.5*kre)*sin(0.5*kur);
    Qo[1]=-cos(0.5*tang)*sin(0.5*kre)*cos(0.5*kur)-sin(0.5*tang)*cos(0.5*kre)*sin(0.5*kur);
    Qo[2]=-sin(0.5*tang)*sin(0.5*kre)*cos(0.5*kur)-cos(0.5*tang)*cos(0.5*kre)*sin(0.5*kur);
    Qo[3]=-sin(0.5*tang)*cos(0.5*kre)*cos(0.5*kur)+cos(0.5*tang)*sin(0.5*kre)*sin(0.5*kur);
}
void makeQuatMat(double *qw,double *a);
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
void quat_for_Az0(double kre,double kur,double tang,double *q);
void quat_for_Az0(double kre,double kur,double tang,double *q)
{
//    q[0]= cos(0.5*tang)*cos(0.5*kre)*cos(0.5*kur);
//    q[1]=-cos(0.5*tang)*sin(0.5*kre)*cos(0.5*kur);
//    q[2]=-sin(0.5*tang)*sin(0.5*kre)*cos(0.5*kur);
//    q[3]=-sin(0.5*tang)*cos(0.5*kre)*cos(0.5*kur);
    q[0]=cos(0.5*tang)*cos(0.5*kre)*cos(0.5*kur)-sin(0.5*tang)*sin(0.5*kre)*sin(0.5*kur);
    q[1]=-cos(0.5*tang)*sin(0.5*kre)*cos(0.5*kur)-sin(0.5*tang)*cos(0.5*kre)*sin(0.5*kur);
    q[2]=-sin(0.5*tang)*sin(0.5*kre)*cos(0.5*kur)-cos(0.5*tang)*cos(0.5*kre)*sin(0.5*kur);
    q[3]=-sin(0.5*tang)*cos(0.5*kre)*cos(0.5*kur)+cos(0.5*tang)*sin(0.5*kre)*sin(0.5*kur);
//    double D[9],ang[3],k,kr,t;
//    ang[0] = kur;
//    ang[1] = kre;
//    ang[2] = tang;
//    printf("%lf %lf %lf\n",ang[0],ang[1],ang[2]);
//    matrixD(ang,D,3);
//    makeMatQuat(D,q);
//    makeQuatMat(q,D);
//    k = -atan2(D[3],D[4]);
//    kr = -atan2(D[2],D[8]);//nen 2*N
//    t = atan2(D[5],sqrt(sqr(D[3])+sqr(D[4])));
//    printf("%lf %lf %lf\n",k,kr,t);
}
void quat_for_B0(double l,double p,double *q);
void quat_for_B0(double l,double p,double *q)
{
    q[0]= -0.5*(-cos(0.5*(l-p))-sin(0.5*(l-p)));
    q[1]= 0.5*(cos(0.5*(l+p))-sin(0.5*(l+p)));
    q[2]= 0.5*(cos(0.5*(l+p))+sin(0.5*(l+p)));
    q[3]= 0.5*(cos(0.5*(l-p))-sin(0.5*(l-p)));
//    double D[9],ang[3],k,kr,t;
//    ang[0] = kur;
//    ang[1] = kre;
//    ang[2] = tang;
//    printf("%lf %lf %lf\n",ang[0],ang[1],ang[2]);
//    matrixD(ang,D,3);
//    makeMatQuat(D,q);
//    makeQuatMat(q,D);
//    k = -atan2(D[3],D[4]);
//    kr = -atan2(D[2],D[8]);//nen 2*N
//    t = atan2(D[5],sqrt(sqr(D[3])+sqr(D[4])));
//    printf("%lf %lf %lf\n",k,kr,t);
}
