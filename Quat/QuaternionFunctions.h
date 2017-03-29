
void QuaternionMultiplication(double *daP, double *daQ, double *daR);
void QuaternionNormalization(double *dQuaternion, double *dQuaternionNormalized);
void QuaternionOfRotation(double* q, double *x, double dt);
void QuaternionNormalization2(double *dQuaternion, double *dQuaternionNormalized);
void QuaternionNormalization2(double *dQuaternion, double *dQuaternionNormalized)
{
    double dC;
    int i;


    for( i = 0; i < 4; i++) dQuaternionNormalized[i] = dQuaternion[i];
    dC = 1;
    for( i = 1; i < 4; i++) dC =dC+(dQuaternion[i]*dQuaternion[i]);
    dC = 1/sqrt(dC);
        for( i = 0; i < 4; i++)
            {
                dQuaternionNormalized[i] = dQuaternion[i]*dC;

            }
    return;
}
void QuaternionNormalization(double *dQuaternion, double *dQuaternionNormalized)
{
    double dC;
    int i;


    for( i = 0; i < 4; i++) dQuaternionNormalized[i] = dQuaternion[i];
    dC = 0;
    for( i = 0; i < 4; i++) dC =dC+(dQuaternion[i]*dQuaternion[i]);
    dC = 1/sqrt(dC);
        for( i = 0; i < 4; i++)
            {
                dQuaternionNormalized[i] = dQuaternion[i]*dC;

            }
    return;
}
void QuaternionOfRotation(double* q, double *x, double dt)
{
        q[0] = 1;
        q[1] = 0.5*dt*x[0];
        q[2] = 0.5*dt*x[1];
        q[3] = 0.5*dt*x[2];
        return;
}
void QuaternionMultiplication(double *daP, double *daQ, double *daR)
{
    daR[0] = daP[0]*daQ[0] - daP[1]*daQ[1] - daP[2]*daQ[2] - daP[3]*daQ[3];
    daR[1] = daQ[0]*daP[1] + daP[0]*daQ[1] - daP[3]*daQ[2] + daP[2]*daQ[3];
    daR[2] = daQ[0]*daP[2] + daP[0]*daQ[2] + daP[3]*daQ[1] - daP[1]*daQ[3];
    daR[3] = daQ[0]*daP[3] + daP[0]*daQ[3] - daP[2]*daQ[1] + daP[1]*daQ[2];

    return;
}
void makeMatQuat(double *a,double *q);
void makeMatQuat(double *a,double *q)
{
    q[0] = 0.5/sqrt(1+a[0]+a[4]+a[8]);
    q[1] = (a[5]-a[7])*0.5/sqrt(1+a[0]+a[4]+a[8]);
    q[2] = (a[6]-a[2])*0.5/sqrt(1+a[0]+a[4]+a[8]);
    q[3] = (a[1]-a[3])*0.5/sqrt(1+a[0]+a[4]+a[8]);
}
