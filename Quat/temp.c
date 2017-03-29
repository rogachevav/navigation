st1[len-1]= (2*sqr(quat_int[0])+2*sqr(quat_int[1])-1)*f[0] + (2*quat_int[1]*quat_int[2]+2*quat_int[0]*quat_int[3])*f[1] + (2*quat_int[1]*quat_int[3]-2*quat_int[1]*quat_int[2])*f[2];
st2[len-1]= (2*sqr(quat_int[0])-2*sqr(quat_int[2])-1)*f[1] + (2*quat_int[1]*quat_int[2]+2*quat_int[0]*quat_int[3])*f[0] + (2*quat_int[2]*quat_int[3]+2*quat_int[1]*quat_int[0])*f[2] - g;
st3[len-1]= (2*sqr(quat_int[0])+2*sqr(quat_int[3])-1)*f[2] + (2*quat_int[3]*quat_int[2]-2*quat_int[0]*quat_int[1])*f[1] + (2*quat_int[1]*quat_int[3]+2*quat_int[0]*quat_int[2])*f[0];
krit_int = sqrt(sqr(norm(st1,len))+sqr(norm(st2,len))+sqr(norm(st3,len)));
printf("krint_int = %lf\n",krit_int);
st1[len-1]= (2*sqr(quat_int[0])+2*sqr(quat_int[1])-1)*f[0] + (2*quat_int[1]*quat_int[2]-2*quat_int[0]*quat_int[3])*f[1] + (2*quat_int[1]*quat_int[3]+2*quat_int[1]*quat_int[2])*f[2];
st2[len-1]= (2*sqr(quat_int[0])+2*sqr(quat_int[2])-1)*f[1] + (2*quat_int[1]*quat_int[2]+2*quat_int[0]*quat_int[3])*f[0] + (2*quat_int[2]*quat_int[3]-2*quat_int[1]*quat_int[0])*f[2] - g;
st3[len-1]= (2*sqr(quat_int[0])+2*sqr(quat_int[3])-1)*f[2] + (2*quat_int[3]*quat_int[2]+2*quat_int[0]*quat_int[1])*f[1] + (2*quat_int[1]*quat_int[3]-2*quat_int[0]*quat_int[2])*f[0];
