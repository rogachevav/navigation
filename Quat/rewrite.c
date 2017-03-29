#include <stdio.h>
int main(void)
{
    FILE *fp;
    FILE *fpp;
    fp=fopen("Data.txt","r");
    fpp=fopen("newData.txt","w");
    double x[3],f[3];
    char ch;
    int i;
   while(1)
   {
        fscanf(fp,"%c",&ch);
        if(ch=='g') break;
        if(ch==',') continue;
        else
        {
            fprintf(fpp,"%c",ch);
            //if(ch=='\n')break;
        }

    }
    return 0;
}
