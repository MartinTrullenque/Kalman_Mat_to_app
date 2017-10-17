/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Martín
 *
 * Created on 7 de octubre de 2017, 10:40
 */

#include <cstdlib>
#include <stdio.h>
#include <math.h>


using namespace std;

typedef struct
{
    float x_axis[700];
    float y_axis[700];
    float x_initial;
    float y_initial;
} tag;
typedef struct{
    float toa_A0[700];
    float toa_A1[700];
    float toa_A2[700];
}ttoa;
typedef struct{
    float toa[];
    float pos_x;
    float pos_y;
}tanchor;
/*
 * 
 */

float determinant(float a[3][3],float k)
{
  float s=1,det=0,b[3][3];
  int i,j,m,n,c;
  if (k==1)
    {
     return (a[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        n=0;
        for (i=0;i<k;i++)
          {
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=a[i][j];
                   if (n<(k-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (a[0][c] * determinant(b,k-1));
          s=-1 * s;
          }
    }
 
    return (det);
}

/*Finding transpose of matrix*/ 
void transpose(float num[3][3],float fac[3][3],float r)
{
  int i,j;
  float b[3][3],inverse[3][3],d;
 
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant(num,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        inverse[i][j]=b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is : \n");
 
   for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         printf("\t%f",inverse[i][j]);
        }
    printf("\n");
     }
}
/*For calculating Determinant of the Matrix */

 
void cofactor(float num[3][3],float f)
{
 float b[3][3],fac[3][3];
 int p,q,m,n,i,j;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     n=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(f-2))
             n++;
            else
             {
               n=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,f-1);
    }
  }
  transpose(num,fac,f);
}
int matrix_operation (float first[3][3], float second[3][3], float difference[3][3], int rows, int columns, int flag)
{
    int c=0;
    int d=0;
    for (c = 0; c < rows; c++) {
     for (d = 0; d < columns; d++) {
         if (flag==1){
         difference[c][d] = first[c][d] + second[c][d];    
         }
         else{
         difference[c][d] = first[c][d] - second[c][d];
         }       
     }    
   }
   return 0;
}
int matrix_operationbis (float first[2][2], float second[2][2], float difference[2][2], int rows, int columns, int flag)
{
    int c=0;
    int d=0;
    for (c = 0; c < rows; c++) {
     for (d = 0; d < columns; d++) {
         if (flag==1){
         difference[c][d] = first[c][d] + second[c][d];    
         }
         else{
         difference[c][d] = first[c][d] - second[c][d];
         }       
     }    
   }
   return 0;
}
int matrixProduct2222(float first[2][2], int m, int n, float second[2][2], int p, int q, float multiply[2][2])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  
  }

int matrixProduct2223(float first[2][2], int m, int n, float second[2][3], int p, int q, float multiply[2][3])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  
  }
int matrixProduct2332(float first[2][3], int m, int n, float second[3][2], int p, int q, float multiply[2][2])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  }
int matrixProduct2333(float first[2][3], int m, int n, float second[3][3], int p, int q, float multiply[2][3])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  
  }
int matrixProduct2331(float first[2][3], int m, int n, float second[3], int p, int q, float multiply[2])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k];
        }
        multiply[c] = sum;
        sum = 0;
      
    }
      
    }
  return 1;
  }
 
  int matrixProduct3222(float first[3][2], int m, int n, float second[2][2], int p, int q, float multiply[3][2])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  
  }
int matrixProduct3223(float first[3][2], int m, int n, float second[2][3], int p, int q, float multiply[3][3])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  }

int matrixProduct3332(float first[3][3], int m, int n, float second[3][2], int p, int q, float multiply[2][3])
{
  int sum = 0;
  if (n != p)
      return 0;
  else
  {
      for (int c = 0; c < m; c++) {
      for (int d = 0; d < q; d++) {
        for (int k = 0; k < p; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
        multiply[c][d] = sum;
        sum = 0;
      }
    }
      return 1;
    }
  
  }
int distanceGeneration(tag paciente, int period){
   
    for ( int i=0; i<700; i++){
         // x_axis generation
            paciente.x_axis[i+1]=paciente.x_axis[i];
            paciente.y_axis[i+1]=paciente.y_axis[i];
            
            if (i>300&&i<350){
               paciente.x_axis[i+1]=paciente.x_axis[i]+period;
               paciente.y_axis[i+1]=paciente.y_axis[i]+period;                
            }                
            else if (i>550 && i<600){
                paciente.y_axis[i+1]=paciente.y_axis[i]-period;
            }
                   
    }
    
    return 1;
}
int distanceComputation(tag paciente,tanchor anc0 ){
    float aux;
    float aux2;
    for(int i=0; i<500; i++){
        aux=(anc0.pos_x-paciente.x_axis[i])*(anc0.pos_x-paciente.x_axis[i]);
        aux2=(anc0.pos_y-paciente.y_axis[i])*(anc0.pos_y-paciente.y_axis[i]);
        anc0.toa[i]=hypot(aux,aux2);
                
    }
    
    
    return 1;
}
int main(int argc, char** argv) {
    int period=0.1;
    int aux=0;
    int a0_y=0;
    int a1_x=2.7;
    int a1_y=0;
    int a2_x=1;
    int a2_y=2.8;
    tag paciente;
    tanchor anc0;
    tanchor anc1;
    tanchor anc2;
    paciente.x_axis[0]=1;
    paciente.y_axis[0]=1;
    anc0.pos_x=0;
    anc0.pos_y=0;
    anc1.pos_x =2.7;
    anc1.pos_y=0;
    anc2.pos_x=1;
    anc2.pos_y=2.8;
    //aux=matrixProduct();
    aux=distanceGeneration(paciente,period);
    aux=sqrt(4);
    aux=hypot(3,4);
    printf("%d",aux);
    
    
    int var1= distanceComputation(paciente,anc0);
    int var2= distanceComputation(paciente,anc1);
    int var3= distanceComputation(paciente,anc2);
    printf("%d",var1);
    for(int i=0;i<7700;i++){
        printf("%f",anc0.toa[i]);
        
    }
        
    return 0;
}
double frand() {
    return 2*((rand()/(double)RAND_MAX) - 0.5);
}

int KalmanFilter(tanchor anchor0,tanchor anchor1, tanchor anchor2, float z_measured[3], float z_var ){
    
    float x_est_last = 0;
    float P_last = 0;
    float inverse [3][3];
    float Q[2][2];
    Q[0][0]=0.0001;
    Q[1][1]=0.0001;
    Q[0][1]=0;
    Q[1][0]=0;
    
    float C[3][3];
    C[0][0]=z_var;
    C[1][1]=z_var;
    C[2][2]=z_var;
    C[0][1]=0;
    C[0][2]=0;
    C[1][0]=0;
    C[1][2]=0;
    C[2][0]=0;
    C[2][1]=0;
    
    
    //float R = 0.617;
    float K[2][3];
    
    
    float P[2][2];
    P[0][0]=10000;
    P[1][1]=10000;
    P[1][0]=0;
    P[0][1]=0;
    
    
    
    float P_temp;
    float x_temp_est[2];
    float x_est[2];
    float A[2][2];
    A[0][0]=1;
    A[1][1]=1;
    A[0][1]=0;
    A[1][0]=0;
    float error_ranging[3];
    float G[3][2];
    float D[3];
    
    // D= sqrt((s(1,n)*ones(1,NBS) - BS_position(1,:)).^2 + (s(2,n)*ones(1,NBS) - BS_position(2,:)).^2);
    float aux=x_est[0]-anchor0.pos_x;
    float aux2=x_est[1]-anchor0.pos_y;
    D[0]=hypot(aux,aux2);
    float auxi=x_est[0]-anchor1.pos_x;
    float auxi2=x_est[1]-anchor1.pos_y;
    D[1]=hypot(auxi,auxi2);
    float auxil=x_est[0]-anchor2.pos_x;
    float auxil2=x_est[1]-anchor2.pos_y;
    D[2]=hypot(auxil,auxil2);
    //  dhat= [(s(1,n)*ones(1,NBS) - BS_position(1,:))' (s(2,n)*ones(1,NBS) - BS_position(2,:))'];
    float dhat[3][2];
    
    dhat[0][0]=x_est[0]-anchor0.pos_x;
    dhat[1][0]=x_est[0]-anchor1.pos_x;
    dhat[2][0]=x_est[0]-anchor2.pos_x;
    
    dhat[0][1]=x_est[1]-anchor0.pos_y;
    dhat[1][1]=x_est[1]-anchor1.pos_y;
    dhat[2][1]=x_est[1]-anchor2.pos_y;
     // G = [dhat./(D'*ones(1,2))];

 
    G[0][0]=dhat[0][0]/D[0];
    G[0][1]=dhat[0][1]/D[0];
    
    G[1][0]=dhat[1][0]/D[1];
    G[1][1]=dhat[1][1]/D[1];
    
    G[2][0]=dhat[2][0]/D[2];
    G[2][1]=dhat[2][1]/D[2];
    
    float Gaux [2][3];
    float aux1[3][2];
    float aux12[3][3];
    float aux13[3][3];
    float aux14[2][3];
    float aux15[2];
    float aux16[2][2];
    float aux17[2][2];
    float aux18[2][2];
    
    Gaux[0][0]=G[0][0];
    Gaux[0][1]=G[1][0];
    Gaux[0][2]=G[2][0];
    
    Gaux[1][0]=G[0][1];
    Gaux[1][1]=G[1][1];
    Gaux[1][2]=G[2][1];
    // er_ranging(:,n)=ranging(:,n)-D'; 	% ranging error from the Kalman predicted value   
    error_ranging[0]=z_measured[0]-D[0];
    error_ranging[1]=z_measured[1]-D[1];
    error_ranging[2]=z_measured[2]-D[2];
    
    // K=P*G'(C+G*P*G')^-1
    int var=matrixProduct3222( G,  3,  2,P, 2,  2,aux1 );
    int var1=matrixProduct3223( aux1,  3,  2,Gaux, 2,  3,aux12 );
    int var2=matrix_operation ( C, aux12,aux13,3,3, 1);
    float d=determinant(aux13,3);
    cofactor(aux13,3);
    int var3=matrixProduct2223( P,  2,  2,Gaux, 2,  3,aux14 ); 
    int var4=matrixProduct2333( aux14,  2,  3,inverse, 3,  3,K );

    // s(:,n)=s(:,n)+K*err;
    int var5=matrixProduct2331( K,2,3,error_ranging, 3,1,aux15 ); 
    x_est[0]=x_est[0]+aux15[0];
    x_est[1]=x_est[1]+aux15[1];

    //P = (eye(2) - K * G) * P;
    int var6=matrixProduct2332( K,2,3,G, 3,2,aux16 );
    int var7=matrix_operationbis ( A, aux16,aux17,2,2, 2);
    int var8=matrixProduct2222( aux17,2,2,P, 2,2,aux18 );
  
    for( int o=0;o<2;o++){
        for(int p=0;p<2;p++){
            P[o][p]=aux18[o][p];
        }
    }
    
    
    //s(:,n+1) = A*s(:,n);
    //Mnn = A*Mnn*A' + Q;
    int var9=matrix_operationbis ( P, Q,aux18,2,2, 1);
    
    for( int o=0;o<2;o++){
        for(int p=0;p<2;p++){
            P[o][p]=aux18[o][p];
        }
    }
 
    /*
    printf("Total error if using raw measured:  %f\n",sum_error_measure);
    printf("Total error if using kalman filter: %f\n",sum_error_kalman);
    printf("Reduction in error: %d%% \n",100-(int)((sum_error_kalman/sum_error_measure)*100));
    */
    
    return 0;
}
