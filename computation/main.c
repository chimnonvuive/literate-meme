#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
// TODO: Student have to implement all required functions here
// All functions is implement before main function
void initializeMatrix(int *matrix, int row, int col);
void printMatrix(int *matrix, int row, int col);

int fightOneByOne(int& HP1, int& HP2, int n) // This function is only declared.
{
    int LTK=999;
    int QQ=777;
    int TT=666;
    if (HP1==LTK && HP2!=QQ && HP2!= TT)
    {   HP2=0;
        return 1;}
    if ((HP1!=LTK && HP2==QQ) || (HP1!=LTK && HP2==TT))
    {   HP1=0;
        return -1;}
    if ((HP1==LTK && HP2==QQ) || (HP1==LTK && HP2==TT))
    {   return 0;}
    if (n==-1) //Quan tong tan cong
    {
        float newHP1, newHP2;
        newHP1=float(HP1)*1.1;
        newHP2=float(HP2)*1.3;
        HP1=ceil(newHP1);
        HP2=ceil(newHP2);
    }
    else if (n==1) //Quan ta tan cong
    {
        float newHP1, newHP2;
        newHP1=float(HP1)*1.5;
        newHP2=float(HP2)/1.3;
        HP1=ceil(newHP1);
        HP2=ceil(newHP2);
    }
    if (HP1<=0)   HP1=0;
    if (HP2<=0)   HP2=0;
    if (HP1>1000) HP1=1000;
    if (HP2>1000) HP2=1000;
    float muy = 2*float(HP1)*HP2/(HP1+HP2);
    int newmuy= ceil(muy);
    HP1=HP1-abs(HP2-newmuy);
    HP2=HP2-abs(HP1-newmuy);
    if (HP1<=0)   HP1=0;
    if (HP2<=0)   HP2=0;
    if (HP1>1000) HP1=1000;
    if (HP2>1000) HP2=1000;
    if (HP1>HP2)      return 1;
    else if (HP1<HP2) return -1;
    else              return 0;
}

int acrossTheRiver(int hp1[],int n, int who, float a)
{
// 1 is LTK, 2 is CHINA
    int i = 0;
    int j = 0;
    int LTK = 999;
    int QQ = 777;
    int TT = 666;
    float temp;
    if (who == 1){
        int tuong=-1;
        for (i = 0;i<n;i++){
            if (hp1[i] == LTK){
                for (j = i;j<n; j++){
                    if (hp1[j] == LTK)
                    tuong=j;
                    }
            }
        }
        if (tuong>=0){
            for (i = 0;i<n;i++){
                if (i!=tuong){
                    temp = float(hp1[i])*(1-a*(i+1)/n);
                    hp1[i]=int(temp);        //lam tron xuong
                }
            }
            return 1;
        }
        else {
            for (i = 0;i<n;i++){
                temp = float(hp1[i])*(1-a*(i+1)/n);
                hp1[i]=int(temp);
            }
        }
//        printf("Hp1 la %d %d %d\n",hp1[0],hp1[1],hp1[2]);
        float muy = 0;
        for (i=0;i<n;i++){
            if (hp1[i]>0)
                muy += float(1)/hp1[i];
        }
        muy = float(n)/muy;
        muy=ceil(muy);
        int num = 0;
        for (i=0;i<n;i++){
            if (hp1[i] <= muy)
                num++;
        }
        float trungbinh;
        trungbinh = float(n)/2;
        int trungbinhtinh;
        trungbinhtinh = ceil(trungbinh);
        for (i=0;i<n;i++)
        {
            hp1[i]=hp1[i];
        }
        if (num >= trungbinhtinh)         return 0;
        else                              return 1;
    }
    else // Quan tong qua song-----------------
    {
        int tuong=-1;
        float temp=0;
        for (i = 0;i<n;i++){
            if (hp1[i] == QQ || hp1[i] == TT){
                for (j = i;j<n; j++){
                    if (hp1[j] == QQ || hp1[j] == TT)
                    tuong=j;
                    }
            }
        }
//        printf("%Tuong la %d\n",tuong);
        if (tuong>=0){
            int sumhp = 0;
            for (i = 0;i<n;i++)   {sumhp += hp1[i];}
            float zeta = 0;
            for (i = 0;i<n;i++)   {zeta += float(i)*hp1[i];}
            zeta = zeta/float(sumhp);
            for (i = 0;i<n;i++){
                if (i!=tuong){
                    temp=float(hp1[i])/zeta;
                    temp=ceil(temp);
                    temp=float(temp)*(1-a);
                    hp1[i]=int(temp);
                }
            }
            return 1;
        }
//        printf("Zeta la %f\n", zeta);
        else {
            int sumhp = 0;
            for (i = 0;i<n;i++)   {sumhp += hp1[i];}
            float zeta = 0;
            for (i = 0;i<n;i++)   {zeta += float(i)*hp1[i];}
            zeta = zeta/float(sumhp);
            for (i=0;i<n;i++){
                temp=float(hp1[i])/zeta;
                temp=ceil(temp);
                temp=float(temp)*(1-a);
                hp1[i]=int(temp);
            }
        }
//        printf("HP2 sau do la %d %d %d\n",hp2[0],hp2[1],hp2[2]);
        float muy = 0;
        for (i=0;i<n;i++){
            if (hp1[i]>0)
                muy += float(1)/hp1[i];
        }
        muy = float(n)/muy;
        muy=ceil(muy);
//        printf("muy la %f\n",muy);
        int num = 0;
        for (i=0;i<n;i++){
            if (hp1[i] <= muy)             //can than, cho int so sanh float
                num++;
        }
        float trungbinh;
        trungbinh = float(n)/2;
        int trungbinhtinh;
        trungbinhtinh = ceil(trungbinh);
//        printf("Trung binh tinh la %d\n",trungbinhtinh);
        if (num >= trungbinhtinh)         return 0;
        else                              return 1;
    }
}

int serialFight(int **matrix1, int **matrix2,int row,int col,int u,float a)
{

    int i=0;int j=0; float temp=0;
    int LTK=999; int QQ =777; int TT = 666;
    int abs_u=abs(u);
    int donvi= abs_u%10;
    int chuc = abs_u/10;
    int siso= col*row-1;

//TIM VI TRI TUONG CUA QUAN TA VA QUAN CHINA
    int tuong1x=-1; int tuong1y= -1; int tuong2x=-1; int tuong2y=-1;
    for (i=0;i<row;i++){
      for (j=0;j<col;j++){
        if (matrix1[i][j]==LTK){
            if (i>=tuong1x) {tuong1x=i; tuong1y=j;                          }
            else if (i==tuong1x && j>tuong1y) {tuong1x=i; tuong1y=j;        }
        }
        if (matrix2[i][j]==QQ || matrix2[i][j]==TT){
            if (i>tuong2x) {tuong2x=i; tuong2y=j;                           }
            else if (i==tuong2x && j>tuong2y) {tuong2x=i; tuong2y=j;        }
        }
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//    printf("Vi tri tuong VN la %d %d\n",tuong1x,tuong1y);
//    printf("Vi tri tuong Tau Khua la %d %d\n\n",tuong2x,tuong2y);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//-----------------------------------------------------------------------------------
// CHINA TAN CONG
    if (u <0){
      int vuotsong[siso+1];
// QUAN TAU KHUA PHAI VUOT SONG NEN CHUYEN MATRIX THANH 1 HANG
      if (chuc == 2) { //Co vuot song
        for (i=0;i<row;i++){
            for (j=0;j<col;j++) {
                vuotsong[i*col+j]=matrix2[i][j];
            }
        }

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//      printf("Quan Tong sau khi xep hang: \n");
//      for (i=0;i<=siso;i++){ if (i==siso) printf("%d\n\n",vuotsong[i]); else printf("%d ",vuotsong[i]);}
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// XET XEM CO TUONG HONG, TINH DIEU KIEN THOI TIET VAO
      if (tuong2x>=0 && u<0){
            int sumhp = 0;
            for (i = 0;i<=siso;i++)   {sumhp += vuotsong[i];}
            float zeta = 0;
            for (i = 0;i<=siso;i++)   {zeta += float(i)*vuotsong[i];}
            zeta = zeta/float(sumhp);                                              //printf("Zeta =%f\n",zeta);
            for (i = 0;i<=siso;i++){
                if (i!=tuong2x*col+tuong2y){ //VI TRI CUA TUONG KHI VUOT SONG
                    temp=float(vuotsong[i])/zeta;
                    temp=ceil(temp);
                    temp=float(temp)*(1-a);
                    vuotsong[i]=int(temp);
                    if (vuotsong[i] <0) vuotsong[i]=0;
                }
            }
        }
      else {
            int sumhp = 0;
            for (i = 0;i<=siso;i++)   {sumhp += vuotsong[i];}
            float zeta = 0;
            for (i = 0;i<=siso;i++)   {zeta += float(i)*vuotsong[i];}
            zeta = zeta/float(sumhp);                                               //printf("Zeta =%f\n",zeta);
            for (i=0;i<=siso;i++){
                temp=float(vuotsong[i])/zeta;
                temp=ceil(temp);                                                   // printf("Sau khi chia zeta lam tron len: %f\n",temp);
                temp=float(temp)*(1-a);
                vuotsong[i]=int(temp);
                if (vuotsong[i] <0) vuotsong[i]=0;
            }
      }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//    printf("Sau khi vuot song la:\n");
//    for (i=0;i<=siso;i++){ if (i==siso) printf("%d\n\n",vuotsong[i]); else printf("%d ",vuotsong[i]);}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// SAU KHI QUA SONG, NEU DOC NAM QUOC SON HA
      if (donvi==1){
          for (i = 0;i<=siso;i++){
            if (i!=tuong2x*col+tuong2y){ // Dieu kien quan khac tuong
              temp=float (vuotsong[i])*10.0/12;                                         //printf("Tong chia zeta %f\n", temp);
              vuotsong[i]=ceil(temp);                                                   //printf("Tong lam tron: %d\n",vuotsong[i]);
              if (vuotsong[i]<0) vuotsong[i]=0; // neu HP lon hon 1000 thi HP=1000
            }
          }
          for (i=0;i<row;i++){
            for (j=0;j<col;j++){
              if ( i==tuong1x && j==tuong1y ) {      // Dieu kien quan khac tuong
                matrix1[i][j]=matrix1[i][j];
                }
              else {
                temp=float (matrix1[i][j])*15/10;
                matrix1[i][j]=ceil(temp);
                if (matrix1[i][j]<0 ) {matrix1[i][j]=0;} // Neu HP nho hon 0 thi HP=0
                if (matrix1[i][j]>1000) matrix1[i][j]=1000; // neu HP lon hon 1000 thi HP=1000
                                                                                        //printf("Quan ta buff 1.5 lan: %d\n",matrix1[i][j]);
              }
            }
          }
        }

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//        printf("Quan Tong chia 1.2:\n");
//        for (i=0;i<=siso;i++){ if (i==siso) printf("%d\n\n",vuotsong[i]); else printf("%d ",vuotsong[i]);}
//        for (i=0;i<row;i++){
//            for (j=0;j<col;j++){
//              printf("Quan ta buff 1.5 lan: %d\n",matrix1[i][j]);
//            }
//        }
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// CHUYEN QUAN TAU KHUA THANH MA TRAN LAI
      for (i=0;i<row;i++){
          for (j=0;j<col;j++){
            int thamso=i*col+j;
            matrix2[i][j]=vuotsong[thamso];
          }
        }
      }

// CHEM GIET NHAU------------------------------------
      for (i=0;i<row;i++){
        for (j=0;j<col;j++){
          int newmuy=0; int HP1=0; int HP2=0;
          HP1=matrix1[i][j];
          HP2=matrix2[i][j];
          float newHP1, newHP2;
          if (i== tuong1x && j ==tuong1y){
                HP1=HP1;
          }
          else{
          newHP1=float(HP1)*11/10;
          HP1=ceil(newHP1);
          if (HP1<0)   HP1=0;
          if (HP1>1000) HP1=1000;
          }                                        //printf("Quan ta la: %d\n",HP1);
          if (i== tuong2x && j ==tuong2y){
                HP2=HP2;
          }
          else{
          newHP2=float(HP2)*13/10;
          HP2=ceil(newHP2);
          if (HP2<0)   HP2=0;
          if (HP2>1000) HP2=1000;
          }                                                 //printf("Quan dich la: %d\n",HP2);
          float muy = 2*float(HP1)*HP2/(HP1+HP2);
          newmuy= ceil(muy);                           //printf("newmuy la: %d\n",newmuy);
// LY THUONG KIET CHEM CHINA
          if ((i== tuong1x) && (j ==tuong1y) && (i!=tuong2x)){
              matrix1[i][j]=HP1;
              matrix2[i][j]=0;
//              printf("LTK chem CHINA\n");
              continue;
          }
// QQ HOAC TT CHEM VN
          else if (i==tuong2x && j == tuong2y){
            if (i!=tuong1x){
              matrix1[i][j]=0;
              matrix2[i][j]=HP2;
//              printf("CHINA chem VN");
              continue;
            }
          }
// CHEM GIET NHU BINH THUONG
          else  {
            int HP1temp=HP1;
            HP1=HP1-abs(HP2-newmuy);
            HP2=HP2-abs(HP1temp-newmuy);
            if (HP1<0)   HP1=0;
            if (HP2<0)   HP2=0;
            if (HP1>1000) HP1=1000;
            if (HP2>1000) HP2=1000;
          }
          matrix1[i][j]=HP1;
          matrix2[i][j]=HP2;
//          printf("tai vi tri %d %d co hai linh %d %d\n",i,j,HP1,HP2);
        }
      }
    }

//-------------------------------------------------------------------------------------//

// QUAN TA TAN CONG
    else {
      int vuotsong[siso+1];
// QUAN VN PHAI VUOT SONG NEN CHUYEN MATRIX THANH 1 HANG
      if (chuc == 2) { //Co vuot song
        for (i=0;i<row;i++){
            for (j=0;j<col;j++) {
                vuotsong[i*col+j]=matrix1[i][j];
            }
        }

// XET XEM CO TUONG HONG, TINH DIEU KIEN THOI TIET VAO
      if (tuong1x>=0 && u<0){
            for (i = 0;i<=siso;i++){
                if (i!=tuong1x*col+tuong1y){ //VI TRI CUA TUONG KHI VUOT SONG
                    temp = float(vuotsong[i])*(1-a*(i+1)/siso);
                    vuotsong[i]=int(temp);        //lam tron xuong
                    if (vuotsong[i] <0) vuotsong[i]=0;
                }
            }
        }
      else {
            for (i = 0;i<=siso;i++){
                temp = float(vuotsong[i])*(1-a*(i+1)/siso);
                vuotsong[i]=ceil(temp);
                if (vuotsong[i] <0) vuotsong[i]=0;
            }
        }

// SAU KHI QUA SONG, NEU DOC NAM QUOC SON HA
      if (donvi==1){
          for (i = 0;i<=siso;i++){
            if (i!=tuong1x*col+tuong1y){
              temp=float (vuotsong[i])*15/10;
              vuotsong[i]=ceil(temp);
            }
          }
          for (i=0;i<row;i++){
            for (j=0;j<col;j++){
              if ( i==tuong2x && j==tuong2y ) {
                matrix2[i][j]=matrix2[i][j];
                }
              else {
                temp=float (matrix2[i][j])*10/12;
                matrix2[i][j]=ceil(temp);
                if (matrix2[i][j]<0 ) {matrix2[i][j]=0;}
                if (matrix2[i][j]>1000) matrix2[i][j]=1000;
              }
            }
          }
        }

// CHUYEN QUAN VN THANH MA TRAN LAI
      for (i=0;i<row;i++){
          for (j=0;j<col;j++){
            int thamso=i*col+j;
            matrix1[i][j]=vuotsong[thamso];
          }
        }
      }

// CHEM GIET NHAU, VN CO LOI----------------------------------------------------
      for (i=0;i<row;i++){
        for (j=0;j<col;j++){
          int newmuy=0; int HP1=0; int HP2=0;
          HP1=matrix1[i][j];
          HP2=matrix2[i][j];
          float newHP1, newHP2;
          if (i== tuong1x && j ==tuong1y){
                HP1=HP1;
          }
          else{
          newHP1=float(HP1)*15/10;
          HP1=ceil(newHP1);
          if (HP1<0)   HP1=0;
          if (HP1>1000) HP1=1000;
          }                                        //printf("Quan ta la: %d\n",HP1);
          if (i== tuong2x && j ==tuong2y){
                HP2=HP2;
          }
          else{
          newHP2=float(HP2)*10/13;
          HP2=ceil(newHP2);
          if (HP2<0)   HP2=0;
          if (HP2>1000) HP2=1000;
          }                                                 //printf("Quan dich la: %d\n",HP2);
          float muy = 2*float(HP1)*HP2/(HP1+HP2);
          newmuy= ceil(muy);                           //printf("newmuy la: %d\n",newmuy);
// LY THUONG KIET CHEM CHINA
          if (i== tuong1x && j ==tuong1y){
            if (i!=tuong2x ){
              matrix1[i][j]=HP1;
              matrix2[i][j]=0;
              continue;
            }
          }
// QQ HOAC TT CHEM VN
          else if (i==tuong2x && j == tuong2y){
            if (i!=tuong1x){
              matrix1[i][j]=0;
              matrix2[i][j]=HP2;
              continue;
            }
          }
// CHEM GIET NHU BINH THUONG
          else {
            int HP1temp=HP1;
            HP1=HP1-abs(HP2-newmuy);
            HP2=HP2-abs(HP1temp-newmuy);
            if (HP1<0)   HP1=0;
            if (HP2<0)   HP2=0;
            if (HP1>1000) HP1=1000;
            if (HP2>1000) HP2=1000;
          }
          matrix1[i][j]=HP1;
          matrix2[i][j]=HP2;
        }
      }
    }

//---------------------------XET AI THANG AI THUA-----------------------------//
    int dem1=0; int dem2 =0;
    for (i=0;i <row; i++){
      for (j=0;j <col;j++){
        if (matrix1[i][j] >0) dem1++;
        if (matrix2[i][j] >0) dem2++;
      }
    }
    if (dem1>dem2) return 1;
    else if (dem1<dem2) return 0;
    else return 0;
}

int main(int argc, char** argv) {

//    For test 1
//    int n1 = 996;
//    int n2 = 997;
//    int n = 1;
//    int result = fightOneByOne(n1, n2, n);
//    printf("Fight result: %d\n", result);
//    printf("HP1: %d\n", n1);
//    printf("HP2: %d\n", n2);

//  For test 3
    int row =3;
    int col =4;
    int u=-10;
    float a=0.7;
    int **matrix1; int **matrix2;
    initializeMatrix(matrix1,row,col);
    initializeMatrix(matrix2,row,col);
//    matrix1[0][0]  = 162; matrix1[0][1]  = 157; matrix1[0][2]  = 139; matrix1[0][3]  = 555;
//    matrix1[1][0]  = 749; matrix1[1][1]  = 999; matrix1[1][2]  = 753; matrix1[1][3]  = 444;
//
//    matrix2[0][0]  = 76; matrix2[0][1]  = 324; matrix2[0][2]  = 567; matrix2[0][3]  = 567;
//    matrix2[1][0]  = 112; matrix2[1][1]  = 987; matrix2[1][2]  = 15; matrix2[1][3]  = 595;

//    matrix1[0][0]  = 157; matrix1[0][1]  = 428; matrix1[0][2]  = 32;
//    matrix1[1][0]  = 759; matrix1[1][1]  = 999; matrix1[1][2]  = 73;
//
//    matrix2[0][0]  = 76; matrix2[0][1]  = 324; matrix2[0][2]  = 567;
//    matrix2[1][0]  = 112; matrix2[1][1]  = 987; matrix2[1][2]  = 15;

    matrix1[0][0]  = 162; matrix1[0][1]  = 157; matrix1[0][2]  =139; matrix1[0][3]  =555;
    matrix1[1][0]  = 749; matrix1[1][1]  = 999; matrix1[1][2]  =753; matrix1[1][3]  =444;
    matrix1[2][0]  = 315; matrix1[2][1]  = 777; matrix1[2][2]  =999; matrix1[2][3]  =112;

    matrix2[0][0]  = 76 ; matrix2[0][1]  = 324; matrix2[0][2]  = 567;matrix1[0][3]  =567;
    matrix2[1][0]  = 112; matrix2[1][1]  = 987; matrix2[1][2]  = 15; matrix1[1][3]  =595;
    matrix2[2][0]  = 112; matrix2[2][1]  = 666; matrix2[2][2]  = 998;matrix1[2][3]  =757;

    printf("HP TRUOC KHI CHEM GIET:\n");
    printMatrix(matrix1,row, col);
    printf("\n");
    printMatrix(matrix2,row, col);
    int ketqua = serialFight(matrix1,matrix2,row,col,u,a);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    printf("ket qua la %d\n",ketqua);
    printf("HP SAU KHI CHEM GIET:\n");
    printMatrix(matrix1,row, col);
    printf("\n");
    printMatrix(matrix2,row, col);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    return 0;
}

void initializeMatrix(int **&matrix, int row, int col)
{
    matrix = new int *[row];
    for (int i = 0; i < row; i++)
    {
        matrix[i] = new int[col];
        for (int j = 0; j < col; j++)
        {
            matrix[i][j] = 0;
        }
    }
}

void printMatrix(int **&matrix, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}
