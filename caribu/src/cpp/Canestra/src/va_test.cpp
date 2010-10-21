
#include <cstdio>
using namespace std;

#include <cstdarg>
#include <cstring>

int f(int first, int second, int third,...){
  va_list ptarg;
  int coeff=1;
  int indice;
  
  
  indice=first+10*second+100*third;
  va_start(ptarg,third);  
  for(register unsigned char i=3;i<4;i++){
    //indice+=coeff*bonind(va_arg(ptarg,int),i);
    indice+=1000*va_arg(ptarg,int);
  }   
  va_end(ptarg);
  return(indice);
}

int f(int first, int second){
  return(first+10*second);
}

main(){

  printf("f(1,2)     = %d\n",f(1,2));
  printf("f(1,2,3)   = %d\n",f(1,2,3));
  printf("f(1,2,3,4) = %d\n",f(1,2,3,4));
  
}//main()
