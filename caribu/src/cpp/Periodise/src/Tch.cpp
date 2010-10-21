// g++ -g -W -Wall -o Tch{,.C} ; Tch

#include<cstdio>
using namespace std;

#include <ctype.h>

char * next_number(char *ch){
  char c;
  int i; 
  
  for (i=0;  (c=ch[i])>0  && isspace(c)!=0; i++) printf(">> ch[%d]=%c\n",i,c); 
  return &(ch[i]); 

}
char * next_nonumber(char *ch){
  char c;
  int i; 
  
  for (i=0; (c=ch[i])>0 && isspace(c)==0; i++) printf(">> ch[%d]=%c\n",i,c); 
  return &(ch[i]); 

}


int main(){
  char line[50], *ch;
  int i; 
  float x; 

  sprintf(line, "3  1.2  \t 3.7 \t\t 1e8\n%c",0);
  printf("==> line=%s\n\n",line);

  ch=line; 
  
  // for (i=0; i<4; i++){
  for (i=0; ch[i]>0; i++){   
    ch=next_number(ch);
    sscanf(ch,"%f",&x);
    printf("==> x=%.2f :: line=%s\n",x,ch);
    ch=next_nonumber(ch);
  }

}


/* Run Time

==> line=3  1.2          3.7             1e8


==> x=3.00 :: line=3  1.2        3.7             1e8

>> ch[0]=3
>> ch[0]= 
>> ch[1]= 
==> x=1.20 :: line=1.2           3.7             1e8

>> ch[0]=1
>> ch[1]=.
>> ch[2]=2
>> ch[0]= 
>> ch[1]= 
>> ch[2]=
>> ch[3]= 
==> x=3.70 :: line=3.7           1e8

>> ch[0]=3
>> ch[1]=.
>> ch[2]=7
>> ch[0]= 
>> ch[1]=
>> ch[2]=
>> ch[3]= 
==> x=100000000.00 :: line=1e8

>> ch[0]=1
>> ch[1]=e
>> ch[2]=8
                  

*/
