#ifndef __T_UTILITAIRES_H__
#define __T_UTILITAIRES_H__

#include <iostream>
using namespace std ;

#include "ferrlog.h"


#ifndef __GNUG__
#include "bool.h"
#include "verbose.h"
#endif

#include "transf.h"
// modif Liste -> liste chainee double - pr QuickSort - Mike 95 

//-************* Tableau Dynamique *** Mike 94 *******************
//typedef __gnuc_va_list va_list;
#include <stdarg.h>
//#include "/usr/local/sparc/lib/gcc-lib/sparc-sun-solaris2.5.1/2.7.2.2.f.2/include/stdarg.h"
//#include "/usr/include/stdarg.h"
#include <stdio.h>
#include <string.h>
/*****************************************************************
           Attention Vorsicht Warning Attenzione 
    L'utilisation des macros de <stdarg.h> empeche le compilateur
    de controler le nombre de parametres des fonctions de cette
    classe / au nombre de dimensions D.
    Le programmeur devra etre tres RIGOUREUX quant aux appels de :
       --> Tabdyn()   : le constructeur
       --> alloue()   : le constructeur a posteriori
       --> operator (): l'operateur d'acces au tableau 
*******************************************************************/
template <class Type, unsigned int D>
class Tabdyn{
  unsigned long taille;
  unsigned int max[D];
  Type *tab;
   bool trie;
public:
  Tabdyn()
    { trie=false; tab=NULL;}
  Tabdyn(int first,...);
  ~Tabdyn(){delete [] tab;}
  void alloue(int first,...);
  void free(){if (tab != NULL) {delete [] tab ;}  tab=NULL;}
  Type& operator()(int first);
  Type& operator()(int first,int second);
  Type& operator()(int first,int second,int third,...);
  Tabdyn<Type,D>& operator =(Tabdyn<Type,D>& rval);
  void maj(Type val);
  unsigned int dim(){return D;}
  unsigned int * maxi(){return max;}
  //Vorsicht : implemente only for 1D - MC
  void  quicksort(int (*clef)(Type *x, Type *y))
    { if(D>1)
      { Ferr <<"Tabdyn[quicksort] implemente que pour les tableaux 1D"<<'\n';
      exit(34);
      }
    if(!trie)
      { trie=true;
      qsort(tab,taille,sizeof(Type),(int (*)(const void *, const void *)) clef);
      }
    }
  
private:  // fonction privee     
   int bonind(int  ind, unsigned char indmax);
   int bonmax(int  ind);
  void lectarg( int& first,va_list & ptarg);    
};

//************ operator = *****************
template <class Type, unsigned int D>   
inline Tabdyn<Type,D>& Tabdyn<Type,D>::operator =(Tabdyn<Type,D>& rval){
  for(register unsigned char i=0;i<D ;i++)
    if(max[i]!=rval.max[i]){
      Ferr<<"\n Tabdyn [operator =] erreur dimension "<<i<<'\n' ;
      exit(35);
    }
  trie=false;
  memcpy(tab,rval.tab,taille*sizeof(Type *));
    
  return *this;
}

//************ maj() *****************
template <class Type, unsigned int D>   
void Tabdyn<Type,D>::maj(Type val){
  if(tab==NULL) {
    Ferr<<"\nTabdyn [maj()]: Tableau dynamique pas alloue: acces impossible!";
    Ferr<<'\n';
    exit(36);
  }
  unsigned long i;
  Type *pt;
  trie=false;
  for(i=0,pt=tab;i<taille;i++,pt++) *pt=val;
}// maj()
//*********** bonind() *****************
template <class Type, unsigned int D>
inline  int  Tabdyn<Type,D>::bonind(int  ind,unsigned  char indmax){
  if(ind<0 || ind>=max[indmax]){
    Ferr<<"\n Tabdyn[operator()]: indice no "<<ind<<" hors-borne!"<<'\n';
    ind=(int)  (1/sin(0.0));
    exit(37);
    return ind; // pour eviter les warning verbeux de SGI
  }
  else  return ind;
}//bonind()
//*********** bonmax() *****************
template <class Type, unsigned int D>
inline   int  Tabdyn<Type,D>::bonmax(int ind)      {
  if(ind<=0) {
    Ferr<<"\n Tabdyn[bornmax()]: dimension de tableau negative ou nulle!"<<'\n';
    exit(38);
    return ind; // pour eviter les warning verbeux de SGI
  }
  else  return ind;
}//bonmax()
//*********** lectarg() *****************
template <class Type, unsigned int D>
inline void  Tabdyn<Type,D>::lectarg ( int &first,va_list & ptarg ){
  // va_start(ptarg,first);
  taille=max[0]=bonmax(first);
  for(register unsigned char  i=1;i<D;i++){
    max[i]=bonmax(va_arg(ptarg,int));
    taille*=max[i];
    // cout<<"taille"<< taille<<"\n";
  }   
  tab=new Type[taille];
  if(tab==NULL) {
    Ferr<<"Plus de mem pour alloc dyn!"<<'\n';
    exit(39);
  }
  va_end(ptarg);
}//lectarg()

//*********** Tabdyn(x1,...) *****************
template <class Type, unsigned int D>
Tabdyn<Type,D>::Tabdyn(int first,... ){
  va_list ptarg;
  va_start(ptarg,first);
  tab=NULL; trie=false; 
  lectarg (first,ptarg);        
}//Tabdyn(x1,...)

//-*********** alloue() *****************
template <class Type, unsigned int D>
void Tabdyn<Type,D>::alloue( int  first,... ){
  if(tab!=NULL) {
    Ferr<< "\nTabdyn [alloue()]: Tableau dynamique deja alloue a la declaration!"<<'\n';
    exit(40);
  }
  va_list ptarg;
  va_start( ptarg,first);
  lectarg (first, ptarg);             
}//alloue()

//-*********** operator (1) *****************
template <class Type,unsigned  int D>
Type & Tabdyn<Type,D>::operator()(int first){
  if(D==1)
    return(tab[first]);
  else {
    Ferr <<"Tabdyn<Type,D>::operator()(int first) incompatbile avec D=";
    Ferr  << D<<'\n' ;
    exit(41);
    return (tab[0]);// pour eviter les wranings verbeux de SGI !!!
  } 
}//operator ()
//-*********** operator (1,2) *****************
template <class Type,unsigned  int D>
Type & Tabdyn<Type,D>::operator()(int first, int second){
  if(D==2)
    return(tab[first+max[0]*second]);
  else  {
    Ferr <<"Tabdyn<Type,D>::operator()(int first,int second) incompatbile";
    Ferr<<" avec D="  << D<<'\n' ;
    exit(42);
    return(tab[0]);// pour eviter les wranings verbeux de SGI !!!
  } 
}//operator ()
//*********** operator (1,2,3,...) *****************
template <class Type,unsigned  int D>
Type & Tabdyn<Type,D>::operator()(int first, int second, int third,...){
  /*  if(tab==NULL) {
    Ferr<< "Tabdyn [operator()]: Tableau dynamique pas alloue : acces impossible!"<<'\n';
    exit(43);
  }
  */
  va_list ptarg;
  int coeff=1;
  int indice;
   
  //indice=bonind(first,0);
  coeff=max[0]*max[1];
  indice=first+max[0]*second+coeff*third;
  va_start(ptarg,third);  
  for(register unsigned char i=3;i<D;i++){
    coeff*=max[i-1];
    //indice+=coeff*bonind(va_arg(ptarg,int),i);
    indice+=coeff*va_arg(ptarg,int);
  }   
  va_end(ptarg);
  return(tab[indice]);
}//operator ()
//*******FIN****ENDE****FINO****END**** Tabdyn****---

enum Booleen {FAUX, VRAI};

// LISTE simple

template <class Type>
class Noeud{
  Type donnee;
  Noeud<Type> *suivant;
  
public:
  //Noeud (Type d=0,Noeud<Type>*n=NULL);
  Noeud (Type d,Noeud<Type>*n);
  Noeud (Type d);
  Noeud ();
  Type& donne(); // renvoie la donnee du noeud
  void donne(Type d); // initialise la donnee du noeud a la valeur d
  inline Noeud<Type>* next(); // renvoie le pointeur sur suivant du noeud
  inline void next(Noeud<Type>* n); //initialise le pointeur sur suivant du noeud a la valeur n
  //show(){cout<<"Noued="<<donnee<<endl;}
};

template <class Type>
class Liste{
private:
  Noeud<Type> *premier;
  Noeud<Type> *courant;
  unsigned int nbe;
  
public:
  Liste() {
    premier=NULL; courant=NULL;  nbe=0;}
  Liste(char a) {
   premier=NULL; courant=NULL;  nbe=0;} 
  ~Liste();     
  void ajoute (Type); // ajoute un element dans la liste
  void init_liste (Type); // initialise la liste
  void ajout_elem (Type); // ajoute un element dans la liste initialisee
  inline Booleen est_fin(); // indique si c'est la fin de la liste
  inline Booleen finito(); // indique si courant=NULL
  inline Booleen est_vide(); // indique si la liste est vide
  inline void debut(); // met le pointeur courant au debut de la liste
  inline void fin(); // met le pointeur courant a la fin de la liste
  inline void suivant(); // met le pointeur courant sur l'element suivant par rapport
  // a la position courante
  inline Type& contenu(); // renvoie le contenu d'un element de la liste
  inline void detruire(); // detruit l'element courant
  inline void detruire_droite(); // detruit l'element a droite de l'element courant
  inline void free_liste(); // libere la liste
  void ranger(Type);
  unsigned int card() // donne le nbre d'elements (cardinal) 
  {return nbe;}
};// class Liste

//template <class Type> inline void detruire_contenu(Liste<Type *>&);

// Liste double [Mike95]

template <class Type>
class NoeudD{
  Type donnee;
  NoeudD<Type> *suivant;
  NoeudD<Type> *prec; //precedent
public:
  NoeudD (Type d);
  NoeudD (Type d,NoeudD<Type>*n);  
  Type& donne(); // renvoie la donnee du NoeudD
  void donne(Type d); // initialise la donnee du NoeudD a la valeur d
  inline NoeudD<Type>* next(); // renvoie le pointeur sur suivant du noeud
  inline void next(NoeudD<Type>* n); //initialise le pointeur sur suivant du noeud a la valeur n
  inline NoeudD<Type>* back(); // renvoie le pointeur sur prec du noeud
  inline void back(NoeudD<Type>* n); //initialise le pointeur sur prec du noeud a la valeur n 
  inline NoeudD& operator =(NoeudD &);
}; //NoeudD

template <class Type>
class ListeD{
  NoeudD<Type> *premier;
  NoeudD<Type> *courant;
  NoeudD<Type> *pause;
  unsigned int nbe; // Nbre d'Elements
  
public:
  ListeD();
  ~ListeD();
  void ajoute (Type); // ajoute un element dans la liste
  void init_liste (Type); // initialise la liste
  void ajout_elem (Type); // ajoute un element dans la liste initialisee
  inline bool est_fin(); // indique si c'est la fin de la liste
  inline bool finito(); // indique si courant=NULL
  inline bool debuto(); // indique si courant=NULL
  inline bool est_debut(); // indique si courant=premier
  inline bool est_vide(); // indique si la liste est vide
  inline void debut(); // met le pointeur courant au debut de la liste
  inline void fin(); // met le pointeur courant a la fin de la liste
  void deb_pause() {pause=courant;}//marche a l'ordre 2
  void fin_pause() {courant=pause;}
  inline void suivant(); // met le pointeur courant sur l'element suivant par rapport  a la position courante
  inline void precedent(); // met le pointeur courant sur l'element precedent par rapport  a la position courante
  inline Type& contenu(); // renvoie le contenu d'un element de la liste
  inline void detruire(); // detruit l'element courant
  inline void detruire_droite(); // detruit l'element a droite de l'element courant
  inline void free_liste(); // libere la liste
  void QuickSort(int (*)(Type,Type));// tri rapide de la liste
  unsigned int card() // donne le nbre d'elements (cardinal) 
  {return nbe;}
private:
  void qsort( NoeudD<Type> *, NoeudD<Type> *,int (*)(Type,Type)); // Fonction recursive utilisee par Quicksort
  inline void invert( NoeudD<Type> *, NoeudD<Type> *);// intervertit deux elements de la liste
};// class ListeD

template <class Type>
inline void detruire_contenu(ListeD<Type *>&);


//***********************************************
// Fin des Declarations - Debut des definitions//
//***********************************************

// Noeud

template <class Type>
//Noeud<Type> :: Noeud (Type d,Noeud<Type>*n=NULL)
Noeud<Type> :: Noeud (Type d,Noeud<Type>*n){
	donnee=d;
	suivant=n;
}
template <class Type>
Noeud<Type> :: Noeud (Type d){
	donnee=d;
	suivant=NULL;
}
template <class Type>
Noeud<Type> :: Noeud (){
	donnee=0;
	suivant=NULL;
}

template <class Type>
inline Type& Noeud<Type> :: donne()
{
	return donnee;
}

template <class Type>
inline void Noeud<Type> :: donne(Type d)
{
	donnee=d;
}
template <class Type>
inline Noeud<Type>* Noeud<Type> :: next()
{
	return suivant;
}

template <class Type>
inline void Noeud<Type> :: next(Noeud<Type>* n)
{
	suivant=n;
}


// LISTE

/*
template <class Type>
Liste<Type>::Liste(){
  premier=NULL;
  courant=NULL;
  nbe=0;
}
*/
template <class Type>
Liste<Type> :: ~Liste()
{
	free_liste();
}

template <class Type>
inline void Liste<Type> :: ajoute (Type contenu)
 { if (est_vide() == VRAI)
      init_liste(contenu);
    else
      ajout_elem(contenu);
 }
template <class Type>
inline void Liste<Type> :: init_liste (Type contenu)
{
	premier=new Noeud<Type>(contenu);
//	premier=new Noeud<Type>;
//	premier->donne(contenu);
//	premier->next(NULL);
	courant=premier;
    nbe++;
}

template <class Type>
inline void Liste<Type> :: ajout_elem (Type contenu)
{
	Noeud<Type>* p=new Noeud<Type>(contenu, courant->next());

	courant->next(p);
	courant=p;
   nbe++;  
}

template <class Type>
inline Booleen Liste<Type> :: est_fin()
{
  if (courant->next() == NULL)
	return VRAI;
    else  return FAUX;
}
template <class Type>
inline Booleen Liste<Type> :: finito()
{
    return((Booleen)(courant==NULL?VRAI:FAUX));
}


template <class Type>
inline Booleen Liste<Type> :: est_vide()
{
   if (premier == NULL) 
         return VRAI;
    else return FAUX;
}

template <class Type>
inline void Liste<Type> :: debut(){
  courant=premier;
}

template <class Type>
inline void Liste<Type> :: fin(){
  if(finito())
    debut();
  while (!est_fin()){
    courant=courant->next();
  }
}

template <class Type>
inline void Liste<Type> :: suivant()
{
//	if(!est_fin())
	if(!finito())
		courant=courant->next();
}


template <class Type>
inline Type& Liste<Type> :: contenu()
{
	return (courant->donne());
}

template <class Type>
inline void Liste<Type>::detruire(){//fucke', revu MC9/10/96
  //cout << "debut de detruire\n"; cout.flush();
  if (!est_vide()){
    // cout << "debut du test dans detruire\n"; cout.flush();
    if (courant==premier){
      // cout << "courant == premier\n"; cout.flush();
      courant=premier->next();
      // cout << "destruction de premier\n" << premier << endl; cout.flush();
      delete premier;
      premier=courant;
      nbe--;
    }
    else{
      //printf("courant!= premier\n");
      //courant->show();
      Noeud<Type>* prec;
      // cout << "declaration du noeud puis boucle\n"; cout.flush();
    
      for (prec=premier; prec->next()!=courant; prec=prec->next());
      prec->next(courant->next());
      // cout << "destruction du courant" << courant << endl; cout.flush();
      delete courant;
      courant=prec;
      nbe--;
    }
  }
  else{
    Ferr << "ERREUR - Impossible de detruire un element car liste deja vide";
    Ferr<<'\n';
    exit (1); // risque de pas voir le msg d'erreus
  }
}

template <class Type>
inline void Liste<Type> :: detruire_droite()
{
	Noeud<Type>* tmp;

	if (!est_vide())
	{
		if (!est_fin())
		{
			tmp=courant->next();
			courant->next(tmp->next());
			delete courant->next();
                        nbe--;
		}
		else
		{
			Ferr << "ERREUR - Impossible de detruire un element car fin de liste"<<'\n';
			exit (44);
		}
	}
	else
	{
		Ferr << "ERREUR - Impossible de detruire un element car liste vide"<<'\n';
		exit (45);
	}
}

template <class Type>
inline void Liste<Type> :: free_liste()
{
	for (debut(); !est_vide(); detruire());
}

template <class Type>
void Liste<Type> :: ranger(Type element)
{
	if (est_vide() == VRAI)
		init_liste(element);
	else
		ajout_elem(element);
    nbe++;
}

template <class Type>
inline void detruire_contenu(Liste<Type*>& liste)
{
	Type *donnee;
	liste.debut();
	while (!liste.est_fin())
	{
		donnee=liste.contenu();
		delete donnee;
		liste.suivant();
	}
}
//-***************** Liste Double *************************************

// NoeudD

template <class Type>
NoeudD<Type>::NoeudD (Type d,NoeudD<Type>*n)
{
	donnee=d;
	suivant=n->next();
	prec=n;
}
template <class Type>
NoeudD<Type>::NoeudD (Type d)
{
	donnee=d;
	suivant=NULL;
	prec=NULL;
}
template <class Type>
inline Type& NoeudD<Type>::donne()
{
	return donnee;
}

template <class Type>
inline void NoeudD<Type>::donne(Type d)
{
	donnee=d;
}
template <class Type>
inline NoeudD<Type>* NoeudD<Type>::next()
{
	return suivant;
}

template <class Type>
inline void NoeudD<Type>::next(NoeudD<Type>* n)
{
	suivant=n;
}

template <class Type>
inline NoeudD<Type>* NoeudD<Type>::back()
{
	return prec;
}

template <class Type>
inline void NoeudD<Type>::back(NoeudD<Type>* n)
{
	prec=n;
}
template <class Type>
inline NoeudD<Type>& NoeudD<Type>::operator =(NoeudD<Type> & rval){
  donnee=rval.donnee;
  suivant=rval.suivant;
  prec=rval.prec;
  
  return *this; 
}
// ListeD

template <class Type>
ListeD<Type>::ListeD()
{
	premier=NULL;
	courant=NULL;
        nbe=0;
}

template <class Type>
ListeD<Type>::~ListeD()
{
	free_liste();
}

template <class Type>
inline void ListeD<Type>::ajoute (Type contenu)
 { if (est_vide() == true)
      init_liste(contenu);
    else
      ajout_elem(contenu);
}
template <class Type>
inline void ListeD<Type>::init_liste (Type contenu)
{
  premier=new NoeudD<Type>(contenu);
  courant=premier;
  nbe++;
}

template <class Type>
inline void ListeD<Type>::ajout_elem (Type contenu)
{
  NoeudD<Type>* p=new NoeudD<Type>(contenu, courant);
  if(courant->next()!=NULL)
    (courant->next())->back(p);
  courant->next(p);
  courant=p;
  nbe++;  
}

template <class Type>
inline bool ListeD<Type>::est_fin()
{
  if (courant->next() == NULL)
    return true;
  else
    return false;
}

template <class Type>
inline  bool ListeD<Type>::finito()
{
    return( (courant==NULL) ?true:false);
}
template <class Type>
inline bool ListeD<Type>::est_debut()
{
    return( (courant==premier));
}
template <class Type>
inline bool ListeD<Type>::debuto()
{
    return( (courant==NULL));
}
template <class Type>
inline bool ListeD<Type>::est_vide()
{
  if (premier == NULL) 
    return true;
  else
    return false;
}

template <class Type>
inline void ListeD<Type>::debut()
{ courant=premier; }

template <class Type>
inline void ListeD<Type>::fin()
{ while(!est_fin())
    courant=courant->next();  
}

template <class Type>
inline void ListeD<Type>::suivant()
{
  if(!finito())
    courant=courant->next();
}
template <class Type>
inline void ListeD<Type>::precedent()
{
  if(!debuto())
    courant=courant->back();
}

template <class Type>
inline Type& ListeD<Type>::contenu()
{ return (courant->donne()); }

template <class Type>
inline void ListeD<Type>::detruire(){
  if (!est_vide())
    {
      if (courant==premier){
	courant=premier->next();
	delete premier;
	premier=courant;
	nbe--;
      }
      else {
	NoeudD <Type>* prec=courant->back();
	prec->next(courant->next());
	(courant->next())->back(prec);
	delete courant;
	courant=prec;
	nbe--;
      }
    }
  else
    Ferr << "ListeD[detruire] Liste deja vide!\n";
}// ListeD::detruire()

template <class Type>
inline void ListeD<Type>::detruire_droite()
{
  NoeudD <Type>* tmp;
  if (!est_vide())
    { if (!est_fin())
	{
	  tmp=courant->next();
	  courant->next(tmp->next());
	  delete courant->next();
	  nbe--;
	}
      else
	  Ferr <<" ListeD[detruire_droite]  Impossible de detruire un element car fin de liste"<<'\n';
    }
  else
      Ferr << "ListeD[detruire_droite] Impossible de detruire un element car liste vide"<<'\n';
}//ListeD::detruire_droite()

template <class Type>
inline void ListeD<Type>::free_liste(){
  if (!est_vide())
    for (debut(); !est_vide(); detruire());
  //else
  //Ferr << "ListeD[detruire] Liste deja vide!\n";
}// ListeD::free_liste()

template <class Type>
inline void ListeD<Type>::invert( NoeudD<Type> *A, NoeudD<Type> *B){
  Type tmp;
  tmp=A->donne();
  A->donne(B->donne());
  B->donne(tmp);
}//ListeD::invert

template <class Type>
void ListeD<Type>::qsort( NoeudD<Type> *g, NoeudD<Type> *d,int (*comp)(Type x, Type y)){
  Type val;
  NoeudD<Type> *i,*j,*d0;

  if(g!=d){
    val=d->donne();
    i=g;
    j=d0=d;
    do {
      while(comp(i->donne(),val)<0 && i!=d){
	i=i->next();
      }//while

      do{
	if(i==j) break;
        j=j->back();
      } while(  comp(j->donne(),val) >0 );
      invert(i,j);
    }while(i!=j);
    
    invert(i,j);
    invert(i,d);
    j=(i!=d0)?i->next():i;
    i=(i!=g)?i->back():i;
    qsort(g,i,comp);
    qsort(j,d,comp);
  }//if
}//ListeD::qsort()
template <class Type>
void ListeD<Type>::QuickSort(int (*comp)(Type x, Type y)){
  NoeudD<Type> *g,*d;

  g=premier;
  for(debut();!est_fin();suivant());
  d=courant;
  qsort(g,d,comp);
}// ListeD::QuikSort()

template <class Type>
inline void detruire_contenu(ListeD<Type *>& liste)
{
	Type *donnee;
	liste.debut();
	while (!liste.est_fin())
	{
		donnee=liste.contenu();
		delete donnee;
		liste.suivant();
	}
}

#endif
