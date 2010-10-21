
#define _Msail

#include <iostream> // définit des namespaces
using namespace std ;

#include "multicou.h"
#include <cmath>
/*class Mlayin{// a allouer pour n couches
public:
  REEL att, sig, ks, sb, sf, ko, uf, ub, w;
   static int i;
  ~Mlayin(){printf("~Mlayin(%d)\n",++i);}
};
static int Mlayin::i=0;
*/
struct Mlayin{// a allouer pour n couches
  REEL att, sig, ks, sb, sf, ko, uf, ub, w;
};

static Mlayin* Tlayin=NULL, layin,*playin;
// prototype
void distrib(double *,double,double);

/**************************************************************************
                         SOUS-PROGRAMME   MSAILAD                           
 
   Calcule les coefficients du systeme de SUITS en fonction des parametres  
   du couvert vegetal (LAI et fonction de distribution des inclinaisons     
   foliaires), des geometries de l'incidence et de la visee                 
   dans le cas d'un couvert multicouches. (Adapte d'apres Verhoef 1984)
   
   les valeurs de la fonction de distribution sont calculees par le         
   sous-programme DISTRIB.                                                  
*************************************************************************/
void msailad(Msailin &msailin){
  register int ia,ic,ili;
  REEL fcum, sbfs, cstl, psir, tsin, sntl, cs2tl, sn2tl;
  REEL cspsi, tanto, tants, t1, btran1, btran2;
  REEL bt1, bt2, bt3, sw1, sw2, sbf;
  REEL bto, sko, bts, sks, rtm, ttl, rtp;
  double u,v,*fbeta;//loi beta
 
  /*
    l   : indice foliaire de la couche                                    
    roo : reflectance des feuilles (lambertiennes) [en pourcents]         
    tau : transmittance des feuilles (lambertiennes) [en pourcents]       
    f   : frequence d'inclinaison foliaire [en fraction]                  
    
    tts : angle zenithal solaire (equivalent hauteur)  [en degres]        
    tto : angle zenithal de visee [en degres]                             
    psi : difference azimutale [en degres]                                
    
    bts : angle de transition beta (solar radiance)                       
    bto : angle de transition beta (canopy radiance)                      
    
    bt1, bt2, bt3 : angles azimutaux auxiliaires                          
    */
  printf("msailad() : Debut, N=%d, nbang=%d\n",N,msailin.nbang);fflush(stdout);
  Tlayin=new Mlayin[N];
  if (msailin.nbang == 13)
    fbeta=new double[13];//Verhoef et Beta(mu,nu)
 
  for (ic = 1; ic < N; ic++) {// boucle globale sur les couches
    playin=&(Tlayin[ic]);
    //  Calculs preliminaires
    rtp = (msailin.roo[ic] + msailin.tau[ic]) * .005;
    if(rtp==0){
      //cas ou s2v met roo=too=0 qund LAI=0
      rtp=rtm=.1;
    }
    else
      rtm = (msailin.roo[ic] - msailin.tau[ic]) * .005;
    
      tants = tan(rd*msailin.tts);
    // printf("%lf, %lf\n",rd,msailin.tto);fflush(stdout);
    tanto = tan(rd*msailin.tto);
    cspsi = cos(rd *msailin.psi);
    psir  = rd*msailin.psi;
    // Appel de la fonction de distribution des angles foliaire, fonction Beta
    if (msailin.nbang == 13) {
      u = msailin.bmu[ic];
      v = msailin.bnu[ic];
      distrib(fbeta,u,v);
      if (ic==1)
	printf("==> Distribution Beta d angle folaire\n");
    }
    // Initialisations
    sbf=sks=sko=sw1=sw2= 0.;
    // Boucle sur les composantes des coefficients
    fcum = 0;
    for (ia = 0; ia<msailin.nbang; ia++) {
      if (msailin.nbang == 13) {
	//13 classe d'angles de Verhoef avec distribution Beta
	msailin.f(ic,ia) = fbeta[ia];
	if (ia < 8)
	  ili = ia*10 + 5;
	else 
	  ili = (ia-8)*2 + 81;
	ttl = (REEL) ili;      
      }
      else {
	// Cas de classes d'angles foliaires mesure'es par s4 

	//msailin.f(ic,ia) *= 100.;//deplacer dans le main a cuase des boucle de direction!
	ttl = 90./(REEL)msailin.nbang *((REEL)ia+.5);
      }
      fcum += msailin.f(ic,ia);
      //printf("tm=%4.1lf - f(%d)=%6.3lf - fcum=%7.3lf\n",ttl,ia,msailin.f(ic,ia),fcum);
      cstl = cos(rd*ttl);
      cs2tl=cstl*cstl;
      sn2tl=1.-cs2tl;
      sntl=sqrt(sn2tl);

      sbf=sbf+cs2tl*msailin.f(ic,ia);

      /*
	Calcul de bts et de bto
	bts et bto sont les angles de transition beta s et beta o pour les
	directions d'incidence et de visee.       
	bts peut etre calcule lorsque ttl+tts > pi/2
	alors si bts < azimut feuille < 2pi - bts , les feuilles sont
	eclairees par leur face inferieure.
	Pour tts+ttl < pi/2 les feuilles sont toujours eclairees par leur
	face superieure; bts = pi pour assurer la continuite des formules
	Rq: bto `a une signification equivalente dans la direction de visee
	*/
	    
      bts = pi;
      if (ttl + msailin.tts > 90) 
	bts = acos(-cstl/(sntl*tants));
      bto = pi;
      if (ttl+msailin.tto > 90){
	printf("Bug\n");
	bto = acos(-cstl/(sntl*tanto));
      }
      // Definition des angles azimutaux auxiliaire bt1, bt2, bt3 utilises
      // dans le calcul du coefficient de diffusion bidirectionnelle w
      btran1 = fabs(bts-bto);
      btran2 = pi*2 - bts - bto;
      if (psir <= btran1) {
	bt1 = psir;
	bt2 = btran1;
	bt3 = btran2;
      }else{
	bt1=btran1;
	if (psir <= btran2) {
	  bt2 = psir;
	  bt3 = btran2;
	}else{
	  bt2 = btran2;
	  bt3 = psir;
	}
      }
      // Calcul des composantes des coefficients
      sks += ((bts - pi*.5)*cstl + sin(bts)*tants*sntl)*msailin.f(ic,ia);
      sko += ((bto - pi*.5)*cstl + sin(bto)*tanto*sntl)*msailin.f(ic,ia);
      tsin = sn2tl*.5 * tants*tanto;
      t1 = cs2tl + tsin*cspsi;
      sw1 += (cs2tl + tsin*cspsi) * msailin.f(ic,ia);
      if (bt2 != 0.) 
	sw2 += (-bt2*t1+ sin(bt2)*(cs2tl/(cos(bts)*cos(bto))
				   + cos(bt1)*cos(bt3)*tsin))*msailin.f(ic,ia);
    }//for ia
    /**********************************************************************
     *                  COEFFICIENTS DU SYSTEME DE SUITS                  *
     *                                                                    *
     *    ATT  coefficient d'attennuation du flux diffus                  *
     *    SIG  coefficient de retrodiffusion du flux diffus               *
     *    KS   coefficient d'extinction du flux direct incident           *
     *    KO   coefficient d'extinction dans la direction de visee        *
     *    SF   coefficient de diffusion du flux direct incident en        *
     *         flux diffus  descendant                                    *
     *    SB   coefficient de diffusion du flux direct incident en        *
     *         flux diffus montant (retrodiffusion)                       *
     *    UF   coefficient de diffusion dans la direction de visee        *
     *         du flux diffus montant                                     *
     *    UB   coefficient de diffusion dans la direction de visee        *
     *         du flux diffus descendant (retrodiffusion)                 *
     *    W    coefficient de diffusion bidirectionnelle                  *
     *                                                                    *
     **********************************************************************/
    sbfs = sbf*rtm * msailin.l[ic]*.01;
    playin->att = msailin.l[ic]*(1-rtp) + sbfs;
    playin->sig = msailin.l[ic]*rtp + sbfs;
     printf("sig[%d]=%g =>rtp=%g, sbfs=%g, sbf=%g, rtm=%g \n",ic,playin->sig,rtp,sbfs,sbf,rtm);
    playin->ks = sks*.02 * msailin.l[ic]/pi;
    printf("<!> sks=%.3g\t msailin.l[%d]=%.3g\t => playin->ks = %.3g\n",sks,ic, msailin.l[ic], playin->ks);
    playin->sb = playin->ks*rtp + sbfs;
    playin->sf = playin->ks*rtp - sbfs;
    playin->ko = sko*.02 * msailin.l[ic]/pi;
    playin->uf = playin->ko*rtp + sbfs;
    playin->ub = playin->ko*rtp - sbfs;
    playin->w=msailin.l[ic]*.01*(msailin.roo[ic]*.01*sw1 + 2*rtp*sw2/pi);
  }//for ic couche
  // mis en comment pour plantage sur sun 11/12/97
  //delete [] fbeta;
  printf("msailad() : Fin\n");fflush(stdout);
}//msailad()

/********************************************************************************
 *                    SOUS-PROGRAMME  LAYER                                     *
 *                                                                              *
 *   Resolution du systeme d'equations differentielles decrivant les            *
 *   les transferts radiatifs dans un couvert de SAIL                           *
 *   VERHOEF 1985                                                               *
 *                                                                              *
 *******************************************************************************/

void mlayer(Msailin &msailin, Mlayout* Tlayout){
  static REEL m, e1, e2, h1, h2;
  static int ic;
  static REEL co, do_, cs, ds, ho, dnd, dno, dns, som;
  Mlayout *layout;
  if(Tlayin==NULL){
    fprintf(stderr,"<!> mlayer() appelee avant msailad() : Tlayin pas alloue\n");
    exit(-1);
  }
  layout=Tlayout;
  //printf("mlayer(): N=%d\n",N);
  for (ic=0; ic <N; layout++,ic++) {
    if (msailin.l[ic] == 0) {
      layout->tss = 1.;
      layout->too = 1.;
      layout->rdd = 0.;
      layout->tdd = 1.;
      layout->rsd = 0.;
      layout->tsd = 0.;
      layout->rdo = 0.;
      layout->tdo = 0.;
      layout->rso = 0.;
    }
    else{
      layin=Tlayin[ic];
      som=layin.att*layin.att - layin.sig*layin.sig;
      if (som<0.0)
	som=0.0;
      m=sqrt(som);
       printf("lai[%d]=%g => layin.att=%g, layin.sig=%g\n",ic,msailin.l[ic],layin.att,layin.sig);
       h1=(layin.att+m)/layin.sig;
      printf("Ok\n");
      h2=1./h1;
      e1=exp(m);
      e2=1./e1;
      dns = layin.ks*layin.ks - m*m;
      dno = layin.ko*layin.ko - m*m;
      cs=(layin.sb*(layin.ks-layin.att)-layin.sf*layin.sig)/dns;
      ds=(-layin.sf*(layin.ks+layin.att)-layin.sb*layin.sig)/dns;
      co=(layin.uf*(layin.ko-layin.att)-layin.ub*layin.sig)/dno;
      do_=(-layin.ub*(layin.ko+layin.att)-layin.uf*layin.sig)/dno;
      // Coefficients permettant le calcul des differents flux sortants
      // `a partir des flux entrants                                           :
      layout->tss=exp(-layin.ks);
      layout->too=exp(-layin.ko);
	if (m != 0.) {
	dnd = h1 * e1 - h2 * e2;
	layout->rdd = (e1 - e2) / dnd;
	layout->tdd = (h1 - h2) / dnd;
	}else{
	  layout->rdd = layin.sig/(1+layin.sig);
	  layout->tdd = 1-layout->rdd;
	}//else
	layout->rsd= cs*(1-layout->tss*layout->tdd) - ds*layout->rdd;
	layout->tsd= ds*(layout->tss-layout->tdd) - cs*layout->tss*layout->rdd;
	layout->rdo= co*(1 - layout->too*layout->tdd) - do_*layout->rdd;
	layout->tdo= do_*(layout->too-layout->tdd) - co*layout->too*layout->rdd;

	ho= (layin.sf*co + layin.sb*do_ + layin.w) / (layin.ko+layin.ks);
	layout->rso = ho * (1 - layout->tss*layout->too)
	  - co*layout->tsd*layout->too - do_*layout->rsd;
    }//else
  }// for ic
  delete [] Tlayin;
}//mlayer()

/****************************************************************************
                          SOUS-PROGRAMME  DISTRIB                            

      Calcul de la valeur de la fonction de distribution des inclinaisons      
      foliaires pour les classes 5,15,25,35,45,55,65,75,81,83,85,87,89         
      Fonction de distribution Beta - GOEL et STREBEL 1984
      Les parametres Mu et Nu ont ete determines par les formules presentees
      par GOEL et STREBEL (1984)  (programme Euler.for)                    
      Les calculs sont effectues par l'approximation de STIRLING           
   *************************************************************************/
void distrib(double *fbeta,double u,double v){
  static double freq, gamm1, gamm2, gamm3, part1, part2, part3;
  register  int ili,i;
  static double w,y1, y2, y3, dif, gam, ang, pie;
  static double pon, sum;

  pie = pi * 2;
  // Calcul des fonctions Gamma (approximation de STIRLING)
  y1=1/(12*u)-1/(360*u*u*u)-u;
  y2=1/(12*v)-1/(360*v*v*v)-v;
  w=u+v;
  y3=1/(12*w)-1/(360*w*w*w)-w;
  y1=exp(y1);
  y2=exp(y2);
  y3=exp(y3);
  gamm1=pow((pie/u),.5) * pow(u,u) * y1;
  gamm2=pow((pie/v),.5) * pow(v,v) * y2;
  gamm3=pow((pie/w),.5) * pow(w,w) * y3;
  gam=gamm3/(gamm1*gamm2);
  //Calcul des valeurs de la fonction de distribution fbeta(ili) 
  sum = 0;
  for (i = 0; i<13; i++) {
    if (i < 8) {
      ili = i * 10 + 5;
      pon = 9;
    }
    else {
      //ili = (i - 9 << 1) + 81;
      ili=(i-8)*2+81;
      pon = 45;
    }
    ang = (double) ili;
    part1 = gam / pon;
    part2=pow((1-ang/90),(u-1));
    part3=pow((ang/90),(v-1));
    freq = part1 * part2 * part3;
    sum += freq;
    fbeta[i] = freq*100;
  }
  dif = 1.-sum;
  for (i = 0; i<13; i++) {
    fbeta[i] += dif*fbeta[i];
  }
}//distrib()

#undef _Msail
