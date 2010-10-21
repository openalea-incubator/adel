######################################################################
##########     mcsail (SAIL multicouches)  ---  MC98     #############
######################################################################


mcsail est une version multicouche en C/C++ du SAIL de Verhoef(1984,1985), 
modifiee par A. Olioso (INRA Bioclimatologie, Avignon) et 
traduite en C et un peu adaptee par moi meme. 

Besoin: Un compilateur C++
======

Installation
============
* detarer l'archive, ce qui cree un repertoire MC-SAIL
 $ gunzip -c mcsail.tar.gz | tar xvf - 
* Aller dans le repertoire src/
* Mettre a jour le Makefile, notamment les variables ROOT, TMP et CC
* Compiler
 $ make

Caracteristiques
================
mcsail calcule les profils de flux Es, E+ et E-, pour un couvert en couches homogenes (indice d'agregation =1), de feuilles tres petites/ hauteur du couvert (pas de traitement du hot spot) et dont les proprietes optiques du sol et des feuilles sont lambertiennes. 

mcsail requiert 4  fichiers d'entree (leur format respecte une compatibilite avec le modele RiRi d'H. Sinoquet (INRA  Bioclimatologie, Clermont-Ferrand)):
* cropchar (Crop Characteristiques)
  Dans ce fichier ne nous interessent que:
  - ligne 1, le premier nombre, qui est le nombre de classes d'angle zenital NBANG, 
    utilisees    pour decrire l'orination foliaire d'une couche (LIDF)
    Rq: Cette LIDF est:
        soit decrite numeriquement via le fichier LEAFAREA,
        soit calculee par la fonction Beta a deux parametre mu et nu, ssi NBANG=13,
             qui correspond au nombre de couches, utilise par Verhoef
  - ligne 2, 
    nombre 1: nombre de couches de vegetation N
    nombre 2: hauteur d'une couche en metre

* leafarea: profil des caracteristiques des couches de vegetation 
  Une ligne decrit une couche, la premiere correspond au bas du couvert.
  Pour chaque ligne:
    nombre 1 et 2: indice de couche (non utilise)
    nombre 3 a` NBANG+3: Proportion de surface foliaire dans la classe d'angle zenital
    nombre  NBANG+4 et NBANG+5: mu et nu utilise quand NBANG=13 pour calculer la LIDF avec la loi Beta
    nombre NBANG+6: LAI de la couche 

* spectral: Proprietes optiques
  ligne 1: nombre de couches NbC ayant des proprietes optiques de feuille differentes
  ligne 2: reflectance du sol
  ligne 3: reflectance et transmittance des feuilles pour la couche pres du sol
  ligne i: reflectance et transmittance des feuilles pour la couche i

  Rq: Si NbC <> N, toutes les couches sont supposees avoir les memes proprietes, qui sont celles de la couche pres du sol (ligne 3). 

* incident: Caracteristiques de l'eclairage par un soleil collimate et un ciel UOC
  Nombre 1: Hauteur du soleil (en degre)
  Nombre 2: Azimuth solaire / rang (inutilise)  
  Nombre 3: Rapport diffus/direct


mcsail genere les profils de flux Es, E+ et E- sous differents formats (se referer a multicou.C)  

Des exemples se trouvent dans data/


Contact
======

Michael CHELLE
U.R. de Bioclimatologie
INRA
F-78850 THIVERVAL-GRIGNON
chelle@bcgn.grignon.inra.fr


References
==========
@article{Verhoef84,
   author = "W. Verhoef",
   journal = rse,
   volume = "16",
   pages = "125--141",
   title = "Light scattering by leaf layers with application to reflectanvce
		  canopy modeling: the {SAIL} model ",
   key = "canopy,SAIL",
   comment = "pionnier",
   year = "1984",
}	

@article{Verhoef85,
   author = "W. Verhoef",
   journal = rse,
   volume = "17",
   pages = "164--178",
   title = "Earth observation modeling based on layer scattering matrices",
   key = "canopy,SAIL",
   comment = "pionnier",
   year = "1985",
}	

@phdthesis{oliosophd,
  author = 	 "Albert Olioso",
  title = 	 "Simulation des échanges d'énergie et de masse d'un couvert
		  végétal, dans le but de relier la transpiration et la
		  photosynthèse aux mesures de réflectance et de température de
		  surface",
  type="Thèse $3^e$ Cycle",		  
  school = 	 "Université de Montpellier II",
  year = 	 "1992",
}

@phdthesis{chellephd,
  author = "Michaël Chelle",
  title  = "Développement d'un modèle de radiosité mixte pour simuler la
            distribution du rayonnement dans les couverts végétaux",
  type   = "Thèse $3^e$ Cycle",	
  school = "Université de Rennes {I} (Informatique)",
  year   = "1997",
  note  = "160 p. --
                  http://www.irisa.fr/EXT\-ERNE/bi\-bli/the\-ses/the\-ses97.html"
}