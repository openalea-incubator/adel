#

# Chaine les appel de s2v, mcsail et canestra, comme le fait cpfg via caribu
# Syntaxe: caribu.csh Ds file.light file.can nz h file.8 file1.opt ... fileN.opt
#
# Example /home/tld/chelle/QHS/Calgary/dev/Canestra/MC-SAIL/Test 
# caribu.csh 0.2 sky2.light Tout.can 3 18 Tout.8 test nir par
# MC98

# Menage
mv *.spec temp/
mv spectral temp/
mv cropchar temp/
mv leafarea temp/
mv *.env temp/
mv *.log temp/

#envt
setenv Sc /tmp/virtualis.can

#
if ( $#argv < 7 ) then 
 echo "Syntax error: caribu.csh Ds file.light file.can nz h file.8 file1.opt [... fileN.opt]"
else
 
# Periodise (caribu le fait il ?)
# 
  echo "periodise"  -m $argv[3] -8 $argv[6] -o $Sc
 periodise -m $argv[3] -8 $argv[6] -o $Sc
# S2V
 echo "s2v" $argv[3-]
 #s2v $argv[3-] >& s2v.log
 s2v $Sc $argv[4-] >& s2v.log
 if( $status != 0) then 
   echo "S2V a plante!"
   exit(-1)
 endif # Loop sur le p.o.
 set i=1
 foreach po ($argv[7-])
   echo "==> "$po, $i
   ## MCSAIL
   cp $po".spec" spectral
   echo "  mcsail" $argv[2]
   mcsail $argv[2]>& sail-$po.log
   if( $status != 0) then 
    echo "Sail a plante!"
    exit(-1)
   endif
   mv mlsail.env $po.env
   ## Canestra
   date
   if ( $i == 1 ) then
      ## 1ere fois: calcul des FF et Bfar
      echo "Appel no. "$i ": canestrad" -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r  $argv[1] -1   
       canestrad -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r $argv[1]  -t /tmp/ -f caribu -v 2 >& nr-$po.log  
  else
      ## ie fois: pas de calcul des FF et Bfar
      echo "Appel no. "$i ": canestrad" -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r  $argv[1] -1   
       canestrad -M $Sc -8 $argv[6] -l $argv[2] -p $po.opt -e $po.env -s -r $argv[1]    -t /tmp/ -w caribu  -v 2>& nr-$po.log 
  # nr-$po'_'$i.log
  endif
  date
  cp E0.dat E_$po.dat
  cp B.dat B_$po.dat
  echo " " 
 @ i = $i  + 1
 end
 # \rm leafarea spectral cropchar
endif
## FIN


### Trucs & Astuces
#
# $status = variable de retour 
#
# Exemple de calcul en variable Cshell
#setenv c `cat $1|wc -w`
#@ c = $c  / $r
