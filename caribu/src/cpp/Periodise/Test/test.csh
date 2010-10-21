#!/bin/csh
# Test de periodise2.0
#  lancer test.csh
#  regarder le fichier periodise' motif.can
# MC98
          
echo \ o--------------------------------oOo--------------------------------o\ 0
periodise  -8 test.8
echo \ o--------------------------------oOo--------------------------------o\ 1
periodise -m rien.can  -8 test.8
echo \ o--------------------------------oOo--------------------------------o\ 2
periodise -m test.can  -8 test.8
echo \ o--------------------------------oOo--------------------------------o\ 3
periodise -m test.can -p rien.nir  -8 test.8
echo \ o--------------------------------oOo--------------------------------o\ 4

