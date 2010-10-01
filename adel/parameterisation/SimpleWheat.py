import csv as csv

class dim_pattern:
	'''
	Class which implements patterns for predicting dimension, Aera_contrib and so on from data.
	'''
	def __init__(self, dimTn):
		'''
		Initialisation of dimension pattern, just need the equivalent of a dimT
		'''
		#Opening dimT file
		csvfile = open(dimT, "rb")
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		reader = csv.reader(csvfile, dialect)
		#Initialisation of the usefull information
		dict {'surface' : [],'bladelength': [], 'sheathlength': [], 'internodelength': [], 'stemdiameter': [] }
		self.emptydict= dict
		self.dictMS
		#dimT construction CAREFULL NO HEADER
		for row in reader:
			self.dictMS['surface'].append(float(row[]))
			self.dictMS['bladelength'].append(float(row[]))
			self.dictMS['sheathlength'].append(float(row[]))
			self.dictMS['internodelength'].append(float(row[]))
			self.dictMS['stemdiameter'].append(float(row[]))
		self.decphy = 1 # decalage des phytomeres
		
	
	def predict_tiller(self, tiller_position):
		'''
		Returns datadict for tiller_position:
		'''
		if tiller_position == 0:
			return self.dictMS
		else:
			dict = self.emptydict
			####
			#Mettre une loi interpolants les valeurs de dict pour le tiller en question
			####
			
			return dict
		
	
	def total_area(dict):
		'''
		Returns total aera of a tiller which is represented by a dict
		'''
		return(sum(dict['surface']))
		
	def PlantAera_contrib(self, nb_tiller):
		surface = [self.total_area()]
		for i in arange(nb_tiller):
			surface.append(self.total_area(self.predict_tiller(i)))
		Stot = sum(array(surface))
		return array(surface)/Stot
		#out = list_poid_0 a nb tiller


soisson = dim_patern(Datasoisson)


def Wheatpop(LAI=3, axe_density = 1000, plant_density= 300, nout = 10, nb_phy=11, 
	dim_pattern = 'soisson' # Qui contient PlantArea_contrib (% de plant area porté par les axes)
	):
	"""
	Generate a population of nout different plants with appropriate tiller number and surface per plant
	To be plug to "agronomic plot" and "SimpleWheat" 
	"""
	
	
	# Tableau en sortie: 
	#nout => nb de ligne du tableau
	# pour chaque ligne est un axe (MB ou talle)
	# boucle pour créer les plantes: loi qui arrondi au supérieur ou inferieur du nombre de talle moyen
	#pour faire les talles: BM: plante i, axe 0
	# out tableau: plante, axe, nbphy mais à terme ce sera un mtg.
	#Pour chaque plante: un numero, un nb d'axe, un nb de phytomere et une surface

#A terme prendra un mtg et cherchera la topo des plantes, + la surface de la plante et nb_thalles
def simpleWheat_param(total_area = 10000 , # aussi dans wheatpop
		total_height = 200,  #Nouveau!!!
		nb_phy = 11, # entrée de Wheatpop
		# leaf area dist function
		dim_pattern = 'soisson', #function des distributions de longueurs d'organes et du passage brin-maitre tiller(gaine entre-noeud et limbe) (3functions)
		geom_pattern = 'soisson', # BD des XYSR de nervures mesurées
		#dred: maximal distance between maintsteem and tiller at flowering
		dred = 7, #Nouveau
		diam_base = 3,#Nouveau
		diam_top = 1,#Nouveau
		nbtalles = 3, #topo: info dans le mtg (mtg.verticies_at_scale(1))
		#No azimuth param
		):
	""" 
	generate simpleMais dictionary from parameter for wheat-like plants
	
	"""

	#creation brin maitre

	nb_phy = int(nb_phy)
	nb_young_phy = int(nb_young_phy)

	# compute the leaf surface
	leaf_area = np.array(bell_shaped_dist(total_area, nb_phy, lad_rmax, lad_skew))
	# derive corresponding length and widths
	lengths = np.sqrt(leaf_area / np.array(shape_factors) / leaf_width_ratio)
	widths = lengths * leaf_width_ratio

	# distances between leaves
	pseudostem = geometric_dist(pseudostem_height, nb_young_phy,pseudostem_dist)
	stem = geometric_dist(total_height - pseudostem_height, nb_phy - nb_young_phy,stem_dist)
	distances = pseudostem + stem
    # sheath and internode true length
    pseudostem_sheath = (np.array(pseudostem)).cumsum()
    pseudostem_internode = [0] * nb_young_phy
    stem_sheath = np.array(stem)
    stem_internode = [0] + (stem_sheath[1:]).tolist() 
    sheaths = pseudostem_sheath.tolist() + stem_sheath.tolist()
    internodes = pseudostem_internode + stem_internode
    # stem diameters
    diameters = [diam_base] * nb_young_phy + np.linspace(diam_base, diam_top, nb_phy - nb_young_phy).tolist()
    # inclinations
    inclinations = np.linspace(basal_insertion, top_insertion, nb_phy)
    # azimuths Attention R
      Azimuth = 
	(n-4 n) 180 + 60 * randn(1),
              (1 n-4)   180 + 20 * randn(1))
         
    axe = [0] * nbphy
    Einc = axe
    #creation des axes

    #Alexis: initialisations des talles
    t_lengths   = talle_expend_length(lengths, nbtalles, nb_phy)
    t_widths    = talle_expend_length(widths, nbtalles, nb_phy)
    t_distances = talle_expend_length(distances, nbtalles, nb_phy)
    t_sheaths   = talle_expend_length(sheaths, nbtalles, nb_phy)
    t_internodes= talle_expend_length(internodes, nbtalles, nb_phy)

    for i in range(1,nbtalles+1):
        t_nb_phy = round(nb_phy - (2.5 + 0.5 * i))
        t_distances = distances[(nb_phy - t_nb_phy -1):]
        distances = np.concatenate((distances,t_distances))
        
        t_diameters = np.linspace(diam_base, diam_top, t_nb_phy)

        inclinations = np.linspace(basal_insertion, top_insertion, t_nb_phy)
        #azimuth : aller chercher regle adel ?
        axe = axe + [i] * t_nb_phy
        #
        t_Einc = [60] + [0] * (t_nb_phy -1)
    # output dict
    dTags = ["plante","axe","phytomere","longFeu","largFeu","inclinaisonFeu","azimuthFeu","longTige","diamTige","longGa","longEn","incEn"] 
    dVals = [ [1] * nb_phy * (nbtalles + 1),
              axe ,
              range(1,nb_phy + 1),
              lengths,
              widths,
              inclinations,
              azimuths,
              distances,
              diameters,
              sheaths,
              internodes,
              Einc
              ]
    dVals = map(np.array,dVals)
    dout = dict(zip(dTags,dVals))

    return dout

#Alexis: propage l'info d'une talle a la suivante 
def talle_expend_length(param,nbtalles,nb_phy):
    exp = param # Initialisation, param appartient uniquement au MB
    for i in range(1,nbtalles+1):
        t_nb_phy = round(nb_phy - (2.5 +0.5 *i))
        t_param = param[(nb_phy - t_nb_phy - 1):]
        exp = np.concatenate((exp, t_param))
    return exp
    