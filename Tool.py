#!/usr/bin/env python3
import math,re

"""
Authors: AGLAVE Marine & ROHMER Coralie
Date: 2018-05
Description:Ceci est un fichier contenant toutes les fonctions du 
projet Python en Dynamique Moléculaire: Etude du Complexe sRNA H/ACA
"""



#============================ Parsers ==================================
def parserConf(infile) :
	"""
    Parse un fichier en format pdb de plusieurs configurations d'une
    protéine en un dictionnaire utilisable par Confo_global.py et Confo_local.py
    :param fichier pdb du type:
    "MODEL
     ATOM      numéro_atome  type_atome type_aa    numéro_aa     x y  z  1.00  température      domaine"  
    :return: un dictionnaire contenant les informations du fichier pdb.
    """
	DicoConf={}
	DicoConf["order"]=[]
	infile=open(infile,"r")
	for i in infile:
		if i[0:5]in "MODEL ":
			i.replace("\n","")
			tab = re.split('( )+',i)
			model=int(tab[2])
			DicoConf[model]={}
			DicoConf["order"] += [model]
			DicoConf[model]["order"]=[]
			
		if i[0:5]in "ATOM  ":
			domaine = i[72:74].replace(" ","")
			if domaine not in DicoConf[model]:
				DicoConf[model][domaine]={}
				DicoConf[model][domaine]["order"]=[]
				DicoConf[model]["order"].append(domaine)
				
			residu = i[23:29].replace(" ","")
			if residu not in DicoConf[model][domaine]:
				DicoConf[model][domaine][residu]={}
				DicoConf[model][domaine][residu]["order"] = []
				DicoConf[model][domaine]["order"].append(residu)
				
			DicoConf[model][domaine][residu]["NomRes"]=i[17:20]
			NomAtom= i[12:16].replace(" ","")
			if NomAtom not in DicoConf[model][domaine][residu]:
				DicoConf[model][domaine][residu][NomAtom]={}
				DicoConf[model][domaine][residu]["order"] += [NomAtom]
			DicoConf[model][domaine][residu][NomAtom]["NumAtom"]=i[6:12].replace(" ","")
			DicoConf[model][domaine][residu][NomAtom]["x"]=float(i[30:38].replace(" ",""))
			DicoConf[model][domaine][residu][NomAtom]["y"]=float(i[38:46].replace(" ",""))
			DicoConf[model][domaine][residu][NomAtom]["z"]=float(i[46:54].replace(" ",""))
	return DicoConf
	
	
def parser(infile) :
	"""
    Parse un fichier en format pdb contennat une configuration de
    référence d'une protéine  en un dictionnaire utilisable par Confo_global.py
    :param fichier pdb du type:
    "ATOM      numéro_atome  type_atome type_aa    numéro_aa     x y  z  1.00  température      domaine"	  
    :return: un dictionnaire contenant les informations du fichier pdb.
    """
	DicoChaine={}
	DicoChaine["order"]=[]
	infile=open(infile,"r")
	for i in infile:
		if i[0:5]in "ATOM  ":
			domaine = i[72:74].replace(" ","")
			if domaine not in DicoChaine:
				DicoChaine[domaine]={}
				DicoChaine[domaine]["order"]=[]
				DicoChaine["order"].append(domaine)
			residu=i[23:29].replace(" ","")
			if residu not in DicoChaine[domaine]:
				DicoChaine[domaine][residu]={}
				DicoChaine[domaine][residu]["order"] = []
				DicoChaine[domaine]["order"].append(residu)
				
			DicoChaine[domaine][residu]["NomRes"]=i[17:20]
			NomAtom= i[12:16].replace(" ","")
			if NomAtom not in DicoChaine[domaine][residu]:
				DicoChaine[domaine][residu][NomAtom]={}
				DicoChaine[domaine][residu]["order"] += [NomAtom]
			DicoChaine[domaine][residu][NomAtom]["NumAtom"]=i[6:12].replace(" ","")
			DicoChaine[domaine][residu][NomAtom]["x"]=float(i[31:38].replace(" ",""))
			DicoChaine[domaine][residu][NomAtom]["y"]=float(i[39:46].replace(" ",""))
			DicoChaine[domaine][residu][NomAtom]["z"]=float(i[47:54].replace(" ",""))
	return DicoChaine

def parserRMSD(infile) :
	"""
    Parse le fichier de sortie du script Confo_global.py
    en un dictionnaire utilisable par graphique_RMSD.py.
    :param fichier pdb du type:
    "Model RMSD[domaine1] RMSD[domaine2] RMSD[domaine3] ...  RMSD_Global"  
    :return: un dictionnaire contenant les informations du fichier.
    """
	dicoVariable = {}
	dicoVariable['order'] = []
	infile=open(infile,"r")
	for i in infile:
		i=i.replace("\n","")
		if i[0] not in "#":
			ligne = i.split("\t")
			for j in range(len(ligne)) :
				dicoVariable[dicoVariable['order'][j]].append(float(ligne[j].replace(" ","")))
		elif i[0:6] in "#Model" :
			ligne = i.split(" ")
			ligne[0] =ligne[0].replace("#","")
			for cle in ligne :
				dicoVariable['order'].append(cle)
				dicoVariable[cle]=[]
	infile.close
	dicoVariable['order'].remove("Model")
	return dicoVariable

#========================= Etude Globale ===============================
def distance(x1,y1,z1,x2,y2,z2):
	"""
    Calcule la distance entre 2 points.
    :param x1
    :param y1
    :param z1
    :param x2
    :param y2
    :param z2
    :return: distance
    """
	return float(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
	
	
def barycentre(dREF,dCONF,model,domaine,residu) :
	"""
    Calcule les coordonnées des barycentres de 2 résidus et calcule la 
    distance entre ces barycentres.
    :param dREF: dictionnaire de la protéine de référence obtenu avec 
    la fonction parser().
    :param dCONF: dictionnaire des conformation de la protéine obtenu 
    avec la fonction parserConf().
    :param model: le model de la conformation d'intrérêt.
    :param domaine: le domaine protéique étudié.
    :param residu: le résidu d'intérêt.
    :return: distance
    """
	# initialise les sommes et le nb atom
	xSomConf=0
	ySomConf=0
	zSomConf=0
	xSomRef=0
	ySomRef=0
	zSomRef=0
	nbAtom=0
	for nomAtom in dCONF[model][domaine][residu]["order"] :
		xSomConf += dCONF[model][domaine][residu][nomAtom]["x"]
		ySomConf += dCONF[model][domaine][residu][nomAtom]["y"]
		zSomConf += dCONF[model][domaine][residu][nomAtom]["z"]
		xSomRef += dREF[domaine][residu][nomAtom]["x"]
		ySomRef += dREF[domaine][residu][nomAtom]["y"]
		zSomRef += dREF[domaine][residu][nomAtom]["z"]
		nbAtom += 1
	# renvoie le calcul de la distance pour le barycentre
	return distance(xSomConf/nbAtom,ySomConf/nbAtom,zSomConf/nbAtom,xSomRef/nbAtom,ySomRef/nbAtom,zSomRef/nbAtom)

#========================== Etude Locale ===============================
def bary_local(dCONF):
	"""
    Conçoit un dictionnaire simplifié des coordonnées des barycentres des résidus
    d[model][domaine][résidus]=[x,y,z].
    :param dCONF: dictionnaire des conformation obtenu avec la fonction parserConf()
    :return: dictionnaire
    """
	dCoordonnes={}
	dCoordonnes["order"]=[]
	for model in dCONF["order"]:
		dCoordonnes["order"].append(model)
		dCoordonnes[model]={}
		dCoordonnes[model]["order"]=[]
		for domaine in dCONF[model]["order"]:
			dCoordonnes[model]["order"].append(domaine)
			dCoordonnes[model][domaine]={}
			dCoordonnes[model][domaine]["order"]=[]
			for residu in dCONF[model][domaine]["order"]:
				sumX=0
				sumY=0
				sumZ=0
				for atome in dCONF[model][domaine][residu]["order"]:
					sumX += dCONF[model][domaine][residu][atome]['x']
					sumY += dCONF[model][domaine][residu][atome]['y']
					sumZ += dCONF[model][domaine][residu][atome]['z']
				nb_atomes=len(dCONF[model][domaine][residu]["order"])
				dCoordonnes[model][domaine][float(residu)]=[sumX/nb_atomes, sumY/nb_atomes, sumZ/nb_atomes]
				dCoordonnes[model][domaine]["order"].append(float(residu))
	return dCoordonnes
	
def atom_local(dCONF):
	"""
    Conçoit un dictionnaire simplifié des coordonnées atomiques des Carbones alpha de chaque résidus
    d[model][domaine][résidus]=[x,y,z] si le domaine est proteique sinon renvoie les coordonnées du barycentre
    :param dCONF: dictionnaire des conformation obtenu avec la fonction parserConf()
    :return: dictionnaire
    """
   
	dCoordonnes={}
	dCoordonnes["order"]=[]
	for model in dCONF["order"]:
		dCoordonnes["order"].append(model)
		dCoordonnes[model]={}
		dCoordonnes[model]["order"]=[] 
		for domaine in dCONF[model]["order"]:
			dCoordonnes[model]["order"].append(domaine)
			dCoordonnes[model][domaine]={}
			dCoordonnes[model][domaine]["order"]=[]
			premierResidu = dCONF[model][domaine]["order"][0] 
			# si le domaine est proteique
			if "CA" in dCONF[model][domaine][premierResidu]:
				for residus in dCONF[model][domaine]["order"]:
					dCoordonnes[model][domaine][float(residus)]=[dCONF[model][domaine][residus]["CA"]['x'],dCONF[model][domaine][residus]["CA"]['y'],dCONF[model][domaine][residus]["CA"]['z']]
					dCoordonnes[model][domaine]["order"].append(float(residus))
			# si le domaine n'est pas proteique
			else :
				for residu in dCONF[model][domaine]["order"]:
					sumX=0
					sumY=0
					sumZ=0
					for atome in dCONF[model][domaine][residu]["order"]:
						sumX += dCONF[model][domaine][residu][atome]['x']
						sumY += dCONF[model][domaine][residu][atome]['y']
						sumZ += dCONF[model][domaine][residu][atome]['z']
					nb_atomes=len(dCONF[model][domaine][residu]["order"])
					dCoordonnes[model][domaine][float(residu)]=[sumX/nb_atomes, sumY/nb_atomes, sumZ/nb_atomes]
					dCoordonnes[model][domaine]["order"].append(float(residu))
	return dCoordonnes
	
def distance_locale(dico,model,dom1,res1,dom2,res2):
	"""
    Calcule la distance entre 2 résidus.
    :param dico: dictionnaire contenant les références des coordonnées des résidus
    :param model: model de la conformation protéique
    :param dom1: domaine protéique du premier résidu
    :param res1: premier résidu
    :param dom2: domaine protéique du deuxième résidu
    :param res2: deuxième résidu
    :return: distance
    """
	x1=dico[model][dom1][res1][0]
	y1=dico[model][dom1][res1][1]
	z1=dico[model][dom1][res1][2]
	x2=dico[model][dom2][res2][0]
	y2=dico[model][dom2][res2][1]
	z2=dico[model][dom2][res2][2]
	return float(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))

def freq_interface(dCoordonnees_totales, tab_groupe_domaine1, tab_groupe_domaine2, seuil_distance):
	"""
    Calcul les fréquences d'appartenance à l'interface de tout les 
    résidus des domaines de tab_groupe_domaine1. L'interface est définie
    par le seuil de distance entre tab_groupe_domaine1 et tab_groupe_domaine2.
    Le résultat est retourné sous la forme d'un dictionnaire.
    :param dCoordonnees_totales: dictionnaire contenant les références des coordonnées des résidus
    :param tab_groupe_domaine1: tableau contenance la liste des domaines1
    :param tab_groupe_domaine2: tableau contenance la liste des domaines2
    :param seuil_distance: seuil définissant l'interface
    :return: dInterface
    """
	dInterface={}                                  # création d'un dictionnaire qui contiendra les fréquences
	dInterface["order"]=[]

	for dom1 in tab_groupe_domaine1:
		dInterface[dom1]={}
		dInterface[dom1]["order"]=[]
		dInterface["order"].append(dom1)
		for model in dCoordonnees_totales["order"]:
			for res1 in dCoordonnees_totales[model][dom1]["order"]:
				if res1 not in dInterface[dom1]["order"]:
					dInterface[dom1][res1] = 0     # initialisation du nombre de fois où le résidu appartient à l'interface à 0
					dInterface[dom1]["order"].append(res1)

				dist_min = 10000000                # initialisation de la distance minimale residu_dom1-residu_dom2
				for dom2 in tab_groupe_domaine2:
					for res2 in dCoordonnees_totales[model][dom2]["order"]:
						d = distance_locale(dCoordonnees_totales,model,dom1,res1,dom2,res2)
						if d < dist_min:           # si la distance entre entre res1 et le res2 le plus proche est plus petite que la valeur minimale
							dist_min = d           # valeur minimale sera donc la distance entre res1 et le res2 le plus proche
				if dist_min <= seuil_distance:     # si res1 appartient à l'interface (suffisement proche d'un res2)
					dInterface[dom1][res1] += 1    # on incrémente de 1 le nombre de fois où le résidu appartient à l'interface
	return dInterface
