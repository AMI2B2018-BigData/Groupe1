#!/usr/bin/env python3
import Tool,sys,math,numpy


"""
		Projet Python en Dynamique Moléculaire:
			Etude du Complexe sRNA H/ACA

Authors: AGLAVE Marine & ROHMER Coralie
Date: 2018-05

Description:
	Analyse des changements conformationnels global:
		- calcule des RMSD par domaines
		- calcul du RMSD global
		- calcul des RMSD croisées par domaines

"""

#Usage du script
def usage():
    print("Programme = Confo_global.py\n",
          "Objectif: Analyse des changements conformationnels globaux\n",
          "To Run: python3 Confo_global.py -ref fichier1.pdb -conf fichier2 -domaines A1,A2,A3\n",
          "Arguments: - obligatoires:  -ref : nom du fichier de reférence .pdb\n",
          "                            -conf: nom du fichier des conformations\n",
		  "                            -domaines: domaine1,domaine2,etc (domaine du calcul)\n",
          "           - optionnels:  -mode: atomique_CA (par défaut), barycentre\n")
          
          
#===================Initialisation des arguments========================
#Vérifie si les deux fichiers sont rentré
try:
	infile_ref=sys.argv[sys.argv.index("-ref")+1]
except:
	print("Le nom du ficher de référence est obligatoire")
	usage()
	sys.exit()
try:
	infile_conf=sys.argv[sys.argv.index("-conf")+1]
except:
	print("Le nom du ficher des conformation est obligatoire")
	usage()
	sys.exit()

#Vérifie que les fichier finissent bien par pdb
nom_ref = infile_ref.split('.')
nom_conf = infile_conf.split('.')
if ( (nom_ref[1] != 'pdb') or (nom_conf[1] != 'pdb')):
	print ("Les fichiers doivent être au format .pdb")
	usage()
	sys.exit()

#Vérifie si un mode est rentré
try:
	mode=sys.argv[sys.argv.index("-mode")+1]
	#vérifie que le mode entré est valide
	if (mode != "atomique_CA" and mode != 'barycentre'):
		print ("Le mode doit etre atomique_CA ou barycentre")
		sys.exit()
except:
	mode="atomique_CA"



#Vérifie que les domaines sont rentrés
try:
	domaines = sys.argv[sys.argv.index("-domaines")+1]
except:
	print ("Les noms des domaines sont obligatoires")
	usage()
	sys.exit()

#Initialise les domaines
domainesUtilisés = domaines.split(',')

#=================================MAIN==================================
#Crée les deux dictionnaires et affiche une erreur si le format ou 
#le nom du fichier est incorrect
try:
	dREF=Tool.parser(infile_ref)
except:
	print ("Erreur sur le format ou le nom pour le fichier de référence")
	sys.exit()
try:
	dCONF=Tool.parserConf(infile_conf)
except:
	print("Erreur sur le format ou le nom pour le fichier de conformation")
	sys.exit()

#Ouvre un fichier de sortie et écrit les paramètres
fichier = open("results/RMSD_"+nom_conf[0], "w")
fichier.write("#REF %s CONF %s Mode %s Domaines %s\n"
%(infile_ref,infile_conf,mode,domaines))

#Ecrit le titre de chaque colonne
fichier.write("#Model ")
for domaine in domainesUtilisés :
	fichier.write("RMSD[%s] "%domaine)
fichier.write("RMSD_Global\n")

#Crée la matrice de distance global
nb_model=len(dCONF["order"])
s=(nb_model,nb_model)
MDistance=numpy.zeros(s)
#Crée un dictionnaire de matrice por stocker les matrices par domaines
listeModel=list(dCONF["order"])
dicoMatrice={}
dicoMatrice["order"]=[]
nb_domaine=len(domainesUtilisés)
for i in range(nb_domaine):
	dicoMatrice[domainesUtilisés[i]]=numpy.zeros(s)


i=0
compteur =1
#Parcours le dictionnaire du fichier des conformations
for model in dCONF["order"] :
	#Initialise le RMSD globale pour ce modèle
	somme=0
	nbPaire=0
	fichier.write("%d\t"%(model))
	for domaine in dCONF[model]["order"] :
		#Passe son tour si ce n'est pas un domaine pris en entrés
		if domaine not in domainesUtilisés :
			continue
		#Initialise le RMSD pour ce domaine
		sommeDomaine=0
		nbPaireDomaine=0
		#Si le domaine est proteique (premier résidu possède un carbone alpha)
		premierResidu = dCONF[model][domaine]["order"][0]
		if "CA" in dCONF[model][domaine][premierResidu]:
			#Pour le mode Barycentre
			if (mode=="barycentre"):
				for residu in dCONF[model][domaine]["order"] :
					sommeDomaine += Tool.barycentre(dREF,dCONF,model,domaine,residu)**2
					nbPaireDomaine +=1
			#Pour le mode Atomique_CA
			else:
				for residu in dCONF[model][domaine]["order"] :
					xConf = dCONF[model][domaine][residu]["CA"]["x"]
					yConf = dCONF[model][domaine][residu]["CA"]["y"]
					zConf = dCONF[model][domaine][residu]["CA"]["z"]
					xRef = dREF[domaine][residu]["CA"]["x"]
					yRef = dREF[domaine][residu]["CA"]["y"]
					zRef = dREF[domaine][residu]["CA"]["z"]
					sommeDomaine += Tool.distance(xConf,yConf,zConf,xRef,yRef,zRef)**2
					nbPaireDomaine +=1
		#Si le domaine n'est pas proteique
		else:
			for residu in dCONF[model][domaine]["order"] :
				sommeDomaine += Tool.barycentre(dREF,dCONF,model,domaine,residu)**2
				nbPaireDomaine +=1
		#Calcul et écrit le RMSD du domaine et ajoute les distance à la somme globale
		fichier.write("%.4f\t"%(math.sqrt(sommeDomaine/nbPaireDomaine)))
		somme += sommeDomaine
		nbPaire += nbPaireDomaine
	#Calcul le RMSD global et l'ecrit
	RMSD= math.sqrt(somme/nbPaire)
	fichier.write("%.4f\n"%(RMSD))
	
	
	#Calcul la matrice RMSD Croisé
	j=i+1
	listeModel.remove(model)
	for model2 in listeModel :
		somme=0
		nbPaire=0
		for domaine in dCONF[model]["order"] :
			if domaine not in domainesUtilisés :
				continue
			sommeDomaine=0
			nbPaireDomaine=0
			premierResidu = dCONF[model][domaine]["order"][0]
			if "CA" in dCONF[model][domaine][premierResidu]:
				for residu in dCONF[model][domaine]["order"] :
						xConf = dCONF[model][domaine][residu]["CA"]["x"]
						yConf = dCONF[model][domaine][residu]["CA"]["y"]
						zConf = dCONF[model][domaine][residu]["CA"]["z"]
						xRef = dCONF[model2][domaine][residu]["CA"]["x"]
						yRef = dCONF[model2][domaine][residu]["CA"]["y"]
						zRef = dCONF[model2][domaine][residu]["CA"]["z"]
						sommeDomaine += Tool.distance(xConf,yConf,zConf,xRef,yRef,zRef)**2
						nbPaireDomaine +=1
			else:
				for residu in dCONF[model][domaine]["order"] :
					sommeDomaine += Tool.barycentre(dCONF[model2],dCONF,model,domaine,residu)**2
					nbPaireDomaine +=1
			dicoMatrice[domaine][i][j]=dicoMatrice[domaine][j][i] = math.sqrt(sommeDomaine/nbPaireDomaine)
			somme += sommeDomaine
			nbPaire += nbPaireDomaine
		MDistance[j][i] = MDistance[i][j] = math.sqrt(somme/nbPaire)
		compteur +=1
		j+=1
	
	i += 1
fichier.close()

fichier2 = open("results/matrice_global", "w")
for i in range(len(dCONF["order"])) :
	for j in  range(len(dCONF["order"])) :
		fichier2.write("%f\t"%(MDistance[i][j]))
	fichier2.write("\n")
fichier2.close()
	
for domaine in domainesUtilisés :
	fichier3 = open("results/matrice_"+domaine, "w")
	for i in range(nb_model) :
		for j in  range(nb_model) :
			fichier3.write("%f\t"%(dicoMatrice[domaine][i][j]))
		fichier3.write("\n")
	fichier3.close()

