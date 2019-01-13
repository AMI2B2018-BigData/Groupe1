### Projet BigData / Classe Inversée
Projet de présentation de Docker

Auteurs: AGLAVE Marine, JELASSi Sarah, ROGIER Alice et ROHMER Coralie

Date: 11/01/2019

******************
Cette présentation a pour but de présenter l'outil Docker à l'ensemble de la classe.
L'exemple d'utilisation de Docker est basé sur la reproductibilité et la portabilité du projet python de l'an dernier.

******************
#Projet python:
Cette étude avait pour but l'identification des fonctions des domaines du complexe sRNP H/ACA

- Confo_global.py permet de réaliser une recherche des domaines flexibles dans la protéine
- Tool.py contient plusieurs outils nécessaires à Confo_global.py.

Usage:
	python3 Confo_global.py -ref pab21_structure_de_ref.pdb -conf pab21_prod_solute_500frames.pdb -domaines A1,A2,A3,A4,B

****************
#Docker usage:

Installation Docker:
https://docs.docker.com/install/

Ouvrir un terminal depuis le dossier du projet.

Construction de l'image:
sudo docker build -t projetpython .

Tag de l'image:
sudo docker tag projetpython ami2b/projet_python:latest

Upload de l'image sur le Hub Docker:
sudo docker push ami2b/projet_python:latest

Lancement de l'image en local:
sudo docker run projetpython

Lancement de l'image à partir du Hub Docker:
sudo docker run ami2b/projet_python:latest

Pour lancer l'image à partir du Hub Docker et récupérer les résultats automatiquement:
Modifier /home/coralie/Bigdata du fichier script_bash.sh par le chemin de votre choix.
Puis dans un terminal, lancer: bash script_bash.sh
NB: Il faudra les droits administrateur pour pouvoir supprimer les résultats.
****************
