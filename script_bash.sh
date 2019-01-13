#lancement de l'image ami2b/projetpython:latest (en provenance du Hub) dans un conteneur
sudo docker run ami2b/projetpython:latest

#récupération du nom du conteneur de l'image ami2b/projetpython:latest
var=`sudo docker container ls -a | cut -d\  -f 1 | sed -n "2 p"`

#récupération du dossier results (du conteneur) sur le bureau (en local)
sudo docker cp $var:./Projet_python_big_data/results /home/coralie/Bigdata

