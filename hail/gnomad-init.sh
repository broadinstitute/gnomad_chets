#!/usr/bin/env bash

sudo apt-get install ipython cmake gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev tmux

sudo easy_install pip
sudo pip install slackclient

wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
sudo gdebi -n rstudio-server-1.0.136-amd64.deb
sudo mkdir /opt/spark
sudo ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

git clone https://github.com/hail-is/hail.git

cd hail
./gradlew shadowJar