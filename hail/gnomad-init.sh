#!/usr/bin/env bash

apt-get install -y ipython cmake gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev tmux

easy_install pip
pip install slackclient

wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
gdebi -n rstudio-server-1.0.136-amd64.deb
mkdir /opt/spark
ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

git clone https://github.com/hail-is/hail.git

cd hail
./gradlew shadowJar