#!/usr/bin/env bash

# iPython and hail master as quickly as possible
apt-get install -y tmux

# R stuff from here
apt-get install -y gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev
export RSTUDIO_DEB=rstudio-server-1.0.143-amd64.deb
wget https://download2.rstudio.org/${RSTUDIO_DEB}
gdebi -n ${RSTUDIO_DEB}
mkdir /opt/spark
ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

# Common R packages, tiered for faster startup
R --vanilla -e "install.packages(c('sparklyr', 'dplyr'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('magrittr', 'ggplot2', 'slackr', 'ggrepel'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('plyr', 'shiny', 'plotly'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('DT', 'tidyverse', 'broom', 'randomForest', 'ROCR', 'shinythemes', 'devtools'), repos='https://cran.rstudio.com')"

# curl http://www.freefontspro.com/d/14454/arial.zip > /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip
# unzip /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip