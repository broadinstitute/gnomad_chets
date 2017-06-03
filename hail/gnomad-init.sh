#!/usr/bin/env bash

# iPython and hail master as quickly as possible
apt-get install -y ipython

pip install slackclient pandas scipy sklearn

export SPARK_HOME=/usr/lib/spark
export HAIL_HOME=/hadoop_gcs_connector_metadata_cache/hail
export HAIL_HASH=$(gsutil cat gs://hail-common/latest-hash.txt)
export HAIL_JAR=hail-hail-is-master-all-spark2.0.2-${HAIL_HASH}.jar
export HAIL_PYTHON_ZIP=pyhail-hail-is-master-${HAIL_HASH}.zip
mkdir $HAIL_HOME
gsutil cp gs://hail-common/${HAIL_JAR} gs://hail-common/${HAIL_PYTHON_ZIP} $HAIL_HOME

# Prepare bashrc and redownload script
cat <<EOT >> /etc/bash.bashrc
export SPARK_HOME=${SPARK_HOME}
export HAIL_HOME=${HAIL_HOME}
export HAIL_HASH=${HAIL_HASH}
export HAIL_JAR=${HAIL_JAR}
export HAIL_PYTHON_ZIP=${HAIL_PYTHON_ZIP}
export _JAVA_OPTIONS='-Xmx8096m'
export PYTHONPATH=${SPARK_HOME}/python:$(ls ${SPARK_HOME}/python/lib/py4j-*-src.zip):${HAIL_HOME}/${HAIL_PYTHON_ZIP}
export SPARK_CLASSPATH=${HAIL_HOME}/${HAIL_JAR}
EOT
cat <<EOT >> /redownload_hail.sh
mkdir -p $HAIL_HOME
gsutil cp gs://hail-common/${HAIL_JAR} gs://hail-common/${HAIL_PYTHON_ZIP} $HAIL_HOME
EOT
chmod +x /redownload_hail.sh

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

# Building Hail. Why not.
apt-get install -y cmake tmux
git clone https://github.com/hail-is/hail.git
cd hail
./gradlew shadowJar
mkdir /hadoop_gcs_connector_metadata_cache/hail_build/
cp -r * /hadoop_gcs_connector_metadata_cache/hail_build/

# curl http://www.freefontspro.com/d/14454/arial.zip > /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip
# unzip /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip