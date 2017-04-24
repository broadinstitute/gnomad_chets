#!/usr/bin/env bash

apt-get install -y ipython cmake gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev tmux

pip install slackclient pandas

export RSTUDIO_DEB=rstudio-server-1.0.143-amd64.deb
wget https://download2.rstudio.org/${RSTUDIO_DEB}
gdebi -n ${RSTUDIO_DEB}
mkdir /opt/spark
ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

export SPARK_HOME=/usr/lib/spark
export HAIL_HOME=/hadoop_gcs_connector_metadata_cache/hail
export HAIL_HASH=`gsutil cat gs://hail-common/latest-hash.txt`
export HAIL_JAR=hail-hail-is-master-all-spark2.0.2-${HAIL_HASH}.jar
export HAIL_PYTHON_ZIP=pyhail-hail-is-master-${HAIL_HASH}.zip

mkdir $HAIL_HOME
gsutil cp gs://hail-common/${HAIL_JAR} gs://hail-common/${HAIL_PYTHON_ZIP} $HAIL_HOME

# Prepare bashrc
cat <<EOT >> /etc/bash.bashrc
export SPARK_HOME=/usr/lib/spark
export HAIL_HOME=/hadoop_gcs_connector_metadata_cache/hail
export _JAVA_OPTIONS='-Xmx8096m'
export PYTHONPATH=${SPARK_HOME}/python:`ls ${SPARK_HOME}/python/lib/py4j-*-src.zip`:${HAIL_HOME}/${HAIL_PYTHON_ZIP}
export SPARK_CLASSPATH=${HAIL_HOME}/${HAIL_JAR}
EOT

# Common R packages, tiered for faster startup
R --vanilla -e "install.packages(c('sparklyr', 'dplyr'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('magrittr', 'ggplot2', 'slackr', 'ggrepel'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('plyr', 'shiny', 'plotly'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('DT', 'tidyverse', 'broom', 'randomForest', 'ROCR', 'shinythemes', 'devtools'), repos='https://cran.rstudio.com')"

# Building Hail. Why not.
git clone https://github.com/hail-is/hail.git
cd hail
./gradlew shadowJar
mkdir /hadoop_gcs_connector_metadata_cache/hail_build/
cp -r * /hadoop_gcs_connector_metadata_cache/hail_build/
