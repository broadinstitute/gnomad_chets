#!/usr/bin/env bash

apt-get install -y ipython cmake gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev tmux

pip install slackclient

wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
gdebi -n rstudio-server-1.0.136-amd64.deb
mkdir /opt/spark
ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

git clone https://github.com/hail-is/hail.git

cd hail
./gradlew shadowJar

cp -r /hail /hadoop_gcs_connector_metadata_cache/
cat <<EOT >> /etc/bash.bashrc
export SPARK_HOME=/usr/lib/spark
export HAIL_HOME=/hadoop_gcs_connector_metadata_cache/hail
export _JAVA_OPTIONS='-Xmx8096m'
export PYTHONPATH=$SPARK_HOME/python:`ls $SPARK_HOME/python/lib/py4j-*-src.zip`:$HAIL_HOME/python:$PYTHONPATH
export SPARK_CLASSPATH=$HAIL_HOME/build/libs/hail-all-spark.jar
#export HAIL_SPARK_PROPERTIES=spark.shuffle.compress=true,spark.kryoserializer.buffer.max=512m,spark.driver.memory=64G,spark.akka.frameSize=1024,spark.driver.maxResultSize=20G,spark.rdd.compress=true,spark.ui.port=54054,spark.executor.memory=200g
EOT

R --vanilla -e "install.packages(c('sparklyr', 'plyr', 'dplyr', 'ggplot2', 'magrittr', 'shiny', 'plotly', 'slackr', 'DT', 'tidyverse', 'broom', 'ggrepel', 'randomForest', 'ROCR', 'shinythemes', 'devtools'), repos='https://cran.rstudio.com')"
