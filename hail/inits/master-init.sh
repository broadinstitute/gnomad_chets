#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    export PATH=/home/anaconda2/bin:${PATH}
    export SPARK_HOME=/usr/lib/spark
    export HAIL_HOME=/home/hail

cat <<EOT >> /etc/bash.bashrc
export PATH=/home/anaconda2/bin:${PATH}
export SPARK_HOME=${SPARK_HOME}
export HAIL_HOME=${HAIL_HOME}
export _JAVA_OPTIONS='-Xmx8096m'
export PYTHONPATH=${SPARK_HOME}/python:$(ls ${SPARK_HOME}/python/lib/py4j-*-src.zip):$(ls ${HAIL_HOME}/pyhail-*zip)
export SPARK_CLASSPATH=${HAIL_HOME}/${HAIL_JAR}
EOT

    pip install slackclient sklearn
    gsutil cp gs://gnomad-public/tools/inits/gnomad-init.sh .
    chmod +x gnomad-init.sh
    ./gnomad-init.sh > startup.log 2>&1 &
fi