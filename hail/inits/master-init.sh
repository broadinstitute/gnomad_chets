#!/bin/bash

export HAIL_HOME=/home/hail
export SPARK_HOME=/usr/lib/spark

HASH=$(gsutil cat gs://hail-common/latest-hash.txt)
JAR=hail-hail-is-master-all-spark2.0.2-${HASH}.jar
ZIP=pyhail-hail-is-master-${HASH}.zip

mkdir -p ${HAIL_HOME}
gsutil cp gs://hail-common/${JAR} gs://hail-common/${ZIP} ${HAIL_HOME}

ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

cat <<EOT >> /etc/bash.bashrc
export HAIL_HOME=/home/hail
export SPARK_HOME=/usr/lib/spark
export PATH=/home/anaconda2/bin:${PATH}
export _JAVA_OPTIONS='-Xmx8096m'
$(ls ${HAIL_HOME})
export PYTHONPATH=${SPARK_HOME}/python:$(ls ${SPARK_HOME}/python/lib/py4j-*-src.zip):$(ls ${HAIL_HOME}/pyhail-*zip)
export SPARK_CLASSPATH=$(ls ${HAIL_HOME}/hail-*jar)
EOT

    pip install slackclient sklearn
    gsutil cp gs://gnomad-public/tools/inits/gnomad-init.sh .
    chmod +x gnomad-init.sh
    ./gnomad-init.sh > startup.log 2>&1 &
fi