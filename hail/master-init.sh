#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then
    gsutil cp gs://gnomad-public/tools/inits/gnomad-init.sh .
    ./gnomad-init.sh &
fi

