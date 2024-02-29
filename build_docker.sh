#!/bin/bash -ex

CONTAINER_PROJECT="broad-qqin"
file="gaftools:0.1.0"

gcloud builds submit "." \
       --tag="us.gcr.io/${CONTAINER_PROJECT}/${file}"
