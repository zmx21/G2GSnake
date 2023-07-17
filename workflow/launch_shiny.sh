#!/bin/bash
docker run --rm --publish 3838:3838 --mount type=bind,source="$(pwd)"/results,target=/home/results -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) zmxu/g2g_shiny_docker
