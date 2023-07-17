#!/bin/bash
docker run --rm --publish 3838:3838 --mount type=bind,source="$(pwd)"/results,target=/home/results zmxu/g2g_shiny_docker
