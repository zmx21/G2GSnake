#!/bin/bash
FILE="$(pwd)"/result/G2G_Results.bmat

if ! [[ $(stat -c "%A" $FILE) =~ "r" ]]; then
  echo "Please add read permission: chmod a+r ./results/G2G_Results.*"
fi

docker run --rm --publish 3838:3838 --mount type=bind,source="$(pwd)"/results,target=/home/results zmxu/g2g_shiny_docker
