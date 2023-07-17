#!/bin/bash
FILE="$(pwd)"/results/G2G_Results.bmat

file_permission=$(stat -f "%A" $FILE)
third_digit=$((file_permission % 10))


if [[ $third_digit -lt 4 ]]; then
  echo "Please add read permission: chmod a+r ./results/G2G_Results.*"
  exit 0
fi

docker run --rm --publish 3838:3838 --mount type=bind,source="$(pwd)"/results,target=/home/results zmxu/g2g_shiny_docker
