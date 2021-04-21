#!/usr/bin/env bash
for dir in $(ls plotVel); do
  ffmpeg -r 5 -pattern_type glob -i "plotVel/$dir/*.jpg" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" "plotVel/vel_$dir.mp4"
done
