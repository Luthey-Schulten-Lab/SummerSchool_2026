#!/bin/bash
file_name=frame
ffmpeg  -i $file_name%04d.tga -b:v 6000k -preset veryslow  -vf "pad=ceil(iw/4)*4:ceil(ih/4)*4" -c:v libx264 -pix_fmt yuv420p replicate_DNA.mp4