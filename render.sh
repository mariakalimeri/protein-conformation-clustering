#!/bin/bash

freqStride=$1
noOfClusters="$(grep 'END' leadersTemp.pdb | wc -l | tr -d '[[:space:]]')"

echo $noOfClusters

start=1

for xAngle in 0 90
do

for yAngle in 0 90
do
# Write customized .vmdrc file
cat > .vmdrc <<EOL
# Set display to orthographic
display projection orthographic
display Depth Cueing off
menu main on
mol delrep 0 top
mol representation NewCartoon
color scale method RGryB
mol color          Timestep
mol material       AOChalky
mol addrep         top
mol drawframes 0 0 0:${freqStride}:${noOfClusters}
rotate y by ${yAngle}
rotate x by ${xAngle}
color Display Background white
display settings Shadows on
display settings Amb. Occl. On
display settings Fec. On
axes location LowerLeft
EOL
cat .vmdrc

# Write a file that contains the rendering command
# Command for mac: 
# render Tachyon scene.dat "/Applications/VMD1.9.3.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 %s -format TGA -res 1024 1024 -o %s.tga"
# Command for linux:
# render Tachyon scene.dat "/usr/local/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TGA -res 1024 1024 -o %s.tga"
cat > renderCommandFile <<EOL
render Tachyon scene.dat "/Applications/VMD1.9.3.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 %s -format TGA -res 1024 1024 -o %s.tga"
quit
EOL
cat renderCommandFile

# Run the whole thing
# Note: The height flag has not been tested for different proteins. There might be cases that don't fit in the screen. 
# For mac
"/Applications/VMD1.9.3.app/Contents/MacOS/startup.command" -dispdev text leadersTemp.pdb -nt -startup .vmdrc -e renderCommandFile -height 5
# For linux
#vmd -dispdev text leadersTemp.pdb -nt -startup .vmdrc -e renderCommandFile -height 4

mv scene.dat.tga renderedImage${start}.tga

rm .vmdrc renderCommandFile scene.dat

((start++))

done
done

mv renderedImage*tga www/
