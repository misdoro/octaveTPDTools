#!/bin/bash
echo "#!/bin/bash" > displayTPD.sh
echo displayTPD $@ >> displayTPD.sh
chmod +x displayTPD.sh
displayTPD.m $@
