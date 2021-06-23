#!/bin/bash
if (( $EUID != 0 )); then
   echo "The installation must be run as root" 1>&2
   exit 1
fi
mkdir -p /opt/local/osra/2.1.1 || { echo "Cannot create /opt/local/osra folder" 1>&2; exit; } 
mkdir -p /usr/local/bin || { echo "Cannot create /usr/local/bin folder" 1>&2; exit; }
cp -f package/* /opt/local/osra/2.1.1/ || { echo "Cannot copy to /opt/local/osra folder" 1>&2; exit; }
echo "Installing binary files in /opt/local/osra"
cp -f osra /usr/local/bin || { echo "Cannot copy to /usr/local/bin" 1>&2; exit; }
echo "Installing osra script in /usr/local/bin"

