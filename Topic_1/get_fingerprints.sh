#!/usr/bin/env bash

# configure ~/.ssh/config with the 3 hosts

hosts=(biol1
       biol2
       biol3
      )
user=ubuntu

keynames=(
    dsa
    ecdsa
    ed25519
    rsa
)

digests=(
    md5
    sha256
)

for h in "${hosts[@]}"; do
    echo "host $h:"
    for k in "${keynames[@]}"; do
	for d in "${digests[@]}"; do
	    while read size hash fqdn algo rest; do
		printf "%10s %5s %s\n" "$algo" "$size" "$hash"
	    done < <(
		ssh -o ForwardAgent=no -o ForwardX11=no "$h" ssh-keygen -E $d "-lf" /etc/ssh/ssh_host_"$k"_key.pub
	    )
	done
    done
    echo
done
