#!/bin/bash

KRAKEN_FILE="$1"
AS_FILE="$2"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <kraken_tsv> <adaptive_sampling_csv>"
    exit 1
fi

awk -F',' '
NR==FNR {
    if (NR > 1) decision[$1]=$2
    next
}
{
    split($0, fields, "\t")
    id=fields[2]
    sub(/^>/, "", id)

    if (id in decision) {
        if (decision[id] == "unblock") unblock++
        else if (decision[id] == "sequence") sequence++
    }
}
END {
    print "unblock:", unblock + 0
    print "sequence:", sequence + 0
    print "total matched:", unblock + sequence
}' "$AS_FILE" "$KRAKEN_FILE"
