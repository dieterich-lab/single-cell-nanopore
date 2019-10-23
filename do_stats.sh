#!/bin/bash

if [ ! $# == 4 ]; then
  echo "Usage: $0 [Primer alignment log] [Barcode alignment left log] [Barcode alignment right log] [log from execution]"
  exit
fi

grep 'Finished extraction' $4

echo "P1 aligned in forward direction on the left piece of the nanopore read (cDNA is located in the middle)"
grep '_Flexbar_removal_cmdline$' $1 | wc -l

echo "P1 aligned in reverse complement direction on the right piece of the nanopore read"
grep '_Flexbar_removal_cmdline_rc$' $1 | wc -l

echo "Assigned Barcode in forward direction"
grep 'Best alignment' $2 | wc -l

echo "Could not assign barcode in forward direction:"
grep 'Unvalid alignment' $2 | wc -l

echo "Assigned Barcode in reverse complement direction"
grep 'Best alignment' $3 | wc -l

echo "Could not assign barcode in reverse complement direction:"
grep 'Unvalid alignment' $3 | wc -l


echo "Unique assigned reads: "
cat $2 $3 | grep '  read id' | awk '{print $3}' | cut -d '_' -f 1 | sort -u | wc -l



