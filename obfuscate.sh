#!/bin/sh
COBF=/usr/local/cobf-1.06
PDA=/home/adriano/empresas/eMediaMusic/library/v1.0-r16
OUTPUT=$PDA/obfuscated/cobf
ORIGINAL=$PDA/original
OTHER=$PDA/additional_files

rm -rf $OUTPUT
mkdir $PDA/obfuscated
mkdir $OUTPUT
cd $ORIGINAL
cobf mitre_pda.c -hi mitre_pda.h -hi _mitre_pda.h -p $COBF/etc/pp_gnu -b -o $OUTPUT -m $COBF/etc/cpp.tok -m $COBF/etc/cansilib.tok -m $PDA/pda_exceptions
mv $OUTPUT/mitre_pda.c $OUTPUT/pda.c
cp $ORIGINAL/mitre_pda.h $OUTPUT/pda.h
cp $OTHER/example.c $OTHER/Makefile $OUTPUT
rm $OUTPUT/cobf.log
