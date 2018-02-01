#!/bin/bash
rm all.tissues.bed
for fls in `ls ./cache/ENCODE/chromHMM.calls/*.bed.gz`
do 
	echo $fls
	tname=`echo $fls | sed "s/.*\///" | sed "s/_.*//"`
	echo $tname
	gunzip -c $fls | sed "s/$/	$tname/" >> all.tissues.bed
done

