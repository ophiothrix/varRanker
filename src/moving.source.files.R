## To make things more manageable, keep "data" directory only for the files that will need to be distributed, i.e. can't be downloaded on-flight. Move all other files to the cache.

sed -i ".bak" 's/data\/ENCODE/cache\/ENCODE/g' *
rm *.bak
sed -i ".bak" 's/data\/conservation.scores/cache\/conservation.scores/g' *
rm *.bak
sed -i ".bak" 's/data\/snps.multicell.bySample/cache\/snps.multicell.bySample/g' *
rm *.bak

sed -i ".bak" 's/data\/ng.3432-S5.txt/cache\/ng.3432-S5.txt/g' *
rm *.bak

grep "ng.3432-S5.txt" *
	snps.multicell.bySample


