for i in *.mp3; do mv "$i" `echo $i | tr ' ' '_'`; done 

for i in *.[Mm][Pp]3; do mv "$i" `echo $i | tr '[A-Z]' '[a-z]'`; done

for i in *.mp3; do lame --decode $i `basename $i .mp3`.wav; done

rm -rf *.mp3
