#If you didn't specify --noSplitSubreads, you'll need to us the sed
#Else, just sort

for i in `ls *.m4`
do
    echo `date` Starting Sorting $i
    #sed 's/\/[0-9]\+_[0-9]\+ / /' $i | sort -k1 -d  > tmp
    sort -k1 -d  $i > tmp
    mv tmp $i
    echo `date` Finished Sorting $i
done
