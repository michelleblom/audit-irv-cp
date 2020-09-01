for d in */ ; do
    bn=`basename $d`
    cd ${bn}
    echo "${bn}"
    python3 ../translate.py "${bn}_statewide.totals" > "${bn}_statewide.raire"
    cd ..
done
