echo "Simulating gene/species trees with SimPhy"
simphy -rs 2 -rl f:10 -sb ln:-15,1 -st u:200000,20000000 -sl f:5 -so f:1 -sp f:100000 -su f:0.00001 -si f:6 -hh ln:1.2,1 -hl ln:1.4,1 -hg f:200 -v 1 -o testwsimphy -cs 6656 -od 1 -op 1 -on 1
if [[ ! -f "INDELIble_wrapper.pl" ]]; then
    echo -e "Downloading INDELIble_wrapper.pl file from repository"
    indelibleWrapperURL="https://raw.githubusercontent.com/adamallo/SimPhy/master/scripts/INDELIble_wrapper.pl"
    if [[ $(uname -s) -eq "Linux" ]]; then
        wget $indelibleWrapperURL
    elif [[ $(uname -s) -eq "Darwin" ]]; then
        curl -0 $indelibleWrapperURL
    fi
fi
echo "Simulating DNA sequences from the gene/species trees above"
perl INDELIble_wrapper.pl testwsimphy/ control.4.txt $RANDOM 1


python ngsphy-refselector/__main__.py  -p test/testwsimphy/ -ip data -op ref -o output -m 0 -n 10 -l DEBUG
