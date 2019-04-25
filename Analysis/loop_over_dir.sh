# This script can loop over the root files in root data/
# The characteristic of TCahin will merge the Tree under all files
# So this is a bypass way for dealing with different data sheet without
# merge them.
# please do "make" before run it

#FILES="CMSSW_output_root/NTU_Inj_Data/V3PCB/unpack_output/*.root"
FILES="CMSSW_output_root/TBHexaboard/module120/unpack_output/*.root"
i = 0
make
for f in $FILES
do
    echo "Processing file $f ... "
    ls $f > data_input.txt
    ./makePlots 0
done
