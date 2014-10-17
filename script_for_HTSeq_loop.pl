#script for running HTSeq in a loop 

for file in `ls -1 *_1.sam`
do
PAIR=$file
echo "$PAIR"
python /home/talwar/.local/bin/htseq-count -s no $PAIR Caenorhabditis_elegans.WBcel235.76.gtf > count_$PAIR.txt
done

