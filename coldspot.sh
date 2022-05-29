#STEP 1 CREAT LIST OF UNIQUE S


function stepone(){
echo ">> COLDSPOT RUNNING STEP ONE"
./analz-spike spikeprot0429/spikeprot0429.fasta 1223 1323 > spikeprot0429/spikeprot-03-05-2022.dat
awk '{print $NF;}' spikeprot0429/spikeprot-03-05-2022.dat  | sort | uniq >  spikeprot0429/spikeprot-uniq.dat


}
 

function blast_uniq(){

while read s; do

echo $s > t.fasta

#SUBJECT IS REF-SPIKE

blastp -subject  data/spike_ref.fasta   -query t.fasta -outfmt 5 > t.ali


#Hsp_QSEQ is QUERY
#Hsp_HSEQ is REF
q=$(grep Hsp_qseq t.ali | awk '{a=substr($1,11); a=substr(a,1,length(a)-11);print a;}')
h=$(grep Hsp_hseq t.ali | awk '{a=substr($1,11); a=substr(a,1,length(a)-11);print a;}')


#HERE ... QUERY QUERY REF-SPIKE
echo $s $q $h | tr -s "* " " "



done < spikeprot0429/spikeprot-uniq.dat > spikeprot0429/alignment.dat

awk '{if(NF==3) print $0; }' spikeprot0429/alignment.dat > tmp
mv tmp alignment.dat

}


#RUN PIPELINE

stepone
blast_uniq

./combine-all spikeprot0429/spikeprot-03-05-2022.dat spikeprot0429/alignment.dat  20200000 20201231 > spike-2020.tsv 
./combine-all spikeprot0429/spikeprot-03-05-2022.dat spikeprot0429/alignment.dat  20200000 20211231 > spike-2021.tsv 
./combine-all spikeprot0429/spikeprot-03-05-2022.dat spikeprot0429/alignment.dat  20200000 20220430 > spike-2022-4.tsv

