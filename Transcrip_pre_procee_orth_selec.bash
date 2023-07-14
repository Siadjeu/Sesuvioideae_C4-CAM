##### Processing initial transcriptome

### Remove all characters after space

sed '/^>/ s/ .*//' Results_Trinity_SSall.Trinity.fasta > Trinity_SSall.fasta
sed '/^>/ s/ .*//' Trinity_cdhit_PP.fasta > cdhit_PP.fasta


### Cluster transcripts with at least 98% similarity using cd-hit-est

cd-hit-est -i Trinity_SSall.fasta -c 0.98 -o Trinity_cdhit_SSall.fasta

nohup cd-hit-est -i Results_Trinity_PP.Trinity.fasta -c 0.98 -o Trinity_cdhit_PP.fasta &

### Select transcripts harboring open reading frames with TransDecoder

# Detect the longest open reading frame

nohup TransDecoder.LongOrfs -t Trinity_cdhit_SSall.fasta --output_dir SSall_transdeco &

nohup TransDecoder.LongOrfs -t cdhit_PP.fasta --output_dir PP_transdeco &

# Predict the open reading frame (Using precomputed blastX alignments to UniProt database to improve prediction)


#!/bin/bash
#SBATCH --partition mpcb.p           # or mpcp.p
#SBATCH --ntasks 1                          # instead of node
#SBATCH --cpus-per-task 16           # use 40 with mpcp.p
##SBATCH --exclusive                      # not needed if you ask for all the cores in a node
#SBATCH --time 1-0:00:00
#SBATCH --mem-per-cpu 30000  # not needed if you use all the cores, you also get all memory per default
#SBATCH --job-name tunnel
#SBATCH --output BlastTrans-log-%J.txt



ml hpc-env/6.4

ml BLAST+/2.10.0-intel-2018a

makeblastdb -in uniprot_sprot.fasta -dbtype prot

blastp -query /jessica/Night/Salmon/transdeco/SSall_transdeco/longest_orfs.pep \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > SSall.outfmt6
	
blastp -query /jessica/Night/Salmon/transdeco/PP_transdeco/longest_orfs.pep \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > PP.outfmt6


TransDecoder.Predict -t Trinity_cdhit_SSall.fasta --retain_blastp_hits SSall.outfmt6 --output_dir SSall_transdeco
TransDecoder.Predict -t cdhit_PP.fasta --retain_blastp_hits SSall.outfmt6 --output_dir PP_transdeco


### Processing the output of TransDecoder

# Remove all characters after space in the peptide sequences

sed '/^>/ s/ .*//' cdhit_PP.fasta.transdecoder.pep > orth_PP.fasta
sed '/^>/ s/ .*//' Trinity_cdhit_SSall.fasta.transdecoder.pep > orth_SS.fasta

### Get IDs of the protein sequences and remove ">" in front of the IDs
## First method (two steps)
# Get IDs
grep -e ">" orth_PP.fasta > PP_ID.txt
grep -e ">" orth_SS.fasta > SS_ID.txt
# Use awk to remove ">" in front of the IDs
awk 'sub(/^>/, "")' PP_ID.txt > PP_IDfinal.txt
awk 'sub(/^>/, "")' SS_ID.txt > SS_IDfinal.txt
## Second method (one step)
# Get IDs and remove ">" (Just for reference, not used in the script)
# Extract IDs header from a fasta file without ">"
grep '^>' Trinity_SS.TransdecoderCdHit.fasta | sed 's,>,,g' > ID_SS.transdecoCdhit.txt

## Use sed to remove everything after "." in the ID file since ".p" was added (This introduces duplicates because a transcript can have more than two proteins, differentiated with p1 or p2,...

sed 's/\..*$//' PP_IDfinal.txt > PP_IDcleaned.txt
sed 's/\..*$//' SS_IDfinal.txt > SS_IDcleaned.txt

## Remove duplicate IDs from the TransDecoder output
# Another way to remove duplicates is to use "sort -u PP_IDcleaned.txt > PP_tes1.txtqq" (Here the file is sorted)
uniq PP_IDcleaned.txt > PP_IDuniq.txt
uniq SS_IDcleaned.txt > SS_IDuniq.txt

## Extract contigs with ORFs from the clustered transcriptome using seqtk

/seqtk/./seqtk subseq cdhit_PP.fasta PP_IDuniq.txt > PP_for_salmonUniq.fasta

/seqtk/./seqtk subseq Trinity_cdhit_SSall.fasta SS_IDuniq.txt > SS_for_salmonUniq.fasta

## Check for duplicate sequences in the fasta file

cat PP_for_salmonUniq.fasta | grep '^>' | sort | uniq -d
cat SS_for_salmonUniq.fasta | grep '^>' | sort | uniq -d

###############################################################################################

##### Search for orthologous transcripts between PP and SS with Orthofinder

### Add species ID at the beginning of the fasta file header for OrthoFinder to ease identification
# Change the file names if necessary
cp PP_for_salmonUniq.fasta PP_for_orthfinUniq.fasta
cp SS_for_salmonUniq.fasta SS_for_orthfinUniq.fasta

# Add species prefix in front of the transcript IDs
perl -pi -e "s/^>/>PP_/g" PP_for_orthfinUniq.fasta
perl -pi -e "s/^>/>SS_/g" SS_for_orthfinUniq.fasta


### Run OrthoFinder to detect orthologous sequences between the two species
ml hpc-env/6.4
ml BLAST+/2.10.0-intel-2018a

nohup orthofinder -d -f /gss/work/wupf2892/jessica/Night/Salmon/transdeco/orthfinder_New/ &

### Process only orthologous sequences

## Extract only orthologous transcripts
# First extract IDs of orthologous sequences
# The "\b" matches a word boundary, i.e., the beginning of a word. Then it has to start with "SS_" followed by any number of other word characters

grep -o '\bPP_\w*' PP_for_orthfinUniq__v__SS_for_orthfinUniq.tsv > Orth_ID-PP.txt
grep -o '\bSS_\w*' SS_for_orthfinUniq__v__PP_for_orthfinUniq.tsv > Orth_ID-SS.txt

# Remove the added species prefix "PP_" and "SS_" in front of the transcript IDs
awk 'sub(/^PP_/, "")' Orth_ID-PP.txt > Orth_IDcl-PP.txt
awk 'sub(/^SS_/, "")' Orth_ID-SS.txt > Orth_IDcl-SS.txt

## Extract orthologous sequences from the clustered transcriptome for Salmon analysis
/gss/work/wupf2892/seqtk/./seqtk subseq /gss/work/wupf2892/jessica/Night/Salmon/transdeco/cdhit_PP.fasta Orth_IDcl-PP.txt > Orth_Sel_forSalm_PP.fasta
/gss/work/wupf2892/seqtk/./seqtk subseq /gss/work/wupf2892/jessica/Night/Salmon/transdeco/Trinity_cdhit_SSall.fasta Orth_IDcl-SS.txt > Orth_Sel_forSalm_SS.fasta


#############################################################

### Run Salmon for counting

# For Sesuvium, we will run two times: one with only orthologous transcripts for comparison with Portulacastrum

# Portulacatrum

#!/bin/bash
#SBATCH --partition mpcb.p           # or mpcp.p
#SBATCH --ntasks 1                          # instead of node
#SBATCH --cpus-per-task 16           # use 40 with mpcp.p
##SBATCH --exclusive                      # not needed if you ask for all the cores in a node
#SBATCH --time 1-0:00:00
##SBATCH --mem-per-cpu 30000  # not needed if you use all the cores, you also get all memory per default
#SBATCH --job-name tunnel
#SBATCH --output Salmon-log-%J.txt

conda activate rna

salmon index -t Orth_Sel_forSalm_PP.fasta -i PP_salmon_index

for sample in *_1.fq.gz*; do
  base=$(basename $sample "_1.fq.gz")
  salmon quant -i PP_salmon_index -l A \
    -1 ${base}_1.fq.gz \
    -2 ${base}_2.fq.gz \
    -p 16 --validateMappings --gcBias --numGibbsSamples 100 -o quant/${base}_quant
done



##### Get the unique ID for differential expression analysis

### Here we use only Trinity IDs: IDs of SS were chosen
## We will blast all orthologs between PP and SS to get 1:1 orthologous relationship

ml hpc-env/6.4
ml BLAST+/2.10.0-intel-2018a

# Create Araport database

makeblastdb -in Orth_Sel_forSalm_SS.fasta -out Orth_Sel_forSalm_SSdb -dbtype nucl

nohup blastn -query Orth_Sel_forSalm_PP.fasta -db Orth_Sel_forSalm_SSdb  -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 20 -out PPvsSS5.out &

## To select the best blast hit

export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr PPvsSS5.out | sort -u -k1,1 --merge > bestHits_PPvsSS5.txt

# Print only the first and second columns of the blast results
awk '{print $1,$2}' bestHits_PPdbvsSS5.txt > PPdbvsSS5_ID_dup.txt

# Remove orthologs with not a 1:1 relationship and keep only the first occurrence (We used SS as the database for blast)
awk '!seen[$2]++' PPdbvsSS5_ID_dup.txt > Test22.txt

# Check if there are no duplicates
# First extract the first column
awk '{print $1}' Test22.txt > Check1.txt
# Check for duplicates and count the number of lines 
sort Check1.txt | uniq -u > Check1_uniq.txt
wc -l Check1_uniq.txt

# For the second column
awk '{print $2}' Test22.txt > Check2.txt
# Check for duplicates and count the number of lines 
sort Check2.txt | uniq -u > Check2_uniq.txt
wc -l Check2_uniq.txt
#

# Extract the unique common SS IDs from the quant.sf file of SS replicates
awk '{print $1}' Test22.txt > SS_comon_II.txt
perl -lanwe 'print if $a{$F[0]}++ == 1;' SS_comon_II.txt quant.sf > quantfi.sf

# Add the headers
echo -e "Name\tLength\tEffectiveLength\tTPM\tNumReads" | cat - quantfi.sf > quantFinal.sf

# Extract PP IDs from the quant.sf file
# Extract the unique IDs of PP to be used for comparison that have a 1:1 relationship with SS
awk '{print $2}' Test22.txt > PP_comon_II.txt

## Extract common PP IDs from the quant.sf file of PP replicates
perl -lanwe 'print if $a{$F[0]}++ == 1;' PP_comon_II.txt quant.sf > quantfi.sf


### Create a replace file with PP IDs in the first column and SS IDs in the second column
# Here we will change PP IDs with SS IDs in all replicate files
# Here, we will show it for one replicate
awk '{print $2,$1}' Test22.txt > replace_ID.txt
# The replace file will be used to change PP IDs in the quant files with the corresponding SS IDs
# Replace PP IDs with SS IDs in PP quant.sf file replicates
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} {if ($1 in a){$1=a[$1]}; print $0}' replace_ID.txt quantfi.sf > quantfin.sf

# Add header
echo -e "Name\tLength\tEffectiveLength\tTPM\tNumReads" | cat - quantfin.sf > quantFinal.sf

# Change all spaces to tab delimited because it may occur that some rows are not in tab delimited format
awk -v OFS="\t" '$1=$1' quantFinal.sf > quantFinalcor.sf
