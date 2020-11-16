#!/bin/bash



# For obtaining the Sequences CH fis
#--------TRAINING:  RUN THIS TRAIN AND SAVE VECTORS ----------
# -i or --speciesfile  <The species input fasta file>   ... training FASTA file
# -d or --dirname      <dir where the training file is> ... training dir 
# -t or --traintype    <Train type>                     ... key value in dictionary

./makeTrainingVecs.py \
    -i train-CH-exon-fish.fasta \
    -d chsfish \
    -t chsfish

#--------RUN THIS WHEN YOU HAVE THE TRAINING VECTORS SAVED ----------
# Analysis of CH in fish
# -i or --speciesfile  <The species fasta file>   ... genome file    (Fasta input genome file)
# -d or --speciesdir   <The species directory>    ... directory for input/output  (where the genome is)
# -T or --trainingdir  <The training directory>   ... directory of stored Training matrices
# -E or --exontype     <Training Type>            ... This is the key value in the dictionary
# -q or --queryfile    <Query File>               ... fasta query  (if Tblast table is to generated)

#./exonFinder.py \
#    -i fTakRub1.pri.cur.20190523.fasta \
#    -d fish/Takifuga_rubripes \
#    -T chsfish \
#    -E chsfish \
#    -q query-IG-fish.fasta

