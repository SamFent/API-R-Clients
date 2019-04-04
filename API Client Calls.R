# API Client Calls

rm(list=ls())

source("ncbiblast.R")

ncbiblast(email= 'test@ebi.ac.uk',
          sequence= 'sp:pak4_human',
          stype= 'protein',
          database= 'uniprotkb_swissprot',
          program= 'blastp',
          outformat= c("accs", "ids"))


source("phobius.R")

phobius(email= 'test@ebi.ac.uk', 
        stype= 'protein', 
        sequence= 'sp:phyb_arath')

source("emboss_backtranseq.R")

emboss_backtranseq(email= "test@ebi.ac.uk",
               sequence= "./Sequence/prot_sequence.fasta",
               outformat= "sequence")


# Order that parameters are put in doesn't matter as long as each wanted parameter is specified in the call
# to specify mulitple outformats,  you need to do outformat= c("accs", "ids")


# Simple Workflow Example 

source("clustalo.R")
source("simple_phylogeny.R")

clustalo(email= 'test@ebi.ac.uk',
         sequence= 'sp:pak4_human,sp:pak2_human,sp:pak5_human',
         stype= 'protein',
         outformat= "aln-clustal_num",
         outfile= "clustal") & 
simple_phylogeny(email= 'test@ebi.ac.uk',
         sequence= "clustal.aln-clustal_num.clustal_num")

# Note: This gives output from both steps 