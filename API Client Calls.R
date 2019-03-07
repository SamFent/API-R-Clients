# API Client Calls

rm(list=ls())

source("ncbiblast.R")

ncbiblast(email= 'test@ebi.ac.uk',
          sequence= 'sp:pak4_human',
          stype= 'protein',
          database= 'uniprotkb_swissprot',
          program= 'blastp')


source("phobius.R")

phobius(email= 'test@ebi.ac.uk', 
        stype= 'protein', 
        sequence= 'sp:phyb_arath')

# Order that parameters are put in doesn't matter if the parameter being given a value is specified 


