# NCBI Blast R Client

# Load Required Libaries 
library(RCurl)
library(readr)
library(stringr)

# ... in a function allows things to be past that aren't specified in the function

ncbiblast <- function(email= NULL, 
                      stype= NULL, 
                      sequence= NULL,
                      program= NULL,
                      database= NULL,
                      outformat= NULL,
                      outfile = NULL,
                      paramDetail= NULL,
                      alignments= 50,
                      matrix= NULL,
                      scores= 50,
                      exp= "10",
                      dropoff= NULL,
                      gapopen= NULL,
                      gapext= NULL,
                      filter= NULL,
                      seqrange= "START-END",
                      gapalign= NULL,
                      compstats= NULL,
                      align= NULL,
                      transltable= NULL,
                      task= NULL,
                      match_scores= NULL,
                      ...){
  # Get list of arguments provided 
  args <- c(as.list(environment()),list(...))
  parameterlist <- c("^sequence$", "^params$", "^stype$", "^email$", "^program$", "^database$", "^outformat$", 
                     "^outfile$", "^paramDetail$", "^help$", "^resultTypes$", "^alignments$", "^matrix$", "^scores$",
                     "^exp$", "^dropoff$", "^gapopen$", "^gapext$", "^filter$", "^seqrange$", "^gapalign$",
                     "^compstats$", "^align$", "^transltable$", "^task$", "^match_scores$")
  paramComp <- grepl(paste(parameterlist, collapse = "|"), names(args))
  paramComp <- grep("FALSE", paramComp)
  # Base URL
  baseURL <- "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
  usage <- "EMBL-EBI NCBI Blast R Client:
    
 Sequence similarity search with NCBI Blast.
  
 [Required (for job submission)]
  email=                E-mail address.
  program=              The BLAST program to be used for the Sequence Similarity
                        Search.
  stype=                Indicates if the sequence is protein or DNA/RNA.
  sequence=             The query sequence can be entered directly into this form.
                        The sequence can be in GCG, FASTA, EMBL (Nucleotide only),
                        GenBank, PIR, NBRF, PHYLIP or UniProtKB/Swiss-Prot (Protein
                        only) format. A partially formatted sequence is not
                        accepted. Adding a return to the end of the sequence may
                        help certain applications understand the input. Note that
                        directly using data from word processors may yield
                        unpredictable results as hidden/control characters may be
                        present.
  database=             Database.
  
 [Optional]
  task=                 Task option (only selectable for blastn).
  matrix=               (Protein searches) The substitution matrix used for scoring
                        alignments when searching the database.
  alignments=           Maximum number of match alignments reported in the result
                        output. 
  scores=               Maximum number of match score summaries reported in the
                        result output.
  exp=                  Limits the number of scores and alignments reported based on
                        the expectation value. This is the maximum number of times
                        the match is expected to occur by chance.
  dropoff=              The amount a score can drop before gapped extension of word
                        hits is halted.
  match_scores=         (Nucleotide searches) The match score is the bonus to the
                        alignment score when matching the same base. The mismatch is
                        the penalty when failing to match.
  gapopen=              Penalty taken away from the score when a gap is created in
                        sequence. Increasing the gap openning penalty will decrease
                        the number of gaps in the final alignment.
  gapext=               Penalty taken away from the score for each base or residue
                        in the gap. Increasing the gap extension penalty favors
                        short gaps in the final alignment, conversly decreasing the
                        gap extension penalty favors long gaps in the final
                        alignment.
  filter=               Filter regions of low sequence complexity. This can avoid
                        issues with low complexity sequences where matches are found
                        due to composition rather than meaningful sequence
                        similarity. However in some cases filtering also masks
                        regions of interest and so should be used with caution.
  seqrange=             Specify a range or section of the input sequence to use in
                        the search. Example: Specifying '34-89' in an input sequence
                        of total length 100, will tell BLAST to only use residues 34
                        to 89, inclusive.
  gapalign=             This is a true/false setting that tells the program the
                        perform optimised alignments within regions involving gaps.
                        If set to true, the program will perform an alignment using
                        gaps. Otherwise, if it is set to false, it will report only
                        individual HSP where two sequence match each other, and thus
                        will not produce alignments with gaps.
  compstats=            Use composition-based statistics.
  align=                Formating for the alignments.
  transltable=          Query Genetic code to use in translation.

 [General]
  help=                 Show this help message and exit.
  resultTypes=          Get available result types for job.
  outfile=              File name for results (default is JobId; for STDOUT).
  outformat=            Result format(s) to retrieve. It accepts comma-separated values.
  params=               List input parameters.
  paramDetail=          Display details for input parameter.
  
  Support/Feedback:
  https://www.ebi.ac.uk/support/"
  # help option
  if("help" %in% names(args)==TRUE){
    cat(usage)
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  else
  # Check arguments given relate to valid input parameters
  if(!length(paramComp)==0){
      cat("Error: Unrecognised parameter(s) detected. Please check your inputs")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()
  } 
  # Get list of parameters 
  if("params" %in% names(args)==TRUE){
    parametersURL <- paste(baseURL, '/parameters', sep="")
    param <- getForm(parametersURL, Accept ='text/plain')
    param <- str_replace_all(param, '[<][^>]*[>]', '\n')
    parameters <- "parameters: "
    cat(parameters, param)  
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  # Get details on specific parameter
  if(!missing(paramDetail)){
    parameter <- paramDetail
    parametersURL <- paste(baseURL, '/parameters', sep="")
    param <- getForm(parametersURL, Accept ='text/plain')
    param <- str_replace_all(param, '[<][^>]*[>]', '\n')
    
    if(grepl(parameter, param)== FALSE){
      stop(paste("parameter",parameter,"not found"))
    }
    if(grepl(parameter, param)== TRUE){
    paramDetURL <- paste(baseURL, '/parameterdetails/', sep="")
    paramDetURL <- paste(paramDetURL, parameter, sep="")
    
    details <- getForm(paramDetURL, Accept = 'text/plain')
    description <- str_replace_all(details, '[<][^>]*[>]', '\n')
    
    cat(description)
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
    }
  }
  # Get Result Types
  if("resultTypes" %in% names(args)==TRUE){
    URL <-  paste(baseURL, '/run', sep="")
    
    
    JobID <- postForm(URL, email= 'test@ebi.ac.uk', 
                      sequence= 'sp:pak4_human', 
                      database= 'uniprotkb_swissprot',
                      stype= 'protein',
                      program= 'blastp')
    
    statusURL <- paste(baseURL, '/status/', sep="")
    statusURL <- paste(statusURL, JobID)
    statusURL <- gsub("\\s+","", statusURL)
    status <- getForm(statusURL, Accept= 'text/plain')
    while(grepl("FINISHED", status) == FALSE){
      statusURL <- paste(baseURL, '/status/', sep="")
      statusURL <- paste(statusURL, JobID)
      statusURL <- gsub("\\s+","", statusURL)
      status <- getForm(statusURL, Accept= 'text/plain')
    }
    if(grepl("FINISHED", status) == TRUE){
      typeURL <- paste(baseURL, '/resulttypes/', sep="")
      typeURL <- paste(typeURL, JobID)
      typeURL <- gsub("\\s+", "", typeURL)
      resultTypes <- getForm(typeURL, Accept= 'text/plain')
      resultTypes <- str_replace_all(resultTypes, '[<][^>]*[>]', '\n')
      cat(resultTypes)
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()
      
    }
  }
  # If arguments are left blank then print usage
  if(is.null(stype)==TRUE && is.null(sequence)==TRUE && is.null(email)==TRUE 
     && is.null(program)==TRUE && is.null(database)==TRUE) {
    
    cat(usage)
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  # Check if required inputs have been entered
  outformats <- c("accs", "ids", "out", "sequence", "xml", "visual-svg", "complete-visual-svg", "ffdp-query-svg",
                  "ffdp-subject-svg")
  if(missing(email)){
    cat("Error: email must be provided")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  if(missing(stype)){
    cat("Error: stype must be provided")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  if(missing(sequence)){
    cat("Error: sequence must be provided")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  if(missing(program)){
    cat("Error: program must be provided")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  if(missing(database)){
    cat("Error: database must be provided")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  if(!missing(outformat)){
    for(format in outformat){
    if(grepl(paste(outformats, collapse = "|"), format)== FALSE){
      cat("Error: outformat invalid. Valid outformats are accs, ids, out, sequence, xml, visual-svg, 
          complete-visual-svg, ffdp-query-svg and ffdp-subject-svg")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()
      }
    }
  }
  # if sequence input is a file then load file
  if(file_test("-f", sequence)== TRUE){
    sequence= read_file(sequence)
  }
  # Check email input relates to a valid email
  if(!missing(email)){
    if(grepl("^.*@.*\\..*$", email)==FALSE){
      cat("Error: Valid email address must be provided")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()
    }
  }
  # Parameter value checker 
  valueCheck <- function(parameter= NULL, ...){
    paramDetURL <- paste(baseURL, '/parameterdetails/', sep="")
    paramDetURL <- paste(paramDetURL, parameter, sep="")
    
    paramdetails <- getForm(paramDetURL, Accept = 'text/plain')
    paramdetails <- str_extract_all(paramdetails, "<value>(-\\w|\\w+)<.value>")
    paramdetails <- unlist(paramdetails)
    paramdetails <- str_remove_all(paramdetails, "<..{4,5}>")
    invisible( list2env(as.list(environment()), parent.frame()) )
    
  } 
  # Check stype input is valid 
  valueCheck(parameter= "stype")
  valueComp <- grepl(stype, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for stype. Check valid inputs using paramDetail= stype")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  # Check program input is valid
  valueCheck(parameter= "program")
  valueComp <- grepl(program, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for program. Check valid inputs using paramDetail= program")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  # Check database input is valid 
  valueCheck(parameter= "database")
  valueComp <- grepl(database, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for database. Check valid inputs using paramDetail= database")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  # Check alignments input is valid - default is 50 
  valueCheck(parameter= "alignments")
  valueComp <- grepl(alignments, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for alignments. Check valid inputs using paramDetail= alignments")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  # Check scores input is valid - default is 50
  valueCheck(parameter= "scores")
  valueComp <- grepl(scores, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for scores. Check valid inputs using paramDetail= scores")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  # Check exp input is valid - default is 10
  valueCheck(parameter= "exp")
  valueComp <- grepl(exp, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for expectation value. Check valid inputs using paramDetail= exp")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  # Check matrix input is valid - default is BLOSUM62
  if(!missing(matrix)){
    valueCheck(parameter= "matrix")
    valueComp <- grepl(matrix, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for matrix. Check valid inputs using paramDetail= matrix")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check dropoff input is valid - default is 0
  if(!missing(dropoff)){
    valueCheck(parameter= "dropoff")
    valueComp <- grepl(dropoff, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for dropoff. Check valid inputs using paramDetail= dropoff")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check gapopen input is valid - default is -1
  if(!missing(gapopen)){
    valueCheck(parameter= "gapopen")
    valueComp <- grepl(gapopen, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for gap open. Check valid inputs using paramDetail= gapopen")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
  }
  # Check gapext input is valid - default is -1
  if(!missing(gapext)){
    valueCheck(parameter= "gapext")
    valueComp <- grepl(gapext, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for gap extend. Check valid inputs using paramDetail= gapext")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
  }
  # Check filter input is valid - default is F
  if(!missing(filter)){
    valueCheck(parameter= "filter")
    valueComp <- grepl(filter, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for filter. Check valid inputs using paramDetail= filter")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check gapalign input is valid - default is true 
  if(!missing(gapalign)){
    valueCheck(parameter= "gapalign")
    valueComp <- grepl(gapalign, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for gapalign. Check valid inputs using paramDetail= gapalign")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check compstats input is valid - default if F
  if(!missing(compstats)){
    valueCheck(parameter= "compstats")
    valueComp <- grepl(compstats, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for composition-based statistics. Check valid inputs using paramDetail= compstats")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check align input is valid - default is 0
  if(!missing(align)){
    valueCheck(parameter= "align")
    valueComp <- grepl(align, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for align. Check valid inputs using paramDetail= align")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check transltable input is valid - default is 1
  if(!missing(transltable)){
    valueCheck(parameter= "transltable")
    valueComp <- grepl(transltable, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for translation table. Check valid inputs using paramDetail= transltable")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check match_scores value if provided
  if(!missing(match_scores)){
    valueCheck(parameter= "match_scores")
    valueComp <- grepl(match_scores, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for match scores. Check valid inputs using paramDetail= match_scores")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }
  # Check task value if provided
  if(!missing(task)){
    valueCheck(parameter= "task")
    valueComp <- grepl(taak, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for task. Check valid inputs using paramDetail= task")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }
  }

  # Check seqrange input is valid - default is START-END
  if(!seqrange== "START-END"){
    if(grepl("^[0-9]+.[0-9]+", seqrange)== FALSE){
      cat("Error: Invalid input for sequence range. Check valid inputs using paramDetail= seqrange")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
  }
  # Check stype and program combination is valid 
  protprograms <- c("blastp", "tblastn")
  dnaprograms <- c("blastx", "blastn", "tblastx")
  
  databaseCollect <- function(sequencetype= NULL){
    parameter <- "database"
    paramDetURL <- paste(baseURL, '/parameterdetails/', sep="")
    paramDetURL <- paste(paramDetURL, parameter, sep="")
    paramdetails <- getForm(paramDetURL, Accept = 'text/plain')
    if(grepl("protein", sequencetype, ignore.case= TRUE)==TRUE){      
    paramdetails <- str_extract_all(paramdetails, "<value>uniprotkb<.*>chembl<.value>")
    paramdetails <- str_extract_all(paramdetails, "<value>(-\\w|\\w+)<.value>")
    paramdetails <- unlist(paramdetails)
    paramdetails <- str_remove_all(paramdetails, "<..{4,5}>")
    paramdetails <- str_remove_all(paramdetails, "protein")
    invisible( list2env(as.list(environment()), parent.frame()) )
    }
    if(grepl("dna", sequencetype, ignore.case= TRUE)==TRUE){
      paramdetails <- str_extract_all(paramdetails, "<value>em_rel<.*>em_ncr_cum_wgs<.value>")
      paramdetails <- str_extract_all(paramdetails, "<value>(-\\w|\\w+)<.value>")
      paramdetails <- unlist(paramdetails)
      paramdetails <- str_remove_all(paramdetails, "<..{4,5}>")
      paramdetails <- str_remove_all(paramdetails, ".*nucleotide")
      invisible( list2env(as.list(environment()), parent.frame()) )
    }
  }


  if(grepl("protein", stype , ignore.case= TRUE)== TRUE){
    if(grepl(paste(protprograms, collapse = "|"), program, ignore.case= TRUE)== FALSE){
      cat("Error: If stype is protein, program must be blastp or tblastn")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
    if(grepl("blastp", program, ignore.case= TRUE)== TRUE){
      databaseCollect(sequencetype= "protein")
      valueComp <- grepl(database, paramdetails)
      valueComp <- grep("TRUE", valueComp)
      if(length(valueComp)==0){
        cat("Error: If program is blastp, please provide a valid protein database")
        opt <- options(show.error.messages=FALSE) 
        on.exit(options(opt)) 
        stop()  
      }
    }
    if(grepl("tblastn", program, ignore.case= TRUE)==TRUE){
      databaseCollect(sequencetype= "dna")
      valueComp <- grepl(database, paramdetails)
      valueComp <- grep("TRUE", valueComp)
      if(length(valueComp)==0){
        cat("Error: If program is tblastn, please provide a valid nucleotide database")
        opt <- options(show.error.messages=FALSE) 
        on.exit(options(opt)) 
        stop()  
      }
    }
  }
  if(grepl("dna", stype, ignore.case= TRUE)== TRUE){
    if(grepl(paste(dnaprograms, collapse = "|"), program, ignore.case= TRUE)== FALSE){
      cat("Error: If stype is DNA, program must be blastx, blastn or tblastx")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
    if(grepl("blastx", program, ignore.case= TRUE)== TRUE){
      databaseCollect(sequencetype= "protein")
      valueComp <- grepl(database, paramdetails)
      valueComp <- grep("TRUE", valueComp)
      if(length(valueComp)==0){
        cat("Error: If program is blastx, please provide a valid protein database")
        opt <- options(show.error.messages=FALSE) 
        on.exit(options(opt)) 
        stop()  
      }
    }
    newdnaprograms <- c("blastn", "tblastx")
    if(grepl(paste(newdnaprograms, collapse= "|"), program, ignore.case= TRUE)== TRUE){
      databaseCollect(sequencetype= "dna")
      valueComp <- grepl(database, paramdetails)
      valueComp <- grep("TRUE", valueComp)
      if(length(valueComp)==0){
        cat("Error: If program is blastn or tblastx, please provide a valid nucleotide database")
        opt <- options(show.error.messages=FALSE) 
        on.exit(options(opt)) 
        stop()  
      }
    }
    
  }
  # Submit Job 
  URL <-  paste(baseURL, '/run', sep="")
  JobID <- postForm(URL, email= email, 
                    sequence= sequence, 
                    database= database,
                    stype= stype,
                    program= program,
                    alignments= alignments,
                    matrix= matrix,
                    scores= scores,
                    exp= exp,
                    dropoff= dropoff,
                    gapopen= gapopen,
                    gapext= gapext,
                    filter= filter,
                    seqrange= seqrange,
                    gapalign= gapalign,
                    compstats= compstats,
                    align= align,
                    transltable= transltable,
                    task= task,
                    match_score= match_scores)
  JobStatus <- "JOB SUBMITTED \n"
  cat(JobStatus)
  cat("JOBID:",JobID, "\n")
  # Job Status
  statusURL <- paste(baseURL, '/status/', sep="")
  statusURL <- paste(statusURL, JobID)
  statusURL <- gsub("\\s+","", statusURL)
  status <- getForm(statusURL, Accept= 'text/plain')
  while(grepl("FINISHED", status) == FALSE){
    statusURL <- paste(baseURL, '/status/', sep="")
    statusURL <- paste(statusURL, JobID)
    statusURL <- gsub("\\s+","", statusURL)
    status <- getForm(statusURL, Accept= 'text/plain')
    JobStatus <- "RUNNING \n"
    cat(JobStatus)
  }
  if(grepl("FINISHED", status) == TRUE){
    JobStatus <- "FINISHED \n"
    cat(JobStatus)
    JobStatus <- "CREATING OUTPUT \n"
    cat(JobStatus)
  # outformat for results - if outformat specified  
    if(!missing(outformat)){
      for(format in outformat){
      resultURL <- paste(baseURL, '/result/', sep="")
      resultURL <- paste(resultURL, JobID)
      resultURL <- paste(resultURL, '/', sep="")
      resultURL <- paste(resultURL, format)
      resultURL <- gsub("\\s+","", resultURL)
      
      results <- getForm(resultURL, Accept= 'text/plain')
      
      if(missing(outfile)){
        name <- JobID
      }
      if(!missing(outfile)){
        name <- outfile

      }
      if(grepl("xml", format)==TRUE){
        sink(paste(name,".", format,".xml", sep=""), append=FALSE)
        cat(results)
        sink()
        output <- paste(name,".", format,".xml\n", sep="")
        cat(output)
        
      }
      txts <- c("accs", "ids", "out", "sequence")
      if(grepl(paste(txts, collapse = "|"), format)==TRUE){
        sink(paste(name,".", format,".txt", sep=""), append=FALSE)
        cat(results)
        sink()
        output <- paste(name,".", format,".txt\n", sep="")
        cat(output)
      }
      svgs <- c("visual-svg", "complete-visual-svg", "ffdp-query-svg", "ffdp-subject-svg")
      if(grepl(paste(svgs, collapse = "|"), format)==TRUE){
        sink(paste(name,".", format,".svg", sep=""), append=FALSE)
        cat(results)
        sink()
        output <- paste(name,".", format,".svg\n", sep="")
        cat(output)
      }
      }
    }
  # outformat for results - if outformat not specified 
    if(missing(outformat)){
      for(outformat in outformats){
        resultURL <- paste(baseURL, '/result/', sep="")
        resultURL <- paste(resultURL, JobID)
        resultURL <- paste(resultURL, '/', sep="")
        resultURL <- paste(resultURL, outformat)
        resultURL <- gsub("\\s+","", resultURL)
        
        results <- getForm(resultURL, Accept= 'text/plain')
        if(missing(outfile)){
          name <- JobID
        }
        
        if(!missing(outfile)){
          name <- outfile
        }
        if(grepl("xml", outformat)==TRUE){
          sink(paste(name,".", outformat,".xml", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".xml\n", sep="")
          cat(output)
          
        }
        txts <- c("accs", "ids", "out", "sequence")
        if(grepl(paste(txts, collapse = "|"), outformat)==TRUE){
          sink(paste(name,".", outformat,".txt", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".txt\n", sep="")
          cat(output)
        }
        svgs <- c("visual-svg", "complete-visual-svg","ffdp-query-svg", "ffdp-subject-svg")
        if(grepl(paste(svgs, collapse = "|"), outformat)==TRUE){
          sink(paste(name,".", outformat,".svg", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".svg\n", sep="")
          cat(output)
        }
      }     
    }       
    
}
  JobStatus <- "OUTPUT CREATED \n"
  cat(JobStatus)
}

