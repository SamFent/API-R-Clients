# Emboss Transeq R Client

# Load Required Libraries 
library(RCurl)
library(readr)
library(stringr)

emboss_transeq <- function(email= NULL,
                           sequence= NULL,
                           outfile= NULL,
                           outformat= NULL,
                           paramDetail= NULL,
                           frame= NULL,
                           codontable= NULL,
                           regions= "START-END",
                           trim= NULL,
                           reverse= NULL,
                           ...){
  # Get list of parameters
  args <- c(as.list(environment()),list(...))
  parameterlist <- c("^sequence$", "^email$", "^outformat$", "^outfile$", "^paramDetail$", 
                     "^help$", "^resultTypes$", "^params$", "^frame$", "^codontable$", "^regions$",
                     "^trim$", "^reverse$")
  results <- grepl(paste(parameterlist, collapse = "|"), names(args))
  results <- grep("FALSE", results)
  baseURL <- "https://www.ebi.ac.uk/Tools/services/rest/emboss_transeq"
  usage <- "EMBL-EBI EMBOSS transeq R Client:

 Sequence translations with transeq.
  
 [Required (for job submission)]
  email=               E-mail address.
  sequence=            Any input formats accepted by EMBOSS can be used, the full
                       list of sequence formats accepted as input by EMBOSS tools
                       can be accessed via the link below. Word processor files may
                       yield unpredictable results as hidden/control characters may
                       be present in the files. It is best to save files with the
                       Unix format option to avoid hidden Windows characters.

 [Optional]
  frame=               The frames to be translated. The order of the frames follows
                       the Staden convention: Frame -1 is the reverse-complement of
                       the sequence having the same codon phase as frame 1. Frame
                       -2 is the same phase as frame 2. Frame -3 is the same phase
                       as frame 3.
  codontable=          Which genetic code table to use. These are kept synchronised
                       with those maintained at the NCBI's Taxonomy Browser.
  regions=             Which regions of the user's DNA molecule are to be
                       translated.
  trim=                Remove '*' and 'X' (stop and ambiguity) symbols from the end
                       of the translation.
  reverse=             Choose this option if you wish to reverse and complement
                       your input sequence before frame translation.
  
 [General]
  help=                Show this help message and exit.
  resultTypes=         Get available result types for job.
  outfile=             File name for results (default is JobId; for STDOUT).
  outformat=           Result format(s) to retrieve. It accepts comma-separated values.
  params=              List input parameters.
  paramDetail=         Display details for input parameter.
  
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
    if(!length(results)==0){
      cat("Error: Unrecognised parameter(s) detected. Please check your inputs")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()
    } 
  # Get Result Types
  if("resultTypes" %in% names(args)==TRUE){
    URL <-  paste(baseURL, '/run', sep="")
    
    testsequence <- "./Sequence/dna_sequence.fasta"
    testsequence= read_file(testsequence)
    
    JobID <- postForm(URL, email= 'test@ebi.ac.uk', 
                      sequence= testsequence)
    
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
  # Get details on specific parameters
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
  # If arguments if left blank then print usage
  if(is.null(sequence)==TRUE && is.null(email)==TRUE) {
    
    cat(usage)
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  # Check if required inputs have been entered
  outformats <- c("out", "sequence")
  if(missing(email)){
    cat("Error: email must be provided")
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
  if(!missing(outformat)){
    for(format in outformat){
      if(grepl(paste(outformats, collapse = "|"), format)== FALSE){
      cat("Error: outformat invalid. Valid outformats are out and sequence")
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
  # Check frame input is valid
  if(!missing(frame)){
  valueCheck(parameter= "frame")
  valueComp <- grepl(frame, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for frame. Check valid inputs using paramDetail= frame")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check codontable input is valid
  if(!missing(codontable)){
  valueCheck(parameter= "codontable")
  valueComp <- grepl(codontable, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for codontable. Check valid inputs using paramDetail= codontable")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check trim input is valid
  if(!missing(trim)){
  valueCheck(parameter= "trim")
  valueComp <- grepl(trim, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for trim. Check valid inputs using paramDetail= trim")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check reverse input is valid
  if(!missing(reverse)){
  valueCheck(parameter= "reverse")
  valueComp <- grepl(reverse, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for reverse. Check valid inputs using paramDetail= reverse")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check regions input is valid - default is START-END
  if(!regions== "START-END"){
    if(grepl("^[0-9]+.[0-9]+", seqrange)== FALSE){
      cat("Error: Invalid input for regions. Check valid inputs using paramDetail= regions")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop() 
    }
  }
  # Submit Job
  URL <-  paste(baseURL, '/run', sep="")
  JobID <- postForm(URL, email= email, 
                    sequence= sequence,
                    frame= frame,
                    codontable= codontable,
                    regions= regions,
                    trim= trim,
                    reverse= reverse)
  
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
        sink(paste(name,".", format,".txt", sep=""), append=FALSE)
        cat(results)
        sink()
        output <- paste(name,".", format,".txt\n", sep="")
        cat(output)
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
        sink(paste(name,".", outformat,".txt", sep=""), append=FALSE)
        cat(results)
        sink()
        output <- paste(name,".", outformat,".txt\n", sep="")
        cat(output)
      }
      
    }    
  }       
  
  JobStatus <- "OUTPUT CREATED \n"
  cat(JobStatus)  
   
}
