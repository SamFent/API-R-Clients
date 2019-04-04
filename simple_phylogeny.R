# Simple Phylogeny R Client 

# Load Required Libaries 
library(RCurl)
library(readr)
library(stringr)

simple_phylogeny  <- function(email= NULL, 
                      sequence= NULL,
                      outfile= NULL,
                      outformat= NULL,
                      paramDetail= NULL,
                      tree= NULL,
                      kimura= NULL,
                      tossgaps= NULL,
                      clustering= NULL,
                      pim= NULL,
                      ...){
  # Get list of parameters
  args <- c(as.list(environment()),list(...))
  parameterlist <- c("^sequence$", "^email$", "^outformat$", "^outfile$", "^paramDetail$", 
                     "^help$", "^resultTypes$", "^params$", "^tree$", "^kimura$", "^tossgaps$", 
                     "^clustering$", "^pim$")
  results <- grepl(paste(parameterlist, collapse = "|"), names(args))
  results <- grep("FALSE", results) 
  baseURL <- "https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny"
  usage <- "EMBL-EBI Simple Phylogeny R Client:

 Generating Phylogenetic Trees with Simple Phylogeny.
  
 [Required (for job submission)]
  email=               E-mail address.
  sequence=            Phylogeny using an alignment directly entered into the input
                       box in a supported format. Alignment formats supported
                       include Clustal, FASTA and MSF. Partially formatted or
                       unaligned sequences are not accepted. Adding a return to the
                       end of the sequence may help the Simple Phylogeny tool
                       understand the input. Note that directly using data from
                       word processors may yield unpredictable results as
                       hidden/control characters may be present. There is currently
                       a limit of 500 sequences and 1MB of data.

 [Optional]
  tree=                Determines the outputs that the Simple Phylogeny tool
                       produces.
  kimura=              Controls whether Simple Phylogeny attempts to correct for
                       multiple substitutions at the same site. This is recommended
                       to be set 'on' for more divergent sequences and has the
                       effect of stretching branch lengths. For very divergent
                       sequences the distances cannot be reliably corrected.
  tossgaps=            With this option enabled columns where any of the sequences
                       in the input have a gap will be excluded, forcing the
                       alignment to use only positions where information can be
                       included from all sequences.
  clustering=          Clustering Methods.
  pim=                 Output the percentage identity matrix.

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
    
    testsequence <- "./Sequence/prot_sequences.aln"
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
  outformats <- c("out", "sequence", "tree")
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
        cat("Error: outformat invalid. Valid outformats are out, sequence and tree")
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
  # Check tree input is valid - default is phylip
  if(!missing(tree)){
  valueCheck(parameter= "tree")
  valueComp <- grepl(tree, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for tree. Check valid inputs using paramDetail= tree")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  } 
  }
  # Check kimura input is valid - default is true
  if(!missing(kimura)){
  valueCheck(parameter= "kimura")
  valueComp <- grepl(kimura, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for kimura. Check valid inputs using paramDetail= kimura")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  } 
  }
  # Check tossgaps input is valid - default is true
  if(!missing(tossgaps)){
  valueCheck(parameter= "tossgaps")
  valueComp <- grepl(tossgaps, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for tossgaps. Check valid inputs using paramDetail= tossgaps")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }   
  }
  # check if clustering input is valid
  if(!missing(clustering)){
    valueCheck(parameter= "clustering")
    valueComp <- grepl(clustering, paramdetails, ignore.case= TRUE)
    valueComp <- grep("TRUE", valueComp)
    if(length(valueComp)==0){
      cat("Error: Invalid input for clustering. Check valid inputs using paramDetail= clustering")
      opt <- options(show.error.messages=FALSE) 
      on.exit(options(opt)) 
      stop()  
    }  
  }
  # Check pim input is valid - default is true
  if(!missing(pim)){
  valueCheck(parameter= "pim")
  valueComp <- grepl(pim, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for pim. Check valid inputs using paramDetail= pim")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  } 
  }
  # Submit Job
  URL <-  paste(baseURL, '/run', sep="")
  JobID <- postForm(URL, email= email, 
                    sequence= sequence,
                    tree= tree,
                    kimura= kimura,
                    tossgaps= tossgaps,
                    clustering= clustering,
                    pim= pim)
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
        txts <- c("out", "sequence")
        if(grepl(paste(txts, collapse = "|"), format)==TRUE){
          sink(paste(name,".", format,".txt", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".txt\n", sep="")
          cat(output)
        }
        if(grepl("tree", format)==TRUE){
          sink(paste(name,".", format,".ph", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".ph\n", sep="")
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
        txts <- c("out", "sequence")
        if(grepl(paste(txts, collapse = "|"), outformat)==TRUE){
          sink(paste(name,".", outformat,".txt", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".txt\n", sep="")
          cat(output)
        }
        if(grepl("tree", outformat)==TRUE){
          sink(paste(name,".", outformat,".ph", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".ph\n", sep="")
          cat(output)
        } 
      }
    }   
  }
  JobStatus <- "OUTPUT CREATED \n"
  cat(JobStatus)    
}
  