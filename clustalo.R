# Clustal Omega R Client 

# Load Required Libaries 
library(RCurl)
library(readr)
library(stringr)

clustalo <- function(email= NULL, 
                     stype= NULL,
                     sequence= NULL,
                     outfile= NULL,
                     outformat= NULL,
                     paramDetail= NULL,
                     guidetreeout= NULL,
                     dismatout= NULL,
                     dealign= NULL,
                     mbed= NULL,
                     mbediteration= NULL,
                     iterations= NULL,
                     gtiterations= NULL,
                     hmmiterations= NULL,
                     outfmt= NULL,
                     order= NULL,
                     ...) {
  
  # Get list of parameters
  args <- c(as.list(environment()),list(...))
  parameterlist <- c("^sequence$", "^email$", "^stype$", "^outfile$", "^outformat$", "^paramDetail$", "^help$", 
                     "^resultTypes$", "^params$", "^guidetreeout$", "^dismatout$", "^dealign$", "^mbed$",
                     "^mbediteration$", "^iterations$", "^gtiterations$", "^hmmiterations$", "^outfmt$",
                     "^order$")
  results <- grepl(paste(parameterlist, collapse = "|"), names(args))
  results <- grep("FALSE", results)
  baseURL <- "https://www.ebi.ac.uk/Tools/services/rest/clustalo"  
  usage <- "EMBL-EBI Clustal Omega Perl Client:

 Multiple sequence alignment with Clustal Omega.
  
 [Required (for job submission)]
  email=               E-mail address.
  stype=               Defines the type of the sequences to be aligned.
  sequence=            Three or more sequences to be aligned can be entered
                       directly into this box. Sequences can be in GCG, FASTA, EMBL
                       (Nucleotide only), GenBank, PIR, NBRF, PHYLIP or
                       UniProtKB/Swiss-Prot (Protein only) format. Partially
                       formatted sequences are not accepted. Adding a return to the
                       end of the sequence may help certain applications understand
                       the input. Note that directly using data from word
                       processors may yield unpredictable results as hidden/control
                       characters may be present. There is currently a sequence
                       input limit of 4000 sequences and 4MB of data.
 [Optional]
  guidetreeout=        Output guide tree.
  dismatout=           Output distance matrix. This is only calculated if the mBed-
                       like clustering guide tree is set to false.
  dealign=             Remove any existing alignment (gaps) from input sequences.
  mbed=                This option uses a sample of the input sequences and then
                       represents all sequences as vectors to these sequences,
                       enabling much more rapid generation of the guide tree,
                       especially when the number of sequences is large.
  mbediteration=       Use mBed-like clustering during subsequent iterations.
  iterations=          Number of (combined guide-tree/HMM) iterations.
  gtiterations=        Having set the number of combined iterations, this parameter
                       can be changed to limit the number of guide tree iterations
                       within the combined iterations.
  hmmiterations=       Having set the number of combined iterations, this parameter
                       can be changed to limit the number of HMM iterations within
                       the combined iterations.
  outfmt=              Format for generated multiple sequence alignment.
  order=               The order in which the sequences appear in the final
                       alignment.

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
    
    JobID <- postForm(URL, email= 'test@ebi.ac.uk',
                      sequence= 'sp:pak4_human,sp:pak2_human,sp:pak5_human',
                      stype= 'protein')
    
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
  if(is.null(stype)==TRUE && is.null(sequence)==TRUE && is.null(email)==TRUE) {
    
    cat(usage)
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()
  }
  # Check if required inputs have been entered
  outformats <- c("out", "sequence", "aln-clustal_num", "tree", "phylotree", "pim")  
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
  if(!missing(outformat)){
    for(format in outformat){
      if(grepl(paste(outformats, collapse = "|"), format)== FALSE){
        cat("Error: outformat invalid. Valid outformats are out, sequence, aln-clustal_num, tree, phylotree and pim")
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
  # Check guidetreeout input is valid - default is true 
  if(!missing(guidetreeout)){
  valueCheck(parameter= "guidetreeout")
  valueComp <- grepl(guidetreeout, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for output guide tree. Check valid inputs using paramDetail= guidetreeout")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
    }
  }  
  # Check dismatout input is valid - default is true 
  if(!missing(dismatout)){
  valueCheck(parameter= "dismatout")
  valueComp <- grepl(dismatout, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for output distance matrix. Check valid inputs using paramDetail= dismatout")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check dealign input is valid - default is true 
  if(!missing(dealign)){
  valueCheck(parameter= "dealign")
  valueComp <- grepl(dealign, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for dealign. Check valid inputs using paramDetail= dealign")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop() 
  }
  }  
  # Check mbed input is valid - default is true 
  if(!missing(mbed)){
  valueCheck(parameter= "mbed")
  valueComp <- grepl(mbed, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for mbed. Check valid inputs using paramDetail= mbed")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }   
  # Check mbediteration input is valid - default is true 
  if(!missing(mbediteration)){
  valueCheck(parameter= "mbediteration")
  valueComp <- grepl(mbediteration, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for mbediteration. Check valid inputs using paramDetail= mbediteration")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  } 
  # Check iterations input is valid - default is 0
  if(!missing(iterations)){
  valueCheck(parameter= "iterations")
  valueComp <- grepl(iterations, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for iterations. Check valid inputs using paramDetail= iterations")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  } 
  }
  # Check gtiterations input is valid - default is -1
  if(!missing(gtiterations)){
  valueCheck(parameter= "gtiterations")
  valueComp <- grepl(gtiterations, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for guide tree iterations. Check valid inputs using paramDetail= gtiterations")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  } 
  }
  # Check hmmiterations input is valid - default is -1
  if(!missing(hmmiterations)){
  valueCheck(parameter= "hmmiterations")
  valueComp <- grepl(hmmiterations, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for HMM iterations. Check valid inputs using paramDetail= hmmiterations")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }   
  # Check outfmt input is valid - default is clustal_num
  if(!missing(outfmt)){
  valueCheck(parameter= "outfmt")
  valueComp <- grepl(outfmt, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for output alignment format. Check valid inputs using paramDetail= outfmt")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Check order input is valid - default is aligned
  if(!missing(order)){
  valueCheck(parameter= "order")
  valueComp <- grepl(order, paramdetails, ignore.case= TRUE)
  valueComp <- grep("TRUE", valueComp)
  if(length(valueComp)==0){
    cat("Error: Invalid input for order. Check valid inputs using paramDetail= order")
    opt <- options(show.error.messages=FALSE) 
    on.exit(options(opt)) 
    stop()  
  }
  }
  # Submit Job
  URL <-  paste(baseURL, '/run', sep="")
  JobID <- postForm(URL, email= email, 
                    sequence= sequence, 
                    stype= stype,
                    guidetreeout= guidetreeout,
                    dismatout= dismatout,
                    dealign= dealign,
                    mbed= mbed,
                    mbediteration= mbediteration,
                    iterations= iterations,
                    gtiterations= gtiterations,
                    hmmiterations= hmmiterations,
                    outfmt= outfmt,
                    order= order)
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
        if(grepl("aln-clustal_num", format)==TRUE){
          sink(paste(name,".", format,".clustal_num", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".clustal_num\n", sep="")
          cat(output)
        }
        if(grepl("^tree$", format)==TRUE){
          sink(paste(name,".", format,".dnd", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".dnd\n", sep="")
          cat(output)
        }
        if(grepl("phylotree", format)==TRUE){
          sink(paste(name,".", format,".ph", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".ph\n", sep="")
          cat(output)
        }    
        if(grepl("pim", format)==TRUE){
          sink(paste(name,".", format,".pim", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", format,".pim\n", sep="")
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
        if(grepl("aln-clustal_num", outformat)==TRUE){
          sink(paste(name,".", outformat,".clustal_num", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".clustal_num\n", sep="")
          cat(output)
        }
        if(grepl("^tree$", outformat)==TRUE){
          sink(paste(name,".", outformat,".dnd", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".dnd\n", sep="")
          cat(output)
        }
        if(grepl("phylotree", outformat)==TRUE){
          sink(paste(name,".", outformat,".ph", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".ph\n", sep="")
          cat(output)
        }    
        if(grepl("pim", outformat)==TRUE){
          sink(paste(name,".", outformat,".pim", sep=""), append=FALSE)
          cat(results)
          sink()
          output <- paste(name,".", outformat,".pim\n", sep="")
          cat(output)
        } 


      }     
    }       
  }
  JobStatus <- "OUTPUT CREATED \n"
  cat(JobStatus)  
  
}
