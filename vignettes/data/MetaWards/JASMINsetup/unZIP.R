## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 2)
    filedir <- args[1]
    hash <- args[2]
} else {
    stop("No arguments")
}

## print progress
print(paste0("Currently evaluating: ", hash))

## extract files
files <- list.files(paste0("../raw_outputs/", hash))
files <- files[grep(glob2rx("age*.db.bz2"), files)]

for(file in files) {

    ## set path
    path <- paste0("../raw_outputs/", hash, "/", file)
    path <- gsub(".bz2", "", path)

    ## unzip DB
    system(paste0("bzip2 -dkf ", path, ".bz2"))
    
    ## write table out
    system(paste0("mkdir -p ", filedir, "raw_outputs/", hash))
    system(paste0("mv ", path, " ", filedir, "raw_outputs/", hash))
    
    print(paste0("Done: ", file))
}
 
print(paste0(hash, ": Finished"))

