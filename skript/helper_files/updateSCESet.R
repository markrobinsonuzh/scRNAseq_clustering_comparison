# Update SCEset data files to SingleCellExperiment and save it 

# function to save SCeset as class SingleCellExperiment under original name
# change the orginial names in import files from res to sceset.....

update_SCEset <- function(files) {
  
  require(scater)
  
data <- vector("list", length(files))
labels <- data
names(data) <- names(labels) <- names(files)
for (i in names(files)) {
  f <- files[[i]]
  load(f)
  print( names(files) )
  stopifnot( class( sceset ) == "SCESet" )
  
  sceset <- updateSCESet(sceset)
  save(sceset, file=files[[i]] )
  }

}

######  source files and update
source("FILES.R")
update_SCEset(files)

source("FILESraw.R")
update_SCEset(files[5])

