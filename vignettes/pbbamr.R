## ---- echo=TRUE, results='asis'------------------------------------------
library(pbbamr)
# For this example, we will load some sample data
# distributed with the package.
bamname = system.file("extdata", "bam_mapping_1.bam", package="pbbamr")
# Load the index table.
ind = loadpbi(bamname)


# Show the first 3 rows (transposed, note kable just formats)
knitr::kable(t(ind[1:3,]))

## ---- echo=TRUE, results='asis'------------------------------------------
library(ggplot2)
# We can calculate the length of the alignment on the reference as 
# the template end minus the template start of that alignment.
ind$alnLength = ind$tend - ind$tstart
ggplot(ind, aes(x=alnLength)) + geom_density(fill="blue") + theme_bw() +
  labs(x="Alignment Length on Reference")

## ---- echo=TRUE, results='asis'------------------------------------------
ind$errorRate = (ind$mismatches + ind$inserts + ind$dels) / ind$alnLength 
ggplot(ind, aes(x=errorRate)) + geom_density(fill="cyan") + theme_bw() +
  labs(x="Error Rate")

## ---- echo=TRUE, results='hold'------------------------------------------
# Get the name of the sample indexed FASTA file
fastaname = system.file("extdata", "lambdaNEB.fa", package="pbbamr")


loadHMMfromBAM(ind$offset, bamname, fastaname)

# let's just grab and plot one alignment.
allAlns = loadDataAtOffsets(ind$offset[1], bamname, fastaname)
# The alignments are returned as a list of data frames
# in this case let's just look at the first
aln = allAlns[[1]]
head(aln[1:6,])

## ---- echo=TRUE, collapse=TRUE, results='hold'---------------------------
# Get all the alignments for all records in the index
allAlns = loadDataAtOffsets(ind$offset, bamname, fastaname)
# let's combine the individual alignments into one big data frame
alns = do.call(rbind, allAlns)
# Now count the number of insertions
insert_cnt = sum(alns$ref=="-")
# Do they match the pbi file?
insert_cnt == sum(ind$inserts)

## ---- echo=TRUE, results='asis'------------------------------------------
# First let's write a function to find single deletion events, and report 
# the deleted base and following base.
findSingleDeletions <- function(aln) {
  # Where are the deletions of size 1?
  # First we will convert it to run length encoding (gap = 5)
  rls = rle(as.integer(aln$read))
  # Now find out what has a deletion of size 1
  del_size1 = which(rls$values == 5 & rls$lengths == 1)
  # and convert back from run length to standard coordinates
  org_positions = cumsum(rls$lengths)
  single_del_locs = org_positions[del_size1]
  # Presumably the alignment should not end with 
  # a deletion, but let's be sure
  single_del_locs = single_del_locs[single_del_locs < nrow(aln)]
  # And now let's grab the deleted base and the one after
  data.frame(AfterDeletion = aln$ref[single_del_locs + 1], 
             Deleted = aln$ref[single_del_locs])
}
# Apply our function to investigate deletions to all alignments
deletions = lapply(allAlns, findSingleDeletions)
# Combine the results
dels = do.call(rbind, deletions)
# And examine the table
knitr::kable(table(dels$AfterDeletion, dels$Deleted))

# Tables can also be visualized as heatmaps
dels$Count = 1
agg = aggregate(Count ~ AfterDeletion + Deleted, dels, length )
ggplot(agg, aes(x=AfterDeletion, y=Deleted, fill=Count)) + geom_tile() +
  theme_bw(base_size=10)  + 
  labs(title="Anyone think G/C homopolymers\nwill be hard in consensus?",
       x="Base After Deletion", y="Deleted Base")

