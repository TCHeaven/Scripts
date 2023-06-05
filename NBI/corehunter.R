options(java.parameters = "-Xmx50G")
library(corehunter)
options(java.parameters = "-Xmx50G")

args <- commandArgs(trailingOnly = TRUE)
OutFile <- args[1]
BASE <- args[2]
BASE <- as.numeric(BASE)

print(OutFile)
print(BASE)

if (BASE > 0) {
	print("Running with set sample number")
	set.seed(1234)
	dist <- distances(file = "p_dis.csv")
	core <- sampleCore(dist, size = BASE)
	write.table(as.list(core$sel), file = OutFile, row.names=FALSE, col.names=FALSE, sep=",") # output

} else if (BASE == 0) {
	print("Running with defaults")
	set.seed(1234)
	dist <- distances(file = "p_dis.csv")
	core <- sampleCore(dist)
	write.table(as.list(core$sel), file = OutFile, row.names=FALSE, col.names=FALSE, sep=",") # output

} else if (BASE < 0) {
	print("error sample number less than zero")
}
