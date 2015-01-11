data_all = read.table("kmer_histo.all.txt")
data_norm = read.table("kmer_histo.NormMaxKCov50.txt")

log_data_all = cbind(data_all[,1], log(data_all[,2]+1))
log_data_norm = cbind(data_norm[,1], log(data_norm[,2]+1))

plot(log_data_norm, col='green', xlim=c(0,200), xlab="kmer occurrence count", ylab="log(number of unique kmers)", main="Kmer composition and \nin silico read normalization")
points(log_data_all, col='red')

