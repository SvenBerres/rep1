library(EpigeneticsCSAMA)
dataDirectory =  system.file("bedfiles", package="EpigeneticsCSAMA")


library(GenomicRanges)
library(rtracklayer)
library(IRanges)

input = import.bed(file.path(dataDirectory, 'ES_input_filtered_ucsc_chr6.bed'))
rep1 = import.bed(file.path(dataDirectory, 'H3K27ac_rep1_filtered_ucsc_chr6.bed'))
rep2 = import.bed(file.path(dataDirectory, 'H3K27ac_rep2_filtered_ucsc_chr6.bed'))
length(input)
length(rep1)
length(rep2)

library(chipseq)

prepareChIPseq = function(reads){
  frag.len = median( estimate.mean.fraglen(reads) )
  cat( paste0( 'Median fragment size for this library is ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  return( trim(reads.extended) )
}

input = prepareChIPseq( input )
rep1 = prepareChIPseq( rep1 )
rep2 = prepareChIPseq( rep2 )

data(si)
si
binsize = 200
bins = tileGenome(si['chr6'], tilewidth=binsize,
                  cut.last.tile.in.chrom=TRUE)
bins

BinChIPseq = function( reads, bins ){
  
  mcols(bins)$score = countOverlaps( bins, reads ) 
  return( bins ) 
}
input.200bins = BinChIPseq( input, bins )
rep1.200bins = BinChIPseq( rep1, bins )
rep2.200bins = BinChIPseq( rep2, bins )

rep1.200bins

plot( 200000:201000, rep1.200bins$score[200000:201000], 
      xlab="chr6", ylab="counts per bin", type="l")
export(input.200bins, 
       con='input_chr6.bedGraph',
       format = "bedGraph")
export(rep1.200bins, 
       con='H3K27ac_rep1_chr6.bedGraph',
       format = "bedGraph")
export(rep2.200bins, 
       con='H3K27ac_rep2_chr6.bedGraph',
       format = "bedGraph")

library(Gviz)
data(bm)
bm
AT = GenomeAxisTrack( )
plotTracks(c( bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto", 
           cex.title=1, fontsize=10 )
input.track = DataTrack(input.200bins, 
                        strand="*", genome="mm9", col.histogram='gray',
                        fill.histogram='black', name="Input", col.axis="black",
                        cex.axis=0.4, ylim=c(0,150))
rep1.track = DataTrack(rep1.200bins, 
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 1", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))

rep2.track = DataTrack(rep2.200bins, 
                       strand="*", genome="mm9", col.histogram='steelblue',
                       fill.histogram='black', name="Rep. 2", col.axis='steelblue',
                       cex.axis=0.4, ylim=c(0,150))

plotTracks(c(input.track, rep1.track, rep2.track, bm, AT),
           from=122530000, to=122900000,
           transcriptAnnotation="symbol", window="auto", 
           type="histogram", cex.title=0.7, fontsize=10 )

#peak finding 
peaks.rep1 = import.bed(file.path(dataDirectory,'Rep1_peaks_ucsc_chr6.bed'))
peaks.rep2 = import.bed(file.path(dataDirectory,'Rep2_peaks_ucsc_chr6.bed'))

peaks1.track = AnnotationTrack(peaks.rep1, 
                               genome="mm9", name='Peaks Rep. 1',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)
peaks2.track = AnnotationTrack(peaks.rep2, 
                               genome="mm9", name='Peaks Rep. 2',
                               chromosome='chr6',
                               shape='box',fill='blue3',size=2)

plotTracks(c(input.track, rep1.track, peaks1.track,
             rep2.track, peaks2.track, bm, AT),
           from=122630000, to=122700000,
           transcriptAnnotation="symbol", window="auto", 
           type="histogram", cex.title=0.7, fontsize=10 )

ovlp = findOverlaps( peaks.rep1, peaks.rep2 )
ovlp