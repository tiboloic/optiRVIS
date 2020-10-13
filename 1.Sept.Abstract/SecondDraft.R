# LT 3/09/2019

# use LIMBR exons to count variants per mutation class

# import BED file

exons = rtracklayer::import.bed("d:/gnomad/intolerance/limbr.ex.bed")

library(rtracklayer)
file = BEDFile("d:/gnomad/intolerance/limbr.ex.bed")
import(file)

# does not work

# import from csv file

limbr = read.csv("D:/gnomad/intolerance/LIMBR.exon.csv")
exons = GRanges(limbr$position, ID=limbr$sub_region, gene= limbr[,1])
