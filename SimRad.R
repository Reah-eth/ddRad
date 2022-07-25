###################################################################################################################################
#### Title: Fragment simulation with "SimRAD"                                                                                 #####      
#### Description: This script can be used to similate the number of fragments to be expected when the ddRad protocal is used. #####
####              Here the genome assemply of Kyus is used for (tape station) peak size of 500 to 610                         #####
#### Author: Reah Gonzales                                                                                                    #####
#### Date: 25-07-2022                                                                                                         #####
#### contact: reah.gonzales@usys.ethz.ch                                                                                      #####
#################################################################################################################################:)

### Since this code in usually implemented on an external server check to ensure you have all the package dependancies for the SimRad package 
### If you do not have all the packages, Start from here 
.libPaths()
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ShortRead")
install.packages("SimRAD")

### If all the packages and the "SimRAD" package is already downloaded start here
library(SimRAD)
### Example3: a double digestion (ddRAD)

### the reference assembly should be downloaded if not used the "sim.DNAseq" 

seq <- ref.DNAseq("ddRad/Kyus.fa",prop.contigs=0.2)


#Restriction Enzyme 1
#TaqI
P1_5 <- "T"
P1_3 <- "CGA"

#Restriction Enzyme 2
#EcoRI NEBU#
P2_5 <- "G"
P2_3 <- "AATTC"

simseq.dig <- insilico.digest(seq, P1_5, P1_3, P2_5 , P2_3, verbose=TRUE) ## the larger the prop.contig the longer this process requires 

#### Result when prop.contig=0.1
## Number of restriction sites for the first enzyme: 53145
## Number of restriction sites for the second enzyme: 2994

#### Result when prop.contig=0.2
### Number of restriction sites for the first enzyme: 2231748
### Number of restriction sites for the second enzyme: 141088


simseq.sel <- adapt.select(simseq.dig, type="AB+BA", P1_5, P1_3, P2_5 , P2_3)


# wide size selection (400-610):
wid.simseq <- size.select(simseq.sel,  min.size = 500, max.size = 610, graph=F, verbose=TRUE)

##### result when prop.contig=0.1
## 683 fragments between 400 and 610 bp

#### Result when prop.contig=0.2
## 15539 fragments between 500 and 610 
### 15539 (fragments between 400 and 610 bp at 20% of halploid genome used) * 10 = 155,390 fragments

# narrow size selection (210-260):
nar.simseq <- size.select(simseq.sel,  min.size = 400, max.size = 610, graph=F, verbose=TRUE)

##### Result when prop.contig=0.1 
### 683 fragments between 400 and 610 bp

##### Result when prop.contig=0.2
### 15539 fragments between 500 and 610 bp

### 15539 (fragments between  500 and 610 bp at 20% of halploid genome used)* 10= 155,390 fragments 


#### the resulting fragment characteristics can be further examined:
boxplot(list(width(simseq.sel), width(wid.simseq), width(nar.simseq)), names=c("All fragments",
                                                                               "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")

