.libPaths()
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ShortRead")
install.packages("SimRAD")
### Example3: a double digestion (ddRAD)

#### If you have already a reference assembly 
simseq <- ref.DNAseq("ddRad/Kyus.fa",prop.contigs=0.2)


#Restriction Enzyme 1
#TaqI
P1_5 <- "T"
P1_3 <- "CGA"

#Restriction Enzyme 2
#EcoRI NEBU#
P2_5 <- "G"
P2_3 <- "AATTC"

simseq.dig <- insilico.digest(simseq, P1_5, P1_3, P2_5 , P2_3, verbose=TRUE)

#### Result 
## Number of restriction sites for the first enzyme: 53145
## Number of restriction sites for the second enzyme: 2994

simseq.sel <- adapt.select(simseq.dig, type="AB+BA", P1_5, P1_3, P2_5 , P2_3)


# wide size selection (200-270):
wid.simseq <- size.select(simseq.sel,  min.size = 400, max.size = 610, graph=F, verbose=TRUE)

##### result 
## 683 fragments between 400 and 610 bp


# narrow size selection (210-260):
nar.simseq <- size.select(simseq.sel,  min.size = 400, max.size = 610, graph=F, verbose=TRUE)

##### Result 
### 683 fragments between 400 and 610 bp


#### the resulting fragment characteristics can be further examined:
boxplot(list(width(simseq.sel), width(wid.simseq), width(nar.simseq)), names=c("All fragments",
                                                                               "Wide size selection", "Narrow size selection"), ylab="Locus size (bp)")

