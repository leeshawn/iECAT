 \name{iECAT}
 \alias{iECAT}
 \alias{iECAT.SSD.OneSet_SetIndex}
 \title{Integrating External Controls to Association Tests}
 \description{
     Test for association of a set of variants with integrating external study samples to improve power  
 }
 \usage{


iECAT(Z, obj, tbl.external.all, weights.beta=c(1,25), weights = NULL, 
r.corr=0, method="davies", missing_cutoff=0.15, MAC.lowlimit=3, 
MAF.limit=1)

iECAT.SSD.OneSet_SetIndex(EC.INFO, SetIndex,  obj, ..., obj.SNPWeight=NULL)

 }
\arguments{
      \item{Z}{a numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. 
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing. A does not need to be a major allele.}
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{tbl.external.all}{a p x 2 matrix of allele count of the matrix. Each row should be matched with the each column in Z.}
       \item{weights.beta}{a numeric vector of parameters for the beta weights for the weighted kernels. 
      If you want to use your own  weights, please use the ``weights'' parameter. It will be ignored if ``weights'' parameter is not null.}
        \item{weights}{a numeric vector of weights for the weighted kernels. 
      It is \eqn{\sqrt{w}} in the SKAT paper. 
      So if you want to use the Madsen and Browning (2009) weight, you should set each element of weights as \eqn{1/ \sqrt{p(1-p)}}, 
      not \eqn{1/ p(1-p)}. When it is NULL, the beta weight with the ``weights.beta'' parameter is used. }    
    \item{r.corr}{the \eqn{\rho} parameter for the compound symmetric correlation structure kernels (default= 0). 
      If you give a vector value, SKAT will conduct the optimal test. It will be ignored if method=``optimal'' or method=``optimal.adj''. See details.}
      \item{method}{a method to compute the p-value (default= "davies"). 
      "davies" represents an exact method that  computes the p-value by inverting the characteristic function of the mixture chisq, 
      "liu" represents an approximation method that matches the first 3 moments, 
      "liu.mod" represents modified "liu" method that matches kurtosis instead of skewness 
      to improve tail probability approximation, "optimal" represents a SKAT-O based on an unified approach.}
        \item{missing_cutoff}{a cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.}
      \item{MAC.lowlimit}{a cutoff for MAC low limit (default=3). Variants with internal study MAC <= MAC.lowlimit will be excluded from the analysis.}
      \item{MAF.limit}{a cutoff for MAF upper limit (default=1). Variants with internal study MAF > MAF.limit will be excluded from the analysis.}
      \item{EC.INFO}{an EC_INFO object returned from Open_SSD_wEC. }
      \item{SetIndex}{a numeric value of Set index. A set index of each set can be found from SetInfo object in EC.INFO$SSD.INFO.}
      \item{\dots}{further arguments to be passed to ``iECAT'' }
      \item{obj.SNPWeight}{ an output object of Read_SNP_WeightFile (default=NULL). 
      If NULL, the beta weight with the ``weights.beta'' parameter will be used.  }
        
}
\value{
	\item{p.value}{p-value of iECAT. }
	\item{p.value.noadj}{p-value without the adjustment of possible batch effects}
	\item{p.value.internal}{SKAT/SKAT-O p-values with only internal study samples}
	\item{param}{estimated parameters of each method.} 
	\item{param$n.marker}{a number of variants in the genotype matrix (Z).}  
	\item{param$n.marker.test}{a number of variants used for the test. } 
	
}

\author{Seunggeun (Shawn) Lee}

\examples{
library(SKAT)

data(Example, package="iECAT")
attach(Example)

# iECAT-O
# test the first gene

obj<-SKAT_Null_Model(Y ~ 1, out_type="D")
Z = Z.list[[1]]
tbl.external.all = tbl.external.all.list[[1]]

iECAT(Z, obj, tbl.external.all, method="optimal")


# test for the first 3 genes in the Example dataset
p.value.all<-rep(0,3)
p.value.internal.all<-rep(0,3)
for(i in 1:3){

	re<-iECAT(Z.list[[i]], obj, tbl.external.all.list[[i]], method="optimal")
	p.value.all[i]<-re$p.value
	p.value.internal.all[i]<-re$p.value.internal

}

# iECAT-O p-values
p.value.all

# SKAT-O p-values
p.value.internal.all

}

