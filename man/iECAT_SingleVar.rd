 \name{iECAT_SingleVar}
 \alias{iECAT_SingleVar}
 \alias{iECAT_SingleVar_Tbl}
 \title{Integrating External Controls to Association Tests}
 \description{
     Test for association of a single variant with integrating external study samples to improve power  
 }
 \usage{

	iECAT_SingleVar(Z, obj, tbl.external, Fisher.test=FALSE, weight.factor=NULL)
	iECAT_SingleVar_Tbl(Tbl, Fisher.test=FALSE, weight.factor=NULL)
 }
\arguments{
      \item{Z}{a numeric genotype vector for the variant.
      Each genotype should be coded as 0, 1, 2, and 9 (or NA) for AA, Aa, aa, and missing. A does not need to be a major allele.}
      \item{obj}{an output object of the SKAT_Null_Model function. }
      \item{tbl.external}{a vector of allele count from the external control. }
      \item{Fisher.test}{an indicator variable (default= "FALSE") to indicate whether to obtain p.value.internal, p.value.noadjust, and p.value.IvsE using Fisher's exact test. If FALSE, 
      asymptotic tests will be used to calculate p-values. }
      \item{weight.factor}{internal use only.}
      \item{Tbl}{a 3 by 2 allele count table. The first column should have allele counts for the minor allele, 
      and the second column should have allele counts for the major allele. 
      The first, second and third row should correspond to case, internal control and external control, respectively. }

            
}
\value{
	\item{OR}{OR estimate from iECAT. }
	\item{SE}{SE estimate of OR.}
	\item{p.value}{p-value of iECAT. }
	\item{p.value.noadj}{p-value without the adjustment of possible batch effects}
	\item{p.value.internal}{p-value with only internal study samples}
	\item{p.value.IvsE}{p-value of comparing internal vs external control samples}
}

\author{Seunggeun (Shawn) Lee}

\examples{
library(SKAT)

data(Example, package="iECAT")
attach(Example)

# iECAT-O
# test the second SNP in the first gene
idx.snp<-2

obj<-SKAT_Null_Model(Y ~ 1, out_type="D")
Z = Z.list[[1]][,idx.snp]
tbl.external= tbl.external.all.list[[1]][idx.snp,]


iECAT_SingleVar(Z, obj, tbl.external)



}

