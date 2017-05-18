 \name{iECAT.SSD.All}
 \alias{iECAT.SSD.All}
 \title{Integrating External Controls to Association Tests}
 \description{
	Iteratively carry out association tests with phenotypes and SNP sets in SSD file. 
 }
 \usage{

	iECAT.SSD.All(EC.INFO, obj, \dots, obj.SNPWeight=NULL)

 }
\arguments{

      \item{EC.INFO}{EC_INFO object returned from Open_SSD_wEC.   }
      \item{obj}{output object from SKAT_Null_Model. }
      \item{\dots}{further arguments to be passed to ``iECAT''. }
      \item{obj.SNPWeight}{output object from Read_SNP_WeightFile (default=NULL). 
      If NULL, the beta weight with the ``weights.beta'' parameter will be used.  }
}
\value{
	\item{results}{dataframe that contains SetID, p-values (P.value, P.value.Noadj, and P.value.Internal), the number of markers in the SNP sets (N.Marker.All), 
	and the number of markers to test for associations (N.Marker.Test).   }
	\item{P.value.Resampling}{currently resampling p-values are not provided.}
}
\details{
Please see iECAT for details.                     
}


\author{Seunggeun Lee}

