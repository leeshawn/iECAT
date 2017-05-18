 \name{Generate_SSD_SetID_wEC}
 \alias{Generate_SSD_SetID_wEC}
 \title{Generate SNP set data file (SSD) with using external control information }
 \description{
	Generate SNP set data (SSD) and external control information files

 }
 \usage{
 
 	Generate_SSD_SetID_wEC(File.Bed, File.Bim, File.Fam, File.SetID
 	, File.EC, File.SSD, File.Info, File.EC.Info)

 }
\arguments{
      \item{File.Bed}{name of the binary ped file (BED).}
      \item{File.Bim}{name of the binary map file (BIM).}
      \item{File.Fam}{name of the FAM file (FAM).}
      \item{File.SetID}{name of the SNP set ID file that defines SNP sets. 
      The first column must be Set ID, and the second column must be SNP ID. There should be no header!! }
      \item{File.EC}{name of the file that has allele count information of the external controls. }
      \item{File.SSD}{name of the SSD file generated. }
      \item{File.Info}{name of the SSD info file generated. }
      \item{File.EC.Info}{name of the external control info file generated. }
      
}

\details{
 The SetID file is a white-space (space or tab) delimitered file with 2 columns:
  SetID and SNP_ID.

 Please keep in mind that there should be no header!
 The SNP_IDs and SetIDs should be less than 50 characters, otherwise, it will return an error message.
    
 File.EC should be a text formated file 
    
 The SSD file is a binary formated file with genotypes. 
 The SSD info file is a text file with general information on data and SNP sets (first 6 rows), 
 and information on each set (after 8th row).
         
}


\author{Seunggeun Lee, Larisa Miropolsky}

