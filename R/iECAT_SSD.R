
#
#	File.EC format: CHR, POS, ID, A1, A2, A1_count, A2_count
#
Generate_SSDOBJ_SetID_wEC<-function(File.Bim, File.SetID, File.EC, File.ECInfo=NULL){

	SetID<-try(read.table(File.SetID, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(SetID)=="try-error"){
		stop("Error in SetID file!") 
	}	
	BIM<-try(read.table(File.Bim, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(BIM)=="try-error"){
		stop("Error in BIM file!") 
	}	
	EC<-try(read.table(File.EC, header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
	if(class(EC)=="try-error"){
		stop("Error in EC file!") 
	}		
	
	colnames(BIM)<-c("CHR", "ID", "D1", "POS", "A1", "A2")
	colnames(EC)<-c("CHR", "POS", "EC_ID", "EC_A1", "EC_A2", "EC_C1", "EC_C2")
	colnames(SetID)<-c("SetID", "ID")

	
	# Get list of unique SNPs in SetID 
	SetID_uniqueSNP.df = data.frame(ID = unique(SetID$ID))
	BIM1<-merge(SetID_uniqueSNP.df, BIM, by.x="ID", by.y="ID", all=FALSE)
	BIM1$POS_ID = sprintf("%d:%d", BIM1$CHR, BIM1$POS)
	EC$POS_ID = sprintf("%d:%d", EC$CHR, EC$POS)	
	
	
	EC_INFO<-merge(BIM1, EC, by.x="POS_ID", by.y="POS_ID", all=FALSE)

	# log
	n.common<-nrow(EC_INFO)
	cat(n.common, "SNPs observed in both SetID (and BIM) and EC files\n")
	
	# Check whether Alleles are correct 
	idx1<-which(EC_INFO$A1 == EC_INFO$EC_A1)
	idx2<-which(EC_INFO$A1 == EC_INFO$EC_A2)
	idx3<-which(EC_INFO$A2 == EC_INFO$EC_A2)
	idx4<-which(EC_INFO$A2 == EC_INFO$EC_A1)

	#include 
	idx.include<-intersect(union(idx1, idx2), union(idx3, idx4))
	
	#matching alleles 
	temp1<-EC_INFO$EC_C1[idx2]
	temp2<-EC_INFO$EC_A1[idx2]
	
	EC_INFO$EC_C1[idx2]<-EC_INFO$EC_C2[idx2]
	EC_INFO$EC_C2[idx2]<-temp1
	EC_INFO$EC_A1[idx2]<-EC_INFO$EC_A2[idx2]
	EC_INFO$EC_A2[idx2]<-temp2
	
	n.common<-nrow(EC_INFO)
	cat(n.common, "After checking alleles,", n.common, " SNPs left\n")

	EC_INFO<-EC_INFO[idx.include,]
	EC_INFO$ID<-as.character(EC_INFO$ID)
		
	# Create hash table
	EC_hash<-new.env()
	for(i in 1:n.common){
		EC_hash[[EC_INFO$ID[i]]]<-i
	}
	
	obj.EC<-list(EC_INFO=EC_INFO, EC_hash=EC_hash)
	if(!is.null(File.ECInfo)){
		save(obj.EC, file=File.ECInfo)
	}
	return(obj.EC)
}


Generate_SSD_SetID_wEC<-function(File.Bed, File.Bim, File.Fam, File.SetID, File.EC, File.SSD, File.Info, File.EC.Info){
	
	#
	obj.EC<-Generate_SSDOBJ_SetID_wEC(File.Bim, File.SetID, File.EC, File.EC.Info)
	SKAT::Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=FALSE)

	#return(obj.EC)
}

Open_SSD_wEC<-function(File.SSD, File.Info, File.EC.Info){

	obj.EC=NULL
	SSD.INFO=SKAT::Open_SSD(File.SSD, File.Info)
	load(File.EC.Info)
	
	EC.INFO<-list(SSD.INFO=SSD.INFO, obj.EC=obj.EC)
	return(EC.INFO)
	
}

Close_SSD_wEC<-function(){

	SKAT::Close_SSD()

}



iECAT.SSD.GetSNP_Weight<-function(SSD.INFO, SetIndex, obj.EC, obj.SNPWeight=NULL){


	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]
	
	try1<-try(SKAT::Get_Genotypes_SSD(SSD.INFO, SetIndex, is_ID=TRUE),silent = TRUE)
	if(class(try1) != "try-error"){
		Z<-try1
		Is.Error<-FALSE	
	} else {
		err.msg<-geterrmessage()
		msg<-sprintf("Error to get genotypes of %s: %s",SetID, err.msg)
		stop(msg)
	}
	
	Is.weights=FALSE
	if(!is.null(obj.SNPWeight)){
		Is.weights=TRUE
	}

	SNP_ID<-colnames(Z)
	p<-ncol(Z)
	weights<-rep(0, p)
	EC_COUNT<-matrix(rep(-1,p*2), ncol=2)
	for(i in 1:p){
		val1<-SNP_ID[i]			
		
		if(Is.weights){
			val2<-obj.SNPWeight$hashset[[val1]]
			
			if(is.null(val2)){
				msg<-sprintf("SNP %s is not found in obj.SNPWeight!", val1)
				stop(msg)
			}

			weights[i]<-val2
		}
		idx<-obj.EC$EC_hash[[val1]]
		if(is.null(idx)){
			msg<-sprintf("SNP %s is not found in obj.SNPWeight!", val1)
			warnings(msg)
		} else {
			EC_COUNT[i,1]<-obj.EC$EC_INFO$EC_C1[idx]
			EC_COUNT[i,2]<-obj.EC$EC_INFO$EC_C2[idx]
		}
	}
	
	# Change FALSE to TRUE
	re=list(Z=Z, Is.weights=Is.weights, weights=weights, EC_COUNT=EC_COUNT)
	return(re)
		

}


iECAT.SSD.OneSet_SetIndex = function(EC.INFO, SetIndex,  obj, ..., obj.SNPWeight=NULL){
	
	#SetIndex<-1; obj.SNPWeight=NULL
	SSD.INFO<-EC.INFO$SSD.INFO
	obj.EC<-EC.INFO$obj.EC
	
	re1 = iECAT.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.EC, obj.SNPWeight=obj.SNPWeight)
	
	weights=NULL
	if(re1$Is.weights){
		weights=re1$weights
	} 
	

	tbl.external.all = re1$EC_COUNT
	# need to check later
	Z<-re1$Z
	Z<-2-Z
	
	# change the coding scheme to the minor allele coding
	Z[Z==9]<-NA
	MAF = colMeans(Z, na.rm=TRUE)/2
	idx.flip<-which(MAF > 0.5)
	if(length(idx.flip) > 0){
		Z[,idx.flip] = 2 - Z[,idx.flip]
		tbl.external.all[idx.flip, c(2,1)]<-tbl.external.all[idx.flip, c(1,2)]
	}
	
	re<-iECAT(Z, obj, tbl.external.all = tbl.external.all, weights = weights, ...)
	return(re)
}



#
# Only SKAT_Null_Model obj can be used
#
iECAT.SSD.All = function(EC.INFO, obj, ..., obj.SNPWeight=NULL){


	SSD.INFO<-EC.INFO$SSD.INFO
	obj.EC<-EC.INFO$obj.EC
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Pvalue.Noadj<-rep(NA,N.Set)
	OUT.Pvalue.Internal<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	# Does not provide resampling 
	Is.Resampling = FALSE
	n.Resampling = 0
	
	if(class(obj) == "SKAT_NULL_Model"){
		if(obj$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$n.Resampling

			#OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	} else if(class(obj) == "SKAT_NULL_Model_ADJ"){
		if(obj$re1$n.Resampling > 0){
			Is.Resampling = TRUE
			n.Resampling = obj$re1$n.Resampling

			#OUT.Pvalue.Resampling<-matrix(rep(0,n.Resampling*N.Set),ncol=n.Resampling)
		}
	}
	Is.Resampling = FALSE
	n.Resampling = 0

	#iECAT.SSD.OneSet_SetIndex = function(EC.INFO, SetIndex,  obj, ..., obj.SNPWeight=NULL){

	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1 = try(iECAT.SSD.OneSet_SetIndex(EC.INFO=EC.INFO, SetIndex=i, obj=obj, ..., obj.SNPWeight=obj.SNPWeight))
		#try1 = try(iECAT.SSD.OneSet_SetIndex(EC.INFO=EC.INFO, SetIndex=i, obj=obj))

		
		if(class(try1) != "try-error"){
			re<-try1
			Is.Error<-FALSE
		} else {

			err.msg<-geterrmessage()
			msg<-sprintf("Error to run iECAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
			
		}
		
		if(!Is.Error){

			OUT.Pvalue[i]<-re$p.value
			OUT.Pvalue.Noadj[i]<-re$p.value.noadj
			OUT.Pvalue.Internal[i]<-re$p.value.internal
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
			if(Is.Resampling){
				OUT.Pvalue.Resampling[i,]<-re$p.value.resampling
			}
		}
		
		if(floor(i/100)*100 == i){
			cat("\r", i, "/", N.Set, "were done");
		}

	}

	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, P.value.Noadj=OUT.Pvalue.Noadj, P.value.Internal=OUT.Pvalue.Internal, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl,P.value.Resampling=OUT.Pvalue.Resampling)
	class(re)<-"SKAT_SSD_ALL"

	return(re)	
}



