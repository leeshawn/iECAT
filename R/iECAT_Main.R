
#
#	Check genotype matrix
#
iECAT_Genotype_Check<-function(Z, obj.res, tbl.external.all, impute.method = "bestguess", missing_cutoff,  MAC.limit, MAF.limit, SetID=NULL, Is.chrX = FALSE, SexVar = NULL){

	n<-dim(Z)[1]
	m<-dim(Z)[2]
	
	id_include = obj.res$id_include
	idx.miss<-which(Z==9)
	if(length(idx.miss) > 0){
		Z[idx.miss]<-NA
	}
	
	MAC<-colSums(Z, na.rm=TRUE)
	MAC1<-colSums(2-Z, na.rm=TRUE)
	MAF<-colMeans(Z, na.rm=TRUE)/2
	idx<-which(MAF > 0.5)
	if(length(idx) > 0){
		MAF[idx]<-1-MAF[idx]
	}
	
	ID_INCLUDE_SNP <- NULL
    for (i in 1:m) {
        missing.ratio <- length(which(is.na(Z[, i])))/n
        sd1 <- sd(Z[id_include, i], na.rm = TRUE)
        if (missing.ratio < missing_cutoff && sd1 > 0) {
            ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP, i)
        } 
    }
    #MAC<<-MAC
    #MAC1<<-MAC1
    ID_INCLUDE_SNP<-intersect(ID_INCLUDE_SNP, intersect(which(MAC > MAC.limit), which(MAC1  > MAC.limit)))
    ID_INCLUDE_SNP<-intersect(ID_INCLUDE_SNP, which(MAF <=MAF.limit))

    
    if (length(ID_INCLUDE_SNP) == 0) {
        if (is.null(SetID)) {
            msg <- sprintf("ALL SNPs have either low MAC, high missing rates or no-variation. P-value=1")
        }
        else {
            msg <- sprintf("In %s, ALL SNPs have either low MAC, high missing rates or no-variation. P-value=1", 
                SetID)
        }
        warning(msg, call. = FALSE)
        re <- list(p.value=1, p.value.noadj=1, p.value.internal = 1, param=list( n.marker= m, n.marker.test = 0), return = 1)
        return(re)
    } else if (m - length(ID_INCLUDE_SNP) > 0) {
        if (is.null(SetID)) {
            msg <- sprintf("%d SNPs with either low MAC, high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
        }
        else {
            msg <- sprintf("In %s, %d SNPs with either low MAC, high missing rates or no-variation are excluded!", SetID, m - length(ID_INCLUDE_SNP))
        }
        warning(msg, call. = FALSE)
    }
	Z1<-cbind(Z[id_include,])

	out.z<-list();
	

	out.z$ID_INCLUDE_SNP=ID_INCLUDE_SNP
	out.z$param<-list(n.marker=m, n.marker.test =length(ID_INCLUDE_SNP))
	out.z$Z1=Z1
	out.z$return=0
		
	return(out.z)
}


iECAT_SingleVar_Tbl<-function(tbl.internal, tbl.external, Fisher.test=FALSE){

	#check dat.tbl
	dat.tbl<-rbind
	n1<-nrow(tbl.internal)
	n2<-ncol(dat.tbl)
	
	if(n1 != 3 || n2 != 2){
		stop("dat.tbl should be 3 x 2 matrix")
	}

	# parameter for internal use
	default.weight.factor=-9
	re.work<-Get_Bayes_Pval(dat.tbl, Y=NULL, X=NULL, G=NULL, weight.factor=default.weight.factor, Fisher.test=Fisher.test)
	
	#		re<-list(OR=OR_new.a[1], Pval.asymptotic=p1, Pval.naive=Pval.naive, Pval.internal=Pval.internal , W=out.b$temp1, 
	#  Pval.internal.f=Pval.internal.f, Pval.naive.f=Pval.naive.f, Pval.ivse=Pval.ivse, SE=OR_new.a[2], Var_Adj.Factor=out.b$Var_Adj.Factor)


	re<-list(OR = re.work$OR, SE = re.work$SE, p.value=re.work$Pval.asymptotic, p.value.noadj=re.work$Pval.naive, 
	p.value.internal = re.work$Pval.internal, p.value.IvsE = re.work$Pval.ivse)
	return(re)

}


iECAT_SingleVar_Tbl<-function(Tbl, Fisher.test=FALSE, weight.factor=NULL){

	#Fisher.test=FALSE
	# parameter for internal use
	default.weight.factor=-9
	
	if(!is.null(weight.factor)){
		default.weight.factor=weight.factor
	}

	re.work<-Get_Bayes_Pval(Tbl, Y=NULL, X=NULL, G=NULL, weight.factor=default.weight.factor, Fisher.test=Fisher.test)
	
	#		re<-list(OR=OR_new.a[1], Pval.asymptotic=p1, Pval.naive=Pval.naive, Pval.internal=Pval.internal , W=out.b$temp1, 
	#  Pval.internal.f=Pval.internal.f, Pval.naive.f=Pval.naive.f, Pval.ivse=Pval.ivse, SE=OR_new.a[2], Var_Adj.Factor=out.b$Var_Adj.Factor)


	re<-list(OR = re.work$OR, SE = re.work$SE, p.value=re.work$Pval.asymptotic, p.value.noadj=re.work$Pval.naive, 
	p.value.internal = re.work$Pval.internal, p.value.IvsE = re.work$Pval.ivse, W=re.work$W[1])
	return(re)


}

iECAT_SingleVar<-function(Z, obj, tbl.external, Fisher.test=FALSE, weight.factor=NULL){

	if(ncol(cbind(Z)) > 1){
		stop("Z should be a vector or a matrix with one column!")
	}

	# obj
	obj.res = SKAT:::SKATExactBin_CheckObj(obj)
	y = round(obj.res$mu + obj.res$res)
	
	# check_Z
	out.Z<-iECAT_Genotype_Check(cbind(Z), obj.res, impute.method = "bestguess", missing_cutoff=1, MAC.limit=0, MAF.limit=1)
	if(out.Z$return == 1){
		return(out.Z)
	}
	Z1<-out.Z$Z1
	
	#
	#	Make tables
	#	
	
	r.case<-sum(Z1[y==1], na.rm = TRUE)
	r.control<-sum(Z1[y==0], na.rm = TRUE)
	r1.case<-sum(2-Z1[y==1], na.rm = TRUE)
	r1.control<-sum(2-Z1[y==0], na.rm = TRUE)
	
	Tbl<-matrix(rep(0, 3*2), ncol=2)
	Tbl[1,]<-c(r.case, r1.case)
	Tbl[2,]<-c(r.control, r1.control)
	Tbl[3,]<-tbl.external

		
	re<-iECAT_SingleVar_Tbl(Tbl=Tbl, Fisher.test=Fisher.test, weight.factor=weight.factor)
	return(re)

}


iECAT<-function(Z, obj, tbl.external.all, weights.beta=c(1,25), weights = NULL, r.corr=0, method="davies", missing_cutoff=0.15, MAC.lowlimit=3, MAF.limit=1){

	re<-iECAT_Work(Z=Z, obj=obj, tbl.external.all= tbl.external.all, weights.beta=weights.beta, weights = weights, r.corr=r.corr, method=method, 
	missing_cutoff=missing_cutoff, MAC.lowlimit=MAC.lowlimit, MAF.limit=MAF.limit)

	n.re<-length(re)
	return(re[-n.re])
	
}

iECAT_Work<-function(Z, obj, tbl.external.all, weights.beta=c(1,25), weights = NULL, r.corr=0, method="davies", missing_cutoff=0.15, MAC.lowlimit=3, MAF.limit=1
, weight.factor=NULL, obj.meta=NULL, idx.error=NULL, Is.Get.Noadj=TRUE, Is.Get.Internal=TRUE, Is.Shrink.Weight=TRUE){

	#SKAT_wExternalControl<-function(Z, Tbl.a.list, is.use.external, weights.beta=c(1,25), r.corr=0, method="davies", weight.factor=-2, obj=NULL, obj.meta=NULL, Is.Shrink.Weight=TRUE,
	#n.int=-1, n.all=-1, Is.Get.Naive=TRUE, Is.Get.Internal=TRUE,  Shrink.Cutoff=1, Mode=1, idx.error=NULL){
	
	# weights.beta=c(1,25); weights = NULL; r.corr=0; method="optimal.adj"; missing_cutoff=0.15
	
	if(method=="optimal.adj"){
		method="optimal"
	}
	
	# parameter for internal use
	default.weight.factor=-9
	
	if(!is.null(weight.factor)){
		default.weight.factor=weight.factor
	}
	
	# obj
	obj.res = SKAT:::SKATExactBin_CheckObj(obj)
	y = round(obj.res$mu + obj.res$res)
	
	out.Z<-iECAT_Genotype_Check(Z, obj.res, impute.method = "bestguess", missing_cutoff=missing_cutoff, MAC.limit=MAC.lowlimit, MAF.limit=MAF.limit)
	if(out.Z$return == 1){
		return(out.Z)
	}
	
	Z1<-out.Z$Z1
	id_include_snp<-out.Z$ID_INCLUDE_SNP
	
	# n.case will be used if weight.factor <= -4. It is the only difference between factor=-3 and -4
	# -5 and -6 test
	n.case=-1
	if( default.weight.factor <= -4  ){
		n.case<-sum(y)
	}
	
	n.int<-length(obj.res$id_include)
	n.all<-n.int + max(rowSums(tbl.external.all)/2)
	
	p<-ncol(Z)
	p.ext<-nrow(tbl.external.all)
	
	if(p != p.ext){
		stop("The number of columns in Z should be the same as the number of rows in tbl.external.all!")
	} 	
	
	
	if(!is.null(weights)){
		weights<-weights[id_include_snp]
	}
	
	#
	#	Make tables
	#	
	Tbl.a.list<-list()
	r.case<-colSums(cbind(Z1[y==1,]), na.rm = TRUE)
	r.control<-colSums(cbind(Z1[y==0,]), na.rm = TRUE)
	r1.case<-colSums(cbind(2-Z1[y==1,]), na.rm = TRUE)
	r1.control<-colSums(cbind(2-Z1[y==0,]), na.rm = TRUE)
	
	is.use.external<-rep(1, length(id_include_snp))
	for(i in 1:length(id_include_snp)){
		j<-id_include_snp[i]
		Tbl<-matrix(rep(0, 3*2), ncol=2)
		Tbl[1,]<-c(r.case[j], r1.case[j])
		Tbl[2,]<-c(r.control[j], r1.control[j])
		Tbl[3,]<-tbl.external.all[j,]
		Tbl.a.list[[i]]<-Tbl
		
		if(tbl.external.all[j,1] < 0){
			is.use.external[i] = 0
		}
	}

	re.work<-SKAT_wExternalControl(cbind(Z1[,id_include_snp]), Tbl.a.list, is.use.external=is.use.external, weights.beta=weights.beta, weights=weights
	, r.corr=r.corr, method=method, weight.factor=default.weight.factor, obj=obj, n.case=n.case, n.int=n.int, 
	n.all=n.all, obj.meta=obj.meta, idx.error=idx.error, Is.Get.Naive = Is.Get.Noadj, Is.Get.Internal=Is.Get.Internal, Is.Shrink.Weight=Is.Shrink.Weight)

	param<-out.Z$param
	re<-list(p.value=re.work$p.value, p.value.noadj=re.work$naive.skat, p.value.internal = re.work$int.skat, param=param, re.work=re.work)
	
	return(re)
	
}









