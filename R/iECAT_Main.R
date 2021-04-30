
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



iECAT_SingleVar_Score<-function(Z, Y, X, internal.indicator, method, MAF.adjust){
	if(ncol(cbind(Z)) > 1){
		stop("Z should be a vector or a matrix with one column!")
	}
	
	if(length(Z) != length(internal.indicator)) {
		stop("Need indicator of internval vs external sources for all samples!")
	}
	
	#--- Choice of method ---#
	# "Internal", "Naive", "iECAT", "iECATminP"
	if (method=="Internal") {
		if (MAF.adjust==TRUE) {print("MAF adjustment not applicable, default to FALSE")}
		re<- SingleVar_Score_internal(Z, Y, X, internal.indicator)
	} else if (method=="Naive") { #"Naive"
		if (MAF.adjust==TRUE) {print("MAF adjustment not applicable, default to FALSE")}
		re<- SingleVar_Score_naive(Z, Y, X, internal.indicator)
	} else {
		re<- SingleVar_Score_iECAT(Z, Y, X, internal.indicator, method, MAF.adjust)
	}
	
	return(re)


}



SingleVar_Score_internal<- function(Z, Y, X, internal.indicator){
	#--- Get index ---#
	idx.internal<- which(internal.indicator==1)

	#--- Create obj ---#
	data<- data.frame(cbind(Y, X))
	G<- Z[idx.internal]
	
	dat.null<-model.frame(Y~., data=data, na.action = na.pass)
	# there is an issue if only one column exists in dat.null
	if(ncol(dat.null)==1){
    	dat.null$add_one_column = 1
	}

	
	Null.internal<-ScoreTest_NULL_Model(Y~., data=as.data.frame(dat.null[idx.internal,]) )
	
	#--- Testing ---#
	A1<- (G[idx.internal]  -  Null.internal$XXVX_inv %*%  (Null.internal$XV %*% G[idx.internal]))[,1]
	S1<- sum(A1 * Null.internal$res)
	Var1 <- sum(A1 *(G[idx.internal] * Null.internal$V) )
	
	#--- SPA and ER calibration and update variance ---#
	S1.q <- sum(A1*(Null.internal$res + Null.internal$mu))/sqrt(sum(G[idx.internal]))
	g1 <- A1/sqrt(sum(G[idx.internal]))
	mu.internal <- sum(Null.internal$mu * g1)
	var.internal <- sum(Null.internal$mu * (1-Null.internal$mu) *g1^2)
	stat_S1.q <- (S1.q-mu.internal)^2/var.internal
	p_S1.q <- pchisq(stat_S1.q, df=1, lower.tail=F)
	zscore_S1.q <- (S1.q-mu.internal)*sqrt(sum(G[idx.internal]))
	pnew_S1.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null[idx.internal,]), G=G[idx.internal], q=S1.q, stat.qtemp=stat_S1.q, mu=Null.internal$mu, g=g1)
	if (pnew_S1.q>0) {Var.spa.S1 <- zscore_S1.q^2/qchisq(pnew_S1.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S1 <- Var1}
	
	p.value<- pchisq(S1^2 / Var.spa.S1, df=1, lower.tail=FALSE)
	return(list=c(method="Internal only", MAF.adjust=NA, Score=S1, VAR=Var.spa.S1, p.value=p.value))
}



SingleVar_Score_naive<- function(Z, Y, X, internal.indicator){
	#--- Create obj ---#
	data<- data.frame(cbind(Y, X))
	G<- Z
	
	dat.null<-model.frame(Y~., data=data, na.action = na.pass)
	# there is an issue if only one column exists in dat.null
	if(ncol(dat.null)==1){
    	dat.null$add_one_column = 1
	}
	
	Null.all<-ScoreTest_NULL_Model(Y~., data=dat.null )
	
	#--- Testing ---#
	A2<- (G  -  Null.all$XXVX_inv %*%  (Null.all$XV %*% G))[,1]
	S2<- sum(A2 * Null.all$res)
	Var2 <- sum(A2 *(G * Null.all$V) )
	
	#--- SPA and ER calibration and update variance ---#
	S2.q <- sum(A2*(Null.all$res + Null.all$mu))/sqrt(sum(G))
	g2 <- A2/sqrt(sum(G))
	mu.all <- sum(Null.all$mu * g2)
	var.all <- sum(Null.all$mu * (1-Null.all$mu) *g2^2)
	stat_S2.q <- (S2.q-mu.all)^2/var.all
	p_S2.q <- pchisq(stat_S2.q, df=1, lower.tail=F)
	zscore_S2.q <- (S2.q-mu.all)*sqrt(sum(G))
	pnew_S2.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null), G=G, q=S2.q, stat.qtemp=stat_S2.q, mu=Null.all$mu, g=g2)
	if (pnew_S2.q>0) {Var.spa.S2 <- zscore_S2.q^2/qchisq(pnew_S2.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S2 <- Sigma[2,2]}
	
	p.value<- pchisq(S2^2 / Var.spa.S2, df=1, lower.tail=FALSE)
	return(list=c(method="Naive method, no adjustment", MAF.adjust=NA, Score=S2, VAR=Var.spa.S2, p.value=p.value))

}


SingleVar_Score_iECAT <- function(Z, Y, X, internal.indicator, method, MAF.adjust){
	#--- Get index ---#
	n<- sum(internal.indicator)
	N<- length(Z)
	idx.internal<- which(internal.indicator==1)
	idx.external<- which(internal.indicator==0)
	ncase<- sum(Y); idx.case<- which(Y==1)
	ncontrol<- N - ncase; idx.control<- which(Y==0)
	idx.internal.case<- intersect(idx.internal, idx.case)
	idx.internal.control<- intersect(idx.internal, idx.control)
	idx.whichinternal.iscontrol<- which(idx.internal %in% idx.control)
	idx.control.internal<- which(idx.control %in% idx.internal)
	idx.control.external<- which(idx.control %in% idx.external)
	n1 <- length(idx.internal.case)
	n2 <- length(idx.internal.control)
	n3 <- length(idx.external)
	a <- n*(n2+n3)/(n2*N)

	#--- Create obj ---#
	data<- data.frame(cbind(Y, X))
	G<- Z
	
	dat.null<-model.frame(Y~., data=data, na.action = na.pass)
	# there is an issue if only one column exists in dat.null
	if(ncol(dat.null)==1){
    	dat.null$add_one_column = 1
	}
	
	Y <- dat.null$Y
	X <- dat.null[,-1]

	Y1<-rep(1,N)
	Y1[idx.external]<-0
	
	dat.null1<-dat.null
	dat.null1[,1]<-Y1

	Null.internal<-ScoreTest_NULL_Model(Y~., data=as.data.frame(dat.null[idx.internal,]) )
	Null.all<-ScoreTest_NULL_Model(Y~., data=dat.null )
	Null.IvE<-ScoreTest_NULL_Model(Y~., data=as.data.frame(dat.null1[idx.control,]) )
	#--- end of creating obj ---#
	
	
	#--- Testing ---#
	A3<- (G[idx.control]  -  Null.IvE$XXVX_inv %*%  (Null.IvE$XV %*% G[idx.control]))[,1]
	
	if (sum(A3)==0){
		stop("Detected monomorphism in controls, unable to use external controls!")
	}else{
		re<- SingleVar_Score_iECAT_kernel(G, Y, X, Null.internal, Null.all, Null.IvE, A3, method, MAF.adjust)
	}#end of monomorphism check
}



SingleVar_Score_iECAT_kernel<- function(G, Y, X, Null.internal, Null.all, Null.IvE, A3, method, MAF.adjust, env = parent.frame()){
	
	#--- Initial calculation of Score and variance ---###
	Score<-rep(0,3)
	Var_EST<-rep(0,5)
	A1<- (G[env$idx.internal]  -  Null.internal$XXVX_inv %*%  (Null.internal$XV %*% G[env$idx.internal]))[,1]
	A2<- (G  -  Null.all$XXVX_inv %*%  (Null.all$XV %*% G))[,1]
	
	Score[1]<-sum(A1 * Null.internal$res)
	Score[2]<-sum(A2 * Null.all$res)
	Score[3]<-sum(A3 * Null.IvE$res)
	
	Var_EST[1] <- sum(A1 *(G[env$idx.internal] * Null.internal$V) )
	Var_EST[2] <- sum(A2 *(G * Null.all$V) )
	Var_EST[3] <- sum(A3 * (G[env$idx.control] * Null.IvE$V))
	Var_EST[4] <- sum((A1 * Null.internal$res)[env$idx.whichinternal.iscontrol] * (A3 * Null.IvE$res)[env$idx.control.internal] )
	Var_EST[5] <- Var_EST[1]
	Var_EST[6] <- sum((A2 * Null.all$res)[env$idx.control] * (A3 * Null.IvE$res) )
	
	Sigma<-matrix(rep(0,9), ncol=3)
	Sigma[1,]<-c(Var_EST[1]/2, Var_EST[5], Var_EST[4])
	Sigma[2,]<-c(0, Var_EST[2]/2, Var_EST[6]/2)
	Sigma[3,]<-c(0, 0, Var_EST[3]/2)
	Sigma<-Sigma +t(Sigma)
	rhorho <- Sigma[1,3]/sqrt(Sigma[1,1]*Sigma[3,3])
	
	S.IvE <-sum(A3 * Null.IvE$res)
	Var.IvE <- sum(A3 * (G[env$idx.control] * Null.IvE$V))
	Z <- S.IvE/sqrt(Var.IvE)
	pval.score.IvE <- pchisq(Z^2, df=1, lower.tail=FALSE)
	
	#--- SPA and ER calibration and update variance ---#
	S1.q <- sum(A1*(Null.internal$res + Null.internal$mu))/sqrt(sum(G[env$idx.internal]))
	S2.q <- sum(A2*(Null.all$res + Null.all$mu))/sqrt(sum(G))
	S3.q <- sum(A3*(Null.IvE$res + Null.IvE$mu))/sqrt(sum(G[env$idx.control]))
	
	g1 <- A1/sqrt(sum(G[env$idx.internal]))
	g2 <- A2/sqrt(sum(G))
	g3 <- A3/sqrt(sum(G[env$idx.control]))
	
	mu.internal <- sum(Null.internal$mu * g1)
	var.internal <- sum(Null.internal$mu * (1-Null.internal$mu) *g1^2)
	mu.all <- sum(Null.all$mu * g2)
	var.all <- sum(Null.all$mu * (1-Null.all$mu) *g2^2)
	mu.IvE <- sum(Null.IvE$mu * g3)
	var.IvE <- sum(Null.IvE$mu * (1-Null.IvE$mu) *g3^2)
	
	
	stat_S1.q <- (S1.q-mu.internal)^2/var.internal
	p_S1.q <- pchisq(stat_S1.q, df=1, lower.tail=F)
	zscore_S1.q <- (S1.q-mu.internal)*sqrt(sum(G[env$idx.internal]))
	pnew_S1.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null[env$idx.internal,]), G=G[env$idx.internal], q=S1.q, stat.qtemp=stat_S1.q, mu=Null.internal$mu, g=g1)
	
	stat_S2.q <- (S2.q-mu.all)^2/var.all
	p_S2.q <- pchisq(stat_S2.q, df=1, lower.tail=F)
	zscore_S2.q <- (S2.q-mu.all)*sqrt(sum(G))
	pnew_S2.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null), G=G, q=S2.q, stat.qtemp=stat_S2.q, mu=Null.all$mu, g=g2)
	
	stat_S3.q <- (S3.q-mu.IvE)^2/var.IvE
	p_S3.q <- pchisq(stat_S3.q, df=1, lower.tail=F)
	zscore_S3.q <- (S3.q-mu.IvE)*sqrt(sum(G[env$idx.control]))
	pnew_S3.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null1[env$idx.control,]), G=G[env$idx.control], q=S3.q, stat.qtemp=stat_S3.q, mu=Null.IvE$mu, g=g3)
	
	if (pnew_S1.q>0) {Var.spa.S1 <- zscore_S1.q^2/qchisq(pnew_S1.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S1 <- Sigma[1,1]}
	if (pnew_S2.q>0) {Var.spa.S2 <- zscore_S2.q^2/qchisq(pnew_S2.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S2 <- Sigma[2,2]}
	if (pnew_S3.q>0) {Var.spa.S3 <- zscore_S3.q^2/qchisq(pnew_S3.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S3 <- Sigma[3,3]}
	
	if (pnew_S3.q==0 | (pnew_S3.q<0 & pval.score.IvE<=0)) {
		tau.spa<-1
	} else {
		tau1.spa<- Score[3]^2 / Var.spa.S3
		tau.spa<- tau1.spa/(1+tau1.spa)
	}
	
	#--- Calculate weight tau and adjustment using MAF ---#
	if (MAF.adjust==TRUE) {
		MAF.case<-mean(G[env$idx.internal.case])/2
		MAF.control<-mean(G[env$idx.internal.control])/2
		MAF.external<-mean(G[env$idx.external])/2
		
		if(MAF.external > MAF.control && MAF.external < MAF.case){
			tau.spa<- 0
		} else if(MAF.external < MAF.control && MAF.external > MAF.case ){
			tau.spa<- 0
		}
	}
	
	#--- Calculate compound score Sw and its variance ---#
	Sigma.spa<- Sigma
	diag(Sigma.spa) <- c(Var.spa.S1, Var.spa.S2, Var.spa.S3)
	Sigma.spa[1,2] <- Sigma.spa[2,1] <- Var.spa.S1
	Sigma.spa[1,3] <- Sigma.spa[3,1] <- rhorho*sqrt(Sigma.spa[1,1]*Sigma.spa[3,3])
	
	d_tau.spa<-(2*Score[3] * Var.spa.S3)/(Score[3]^2 + Var.spa.S3)^2
	if(tau.spa==1 || tau.spa==0){d_tau.spa=0}
	
	G_d.spa<-c(a*tau.spa, 1- tau.spa, d_tau.spa * (a*Score[1] - Score[2] ))
	
	
	Sw.spa <- tau.spa*a* Score[1] + (1-tau.spa)*Score[2] #new
	VARw.spa <- t(G_d.spa) %*% Sigma.spa %*% G_d.spa #new

	#--- Calculate p-value ---#	
	p.value<- pchisq(Sw.spa^2 / VARw.spa, df=1, lower.tail=FALSE)
	p.value.tau1<- pchisq(Score[1]^2 / Var.spa.S1, df=1, lower.tail=FALSE)
	re<- list(c(method="iECAT", MAF.adjust=MAF.adjust, Score=Sw.spa, VAR=VARw.spa, p.value=p.value))
	
	if (method=="iECATminP") {
		rhorhorho <- sqrt((env$a*tau.spa+1-tau.spa)*Var.spa.S1/VARw.spa)
		if (abs(rhorhorho)<=1) {cmat <- matrix(c(1,rhorhorho,rhorhorho,1), nrow=2)} else{cmat <- matrix(c(1,sign(rhorhorho),sign(rhorhorho),1), nrow=2)}
		Z2 <- qnorm(min(p.value, p.value.tau1)/2)
		p.value <- pmvnorm(lower=c(-Inf,-abs(Z2)), upper=c(-abs(Z2),Inf), mean=c(0,0), corr=cmat)[1] + pmvnorm(lower=c(-abs(Z2),abs(Z2)), upper=c(Inf,Inf), mean=c(0,0), corr=cmat) + pmvnorm(lower=c(abs(Z2),-Inf), upper=c(Inf,abs(Z2)), mean=c(0,0), corr=cmat)[1] + pmvnorm(lower=c(-Inf,-Inf), upper=c(abs(Z2),-abs(Z2)), mean=c(0,0), corr=cmat)[1]
		re<- list(c(method="iECAT minP", MAF.adjust=MAF.adjust, Score=NA, VAR=NA, p.value=p.value))
	}
	
	return(re)
}







