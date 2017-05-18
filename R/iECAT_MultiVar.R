SKAT_wExternalControl.Work<-function(nSNP, Tbl.a.list, W, is.use.external, Is.Shrink.Weight, n.case, n.int, n.all, weight.factor, Shrink.Cutoff){



	OR<-matrix(rep(0, nSNP *6),ncol=6)
	MAC<-matrix(rep(0, nSNP*2), ncol=2)
	Shrink.W<-rep(0, nSNP)
	for(i in 1:nSNP){
	
		Tbl.a<-Tbl.a.list[[i]]
		MAC[i,1]<-sum(Tbl.a[1:2,1])
		MAC[i,2]<-sum(Tbl.a[,1], na.rm=TRUE)
		
		if(is.use.external[i]==1){
			
			if(length(which(Tbl.a==0)) > 0){
				Tbl.a<-Tbl.a + 0.5
			}
		
			re<-Get_Bayes_Pval(Tbl.a=Tbl.a, weight.factor=weight.factor)
			OR[i,1]<-re$OR 
			OR[i,2]<-re$SE
			Shrink.W[i]<-re$W
		
			Tbl1<-Tbl.a[1:2,]
			OR.in<-GetOR(Tbl1)
			OR[i,3]<-OR.in$OR 
			OR[i,4]<-OR.in$se
			
			Tbl2<-Tbl.a[1:2,] 
			Tbl2[2,] <- Tbl.a[2,] + Tbl.a[3,]
			OR.naive<-GetOR(Tbl2)
			OR[i,5]<-OR.naive$OR 
			OR[i,6]<-OR.naive$se
			
			
		} else {
			Tbl1<-Tbl.a[1:2,]
			if(length(which(Tbl1==0)) > 0){
				Tbl1<-Tbl1 + 0.5
			}
			
			OR.in<-GetOR(Tbl1)
			OR[i,1]<-OR.in$OR 
			OR[i,2]<-OR.in$se
			OR[i,3]<-OR.in$OR 
			OR[i,4]<-OR.in$se
			OR[i,5]<-OR.in$OR 
			OR[i,6]<-OR.in$se
		}
	}
	
	p1<-1-pchisq((OR[,1]/OR[,2])^2, df=1)
	p2<-1-pchisq((OR[,3]/OR[,4])^2, df=1)
	OR1<-OR
	#OR1[,1:2]<-OR1[,1:2]/OR1[,2]^2
	#OR1[,3:4]<-OR1[,3:4]/OR1[,2]^2
	OR1<- OR1 * W
	#OR1<-OR
	
	#Is.Shrink.Weight=FALSE
	Shrink.weight<-rep(1, nSNP)
	if(Is.Shrink.Weight){
		n1<-n.int
		n2<-n.all
		Shrink.weight<-(n2-n1)*Shrink.W + n1
		
		# n.case is not -1
		if(n.case > 0){
			y_bar1<-n.case/n.int
			y_bar2<-n.case/n.all			
			
			n1<-n.int * y_bar1 * (1-y_bar1)
			n2<-n.all * y_bar2 * (1-y_bar2)
			Shrink.weight<-(n2-n1)*Shrink.W + n1			
		}
		
		#if(Shrink.Cutoff==2){
		#	idx.w<-which(Shrink.W < 0.5)
		#} else if(Shrink.Cutoff==3){
		#	idx.w<-which(OR[,2] - OR[,4] > 0) 
		#} else { 
		#	idx.w<-which(Shrink.W < 0.2)
		#}
		
		#if(length(idx.w) > 0){
		#	OR1[idx.w,1:2]<-OR1[idx.w,3:4]
		#	Shrink.weight[idx.w]<-n1
		#}

		OR1[,1:2]<-rbind(OR1[,1:2]) * Shrink.weight
	}


	return(list(OR=OR, OR1=OR1, Shrink.W=Shrink.W, p1=p1, p2=p2, MAC=MAC, Shrink.weight=Shrink.weight))

}



SKAT_wExternalControl<-function(Z, Tbl.a.list, is.use.external, weights.beta=c(1,25), weights = NULL, r.corr=0, method="davies", weight.factor=-2, obj=NULL, obj.meta=NULL, Is.Shrink.Weight=TRUE,
n.case=-1, n.int=-1, n.all=-1, Is.Get.Naive=TRUE, Is.Get.Internal=TRUE,  Shrink.Cutoff=1, idx.error=NULL){


	#r.corr=0; method="optimal"; weight.factor=-2; weights.beta=c(1,25); Is.Shrink.Weight=TRUE;obj=obj; n.int=N.Sample; n.all=N.Sample.Total; Is.Get.Naive=TRUE;Is.Get.Internal=TRUE
	nSNP<-length(Tbl.a.list)
	IDX_MISS<-union(which(is.na(Z)),which(Z == 9))
	if(length(IDX_MISS) > 0){
		Z[IDX_MISS]<-NA
	} 
	if(length(IDX_MISS) > 0){

		msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )

		warning(msg,call.=FALSE)
		Z<-SKAT:::Impute(Z,impute.method="fixed")
	} 

	MAF<-colMeans(Z)/2
	idx<-which(MAF > 0.5)
	if(length(idx) > 0){
		MAF[idx]<-1-MAF[idx]
		Z[,idx]<-2-Z[,idx]
	}
	
	W1<-MAF * (1-MAF)
	#W1<-MAF
	if(is.null(weights)){
		W2<-dbeta(MAF, weights.beta[1], weights.beta[2])
	} else {
		W2<-weights
	}
	#W3<-rep(1, length(MAF)); W3[is.use.external==1]<-2
	W<- W1* W2

	COR<-cor(Z)
	
	##############################################
	# Run
	
	re.work<-SKAT_wExternalControl.Work(nSNP=nSNP, Tbl.a.list=Tbl.a.list, W=W, is.use.external=is.use.external, Is.Shrink.Weight=Is.Shrink.Weight, n.case=n.case, n.int=n.int, n.all=n.all, 
	weight.factor=weight.factor, Shrink.Cutoff=Shrink.Cutoff)
	
	Phi1<-t(t(COR * re.work$OR1[,2]) * re.work$OR1[,2])
	Phi2<-t(t(COR * re.work$OR1[,4]) * re.work$OR1[,4])
	Score1<-re.work$OR1[,1]
	Score2<-re.work$OR1[,3]	

	# For debug
	#OR1<<-OR1
	#Phi1<<-Phi1
	#COR<<-COR
	
	if(method=="optimal"){
		method="optimal.mod"
	}

	re.method<-SKAT:::SKAT_Check_Method(method,r.corr)

	
	re<-MetaSKAT:::Met_SKAT_Get_Pvalue(Score1, Phi1, re.method$r.corr, re.method$method, Score.Resampling=NULL)
	re$OR = re.work$OR
	re$OR1 = re.work$OR1
	re$p1 = re.work$p1
	re$p2 =re.work$p2
	re$n.SNP = length(MAF)
	re$MAC.single= cbind(re.work$MAC, MAF)
	re$Shrink.W = re.work$Shrink.W
	re$Shrink.weight = re.work$Shrink.weight
	
	re$int.p= 2
	re$int.skat =2
	re$naive.skat=2
	
	if(Is.Get.Internal){
		re2<-MetaSKAT:::Met_SKAT_Get_Pvalue(Score2, Phi2, re.method$r.corr, re.method$method, Score.Resampling=NULL)
		re$int.p = re2$p.value
		
		if(!is.null(obj.meta)){
			
			re$int.skat<-MetaSKAT_wZ(Z, obj.meta, method=method, r.corr=r.corr)$p.value
	
	
		} else if(!is.null(obj)){
			re$int.skat = SKAT(Z, obj, method=method, r.corr=r.corr)$p.value
		}
	}
	
	if(Is.Get.Naive){
	
		if(!is.null(idx.error)){
			n1<-n.int
			n2<-n.all
			Shrink.W<-rep(1,nSNP)
			Shrink.W[idx.error]<-0
			Shrink.weight<-sqrt((n2-n1)*Shrink.W + n1)
		
			OR.naive<-rbind(re.work$OR1[,5:6]) * Shrink.weight
			OR.2<-rbind(re.work$OR1[,3:4]) * Shrink.weight
			OR.naive[idx.error,]<-OR.2[idx.error,]
			
			Score3<-OR.naive[,1]	
			Phi3<-t(t(COR * OR.naive[,2]) * OR.naive[,2])	
		
		} else {
		
			Phi3<-t(t(COR * re.work$OR1[,6]) * re.work$OR1[,6])
			Score3<-re.work$OR1[,5]	
		}
		re3<-MetaSKAT:::Met_SKAT_Get_Pvalue(Score3, Phi3, re.method$r.corr, re.method$method, Score.Resampling=NULL)
		re$naive.skat<-re3$p.value 
	}
	
	return(re)

}




Get_Tbl_List<-function(dat1, n2, n2.case){

	Tbl.a.list<-list()
	n.SNP<-nrow(dat1)
	for(i in 1:n.SNP){
	
		p2<-dat1$MAF[i] * n2 *2
		p2.case<-dat1$MAF.Case[i] * n2.case *2	
		c1<-dat1$C1.minor[i]
		c2<-dat1$C2.minor[i]
	
		tbl1<-c(c1, c2)
		tbl2<-c(p2, n2*2-p2) 
		tbl.case<-c(p2.case, n2.case*2-p2.case)
		Tbl.a<-rbind(tbl.case, tbl2, tbl1)
		Tbl.a.list[[i]]<-Tbl.a
	}
	return(Tbl.a.list)
}

Get_Tbl_List_wG<-function(dat1, Z, idx.case, idx.control, n2, n2.case){

	Tbl.a.list<-list()
	n.SNP<-nrow(dat1)
	for(i in 1:n.SNP){
		if(mean(Z[,i], na.rm=TRUE) > 1){
			Z[,i]<-2 - Z[,i]
		}
	
		p2<-sum(Z[idx.control,i], na.rm = TRUE)
		p2.case<-sum(Z[idx.case,i], na.rm = TRUE)	
		
		
		c1<-dat1$C1.minor[i]
		c2<-dat1$C2.minor[i]
	
		tbl1<-c(c1, c2)
		tbl2<-c(p2, n2*2-p2) 
		tbl.case<-c(p2.case, n2.case*2-p2.case)
		Tbl.a<-rbind(tbl.case, tbl2, tbl1)
		Tbl.a.list[[i]]<-Tbl.a
	}
	return(Tbl.a.list)
}
