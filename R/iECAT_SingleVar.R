##########################################
# Get OR

Conti_Correction_Tbl.a<-function(Tbl.a){
	
	Tbl.a[Tbl.a < 0]<-0
	
	Tbl1<-Tbl.a[1:2,]
	if(length(which(Tbl1==0)) > 0){
		Tbl1 = Tbl1 + 0.5
		Tbl.a[1:2,]<-Tbl1
	}	
	
	if(length(which(Tbl.a[3,]==0)) > 0){
		Tbl.a[3,] = Tbl.a[3,]  + 0.5
	}	
	
	return(Tbl.a)

}

GetOR<-function(Tbl){
	
	Tbl[Tbl < 0]<-0
	if(length(which(Tbl==0)) > 0){
		Tbl = Tbl + 0.5
	}
	
	OR1<-log(Tbl[1,1]*Tbl[2,2]/(Tbl[1,2]*Tbl[2,1]))
	se1<-sqrt(sum(1/Tbl))
	pval<-pchisq(OR1^2/se1^2, df=1, lower.tail = FALSE)
	return(list(OR=OR1, se=se1, pval=pval))

}


GetOR_Tbl3<-function(Tbl){

	if(length(which(Tbl==0)) > 0){
		Tbl = Tbl + 0.5
	}

	a1<-Tbl[1,1]
	b1<-Tbl[2,1]
	n1<-sum(Tbl[1,])
	n2<-sum(Tbl[2,])

	#a1<-a[1]; b1<-b[1]
	p1<-a1/n1
	p2<-b1/n2
	r<-n2/n1
	
	OR1<-log(a1) + log(n1+n2 -a1 -b1) -log(n1-a1) - log(a1 +b1)
	
	S1<-p1 * (1-p1) / n1
	S2<-p2 * (1-p2) / n2
	Sigma<-diag(c(S1,S2))

	f1<-1/p1 + 1/(1-p1) -1/(1+r -p1 -r*p2) -1/(p1+r*p2)
	f2<- -1/(1+1/r -p1/r - p2) -1/(p1/r + p2)
	
	f.all<-c(f1, f2)
	
	re.currect<-t(f.all) %*% Sigma %*% f.all
	#re.misspecified<-1/a1 + 1/(a1 +b1) + 1/(n1-a1) + 1/(n1+n2 -a1-b1)
	
	se1= sqrt(re.currect)
	pval<-pchisq(OR1^2/se1^2, df=1, lower.tail = FALSE)
	
	re<-list(OR=OR1, se = se1, pval=pval)
	return(re)
}


Get_Corr1<-function(ncase, ncontrol, n1){
	
	ncasek = ncase
	ncontrolk = ncontrol + n1
	nk = ncasek+ncontrolk
	
	ncasel = ncase
	ncontroll = ncontrol
	nl = ncasel + ncontroll
	
	ncasekl = ncase
	ncontrolkl = ncontrol
	
	Corr<-(ncontrolkl * sqrt(ncasek * ncasel / (ncontrolk* ncontroll)) +  ncasekl * sqrt(ncontrolk* ncontroll  / (ncasek * ncasel)))/ sqrt(nk * nl)
	
	return(Corr)

}



Get_Corr2<-function(ncase, ncontrol, n1){
	
	ncasek = ncase
	ncontrolk = ncontrol
	nk = ncasek+ncontrolk
	
	ncasel = n1
	ncontroll = ncontrol
	nl = ncasel + ncontroll
	
	ncasekl = 0
	ncontrolkl = ncontrol
	
	Corr<-(ncontrolkl * sqrt(ncasek * ncasel / (ncontrolk* ncontroll)) +  ncasekl * sqrt(ncontrolk* ncontroll  / (ncasek * ncasel)))/ sqrt(nk * nl)
	
	return(Corr)

}

Get_ORW<-function(gamma, b1, b2, b3, SE1, SE2, SE3, ncase, ncontrol, nexternal){

	Rmat<-diag(1, 3)
	Rmat[1,2]<-Rmat[2,1]<-Get_Corr1(ncase, ncontrol, nexternal)
	Rmat[2,3]<-Rmat[3,2]<-Get_Corr2(ncase, ncontrol, nexternal)
	
	SE<-c(SE1, SE2, SE3)
	
	Sigma<-t(t(Rmat * SE) * SE) 
	tau<-gamma * b3^2 / (gamma *b3^2 + SE3^2)
	tau_d<-(2* gamma * b3 * SE3^2)/(gamma *b3^2 + SE3^2)^2
	
	g1<-(1-tau)
	g2<-tau
	g3<-(1-tau_d) * b1 + tau_d * b2

	g.a<-c(g1, g2, g3)
	
	b.W<-(1-tau) * b1 + tau*b2 
	Var.W<-t(g.a) %*% Sigma %*% (g.a)
	
	re = list(OR=b.W, SE=sqrt(Var.W), tau=tau )
	return(re)

}


Estimate_OROnly<-function(n1, n2, n3, a1, b1, c1, gamma){
	
	#a1<-a[1]; b1<-b[1]; c1<-c[1]; gamma<-1/4
	out_bd<-EstimateVar_B1(n2, n3, b1, c1) # Get Corr between b and d
	out_cd<-EstimateVar_B1(n3, n2, c1, b1) # Get Corr between c and d
	
	ahat<-log(a1) - log(n1-a1)
	bhat<-log(b1) - log(n2-b1)
	chat<-log(c1) - log(n3-c1)
	dhat<-log(b1 + c1) - log(n2+n3-b1-c1)
		
	Var_a<-1/a1 + 1/(n1-a1)
	Var_b<-1/b1 + 1/(n2-b1)
	Var_c<-1/c1 + 1/(n3-c1)
	Var_d<-1/(b1+c1) + 1/(n2+n3-b1 -c1)
	
	Var_Be<-Var_a + Var_d
	Var_Bint<-Var_a + Var_b
	Var_Bi<-Var_b + Var_c
	
	Cov_Be_Bint<-Var_a + out_bd$Cov_bd
	Cov_Be_Bi<- - out_bd$Cov_bd + out_cd$Cov_bd
	Cov_Bint_Bi<-- Var_b 
	
	Sigma<-diag(c(Var_Be, Var_Bint, Var_Bi))
	Sigma[1,2] = Sigma[2,1] = Cov_Be_Bint
	Sigma[1,3] = Sigma[3,1] = Cov_Be_Bi
	Sigma[2,3] = Sigma[3,2] = Cov_Bint_Bi
	SE3<-sqrt(Var_Bi)
	
	
	
	b1<-ahat -dhat
	b2<-ahat - bhat
	b3<-bhat - chat

	
	tau<-gamma * b3^2 / (gamma *b3^2 + SE3^2)
	tau_d<-(2* gamma * b3 * SE3^2)/(gamma *b3^2 + SE3^2)^2
	
	g1<-(1-tau)
	g2<-tau
	g3<-(1-tau_d) * b1 + tau_d * b2

	g.a<-c(g1, g2, g3)
	
	b.W<-(1-tau) * b1 + tau*b2 
	Var.W<-t(g.a) %*% Sigma %*% (g.a)
	
	re = list(OR=b.W, SE=sqrt(Var.W), Sigma=Sigma, tau=tau)

	return(re)
	
}


WithCorr_Get_Var_Overlap<-function(out, X, G, Y, idx.share1, idx.share2){

	# Y<-Y4; G<-G4; X<-X4; out<-glm(Y4 	~ X4+G4, family=binomial); 
	#idx.share1<-1:ncontrol; idx.share2<-1:ncontrol + ncontrol
	
	n.y<-length(Y)
	X1<-cbind(rep(1, n.y), X, G)
	mu<-out$fitted.values
	I<- t(X1) %*% (X1 * mu *(1-mu))
	U<-X1 *(Y - mu) 
	
	idx.share<-c(idx.share1, idx.share2)
	
	D<-matrix(rep(0, ncol(X1)*ncol(X1)), ncol=ncol(X1)) 
	for(i in 1:ncol(X1)){
		for(j in 1:ncol(X1)){
			
			D[i,j]<-sum(U[idx.share1,i] * U[idx.share2,j])*2
		}
	}

	I1<-I + D
	V<- solve(I) %*% I1 %*% solve(I)
	V1<- solve(I) 
	return(list(Var=V, V1=V1, D=D, I=I))
	
}


WithCorr_Get_Cov_Overlap<-function(out1, X1, G1, Y1, out4, X4, G4, Y4, ncase, ncontrol, n1){


	# out1<-glm(Y1 	~ X1+G1, family=binomial); out4<-glm(Y4 	~ X4+G4, family=binomial); 
	 idx.1A<-(ncase+1):(ncase+ncontrol);idx.1B<-(ncase+ncontrol+1):(ncase+ncontrol+n1)
	 idx.2A<-1:ncontrol; idx.2A1<-(ncontrol+1):(2*ncontrol);idx.2B<-(2*ncontrol+1):(2*ncontrol + n1) 
	
	Y<-Y1
	n.y<-length(Y)
	X1.1<-cbind(rep(1, n.y), X1, G1)
	mu<-out1$fitted.values
	I1<- t(X1.1) %*% (X1.1 * mu *(1-mu))
	U1<-X1.1 *(Y1 - mu) 

	Y<-Y4
	n.y<-length(Y)
	X4.1<-cbind(rep(1, n.y), X4, G4)
	mu<-out4$fitted.values
	I2<- t(X4.1) %*% (X4.1 * mu *(1-mu))
	U2<-X4.1 *(Y - mu) 
	
	D<-matrix(rep(0, ncol(X1.1)*ncol(X1.1)), ncol=ncol(X1.1)) 
	for(i in 1:ncol(X1.1)){
		for(j in 1:ncol(X1.1)){
			
			temp1<-sum(U1[idx.1A,i] * U2[idx.2A,j])
			temp2<-sum(U1[idx.1A,i] * U2[idx.2A1,j])
			temp3<-sum(U1[idx.1B,i] * U2[idx.2B,j])
						
			D[i,j]<-temp1+temp2+temp3
		}
	}

	Cov1<-solve(I1) %*% D %*% solve(I2)

	return(list(Cov=Cov1, I1=I1, I2=I2, D=D))
	
}





EstimateVar_B1<-function(n2, n3, b1, c1){

	#b1<-b[1]; c1<-c[1]
	Var_b<-1/b1 + 1/(n2-b1)
	Var_d<-1/(b1+c1) + 1/(n2+n3-b1 -c1)
	
	
	p1<-b1/n2
	p2<-c1/n3
	r<-n3/n2
	
	S1<-p1 * (1-p1) / n2
	S2<-p2 * (1-p2) / n3
	Sigma<-diag(c(S1,S2))

	f1<-1/p1 + 1/(1-p1) -1/(1+r -p1 -r*p2) -1/(p1+r*p2)
	f2<- -1/(1+1/r -p1/r - p2) -1/(p1/r + p2)
	
	f.all<-c(f1, f2)
	
	Var_bd<-t(f.all) %*% Sigma %*% f.all
	Cov_bd<-(Var_b + Var_d - Var_bd)/2
	re<-list(Var_bd=Var_bd, Cov_bd=Cov_bd)
	return(re)
}

EstimateVar<-function(n.a, a.a, va1 =NULL){
	
	a1<-a.a[1]; b1<-a.a[2]; c1<-a.a[3]
	n1<-n.a[1]; n2<-n.a[2]; n3<-n.a[3]
	
	out_b1<-EstimateVar_B1(n2, n3, b1, c1)

	Var_a<-1/a1 + 1/(n1-a1)
	Var_d<-1/(b1+c1) + 1/(n2+n3-b1 -c1)
	
	Var_Be<-Var_a + Var_d
	Var_B1<-out_b1$Var_bd
	Cov_Be1<-Var_d - out_b1$Cov_bd
	
	if(is.null(va1)){
		Var_Bw=NULL
	} else {
		Var_Bw<-Var_Be + va1^2 *Var_B1 - 2 * va1 * Cov_Be1	
	}
	re<-list(Var_Bw=Var_Bw, Var_Be=Var_Be, Var_B1=Var_B1, Cov_Be1=Cov_Be1, Var_bd = out_b1$Var_bd, Cov_bd = out_b1$Cov_bd)
	return(re)
	
}

Shrinkage_TestStat_Work<-function(gamma_param, b1, b2, b3, SE1, SE2, SE3, Pval3, Cov_b1_b3, Var_Adj.Factor, Cutoff, Pval.Cutoff){


	if( Pval3 > Pval.Cutoff){
		tau2=0
	} else {
		tau2 = b3^2 * gamma_param
	}
	
	# Get Tau
	tau= tau2 /(SE3^2  + tau2)
	
	# Variance calculation
	va1<-tau2*(tau2 + 3* SE3^2)/(SE3^2 + tau2)^2


	
	bw1 =   b1*(1-tau) + b2* tau
	bw =  b1 - b3^3 * gamma_param / (b3^2 * gamma_param + SE3^2)
	
	Var_Bw<-SE1^2 + va1^2 *SE3^2 - 2 * va1 * Cov_b1_b3	
	Va<-Var_Bw * Var_Adj.Factor

	# Use cutoff for numerical reason
	if(1-tau < Cutoff || Pval3 < 0.05){
		Va = SE2^2
		bw = b2
		tau = 1
	}		

	#out_print1<-c(tau2,tau, b3^2, SE3^2, va1, bw, bw1, Var_Bw, Va, gamma_param)
	#cat("[*3][", out_print1, "]\n")
	
	return(list(OR_new=bw, Se_new=sqrt(Va),  tau=tau, temp1=1-tau))

}

Shrinkage_TestStat_WithCorr<-function(Y, X, G, gamma_param, ncase, ncontrol, nexternal, Cutoff=0, Pval.Cutoff=0.5, Var_Adj.Factor=1){


	n1<-nexternal
	n<-ncase+ncontrol
	idx1<-1:(n+n1)
	idx2<-1:n
	idx3<-c((ncase+1):n, c((ncase+1):(n+n1)))

	idx.share1<-1:ncontrol; idx.share2<-1:ncontrol + ncontrol

	Y1<-Y; G1<-G;
	Y2<-Y[idx2]; G2<-G[idx2];
	Y3<-c(rep(1,ncontrol), rep(0, n1+ncontrol)); G3<-G[idx3]; X3<-X[idx3]

	if(is.null(X)){

		out1<-glm(Y1 	~ G1, family=binomial)
		out2<-glm(Y2 	~ G2, family=binomial)
		out3<-glm(Y3 	~ G3, family=binomial)	
		X1<-NULL; X2<-NULL; X3<-NULL;


		
	} else {
		X1<-X
		X2<-X[idx2,]
		X3<-X[idx3,]
	
		out1<-glm(Y1 	~ G1+X1, family=binomial)
		out2<-glm(Y2 	~ G2+X2, family=binomial)
		out3<-glm(Y3 	~ G3+X3, family=binomial)
	}
	
	Var3<-WithCorr_Get_Var_Overlap(out3, NULL, G3, Y3, idx.share1, idx.share2)
	Cov<-WithCorr_Get_Cov_Overlap(out1, X1, G1, Y1, out3, X3, G3, Y3, ncase, ncontrol, n1)

	Var_b3<-Var3$Var[2,2]
	Cov_b1_b3<-Cov$Cov[2,2]
	
	b1<-summary(out1)$coefficients[2,1]
	b2<-summary(out2)$coefficients[2,1]
	b3<-summary(out3)$coefficients[2,1]

	SE1<-summary(out1)$coefficients[2,2]
	SE2<-summary(out2)$coefficients[2,2]
	SE3<-sqrt(Var_b3)	
	Pval3<-pchisq((b3/SE3)^2, df=1, lower.tail=FALSE)
	
	#out1.1<<-out1; out2.1<<-out2; out3.1<<-out3
	#out_print1<-c(b1, b2, b3, SE1, SE2, SE3, Pval3, Var_b3, Cov_b1_b3)
	#cat("[*2][", out_print1, "]\n")
	
	re<-Shrinkage_TestStat_Work(gamma_param=gamma_param, b1=b1, b2=b2, b3=b3, SE1=SE1, SE2=SE2, SE3=SE3, Pval3=Pval3
	, Cov_b1_b3=Cov_b1_b3, Var_Adj.Factor=Var_Adj.Factor, Cutoff=Cutoff, Pval.Cutoff=Pval.Cutoff)
	
	
	return(re)

}



Shrinkage_TestStat<-function(Tbl.a, gamma_param, Cutoff=0, Pval.Cutoff, Var_Adj.Factor=1){
	
	# Tbl1 : Internal control only
	# Tbl2 : Internal + External
	# Tbl3 : Interal control vs External control (or Internal vs Augmented (Internal + External))
	
	#weight.factor=1; Cutoff=0.2; Pval.Cutoff=0.05
	Tbl.a=Conti_Correction_Tbl.a(Tbl.a)

	Tbl1<-Tbl.a[1:2,]
	Tbl2<-Tbl1
	Tbl2[2,]<-Tbl.a[2,] + Tbl.a[3,]
	Tbl3<-Tbl.a[2:3,]

	a.a<-Tbl.a[,1]
	n.a<-rowSums(Tbl.a)
	
	OR1<-GetOR(Tbl1)
	OR2<-GetOR(Tbl2)	
	OR3<-GetOR_Tbl3(Tbl3)	# For GetOR_Tbl3 function, do not need to collapse 2 and 3 column
	
	b1<-OR2$OR
	b2<-OR1$OR
	b3<-OR3$OR
	SE1<-OR2$se
	SE2<-OR1$se
	SE3<-OR3$se	

	Pval3=OR3$pval
	
	Var_out=EstimateVar(n.a, a.a, va1 =NULL)
	Var_b3 = Var_out$Var_B1
	Cov_b1_b3 = Var_out$Cov_Be1
	
	out_print1<-c(b1, b2, b3, SE1, SE2, SE3, Pval3, Var_b3, Cov_b1_b3)
	#cat("[*1][", out_print1, "]\n")

	re<-Shrinkage_TestStat_Work(gamma_param=gamma_param, b1=b1, b2=b2, b3=b3, SE1=SE1, SE2=SE2, SE3=SE3, Pval3=Pval3,
	Cov_b1_b3=Cov_b1_b3, Var_Adj.Factor=Var_Adj.Factor, Cutoff=Cutoff, Pval.Cutoff=Pval.Cutoff)
	
	
	return(re)
	
}




Get_Gamma_Var_Param_OLD<-function(weight.factor, Tbl.a){

	Var_Adj.Factor=1
	Pval.Cutoff=1
	if(weight.factor==-2){
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param=1/8
		} else {
			gamma_param=1/2
		}
		Var_Adj.Factor=1.02^2
	} else if(weight.factor==-3 || weight.factor==-4){
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param=1/5
		} else {
			gamma_param=1
		}
	} else if(weight.factor== -5 ){
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param=1/2
		} else {
			gamma_param=1
		}
	} else if(weight.factor== -6 ){
		gamma_param = 1
	} else if(weight.factor== -7 ){
		gamma_param = 1
		Pval.Cutoff=0.317 #1-pchisq(1, df=1)
	} else if(weight.factor== -8 ){		# hybrid of -7 and -9
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param= 0
		} else {
			gamma_param=1
			Pval.Cutoff=0.317
		}
	} else if(weight.factor== -9 ){
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param=0
		} else {
			gamma_param=1
		}
	}  	else {
		gamma_param=1/weight.factor
	}
	
	
	re<-list(gamma_param=gamma_param, Var_Adj.Factor=Var_Adj.Factor, Pval.Cutoff=Pval.Cutoff)
	return(re)
}



Get_Gamma_Var_Param<-function(weight.factor, Tbl.a){

	Var_Adj.Factor=1
	Pval.Cutoff=1
	if(weight.factor== -9 ){
		MAF1<-Tbl.a[,1] / rowSums(Tbl.a)
		if(abs(MAF1[1] - MAF1[2]) > abs(MAF1[1] - MAF1[3]) && sign(MAF1[1] - MAF1[2]) == sign(MAF1[1] - MAF1[3])){
			gamma_param=0
		} else {
			gamma_param=1
		}
	}  	else {
		gamma_param=1/weight.factor
	}
	
	
	re<-list(gamma_param=gamma_param, Var_Adj.Factor=Var_Adj.Factor, Pval.Cutoff=Pval.Cutoff)
	return(re)
}


Get_Bayes_Pval<-function(Tbl.a, Y=NULL, X=NULL, G=NULL, weight.factor=1, Cutoff=0.4, Fisher.test=FALSE){

	#weight.factor=1; CutOff=0.4;Fisher.test=FALSE; 
	#
	
	# change later
	ncase<-sum(Tbl.a[1,]) /2
	ncontrol<-sum(Tbl.a[2,]) /2 
	nexternal<-sum(Tbl.a[3,]) /2
	
	Tbl.a=Conti_Correction_Tbl.a(Tbl.a)
	
	Tbl1<-Tbl.a[1:2,]
	Tbl2<-Tbl1
	Tbl2[2,]<-Tbl1[2,] + Tbl.a[3,]
	Tbl3<-Tbl.a[2:3,]
	
	###################################
	# Get OR and fisher test
	

	if(! Fisher.test){	
		
		OR1<-GetOR(Tbl1)
		OR2<-GetOR(Tbl2)
		OR3<-GetOR(Tbl3)
	
		Pval.internal=OR1$pval
		Pval.naive=OR2$pval
		Pval.ivse=OR3$pval
	
	} else {  
		Pval.internal=fisher.test(Tbl1)$p.value
		Pval.naive=fisher.test(Tbl2)$p.value	
		Pval.ivse=fisher.test(Tbl3)$p.value
	} 
	
	################################
	# Use
	re.gamma = Get_Gamma_Var_Param(weight.factor=weight.factor, Tbl.a=Tbl.a)
 	
 	if(is.null(Y)){
 		out.b<-Shrinkage_TestStat(Tbl.a, gamma_param= re.gamma$gamma_param, Cutoff=Cutoff, Pval.Cutoff=re.gamma$Pval.Cutoff, Var_Adj.Factor=re.gamma$Var_Adj.Factor)
 	
 	} else {
		
		out.b<-Shrinkage_TestStat_WithCorr(Y, X, G, gamma_param=re.gamma$gamma_param, ncase=ncase, ncontrol=ncontrol, nexternal=nexternal
		, Cutoff=Cutoff, Pval.Cutoff=re.gamma$Pval.Cutoff, Var_Adj.Factor=re.gamma$Var_Adj.Factor)
 	
 	}
 	
	
	OR_new.a<-c(out.b$OR_new, out.b$Se_new)
	OR_new.std1<-OR_new.a[1] / OR_new.a[2]
	
	p1<-pchisq(OR_new.std1^2, df=1, lower.tail=FALSE)
	

	re<-list(OR=OR_new.a[1], Pval.asymptotic=p1, Pval.naive=Pval.naive, Pval.internal=Pval.internal , Pval.ivse=Pval.ivse,  W=out.b$temp1, SE=OR_new.a[2], Var_Adj.Factor=out.b$Var_Adj.Factor)
	return(re) 
	
	
}


#
#	For internal use only 
#
QQPlot<-function(pval, main, MAC=NULL, MAC.cut=10, xlab="-log10 expected", ylab="-log10 observed"){

	#pval=OUT[,2]; main="New.t";MAC=NULL;MAC.cut=10
	idx<-which(pval < 10^-10)
	pval[idx] =10^-10
	
	idx1<-which(pval <= 1)
	if(!is.null(MAC)){
		idx2<-which(MAC > MAC.cut)
		idx1<-intersect(idx1, idx2)
	}
	
	pval<-pval[idx1]
	n1<-length(pval)
	max1<-max(c(-log10((1:n1)/(n1+1))),-log10(pval))
	
	qqplot( -log10((1:n1)/(n1+1)),-log10(pval), xlab=xlab, ylab=ylab, main=main, xlim=c(0,max1), ylim=c(0, max1))
	abline(0,1)
	
}

#
#	For internal use only 
#

PTable<-function(PVal, alpha=c(0.05,0.01, 0.005, 0.001,0.0005, 0.0001,0.00005, 0.00001 ,0.000005, 0.0000025, 0.000001), onlyCount=FALSE, noInfo=FALSE){
	
	#PVal<-out1; alpha=c(0.05, 10^-3, 2.5*10^-6);onlyCount=FALSE
	#PVal<-t(OUT_Single[,1:n1*3 -2]); alpha=alpha.a; onlyCount=TRUE;noInfo=TRUE
	p<-dim(PVal)[1]
	n<-dim(PVal)[2]
	#
	out<-matrix(rep(0,n*length(alpha)),ncol=n)

	for(i in 1:length(alpha)){
		for(j in 1:n){
			if(!onlyCount){
				out[i,j]<-length(which(PVal[,j] < alpha[i])) / p 
			} else {
				out[i,j]<-length(which(PVal[,j] < alpha[i]))
			}

		}
	}
	if(!is.null(colnames(PVal))){
		colnames(out)<-colnames(PVal)
	}
	if(onlyCount && noInfo==FALSE){
		out<-cbind(out, p)
	}
	if(noInfo){
		return(out)
	} 
	
	return(cbind(alpha, out))
	

}


	
	
ScoreTest_Get_X1 = function(X1){
  
  qr1<-qr(X1)
  q1<-ncol(X1)
  if(qr1$rank < q1){
    
    X1.svd<-svd(X1)
    X1 = X1.svd$u	#Ques: Why?
  } 
  
  return(X1)
  
}


ScoreTest_NULL_Model = function(formula, data ){
  
  mod = lm(formula, data=data)
  X1<-model.matrix(formula,data=data)
  X1<-ScoreTest_Get_X1(X1)
  
  glmfit= glm(formula, data=data, family = "binomial")
  mu = glmfit$fitted.values
  
  V = mu*(1-mu)
  res = glmfit$y- mu
  n1<-length(res)
  
  #
  XV = t(X1 * V)
  XVX_inv= solve(t(X1)%*%(X1 * V))
  XXVX_inv= X1 %*% XVX_inv   
  
  re<-list(y=glmfit$y, mu=mu, res=res, V=V, X1=X1, XV=XV, XXVX_inv =XXVX_inv, glmfit=glmfit)
  return(re)
  
}


SPA_ER_pval <- function(tempdat, G, q, stat.qtemp, mu, g, Cutoff=2) {
	if (stat.qtemp > Cutoff^2) {
		if (sum(G)<10){ #ER
			obj <- SKAT_Null_Model(y~., out_type="D", Adjustment=FALSE, data=tempdat)
			temp_binary=SKATBinary(as.matrix(G),obj, method.bin="Hybrid")
			#p_temp=c(temp_binary$p.value,temp_binary$p.value.resampling)
			p_temp <- temp_binary$p.value
		} else { #MAC>=10, no resampling
			p_temp0 = SPAtest:::Saddle_Prob(q, mu=mu, g=g, Cutoff=2,alpha=5*10^-8)$p.value
			if (p_temp0!=0){p_temp<-p_temp0} else{p_temp<- -2}
		}
	updated.pval <- p_temp} else {
	updated.pval <- -1	
	} #end of if (stat.qtemp > Cutoff^2)
	return(updated.pval)
}
