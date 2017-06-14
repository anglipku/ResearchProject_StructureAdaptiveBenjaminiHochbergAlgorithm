#####################################################################################
## download & process the gene/drug data
#####################################################################################

gene_drug_get_data = function(data_downloaded = FALSE){
	if(data_downloaded){Data = read.table('gene_drug_data.txt')}else{
		if(!require(GEOquery)){
			install.packages("GEOquery")
		}
		library("GEOquery")
		gds2324 = getGEO('GDS2324')
		eset = GDS2eSet(gds2324, do.log2=TRUE)
		Data = exprs(eset)
		write.table(Data,'gene_drug_data.txt')
	}
	Data
}

gene_drug_get_pvals = function(Data){
	
	# code for processing this data set & for accumulation tests
	# the code in this function is from the paper:
	# Ang Li and Rina Foygel Barber, "Accumulation tests for FDR control in ordered hypothesis testing" (2015). Available from http://arxiv.org/abs/1505.07352
	
	# the following code is copied (with some edits) from:
	# http://www.stat.uchicago.edu/~rina/accumulationtests/gene_dosage_experiment.R
	
	##### Two-sample t-tests on a data matrix
	# This function inputs a data matrix X with n columns. The columns indexed by g_1 and g_2 belong to sample 1 and sample 2, respectively (e.g. cases and controls). For each row of X, the function runs a two-sample t-test for that row.
	ttest_mat=function(X,g1,g2){
		# g1 & g2 give columns for groups 1 & 2
		n1=length(g1);n2=length(g2)
		means1=rowMeans(X[,g1]);means2=rowMeans(X[,g2])
		vars1=rowSums((X[,g1]-means1%*%t(rep(1,n1)))^2)/(n1-1);vars2=rowSums((X[,g2]-means2%*%t(rep(1,n2)))^2)/(n2-1)
		sds_diff = sqrt(vars1/n1 + vars2/n2)
		tstats = (means1-means2)/sds_diff
		dfs = (vars1/n1+vars2/n2)^2 / ((vars1/n1)^2/(n1-1) + (vars2/n2)^2/(n2-1))
		pvals=2*(1-pt(abs(tstats),dfs))
		output=list()
		output$tstats=tstats
		output$dfs= dfs
		output$pvals= pvals
		output
	}

	# The next function is the same, except that instead of performing two-sided t-tests, we perform a one-sided t-test for each row of X. The signs s_i specify, for the i-th row of X, whether we are testing for a positive effect or a negative effect.
	signed_ttest_mat=function(X,g1,g2,s){
		n1=length(g1);n2=length(g2)
		means1=rowMeans(X[,g1]);means2=rowMeans(X[,g2])
		vars1=rowSums((X[,g1]-means1%*%t(rep(1,n1)))^2)/(n1-1);vars2=rowSums((X[,g2]-means2%*%t(rep(1,n2)))^2)/(n2-1)
		sds_diff = sqrt(vars1/n1 + vars2/n2)
		tstats = s*(means1 - means2)/sds_diff
		dfs = (vars1/n1+vars2/n2)^2 / ((vars1/n1)^2/(n1-1) + (vars2/n2)^2/(n2-1))
		pvals=(1-pt(tstats,dfs))
		output=list()
		output$tstats=tstats
		output$dfs= dfs
		output$pvals= pvals
		output
	}

	# Each row of "Data" is a gene (total: n=22283 genes). 
	# The 25 columns of "Data" correspond to 5 trials each at 5 different dosages: columns 1-5 are zero dose (control group), followed by columns 6-10 at the lowest dose, etc. The entries of "Data" are log-transformed gene expression levels.
	n=dim(Data)[1]
	highdose=21:25;lowdose=6:10;control=1:5
	
	### Computing p-values
	# Next, we will use the highest-dosage data to produce an ordering on the low-dosage vs control-group p-values, resulting in an ordered sequence of p-values which we will use for the accumulation tests.
	# First, for each gene, produce a pval for highest dose compared to mean of lowest dose + control. Record also the sign of the effect (increased or reduced gene expression at the highest dose, compared to pooled low-dose and control-group data).
	ttest_highdose=ttest_mat(Data,highdose,c(lowdose, control))
	pvals_highdose=ttest_highdose$pvals
	test_signs=sign(ttest_highdose$tstats)
	
	pvals_highdose_small=ttest_mat(Data,highdose[1:2],c(lowdose,control))$pvals
	set.seed(1);pvals_random=runif(n)
	
	# Next, for each gene we will perform a one-sided t-test (using the sign of the estimated high-dose effect), and then use a permutation test to get a final p-value for this gene. These p-values will then be reordered according to the high-dose results above.
	signed_pvals_lowdose_ttest=signed_ttest_mat(Data,lowdose,control,test_signs)$pvals
	
	signed_pvals_lowdose_ttest_permuted=matrix(0,choose(10,5),n)
	nchoosek=combn(1:10,5)
	for(i in 1:choose(10,5)){
		permord=c(nchoosek[,i],setdiff(1:10,nchoosek[,i]))
		signed_pvals_lowdose_ttest_permuted[i,]=signed_ttest_mat(Data[,permord],lowdose,control,test_signs)$pvals
	}
	signed_pvals_permutation_test=colSums(abs(signed_pvals_lowdose_ttest_permuted)-rep(1,choose(10,5))%*%t(abs(signed_pvals_lowdose_ttest))<=0)/choose(10,5)
	
	output = list()
	output$pvals = signed_pvals_permutation_test*(1-1/(1+choose(10,5))) # multiplying to avoid pvalues exactly equal to 1 due to discretization
	output$ord = order(abs(pvals_highdose))
	output$ord_small = order(abs(pvals_highdose_small))
	output$ord_random = order(abs(pvals_random))
	output
	
}


#####################################################################################
## set up plotting functions for gene/drug data example
#####################################################################################

gene_drug_plot_results = function(NumRej,filename){
	### Plot results for each method
	cols=c(rep('black',4), 'blue', 'red','orange')
	ltys=c(1,5,3,4,rep(1,2),1)
	pchs=c(20,15,17,19,4,1,15)
	methods=c('SeqStep', 'Accum. Test (HingeExp)', 'ForwardStop', 'Adaptive SeqStep', 'BH', 'Storey', 'SABHA')

	pdf(filename,width=12.5,height=4)
	par(mfrow=c(1,3),oma=c(0,0,3,0))
	titles=c('Highly informative ordering', 'Moderately informative ordering', 'Random ordering')
	for(j in 1:3){
		plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,8000),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries', main = titles[j], axes=FALSE, font.main=1, cex.main = 1.5, cex.lab = 1.5)
		axis(side=1,at=0:5/10, cex.axis = 1.5)
		axis(side=2, cex.axis = 1.5)
		alpha_pt=(1:(10*max_alpha))*10
		for(i in length(methods):1){
			points(alphalist,NumRej[i,,j],type='l',col=cols[i],lty=ltys[i])
			points(alphalist[alpha_pt],NumRej[i,alpha_pt,j],col=cols[i],pch=pchs[i])
		}
		if(j==1){
			legend(0,8000,methods,col=cols,lty=ltys,pch=pchs,seg.len=3,cex=1.2)
	  }
  }   
  mtext("Gene/drug response data: results",outer=TRUE,cex=1.3,font=2)
	box(); dev.off()
}

#####################################################################################
## run gene/drug data example
#####################################################################################

alphalist = seq(0.01,0.5,by=0.01) # target FDR level
tau = 0.5; eps = 0.1 # parameters for SABHA
thr = 0.5 # parameter for Storey-BH
thr1 = 0.1; thr2 = 0.5 # parameters for adaptive SeqStep

Data = gene_drug_get_data()
# if data already downloaded & saved:
# Data = gene_drug_get_data(data_downloaded = TRUE)



output = gene_drug_get_pvals(Data)

### set up all methods
source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')

BH_method = function(pvals,alpha){
	khat=max(c(0,which(sort(pvals)<=alpha*(1:length(pvals))/length(pvals))))
	which(pvals<=alpha*khat/length(pvals))
}

Storey_method = function(pvals,alpha,thr){
	est_proportion_nulls=min(1,mean(pvals>thr)/(1-thr))
	pvals[pvals>thr] = Inf
	khat=max(c(0,which(sort(pvals)<=alpha/est_proportion_nulls*(1:length(pvals))/length(pvals))))
	which(pvals<=alpha/est_proportion_nulls*khat/length(pvals))
}
	
	
SABHA_method = function(pvals, qhat, alpha, tau){
	# Use the original, or estimated q as input
	pvals[pvals>tau] = Inf
	khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
    which(qhat*pvals<=alpha*khat/length(pvals))
}

Adaptive_SeqStep_method = function(pvals, alpha, thr1, thr2){ # Lei & Fithian 2016's method
	# thr1 & thr2 correspond to s & lambda (with s<=lambda) in their paper
	fdphat = thr1 / (1-thr2) * (1 + cumsum(pvals>thr2)) / pmax(1, cumsum(pvals<=thr1))
	if(any(fdphat<=alpha)){
		khat = max(which(fdphat<=alpha))
		return(which(pvals[1:khat]<=thr1))
	}else{
		return(NULL)
	}
	
}

source('All_q_est_functions.R')
num_alpha=length(alphalist); max_alpha = max(alphalist)

# gather results
NumRej = array(0,c(7,num_alpha,3))
# methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 SABHA
for(j in 1:3){
	if(j==1){pvals=output$pvals[output$ord]}else{if(j==2){pvals=output$pvals[output$ord_small]}else{pvals=output$pvals[output$ord_random]}}
	qhat = Solve_q_step(pvals,tau,eps)
	for(i in 1:num_alpha){
		NumRej[1,i,j] = SeqStep(pvals,alpha=alphalist[i],C=2)
		NumRej[2,i,j] = HingeExp(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i],C=2)
		NumRej[3,i,j] = ForwardStop(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i])
		NumRej[4,i,j] = length(Adaptive_SeqStep_method(pvals, alphalist[i], thr1, thr2))
		NumRej[5,i,j] = length(BH_method(pvals, alphalist[i]))
		NumRej[6,i,j] = length(Storey_method(pvals, alphalist[i], thr))
		NumRej[7,i,j] = length(SABHA_method(pvals, qhat, alphalist[i], tau))
	}
}

gene_drug_plot_results(NumRej,'gene_drug_results.pdf')

