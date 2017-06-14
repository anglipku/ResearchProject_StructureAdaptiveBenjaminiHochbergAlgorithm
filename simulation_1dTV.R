### Setup
# We first define some functions for running the simulation and for plotting the results.

Setup_simulation_functions = function(){
# PSeq_SABHA: returns the simulated hypotheses, with same dim as true_q
# column 1 contains the vector of p-values
# column 2 contains the labels of nonnulls (1) and nulls (0)
# inputs: true_q (prob to be null for each hypothesis i), mu (signal strength, a scalar)
PSeq_SABHA <<- function(true_q, mu, seed){
  
    set.seed(seed)
    n = length(true_q)
    labels = (runif(n) > true_q) # 0 means null, 1 means nonnull
    Pvals = rep(0,n)
    Pvals[which(labels == 0)] = runif(sum(labels == 0))
    z_nonnull = rnorm(sum(labels == 1), mean = mu, sd = 1)
    Pvals[which(labels == 1)] = 2*(1-pnorm(abs(z_nonnull)))
    
    list(Pvals, labels)
}

# CompEffectFunc_SABHA_mu: returns the average power and observed FDR for BH, Storey, SABHA with true q, SABHA with estimated q, under different signal strengths (mu)
# thr: threshold in Storey's method
# tau, eps: parameters for SABHA
# alpha: target FDR level
CompEffectFunc_SABHA_mu <<- function(true_q, mu_seq, mu_ind, seeds, alpha, ADMM_params, tau, eps, thr, TV_bd_seq){
    numSimu=length(seeds)
    Power = FDR = matrix(0,length(mu_seq),6) 
	
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
		pvals[pvals>tau] = Inf
		khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
	  which(qhat*pvals<=alpha*khat/length(pvals))
	}
	
	SABHA_priors = matrix(0, length(true_q), length(TV_bd_seq))

    
    for (i in 1:numSimu){
        for(k in 1:length(mu_seq)){
            Pvals = PSeq_SABHA(true_q, mu_seq[k], seeds[i])          
            
            for(j in 1:6){
            	# methods: 1 BH, 2 Storey-BH, 3,4,5 SABHA, 6 oracle SABHA
                if(j == 1){ # BH
                    rej = BH_method(Pvals[[1]], alpha)
                } else if (j == 2) { # Storey-BH
                    rej = Storey_method(Pvals[[1]], alpha, thr)
                    if(i==1 & k==mu_ind){est_proportion_nulls = min(1,mean(Pvals[[1]]>thr)/(1-thr))}
                } else if (j == 3) { # SABHA with estimated q
					est_q = Solve_q_TV_1dim(Pvals[[1]], tau, eps, TV_bd_seq[1], ADMM_params)
                    rej = SABHA_method(Pvals[[1]], est_q, alpha, tau)
                    if(i==1 & k==mu_ind){SABHA_priors[,1] = est_q}
                } else if (j == 4) { # SABHA with estimated q
					est_q = Solve_q_TV_1dim(Pvals[[1]], tau, eps, TV_bd_seq[2], ADMM_params)
                    rej = SABHA_method(Pvals[[1]], est_q, alpha, tau)
                    if(i==1 & k==mu_ind){SABHA_priors[,2] = est_q}
                } else if (j == 5) { # SABHA with estimated q
					est_q = Solve_q_TV_1dim(Pvals[[1]], tau, eps, TV_bd_seq[3], ADMM_params)
                    rej = SABHA_method(Pvals[[1]], est_q, alpha, tau)
                    if(i==1 & k==mu_ind){SABHA_priors[,3] = est_q}
                } else { # SABHA with true q
                    rej = SABHA_method(Pvals[[1]], true_q, alpha, tau)
                }
                
                Power[k, j] = Power[k, j] + length(intersect(rej, which(Pvals[[2]]==1)))/sum(Pvals[[2]]==1)
                FDR[k, j] = FDR[k, j] + length(setdiff(rej, which(Pvals[[2]]==1)))/max(length(rej), 1)
        
            }
        }
        print(paste0("Trial ", i, " finished"))
    }
    Power=Power/numSimu
    FDR=FDR/numSimu
    return(list(Power, FDR, est_proportion_nulls, SABHA_priors))
    
    # Power & FDR are averaged over all trials
    # est_proportion_nulls & SABHA_priors are for trial 1 only
}

print('Simulation functions defined.')
}


Setup_simulation_plotting_functions = function(){
	
# plot FDR and power results for each method across range of signal strengths mu
FDR_and_power_plot <<- function(FDRmat, powermat, alpha, mu_seq, filename){
	
	pdf(filename, 10, 12.6)
	par(mar=c(5,0,2,1))
	
	names = c('BH', 'Storey-BH', paste0('SABHA (m=', TV_bd_seq[1], ')'), paste0('SABHA (m=', TV_bd_seq[2], ')'), paste0('SABHA (m=', TV_bd_seq[3], ')'), 'oracle SABHA')
	cols = c('blue', 'green', 'orange', 'orange', 'orange', 'red'); ltys = c(1, 1, 1, 2, 3, 1);  pchs = c(15, 17,  19, 19, 19,  18)
	
	xzero=min(mu_seq)-.05*(max(mu_seq)-min(mu_seq))
	xzero1=min(mu_seq)-.17*(max(mu_seq)-min(mu_seq))
	xzero2=min(mu_seq)-.07*(max(mu_seq)-min(mu_seq))
	xzero3=min(mu_seq)-.1*(max(mu_seq)-min(mu_seq))
	plot(0:1,0:1,type='n',xlab='', ylab='',xlim=c(xzero1,max(mu_seq)),ylim=c(-.25,1.05),axes=FALSE,main='Simulated data: results',cex.main=2.5)
	mtext(expression(paste('Signal Strength ', mu[sig])), side = 1, line = 3.5, cex=2, at = 2)
  segments(xzero2,0,max(mu_seq),0)
	axis(side=1,at=mu_seq,cex.axis=2)
  segments(xzero,-.25,xzero,-.05)
	segments(xzero,.05,xzero,1.05)
	for(i in 0:2){ # set up y axis for FDR
		  segments(xzero,-.25+i/10,xzero2,-.25+i/10)
		  text(xzero3,-.25+i/10,i/10,cex=2)
	}
	for(i in 0:10){ # set up y axis for power
		  segments(xzero,.05+i/10,xzero2,0.05+i/10)
		  text(xzero3,.05+i/10,i/10,cex=2)
	}
	text(xzero1,-.125,'FDR',srt=90,cex=2)
	text(xzero1,.55,'Power',srt=90,cex=2)
	for(i in 1:length(names)){
      lines(mu_seq,FDRmat[,i]-.25, lty = ltys[i], col=cols[i], lwd = 1.5)
		  lines(mu_seq,powermat[,i]+.05, lty = ltys[i], col=cols[i], lwd = 1.5)
  		points(mu_seq,FDRmat[,i]-.25,pch=pchs[i],col=cols[i],cex=1.5)
		  points(mu_seq,powermat[,i]+.05,pch=pchs[i],col=cols[i],cex=1.5)		
	}
    
  points(mu_seq,rep(alpha-.25, length(mu_seq)),type='l',lty='dotted',col='gray50',lwd=2) # target FDR
  legend(2,0.38,legend=names,col=cols, pch=pchs,lty=ltys,cex=2,pt.cex=1.5,lwd=1.5)
  
  dev.off()
}

# Plot the estimated priors (i.e. q-hat) for each method
plot_Estimated_priors <<- function(true_q, est_proportion_nulls, SABHA_priors, filename){
	
	pdf(filename, 7.14, 9)
	par(mfcol=c(5,1),mar=c(5,4,2,1))
	
	names = c('BH', 'Storey-BH', paste0('SABHA (m=', TV_bd_seq[1], ')'), paste0('SABHA (m=', TV_bd_seq[2], ')'), paste0('SABHA (m=', TV_bd_seq[3], ')'))
	
    n = length(true_q)
    Est_q = list()
    Est_q[[1]] = rep(1,n)
    Est_q[[2]] = est_proportion_nulls*rep(1,n)
    Est_q[[3]] = SABHA_priors[,1]
    Est_q[[4]] = SABHA_priors[,2]
    Est_q[[5]] = SABHA_priors[,3]
    
    title1 = 'True vs estimated probability of a null'
    xlab1 = 'Index i'
    for(i in 1:5){
    	  if(i==1){title=title1}else{title=''}
    	  if(i==5){xlab=xlab1}else{xlab=''}
    
	      plot(0:1,0:1,type='n',xlab = xlab, ylab = '', main = title, axes=FALSE, xlim = c(1, n), ylim = c(0, 1), cex.main=2.5, cex.lab=2)
    	  axis(side=1, cex.axis=2)
	      axis(side=2, cex.axis=2, las=1, at=c(0,0.5,1))
	      lines(1:n, true_q, lty = 3) # true q
    	  lines(1:n, Est_q[[i]], col='blue')
		    if(i==1){
			      legend(1,.95,legend=c('True probability',names[i]),col=c('black','blue'),lty=c(3,1),cex=2)
		    }else{
			      legend(1,.95,legend=names[i],col='blue',lty=1,cex=2)
		    }
	  }  
	  dev.off()
}


print('Simulation plotting functions defined.')
}


Setup_simulation_functions()
Setup_simulation_plotting_functions()
source('All_q_est_functions.R')


### Simulation

n = 500 # number of hypotheses
true_q = c(rep(0.1, n/2), rep(0.9, n/2)) # true_q[i] = prob that hypothesis i is null
TV_bd_seq = c(0.25, 2, 20) # value of "m" for the SABHA total variation norm bound

ntrials = 10; seeds = 1:ntrials;
mu_seq = seq(0.5, 3.5, 0.5) # signal strength: range
mu_ind = 5 # signal strength mu_seq[mu_ind]: choose one for comparing q-hat's
alpha = 0.1 # target FDR level
alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
tau = 0.5; eps = 0.1 # parameters for SABHA total variation method
thr = 0.5 # tau parameter for Storey-BH method

# run the full experiment over the range of signal strength values
SABHA_simulation_Res = CompEffectFunc_SABHA_mu(true_q, mu_seq, mu_ind, seeds, alpha, ADMM_params, tau, eps, thr, TV_bd_seq)

### Save the results
save(SABHA_simulation_Res, file = "simulation_1dTV.RData")

# If previously computed result is available
# load("simulation_1dTV.RData")

### Plots
# First figure: compare methods among signal strength (mu). For each method and each mu, we record the actual FDR attained, and the power to detect signals, averaged over all trials.
# Second figure: how the methods perform as TV_bd varies

FDR_and_power_plot(SABHA_simulation_Res[[2]], SABHA_simulation_Res[[1]], alpha, mu_seq, 'simulation_1dTV_FDR_power.pdf')

plot_Estimated_priors(true_q, SABHA_simulation_Res[[3]], SABHA_simulation_Res[[4]], 'simulation_1dTV_est_q.pdf')


