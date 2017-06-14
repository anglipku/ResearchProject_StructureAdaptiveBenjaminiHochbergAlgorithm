

#####################################################################################
## download & process the fMRI data
#####################################################################################

paired_t_test = function(X1, X2){
    n = dim(X1)[1] # X1 & X2 are vectors of length n
    Diff = X1 - X2
    avg_D = colMeans(Diff)
    var_D = colSums((Diff - matrix(1, n, 1)%*%avg_D)^2) / (n-1)
    t_stat = abs(avg_D)/sqrt(var_D/n)
    p_val = (1-pt(t_stat, df = n-1))*2    
    p_val
}

fMRI_get_data_and_pvals = function(data_downloaded = FALSE){
	
	if(!require(R.matlab)){
		install.packages("R.matlab")
	}
	library("R.matlab")
	
	### Download CMU StarPlus fMRI data (subject 04847) from web
	if (!data_downloaded){
    		download.file(url = "http://www.cs.cmu.edu/afs/cs.cmu.edu/project/theo-81/www/data-starplus-04847-v7.mat", destfile = "./fMRI_data.mat")
	}
	
	### Read "meta" and "data" from data-starplus-04847-v7.mat:
	Subject = readMat("./fMRI_data.mat")
	n = Subject$meta[[3]] # number of voxels = 4698
	voxel_coords = Subject$meta[[1]] # 3D coordinates of each voxel = n-by-3 matrix
	
	### organize labels of each voxel into the Regions Of Interest (ROIs)
	ROI_names = rep("0", 24)  
	for (i in 1:24){
		ind = c(1:3,5:25)[i] # ROI #4 , called 'LIPG', is not assigned to any voxel
		ROI_names[i] = Subject$meta[[16]][[1+3*(ind-1)]]
	}
	
	ROI_voxels = list() # ROI_voxels[[i]] contains voxel numbers (1 - 4698) of each voxel in ROI #i
	ROI_labels = rep(0,n) # ROI_labels[i] is the ROI label of voxel #i
	for(i in 1:n){
		which_ROI = which(ROI_names == toString(Subject$meta[[17]][[i]][[1]]))
		if(length(which_ROI)>0){
			ROI_labels[i] = which_ROI
		}
	}

	# remove voxels which are not assigned to an ROI
	remove_voxels = which(ROI_labels==0)
	n = n - length(remove_voxels)
	ROI_labels = ROI_labels[-remove_voxels]
	voxel_coords = voxel_coords[-remove_voxels,]
	
	### Process the activity data
	# Find which trials have Picture (P) phase first, Sentence (S) phase second
	trials = which((Subject$info[1,,] > 1) & (Subject$info[14,,] == "P"))
	
	# Get the average activity levels recorded in each of the selected trials, for the P & S phases
	P_act = S_act = matrix(0, length(trials), n)
	for (i in 1:length(trials)){
		  P_act[i,] = colMeans(Subject$data[[trials[i]]][[1]][1:16,-remove_voxels]) # time 1-16 = P phase
	    S_act[i,] = colMeans(Subject$data[[trials[i]]][[1]][17:32,-remove_voxels]) # time 17-32 = S phase
	}
	### Compute p-values via a paired t-test
	pvals = matrix(paired_t_test(P_act, S_act), n, 1)  
	
	output = list()
	output$P_act = P_act
	output$S_act = S_act
	output$pvals = pvals
	output$ROI_names = ROI_names
	output$ROI_labels = ROI_labels
	output$voxel_coords = voxel_coords
	output
}


#####################################################################################
## set up plotting functions for fMRI data example
#####################################################################################

matrix_for_slice_RGB = function(pixel_coords, pixel_RGB,background_black=FALSE){
	# returns matrix to plot 64-by-64 slices of brain with RGB color specification
	# pixel_coords is n-by-2 coordinates, pixel_coords[i,] is in {1..64}-by-{1..64}
	# pixel_RGB is n-by-4 with row i specifying the RGB color for pixel i (normalized to [0,1])
	# returns in RGB format (normalized to [0,1])
	if(background_black){output = array(0,c(64,64,3))}else{output = array(1,c(64,64,3))}
	for(i in 1:dim(pixel_coords)[1]){
		output[pixel_coords[i,1],pixel_coords[i,2],] = pixel_RGB[i,]
	}
	output
}

matrix_for_slice_BW = function(pixel_coords,pixel_BW,background_black=FALSE){
	# returns matrix to plot 64-by-64 slices of brain with CMYK color specification
	# pixel_coords is n-by-2 coordinates, pixel_coords[i,] is in {1..64}-by-{1..64}
	# pixel_BW is n-by-1 with entry i specifying grayscale value for pixel i in [0,1] (0 = black)
	# returns in grayscale (black = 0)
	if(background_black){output = matrix(0,64,64)}else{output = matrix(1,64,64)}
	for(i in 1:dim(pixel_coords)[1]){
		output[pixel_coords[i,1],pixel_coords[i,2]] = pixel_BW[i]
	}
	output
}

create_fMRI_results_image = function(output,qhat,BH_Res,Storey_Res,SABHA_Res,filename){
	gap_h = 10; gap_v = 0 # horizontal and vertical gaps around each slice
	image_matrix = array(1,c(8*(64+2*gap_v),8*(64+2*gap_h),3))
	# from top to bottom: 8 slices of the brain
	# from left to right: ROIs; P data; S data; p-values; q-hat; BH; Storey; SABHA
	
	# ROIs
	set.seed(1)
	RGB_ROIs = matrix(runif(24*3),24,3)
	for(slice in 1:8){
		which_voxels = which(output$voxel_coords[,3]==slice)
		x0 = gap_h; y0 = gap_v + (8-slice)*(64+2*gap_v)
		pixel_RGB = matrix(0,length(which_voxels),3)
		for(i in 1:length(which_voxels)){
			pixel_RGB[i,] = RGB_ROIs[output$ROI_labels[which_voxels[i]],]
		}
		image_matrix_temp = matrix_for_slice_RGB(output$voxel_coords[which_voxels,1:2],pixel_RGB)
		for(i in 1:3){
			image_matrix[y0+(64:1),x0+(1:64),i] = t(image_matrix_temp[,,i])
		}
	}
	
	# P data & S data
	image_matrix[,64+2*gap_h+(1:(2*(64+2*gap_h))),] = 0 # black background
	for(slice in 1:8){
		which_voxels = which(output$voxel_coords[,3]==slice)
		x0_P = 64+2*gap_h+ gap_h; x0_S = 2*(64+2*gap_h)+ gap_h; y0 = gap_v + (8-slice)*(64+2*gap_v)
		Prange = range(colMeans(output$P_act)); Srange = range(colMeans(output$S_act))
		pixel_BW_P = (colMeans(output$P_act[,which_voxels])-Prange[1])/(Prange[2]-Prange[1])
		pixel_BW_S = (colMeans(output$S_act[,which_voxels])-Srange[1])/(Srange[2]-Srange[1])
		for(i in 1:3){
			image_matrix[y0+(64:1),x0_P+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],pixel_BW_P,background_black=TRUE))
			image_matrix[y0+(64:1),x0_S+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],pixel_BW_S,background_black=TRUE))
		}
	}
	
	# p-values and qhat
	for(slice in 1:8){
		which_voxels = which(output$voxel_coords[,3]==slice)
		x0_pvals = 3*(64+2*gap_h)+ gap_h; x0_qhat = 4*(64+2*gap_h)+ gap_h; y0 = gap_v + (8-slice)*(64+2*gap_v)
		for(i in 1:3){
			image_matrix[y0+(64:1),x0_pvals+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],output$pvals[which_voxels]))
			image_matrix[y0+(64:1),x0_qhat+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],qhat[which_voxels]))
		}
	}
	
	image_matrix[,x0_qhat+(64+gap_h)+(0:1),] = 0
	
	# BH, Storey-BH, SABHA results
	for(slice in 1:8){
		which_voxels = which(output$voxel_coords[,3]==slice)
		x0_BH = 5*(64+2*gap_h)+ gap_h; x0_Storey = 6*(64+2*gap_h)+ gap_h; x0_SABHA = 7*(64+2*gap_h)+ gap_h; y0 = gap_v + (8-slice)*(64+2*gap_v)
		for(i in 1:3){
			image_matrix[y0+(64:1),x0_BH+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],1 - BH_Res[which_voxels]))
			image_matrix[y0+(64:1),x0_Storey+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],1 - Storey_Res[which_voxels]))
			image_matrix[y0+(64:1),x0_SABHA+(1:64),i] = t(matrix_for_slice_BW(output$voxel_coords[which_voxels,1:2],1 - SABHA_Res[which_voxels]))
		}
	}
	
	col = rgb(image_matrix[,,1],image_matrix[,,2],image_matrix[,,3])
	dim(col) = dim(image_matrix[,,1])
	if(!require(grid)){
		install.packages("grid")
	}
	library(grid)
	
	pdf(filename,8,6)
	
	par(mar=c(0,0,2,0))
	plot(0:1,0:1,type='n',axes=FALSE,xlab='',ylab='',main='fMRI data: results')
	grid.raster(col, interpolate=FALSE, x=0.5, y=0.52, width=0.92, height = 0.84)
	text(0.5/8,0.01,'(a) ROIs')
	text(2/8,0.01,'(b) Data')
	text(1.5/8,0.05,'Picture')
	text(2.5/8,0.05,'Sentence')
	text(3.5/8,0.01,'(c) p-values')
	text(4.5/8,0.01,expression(paste('(d) ',widehat(q))))
	text(6.5/8,0.01,'(e) Results')
	text(5.5/8,0.05,'BH')
	text(6.5/8,0.05,'Storey-BH')
	text(7.5/8,0.05,'SABHA')
	
	dev.off()
}

create_fMRI_ROI_barplot = function(ROI_names,avg_per_ROI,filename){
	pdf(filename,10,4)
	ord = order(avg_per_ROI[,1],decreasing=TRUE)
  cols = c('grey','blue','green','orange')	
  par(mar = c (5,5,2,2))
	barplot(t(avg_per_ROI[ord,]),beside=TRUE,col=cols,names=ROI_names[ord],las=2,ylim=c(0,1),main='fMRI data: discoveries by ROI')
	legend('topright', c(expression(paste(widehat(q), " for the ROI")),"BH", "Storey-BH", "SABHA"), cex=1, bty = 'n', fill=cols)
	dev.off()
}

	
	
#####################################################################################
## run fMRI data example
#####################################################################################

alpha = 0.2 # target FDR level
tau = 0.5; eps = 0.1 # parameters for SABHA
thr = 0.5 # parameter for Storey-BH
ADMM_params = c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr

output = fMRI_get_data_and_pvals()
# if data is already downloaded:
# output = fMRI_get_data_and_pvals(data_downloaded=TRUE)
n = length(output$pvals)

### set up the three methods
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

source('All_q_est_functions.R')

# gather results
BH_Res = Storey_Res = SABHA_Res = rep(0,n)
BH_Res[BH_method(output$pvals, alpha)] = 1
Storey_Res[Storey_method(output$pvals, alpha, thr)] = 1
qhat = Solve_q_block(output$pvals,tau,eps,output$ROI_labels,ADMM_params)
SABHA_Res[SABHA_method(output$pvals,qhat,alpha,tau)] = 1

num_discoveries = c(sum(BH_Res),sum(Storey_Res),sum(SABHA_Res))
names(num_discoveries) = c('BH','Storey-BH','SABHA')
print(num_discoveries)

# plot all results
create_fMRI_results_image(output,qhat,BH_Res,Storey_Res,SABHA_Res,'fMRI_results.pdf')


avg_per_ROI = matrix(0,24,4) 
# avg_per_ROI[i,1] = avg qhat value in ROI i
# avg_per_ROI[i,2:4] = proportion of discoveries made in ROI i by each method
# (2 = BH, 3 = Storey-BH, 4 = SABHA)
for(i in 1:24){
	avg_per_ROI[i,1] = mean(qhat[which(output$ROI_labels==i)])
	avg_per_ROI[i,2] = mean(BH_Res[which(output$ROI_labels==i)])
	avg_per_ROI[i,3] = mean(Storey_Res[which(output$ROI_labels==i)])
	avg_per_ROI[i,4] = mean(SABHA_Res[which(output$ROI_labels==i)])
}

# plot number of discoveries per ROI
create_fMRI_ROI_barplot(output$ROI_names,avg_per_ROI,'fMRI_ROI_barplot.pdf')
