model_run <- function(param_list,response,model_statement,n_vars,n_mcmc_iter){
  #structure for plugging in the params to generate a full prior -- this is what MCMCglmm takes
  G <- list()
  for(i in 1:n_vars){
    varname <- paste("G",i,sep="")
    G[[varname]] <- param_list
  }
  if(length(param_list)==2) R=param_list else R=list(V = 1,nu = .002) #set prior for residual
  prior <- list(R=R, #variance prior for residual
                G=G)  #variance priors for random effects
  
  #run the model
  model <- MCMCglmm(formula(paste0(response,'~','1+age')),
                    random=formula(paste0('~',model_statement)),
                    prior=prior,
                    data=phens,nitt=n_mcmc_iter,burnin=round(.01*n_mcmc_iter),thin=10,verbose=TRUE)
  #I think the "missing values in random predictors" warning is just because the parents in phens don't have predictor values
  return(model)
}


#function that takes an mcmc model object (output of MCMCglmm) and looks at the results
model_diagnostics <- function(model){
  #trace and density plots
  plot(model)
  
  #proportional variances as new mcmc objects
  scaled_vars <- list()
  for(i in 1:ncol(model$VCV)){
    varname <- paste('scale',colnames(model$VCV)[i],sep='_')
    scaled_vars[[varname]] <- model$VCV[,i]/rowSums(as.data.frame(model$VCV)) #store them in a list
  }
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  
  #store convergence info by variable (original and scaled)
  diagnostics=data.frame()
  #original variables
  for(i in 1:ncol(model$VCV)){
    v_name <- colnames(model$VCV)[i]
    v_mcmc <- model$VCV[,i]
    v_diagnostics <- data.frame(variable=v_name,
                                effSamp=effectiveSize(v_mcmc),
                                firstAutocorr=autocorr(v_mcmc)[2],lastAutocorr=autocorr(v_mcmc)[nrow(autocorr(v_mcmc))])
    diagnostics=rbind(diagnostics,v_diagnostics)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- scaled_vars[[i]]
    v_diagnostics <- data.frame(variable=v_name,
                                effSamp=effectiveSize(v_mcmc),
                                firstAutocorr=autocorr(v_mcmc)[2],lastAutocorr=autocorr(v_mcmc)[nrow(autocorr(v_mcmc))])
    diagnostics=rbind(diagnostics,v_diagnostics)
  }
  rownames(diagnostics)=c(colnames(model$VCV),colnames(matrix_scaled_vars))
  print(diagnostics)
  return(diagnostics)
}


model_results <- function(model){
  #proportional variances as new mcmc objects
  scaled_vars <- list()
  for(i in 1:ncol(model$VCV)){
    varname <- paste('scale',colnames(model$VCV)[i],sep='_')
    scaled_vars[[varname]] <- model$VCV[,i]/rowSums(as.data.frame(model$VCV)) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  plot(mcmc_areas(matrix_scaled_vars,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:ncol(model$VCV)){
    v_name <- colnames(model$VCV)[i]
    v_mcmc <- model$VCV[,i]
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- scaled_vars[[i]]
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model$VCV),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
  return(estimates)
  
}