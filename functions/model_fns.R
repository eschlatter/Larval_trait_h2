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
                    prior=prior, ginverse=list(animal=Ainv),
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

#same as above function, but takes a dataframe as input; last column should be chain number
model_results_df <- function(model_df){
  #proportional variances as new mcmc objects
  num_vars <- ncol(model_df)-1
  scaled_vars <- list()
  for(i in 1:num_vars){
    varname <- paste('scale',colnames(model_df[i]),sep='_')
    scaled_vars[[varname]] <- model_df[,i]/rowSums(as.data.frame(model_df[1:num_vars])) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  plot(mcmc_areas(matrix_scaled_vars,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:num_vars){
    v_name <- colnames(model_df)[i]
    v_mcmc <- as.mcmc(model_df[,i])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- as.mcmc(scaled_vars[[i]])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model_df[1:num_vars]),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
  return(estimates)
  
}


#same as above function, but takes a list of 3 model objects as input; should be 3 repeated runs (chains) of the same model
model_results_chains <- function(model_chains){
  chain1 <- as.data.frame(model_chains[[1]]$VCV)
  chain1$chain <- 1
  chain2 <- as.data.frame(model_chains[[2]]$VCV)
  chain2$chain <- 2
  chain3 <- as.data.frame(model_chains[[3]]$VCV)
  chain3$chain <- 3
  
  model_df <- rbind(chain1,chain2,chain3)
  
  #proportional variances as new mcmc objects
  num_vars <- ncol(model_df)-1
  scaled_vars <- list()
  for(i in 1:num_vars){
    varname <- paste('scale',colnames(model_df[i]),sep='_')
    scaled_vars[[varname]] <- model_df[,i]/rowSums(as.data.frame(model_df[1:num_vars])) #store them in a list
  }
  
  #plot histograms of variance estimates
  matrix_scaled_vars <- do.call(cbind,scaled_vars) #put all variables together into a matrix for plotting
  
  #add dummy columns for dam and sire, if necessary
  if(num_vars==3){
    matrix_scaled_vars_plot <- cbind(matrix_scaled_vars,rep(0,nrow(matrix_scaled_vars)),rep(0,nrow(matrix_scaled_vars)))
    matrix_scaled_vars_plot <- matrix_scaled_vars_plot[,c(1,2,4,5,3)]
    colnames(matrix_scaled_vars_plot)[c(3,4)] <- c('scale_dam',"scale_sire")
  }
  else matrix_scaled_vars_plot <- matrix_scaled_vars
  
  plot(mcmc_areas(matrix_scaled_vars_plot,prob=0.95,area_method='scaled height')) #scaled to 1
  
  #store estimates (mean, median, upper+lower95) for each variable (original and scaled)
  estimates=data.frame()
  #original variables
  for(i in 1:num_vars){
    v_name <- colnames(model_df)[i]
    v_mcmc <- as.mcmc(model_df[,i])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  #scaled variables  
  for(i in 1:ncol(matrix_scaled_vars)){
    v_name <- colnames(matrix_scaled_vars)[i]
    v_mcmc <- as.mcmc(scaled_vars[[i]])
    v_estimates <- data.frame(variable=v_name,
                              mean=mean(v_mcmc),median=median(v_mcmc),upper95=HPDinterval(v_mcmc)[,2],lower95=HPDinterval(v_mcmc)[,1])
    estimates=rbind(estimates,v_estimates)
  }
  rownames(estimates)=c(colnames(model_df[1:num_vars]),colnames(matrix_scaled_vars))
  
  #estimates_long <- pivot_longer(estimates,mean:lower95,names_to='SummStat',values_to='Value') #in "long" format -- not sure yet which will be most useful
  print(estimates)
}

#function to run and plot the gelman-rubin diagnostic (using gelman-rubin functions from the coda package)
#input: a list of 3 model objects (should be 3 repeated runs/chains of the same model)
run_gelman_diagnostics <- function(model_chains){
  chain1 <- as.data.frame(model_chains[[1]]$VCV)
  chain1$chain <- 1
  chain2 <- as.data.frame(model_chains[[2]]$VCV)
  chain2$chain <- 2
  chain3 <- as.data.frame(model_chains[[3]]$VCV)
  chain3$chain <- 3
  
  model_chains_df <- rbind(chain1,chain2,chain3)
  mcmc_trace(model_chains_df)
  
  ## Gelman-Rubin diagnostic: looks good; upper limit is basically 1
  
  intercept_gel=mcmc.list(model_chains[[1]]$Sol,model_chains[[2]]$Sol,model_chains[[3]]$Sol)
  clutch_gel=mcmc.list(model_chains[[1]]$VCV[,'clutch'],model_chains[[2]]$VCV[,'clutch'],model_chains[[3]]$VCV[,'clutch'])
  animal_gel=mcmc.list(model_chains[[1]]$VCV[,'animal'],model_chains[[2]]$VCV[,'animal'],model_chains[[3]]$VCV[,'animal'])
  units_gel=mcmc.list(model_chains[[1]]$VCV[,'units'],model_chains[[2]]$VCV[,'units'],model_chains[[3]]$VCV[,'units'])
  if(ncol(model_chains[[1]]$VCV)==5){
    dam_gel=mcmc.list(model_chains[[1]]$VCV[,'dam'],model_chains[[2]]$VCV[,'dam'],model_chains[[3]]$VCV[,'dam'])
    sire_gel=mcmc.list(model_chains[[1]]$VCV[,'sire'],model_chains[[2]]$VCV[,'sire'],model_chains[[3]]$VCV[,'sire'])
  }
  
  #convergence diagnostic values (potential scale reduction factor):
  #"Approximate convergence is diagnosed when the upper limit is close to 1"
  #https://rdrr.io/cran/coda/man/gelman.diag.html
  print(gelman.diag(intercept_gel))
  print(gelman.diag(clutch_gel))
  print(gelman.diag(animal_gel))
  print(gelman.diag(units_gel))
  if(ncol(model_chains[[1]]$VCV)==5){
    print(gelman.diag(dam_gel))
    print(gelman.diag(sire_gel))
  }
  
  #accompanying plots
  gelman.plot(intercept_gel)
  gelman.plot(clutch_gel,main="clutch")
  gelman.plot(animal_gel,main="animal")
  gelman.plot(units_gel,main="units")
  if(ncol(model_chains[[1]]$VCV)==5){
    gelman.plot(dam_gel,main='dam')
    gelman.plot(sire_gel,main='sire')
  }
}
