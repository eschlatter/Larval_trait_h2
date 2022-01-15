#parameter-expanded prior
prior_expand_params <- list(V = 1,nu = 1000, alpha.mu=0, alpha.V=1)

#also parameter-expanded, from tuto_en; with residual edited (nu=0.002 instead of nu=1) to be proper
prior_ext_params <- list(V=1,nu=1,alpha.mu=0,alpha.V=1000)

#inverse gamma prior
prior_invgamma_params <- list(V = 1,nu = .002)

#gentler inverse gamma prior
prior_invgamma2_params <- list(V = 1,nu = .02)

#flat prior
prior_flat_params <- list(V = 1e-10,nu = -1) #trusting based on course notes that this is actually flat, because I can't tell what a negative nu means