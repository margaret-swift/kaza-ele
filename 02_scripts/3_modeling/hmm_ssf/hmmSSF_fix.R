.hmmSSF <- function (
    ssf_formula, tpm_formula = ~1, n_states, data, ssf_par0, 
    tpm_par0 = NULL, optim_opts = list(trace = 0, maxit = 50000), 
    method = "Nelder-Mead") 
{
  message('hmmSSF function [Maggie]')
  data <- data[order(data$ID, data$stratum, -data$obs), ]
  
  obs <- subset(data, obs == 1)
  options(na.action = "na.pass")
  ssf_MM <- model.matrix(ssf_formula, data)
  ssf_MM <- ssf_MM[, !colnames(ssf_MM) == "(Intercept)"]
  options(na.action = "na.pass")
  tpm_MM <- model.matrix(tpm_formula, obs)
  if (is.null(tpm_par0)) {
    tpm_par0 <- matrix(0, nrow = ncol(tpm_MM), ncol = n_states * 
                         (n_states - 1))
    tpm_par0[1, ] <- -2
  }
  par <- c(ssf_par0, tpm_par0)
  message(par)
  sampling_densities <- data$w
  fit <- optim(par = par, fn = nllk, ssf_MM = ssf_MM, tpm_MM = tpm_MM, 
               sampling_densities = sampling_densities, stratum = data$stratum, 
               ID = data$ID, n_states = n_states, n_obs = nrow(obs), 
               control = optim_opts, hessian = TRUE, method = method)
  par <- format_par(par = fit$par, n_states = n_states, ssf_cov = colnames(ssf_MM), 
                    tpm_cov = colnames(tpm_MM))
  args <- list(tpm_formula = tpm_formula, ssf_formula = ssf_formula, 
               data = data, n_states = n_states, ssf_cov = colnames(ssf_MM), 
               tpm_cov = colnames(tpm_MM))
  mod <- list(par = par, fit = fit, args = args)
  class(mod) <- append("hmmSSF", class(mod))
  return(mod)
}
