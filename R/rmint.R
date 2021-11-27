#' Measurement invariance test
#'
#' rmint examines across-groups equivalence of confirmatory factor
#' analysis (CFA) measurement model parameters as well as testing
#' the equality of factor means among groups. The sequence of
#' measurement invariance followed by rmint is:

#' 1) Equal form/configural invariance solution has the same form/
#' structure with the same indicators loading on the latent variables
#' for each group.

#' 2) Equal loadings/weak invariance solution constrains the loadings
#' to be equal for each group.

#' 3) Equal intercepts/strong invariance solution constrains the intercepts
#' (+ loadings) to be equal for each group.

#' 4) Equal error variances/strict invariance solution constrains the error
#' variances (+ loadings and intercepts) to be equal for each group.

#' 5) Equal factor means solution constrains the latent means (+ loadings
#' intercepts, and error variances) to be equal for each group.
#'
#' If the Chi-square difference test prøves to be significant, the equality of
#' the tested parameter is not established for the specified model.
#'
#' For the examples, we use the following model:
#' library(psych)
#' data(bfi, package = "psych")
#' library(lavaan)
#' PERS.model <- '
#' Agree =~ A1 + A2 + A3 + A4 + A5
#' Consc =~ C1 + C2 + C3 + C4 + C5
#' Extra =~ E1 + E2 + E3 + E4 + E5
#' Neuro =~ N1 + N2 + N3 + N4 + N5
#' Openn =~ O1 + O2 + O3 + O4 + O5
#' @examples
#' rmint(meas.model=PERS.model, group=gender, data=bfi)
#' rmint(meas.model=PERS.model, group=gender, data=bfi, fitmeas = FALSE)
#' rmint(meas.model=PERS.model, group=gender, data=bfi, fitmeas = TRUE)
#' @export
rmint <- function(meas.model, dataset, gvar, fitmeas = FALSE){

 library(lavaan)

 #equal form
 eqform <- sem(meas.model, dataset, group=gvar)
 eqform.fit <- fitMeasures(eqform, c("chisq", "df", "pvalue", "aic", "bic"))
 eqform.fit2 <- fitMeasures(eqform, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli"))

 #equal loadings
 eqload <- sem(meas.model, dataset, group=gvar, group.equal = c("loadings"))
 eqload.fit <- fitMeasures(eqload, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli"))
 #equal intercepts
 eqint <- sem(meas.model, dataset, group=gvar, group.equal = c("loadings", "intercepts"))
 eqint.fit <- fitMeasures(eqint, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli"))
 #equal error variance
 eqerrvar <- sem(meas.model, dataset, group=gvar, group.equal = c("loadings", "intercepts", "residuals"))
 eqerrvar.fit <- fitMeasures(eqerrvar, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli"))
 #equal factor means
 eqfmeans <- sem(meas.model, dataset, group=gvar, group.equal = c("loadings", "intercepts", "residuals", "means"))
 eqfmeans.fit <- fitMeasures(eqfmeans, c("chisq", "df", "pvalue", "rmsea", "cfi", "tli"))

 #comparisons
 # 1
 eqform.fit
 # 2 vs 1
 eqloadVSeqform <- anova(eqload, eqform)
 eqloadVSeqform
 # 3 vs 2
 eqintVSeqload <- anova(eqint, eqload)
 eqintVSeqload
 # 4 vs 3
 eqerrvarVSeqint <- anova(eqerrvar, eqint)
 eqerrvarVSeqint
 # 5 vs 4
 eqfmeansVSeqerrvar <- anova(eqfmeans, eqerrvar)
 eqfmeansVSeqerrvar
 out <- list(paste("Tests of measurement invariance and factor mean difference"),
             "Equal.Form?" = eqform.fit,
             "Equal.Loadings?" = eqloadVSeqform,
             "Equal.Intercepts?" = eqintVSeqload,
             "Equal.Error.Variances?" = eqerrvarVSeqint,
             "Equal.Factor.Means?" = eqfmeansVSeqerrvar)

 if(fitmeas == TRUE){
  return(list(out,
              paste("Fit Indices for all the models estimated:"),
              Equal.form = eqform.fit2,
              Equal.loadings = eqload.fit,
              Equal.Intercepts = eqint.fit,
              Equal.Error.Variances = eqerrvar.fit,
              Equal.Factor.Means = eqfmeans.fit))
 }
 else {
  return(out)
 }
}









