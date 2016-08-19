if (!isGeneric("plot"))
      setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

if (!isGeneric("coef"))
      setGeneric("coef", function(object, ...) standardGeneric("coef"))

if (!isGeneric("vcov"))
      setGeneric("vcov", function(object, ...) standardGeneric("vcov"))

if (!isGeneric("nobs"))
      setGeneric("nobs", function(object, ...) standardGeneric("nobs"))

if (!isGeneric("logLik"))
      setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

if (!isGeneric("deviance"))
      setGeneric("deviance", function(object, ...) standardGeneric("deviance"))

if (!isGeneric("AIC"))
      setGeneric("AIC", function(object, ..., k=2) standardGeneric("AIC"))

if (!isGeneric("pairs"))
      setGeneric("pairs", function(x, ...) standardGeneric("pairs"))