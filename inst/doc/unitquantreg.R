## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----distr-families-----------------------------------------------------------
lt_families <- list("unit-Weibull" = "uweibull",
                    "Kumaraswamy" = "kum",
                    "unit-Logistic" = "ulogistic",
                    "unit-Birnbaum-Saunders" = "ubs",
                    "log-extended Exponential-Geometric" = "leeg",
                    "unit-Chen" = "uchen",
                    "unit-Generalized Half-Normal-E" = "ughne",
                    "unit-Generalized Half-Normal-X" = "ughnx",
                    "unit-Gompertz" = "ugompertz",
                    "Johnson-SB" = "johnsonsb",
                    "unit-Burr-XII" = "uburrxii",
                    "arc-secant hyperbolic Weibull" = "ashw",
                    "unit-Gumbel" = "ugumbel")

## ----example------------------------------------------------------------------
library(unitquantreg)
data(water)
lt_fits <- lapply(lt_families, function(fam) {
  unitquantreg(
    formula = phpws ~ mhdi + incpc + region + log(pop), data = water, tau = 0.5,
    family = fam, link = "logit", link.theta = "log")
})
t(sapply(lt_fits, coef))

## ----likelihood-stats---------------------------------------------------------
likelihood_stats(lt = lt_fits)

## ----pairwise-vuong-----------------------------------------------------------
# Select just a few model to not mess the output
lt_chosen <- lt_fits[c("unit-Logistic", "Johnson-SB", "unit-Burr-XII", "unit-Weibull")]
pairwise.vuong.test(lt = lt_chosen)

## ----methods------------------------------------------------------------------
methods(class = "unitquantreg")

## ----various-quantiles--------------------------------------------------------
fits <- unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                     data = water, tau = 1:49/50, family = "uweibull",
                     link = "logit", link.theta = "log")
class(fits)

## ----plot-unitquantregs-------------------------------------------------------
plot(fits, which = "coef", mean_effect = FALSE)

## -----------------------------------------------------------------------------
methods(class = "unitquantregs")

