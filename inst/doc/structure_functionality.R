## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 12,
  fig.height = 6)

## ---- echo=FALSE--------------------------------------------------------------
fam_name <- c("unit-Weibull", "Kumaraswamy", "unit-Logistic", "unit-Chen",
              "unit-Birnbaum-Saunders", "log-extended Exponential-Geometric",
              "unit-Generalized Half-Normal-E", "unit-Generalized Half-Normal-X",
              "unit-Gompertz", "unit-Burr-XII", "Johnson-SB",
              "arc-secant hyperbolic Weibull", "unit-Gumbel")
abbrev_name <- c("uweibull", "kum", "ulogistic", "uchen", "ubs", "leeg", "ughne",
                 "ughnx", "ugompertz", "uburrxii", "johnsonsb", "ashw", "ugumbel")
refs <- c("[Mazucheli, et al. (2018)](http://japs.isoss.net/13(2)1%2011046.pdf)",
          "[Kumaraswamy, (1980)](https://www.sciencedirect.com/science/article/abs/pii/0022169480900360)",
          "[Tadikamalla and Johnson (1982)](https://doi.org/10.2307/2335422)",
          "[Korkmaz, et al. (2020)](https://doi.org/10.1515/ms-2022-0052)",
          "[Mazucheli, et al. (2021)](https://www.mdpi.com/2073-8994/13/4/682)",
          "[Jodrá and Jiménez-Gamero (2020)](https://doi.org/10.57805/revstat.v18i4.309)",
          "[Korkmaz MÇ (2020)](https://www.scientificbulletin.upb.ro/rev_docs_arhiva/full6b9_464742.pdf)",
          "New",
          "[Mazucheli et al. (2019)](https://rivista-statistica.unibo.it/article/view/8497)",
          "[Korkmaz and Chesneau (2021)](https://link.springer.com/article/10.1007/s40314-021-01418-5)",
          "[Johnson (1949)](https://doi.org/10.2307/2332539)",
          "[Korkmaz et al. (2021)](https://www.tandfonline.com/doi/full/10.1080/02664763.2021.1981834)",
          "New")
tab <- data.frame(fam_name, abbrev_name, refs)
knitr::kable(tab[order(tab$abbrev_name), ], col.names = c("Family", "Abbreviation", "Reference"),
             caption = "Available families of distributions their abbreviations and reference.",
             label = "distributions", row.names = FALSE)

## ----structure, echo=FALSE----------------------------------------------------
library(unitquantreg)
args(unitquantreg)

## -----------------------------------------------------------------------------
unlist(unitquantreg.control())

## ----methods-unitquantreg-----------------------------------------------------
methods(class = "unitquantreg")

## ----methods-unitquantregs----------------------------------------------------
methods(class = "unitquantregs")

## ----water-data---------------------------------------------------------------
data(water)
head(water)

## ----fitting------------------------------------------------------------------
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
lt_fits <- lapply(lt_families, function(fam) {
  unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
               data = water, tau = 0.5, family = fam, link = "logit",
               link.theta = "log")
})
t(sapply(lt_fits, coef))

## ----like-stats---------------------------------------------------------------
likelihood_stats(lt = lt_fits)

## ---- vuong-tests-------------------------------------------------------------
lt_chosen <- lt_fits[c("unit-Logistic", "Johnson-SB", "unit-Burr-XII", "unit-Weibull")]
pairwise.vuong.test(lt = lt_chosen)

## ----plots-diagnostic, cache=TRUE---------------------------------------------
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
plot(lt_fits[["unit-Logistic"]])
par(oldpar)

## ----fits-ulogistic, cache=TRUE-----------------------------------------------
system.time(
  fits_ulogistic <- unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                                 data = water, tau = 1:49/50,
                                 family = "ulogistic", link = "logit",
                                 link.theta = "log"))

## ----plots-hnp, cache=TRUE----------------------------------------------------
library(ggplot2)
get_data <- function(obj) {
  tmp <- hnp(obj, halfnormal = FALSE, plot = FALSE, nsim = 10)
  tmp <- as.data.frame(do.call("cbind", tmp))
  tmp$tau <- as.character(obj$tau)
  tmp
}
chosen_taus <- c("0.02", "0.5", "0.98")
df_plot <- do.call("rbind", lapply(fits_ulogistic[chosen_taus], get_data))
df_plot$tau <- paste0(expression(tau), " == ", df_plot$tau)

ggplot(df_plot, aes(x = teo, y = obs)) +
  facet_wrap(~tau, labeller = label_parsed) +
  geom_point(shape = 3, size = 1.4) +
  geom_line(aes(y = median), linetype = "dashed") +
  geom_line(aes(y = lower), col = "#0080ff") +
  geom_line(aes(y = upper), col = "#0080ff") +
  theme_bw() +
  labs(x = "Theoretical quantiles", y = "Randomized quantile residuals") +
  scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
  scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
  theme_bw() +
  theme(text = element_text(size = 16, family = "Palatino"), 
        panel.grid.minor = element_blank())

## ----summary-fits-------------------------------------------------------------
summary(lt_fits[["unit-Logistic"]])

## ----plot-ulogistic-----------------------------------------------------------
plot(fits_ulogistic, which = "coef")

## ----fits-uweibull, cache=TRUE------------------------------------------------
system.time(
  fits_uweibull <- unitquantreg(formula = phpws ~ mhdi + incpc + region + log(pop),
                                 data = water, tau = 1:49/50,
                                 family = "uweibull", link = "logit",
                                 link.theta = "log"))
plot(fits_uweibull, which = "coef")

## ----plot-conddis-------------------------------------------------------------
lt_data <- list(mhdi = c(0.5, 0.7), incpc = round(mean(water$incpc)),
                region = c(1, 0), pop = round(mean(water$pop)))
plot(fits_ulogistic, which = "conddist", at_obs = lt_data, at_avg = FALSE,
     dist_type = "density")
plot(fits_ulogistic, which = "conddist", at_obs = lt_data, at_avg = FALSE,
     dist_type = "cdf")

## ----seesion-info-------------------------------------------------------------
sessionInfo()

