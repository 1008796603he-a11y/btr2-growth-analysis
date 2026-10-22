## ============================================================================
## FINAL growth-rate analysis — Btr2 nine-construct experiment
## ----------------------------------------------------------------------------
## This is the consolidated version of the analysis that produced the final
## nine-construct boxplot and the Tukey/CLD significance letters.
##
## Pipeline (what was actually done):
##   1. Two independent runs were combined:  Aug 27  +  Sep 19.
##   2. Per-well specific growth rate (mu, h^-1) was estimated by fitting a
##      linear regression to ln(OD600) vs time over an automatically selected
##      exponential-phase window (>=12 consecutive points, slope > 0,
##      R^2 >= 0.95; the eligible window with the HIGHEST R^2 was chosen).
##   3. OD600 was NOT blank-subtracted before fitting (unblanked OD600).
##   4. One-way ANOVA  aov(mu ~ Construct)  --> TukeyHSD --> compact-letter
##      display (multcompLetters4).  Run was NOT included as a model term.
##
## Two data-entry points are provided:
##   SECTION A  rebuilds mu from raw per-well OD time series (needs raw data).
##   SECTION B  reproduces the final stats + figures directly from the saved
##              per-well rate table (comb2_growth_rates.csv) — guaranteed to
##              match the final published numbers/letters.
## ============================================================================

library(data.table)
library(ggplot2)
library(multcompView)

## ===========================================================================
## SECTION A — estimate per-well mu from raw OD600 time series
## ---------------------------------------------------------------------------
## Run this only if starting from raw plate-reader data. It expects a long
## data.table `raw_dt` with columns:
##     Construct   (character; "blank" for blank wells)
##     Run         (character; e.g. "Aug27", "Sep19")
##     Well        (character; unique well id within a run)
##     Time_hours  (numeric; hours from start)
##     OD600       (numeric; raw, UNBLANKED)
##
## Combine rules that were applied when building `raw_dt`:
##   * B2P3 wells were taken from BOTH retained runs (Aug 27 and Sep 19; 8 wells
##     each). The Aug 15 B2P3 was a suspected wrong-plasmid clone and was
##     excluded entirely.
##   * No blank subtraction (unblanked OD600).
##   * blank wells are dropped before rate fitting.
## ===========================================================================

## --- exponential-phase window finder (as used) -----------------------------
find_exp_phase <- function(dt, min_points = 12, min_r2 = 0.95) {
  dt <- dt[order(Time_hours)]
  n  <- nrow(dt)
  if (n < min_points) return(NULL)

  res <- list(); k <- 1L
  for (start in 1:(n - min_points + 1)) {
    for (end in (start + min_points - 1):n) {
      x   <- dt[start:end]
      fit <- lm(log_OD ~ Time_hours, data = x)
      res[[k]] <- data.table(
        start_time = min(x$Time_hours),
        end_time   = max(x$Time_hours),
        n_points   = nrow(x),
        slope      = coef(fit)[["Time_hours"]],
        r_squared  = summary(fit)$r.squared
      )
      k <- k + 1L
    }
  }
  res <- rbindlist(res)
  res <- res[slope > 0 & r_squared >= min_r2]      # eligible windows only
  if (nrow(res) == 0) return(NULL)
  res[which.max(r_squared)]                         # highest-R^2 window
}

## --- fit mu per well (uncomment when raw_dt is available) -------------------
# growth_dt <- raw_dt[OD600 > 0 & !is.na(OD600) & Construct != "blank"]
# growth_dt[, log_OD := log(OD600)]
# rc <- growth_dt[, find_exp_phase(.SD), by = .(Construct, Run, Well)]
# rc <- rc[, .(Construct, Run, Well, mu = slope)]
# fwrite(rc, "comb2_growth_rates.csv")

## ===========================================================================
## SECTION B — final statistics + figures from the per-well rate table
## ===========================================================================

rc <- fread("comb2_growth_rates.csv")          # columns: Construct, Run, Well, mu

## --- display order, colours, display names ---------------------------------
ord <- c("CO7","CO13","CO13-8bpdel","B2P3","B2P3-IS10",
         "B2P3-8bpdel","B2P3-AGA","Btr2-8bpdel","Btr1")
## NOTE: names on the left (CO7, Btr1, B2P3-AGA) are the data keys used in the
## rate table; the right-hand side is the label shown on the figure.
disp <- c(CO7="2.5 kb insert control", CO13="CO13", `CO13-8bpdel`="CO13-\u03948", B2P3="B2P3",
          `B2P3-IS10`="B2P3-IS10", `B2P3-8bpdel`="B2P3-\u03948",
          `B2P3-AGA`="B2P3-T>A", `Btr2-8bpdel`="Btr2-\u03948", Btr1="5 kb insert control")
cols <- c(CO7="#848484", CO13="#3F78AB", `CO13-8bpdel`="#5DC7F7", B2P3="#9753B3",
          `B2P3-IS10`="#F58614", `B2P3-8bpdel`="#E9B914", `B2P3-AGA`="#E62C2E",
          `Btr2-8bpdel`="#94673C", Btr1="#5BB558")

## --- ANOVA + Tukey + compact letter display --------------------------------
rc[, Ctest := gsub("-", "_", Construct)]
fit   <- aov(mu ~ Ctest, data = rc)
tukey <- TukeyHSD(fit)
lett  <- multcompLetters4(fit, tukey)

letters_df <- data.frame(
  Construct = gsub("_", "-", names(lett$Ctest$Letters)),
  Letters   = unname(lett$Ctest$Letters),
  stringsAsFactors = FALSE
)

## --- summary table (mean, sd, n, letter) -----------------------------------
summ <- rc[, .(mean_mu = round(mean(mu), 4),
               sd      = round(sd(mu),   4),
               n       = .N), by = Construct]
summ <- merge(as.data.frame(summ), letters_df, by = "Construct")
summ <- summ[match(ord, summ$Construct), ]
print(summ, row.names = FALSE)
print(summary(fit))
fwrite(summ, "comb2_growth_rate_summary.csv")

## --- boxplot (ggplot; full-colour, significance letters) -------------------
rc[, disp := factor(disp[Construct], levels = disp[ord])]
rc[, fillcol := cols[Construct]]
lab <- data.frame(disp = disp[ord],
                  Letters = letters_df$Letters[match(ord, letters_df$Construct)],
                  y = 0.43)

set.seed(7)
p <- ggplot(rc, aes(disp, mu)) +
  stat_boxplot(geom = "errorbar", width = 0.28, linewidth = 0.6, colour = "black") +
  geom_boxplot(aes(fill = fillcol), colour = "black",
               outlier.shape = NA, width = 0.62, linewidth = 0.6) +
  geom_jitter(width = 0.14, height = 0, size = 1.7, alpha = 0.55, colour = "#333333") +
  geom_text(data = lab, aes(x = disp, y = y, label = Letters),
            inherit.aes = FALSE, fontface = "bold", size = 6) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 0.46), breaks = seq(0, 0.4, 0.1), expand = c(0, 0)) +
  labs(x = "Construct",
       y = expression(paste("Specific growth rate (", mu, " h"^-1, ")")),
       title = "Specific growth rates of all nine constructs") +
  theme_classic(base_size = 15) +
  theme(plot.title  = element_text(size = 17, hjust = 0, margin = margin(b = 12)),
        axis.title  = element_text(size = 15),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.line   = element_line(linewidth = 0.6, colour = "black"),
        legend.position = "none",
        plot.margin = margin(12, 16, 8, 12))
ggsave("growth_rate_boxplot_final.png", plot = p,
       width = 11.2, height = 6.0, dpi = 300, bg = "white")

cat("\nDone. Wrote comb2_growth_rate_summary.csv and growth_rate_boxplot_final.png\n")
