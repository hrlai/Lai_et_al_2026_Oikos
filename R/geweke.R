# Script to do Geweke diagnostics

geweke_greta <- function(draws, LV = FALSE) {
    geweke <- coda::geweke.diag(draws)

    # Ignoring LVs and loadings when looking at Geweke diagnostic
    sub_geweke <- geweke$`1`$z
    if (LV) {
        sub_geweke <- sub_geweke[-grep("loadings\\[", names(sub_geweke))]
        sub_geweke <- sub_geweke[-grep("lvs\\[", names(sub_geweke))]
    }

    table(2 * pnorm(abs(sub_geweke), lower.tail = FALSE) < 0.05)
}
