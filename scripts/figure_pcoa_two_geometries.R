## figure_pcoa_two_geometries.R
##
## Ordination of the nine matched pairs under two distance measures --
## Bray-Curtis and Aitchison (CLR-based) -- with the PERMANOVA result for
## each printed on the panel. If the site effect only showed up under one
## geometry it would suggest an artefact of the distance choice; here it
## doesn't. Computed directly from the spot table.

library(data.table)
library(ggplot2)
library(patchwork)

CSV <- "/Users/akhaliq/Desktop/rebuttle_natgen/new/pdac_spot_level.csv.gz"
md <- as.data.table(read.csv(gzfile(CSV), stringsAsFactors = FALSE))
setnames(md, "ecotype", "cc")
md[, site := c("Pancreas" = "PT", "Liver" = "HM")[site]]

patients <- c("PT_2", "PT_3", "PT_4", "PT_6", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12")
classes  <- c("CC1", "CC2", "CC3", "CC4", "CC5", "CC7")

s <- md[patient %in% patients & site %in% c("PT","HM") & cc %in% classes]
counts <- s[, .N, by = .(patient, site, cc)]
counts[, prop := N / sum(N), by = .(patient, site)]
wide <- dcast(counts, patient + site ~ cc, value.var = "prop", fill = 0)
wide <- wide[order(match(patient, patients), factor(site, levels = c("HM","PT")))]
X <- as.matrix(wide[, ..classes])
grp <- wide$site

clr <- function(X) {
  Z <- X; Z[Z == 0] <- min(Z[Z > 0]) / 2; Z <- Z / rowSums(Z)
  t(apply(log(Z), 1, function(r) r - mean(r)))
}

exact_permanova <- function(X, method) {
  D <- if (method == "bray") as.matrix(dist(X, method = "manhattan")) * 0 else NULL # placeholder unused
  D <- if (method == "bray") {
    n <- nrow(X)
    d <- matrix(0, n, n)
    for (i in 1:n) for (j in 1:n) d[i,j] <- sum(abs(X[i,]-X[j,])) / sum(X[i,]+X[j,])
    d
  } else as.matrix(dist(clr(X)))
  n <- nrow(D); A <- -0.5*D^2; C <- diag(n)-1/n; G <- C %*% A %*% C
  pseudoF <- function(g) {
    SST <- sum(diag(G)); tot <- 0
    for (l in unique(g)) { i <- which(g==l); tot <- tot + sum(G[i,i])/length(i) }
    SSA <- tot - sum(G)/n
    c(F=(SSA/1)/((SST-SSA)/(n-2)), R2=SSA/SST)
  }
  grp0 <- grp; obs <- pseudoF(grp0); hits <- 0
  for (i in 0:511) {
    bits <- as.integer(intToBits(i))[1:9]; g <- grp0
    for (k in which(bits==1)) { a<-2*k-1; b<-2*k; tmp<-g[a]; g[a]<-g[b]; g[b]<-tmp }
    if (pseudoF(g)["F"] >= obs["F"] - 1e-8) hits <- hits + 1
  }
  c(obs, P = hits/512)
}

res_bray <- exact_permanova(X, "bray")
res_ait  <- exact_permanova(X, "aitchison")

pcoa_plot <- function(X, method, res, title) {
  D <- if (method == "bray") {
    n <- nrow(X); d <- matrix(0,n,n)
    for (i in 1:n) for (j in 1:n) d[i,j] <- sum(abs(X[i,]-X[j,])) / sum(X[i,]+X[j,])
    as.dist(d)
  } else dist(clr(X))
  mds <- cmdscale(D, k = 2, eig = TRUE)
  var_explained <- round(100 * mds$eig[1:2] / sum(abs(mds$eig)), 1)
  df <- data.table(x = mds$points[,1], y = mds$points[,2], site = grp, patient = wide$patient)
  ggplot(df, aes(x, y, colour = site)) +
    geom_line(aes(group = patient), colour = "grey70") +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c("PT" = "#1F3864", "HM" = "#C0504D"),
                         labels = c("PT" = "primary", "HM" = "hepatic metastasis"),
                         name = NULL) +
    labs(title = title,
         subtitle = sprintf("F = %.2f, R\u00b2 = %.3f, P = %.4f", res["F"], res["R2"], res["P"]),
         x = sprintf("axis 1 (%.1f%%)", var_explained[1]),
         y = sprintf("axis 2 (%.1f%%)", var_explained[2])) +
    theme_minimal(base_size = 10)
}

p1 <- pcoa_plot(X, "bray", res_bray, "Bray-Curtis")
p2 <- pcoa_plot(X, "aitchison", res_ait, "Aitchison")

combined <- (p1 | p2) + plot_layout(guides = "collect") +
  plot_annotation(title = "The site effect doesn't depend on the distance measure") &
  theme(legend.position = "bottom")

ggsave("figure_pcoa_two_geometries.png", combined, width = 9.5, height = 4.6, dpi = 300)

cat("Bray-Curtis: F =", res_bray["F"], "R2 =", res_bray["R2"], "P =", res_bray["P"], "\n")
cat("Aitchison:   F =", res_ait["F"],  "R2 =", res_ait["R2"],  "P =", res_ait["P"], "\n")
