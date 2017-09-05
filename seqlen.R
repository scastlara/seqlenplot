#!/usr/bin/R
library(ggplot2);
library(scales);

arguments <- commandArgs(TRUE);
TRLEN     <- arguments[1];
PROTLEN   <- arguments[2];
OUTPUT    <- arguments[3]

trlength  <- read.table(file=TRLEN,   header=T, comment.char = "");
orflength <- read.table(file=PROTLEN, header=T, comment.char = "");

merged <- merge(trlength, orflength, by="Name");
colnames(merged) <- c("NAME", "TRANSCRIPT", "PROTEIN");
merged$TRANSCRIPT <- merged$TRANSCRIPT / 3;

# BREAKS
all_lengths <- c(merged$TRANSCRIPT, merged$PROTEIN);
breaks <- pretty(log10(all_lengths), n=8);
breaks2 <- as.integer(10^breaks);
RoundUp <- function(from,to) ceiling(from/to)*to
breaks2 <- RoundUp(breaks2, 10);

# Compute label placement
xlabel <- max(merged$TRANSCRIPT) + 0.05*max(merged$TRANSCRIPT)
y100 <- xlabel
y75  <- xlabel*0.75
y50  <- xlabel*0.50
y25  <- xlabel*0.25
y10  <- xlabel*0.10
minx <- min(merged$TRANSCRIPT)
miny <- min(merged$PROTEIN)
mintotal <- min(c(minx, miny)) - 0.05*min(c(minx, miny))

# Plot the thing
plot <- ggplot(merged, aes(x=TRANSCRIPT, y=PROTEIN)) +
    geom_point(alpha = 0.5, color="#EBD170") +
    xlab("\nTranscript Length (number of codons)") +
    ylab("Protein Length (number of aa)\n") +
    stat_density2d(color="#002640", alpha=0.5) +
    scale_x_log10(breaks=breaks2, limits=c(mintotal,xlabel)) + scale_y_log10(breaks=breaks2, limits=c(mintotal,xlabel)) +
    stat_function(fun=function(x)x+log10(1),    geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.75), geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.5),  geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.25), geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    stat_function(fun=function(x)x+log10(0.1),  geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
    annotate("text", label=c("100%","75%","50%","25%", "10%"), x=xlabel, y=c(y100, y75, y50, y25, y10),
             angle=45, hjust=c(1.2,1.2,1.2,1.2,1.2), vjust=c(1.1,1.1,1.1,1.1,1.1), size=3.5) +
    theme_bw()

ggsave(file=OUTPUT, plot);
