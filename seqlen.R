#!/usr/bin/R
library(ggplot2);
library(scales);

read_cmdlineargs <- function() {
    arguments <- commandArgs(TRUE);
    TRLEN     <- arguments[1];
    PROTLEN   <- arguments[2];
    OUTPUT    <- arguments[3];
    XLAB      <- "Transcript Length (bp)";
    YLAB      <- "Protein Length (codons)";
    if (length(arguments) > 3) {
        XLAB <- arguments[4];
    }
    if (length(arguments) > 4) {
        YLAB <- arguments[5];
    }
    trlength  <- read.table(file=TRLEN,   header=T, comment.char = "");
    protlength <- read.table(file=PROTLEN, header=T, comment.char = "");
    return(list(trlength, protlength, OUTPUT, XLAB, YLAB));
}

merge_df <- function(trlength, protlength) {
    merged <- merge(trlength, protlength, by="Name");
    colnames(merged) <- c("NAME", "TRANSCRIPT", "PROTEIN");
    merged$TRANSCRIPT <- merged$TRANSCRIPT / 3;
    return(merged);
}

compute_breaks <- function(df) {
    all_lengths <- c(df$TRANSCRIPT, df$PROTEIN);
    breaks <- pretty(log10(all_lengths), n=8);
    breaks2 <- as.integer(10^breaks);
    RoundUp <- function(from,to) ceiling(from/to)*to
    breaks2 <- RoundUp(breaks2, 10);
    return(breaks2);
}

compute_perc_position <- function(df) {
    return(max(df$TRANSCRIPT) + 0.05*max(df$TRANSCRIPT));
}

compute_mintotal <- function(df) {
    minx <- min(df$TRANSCRIPT);
    miny <- min(df$PROTEIN);
    mintotal <- min(c(minx, miny)) - 0.05*min(c(minx, miny));
    return(mintotal);
}

seqplot <- function(df, mintotal, perc_position, breaks, xlab, ylab) {
    scatt <- ggplot(df, aes(x=TRANSCRIPT, y=PROTEIN)) +
             geom_point(alpha = 0.3, color="#F7B043") +
                coord_fixed() +
                xlab(paste0("\n", xlab)) +
                ylab(paste0("\n", ylab)) +
                theme_bw() +
                theme(legend.position = "none",
                      plot.title=element_text(vjust=3)) +
                stat_density2d(color="#002640") +
                scale_x_log10(breaks=breaks, limits=c(mintotal, perc_position)) +
                scale_y_log10(breaks=breaks, limits=c(mintotal, perc_position)) +
                stat_function(fun=function(x) log10(x),
                              geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
                stat_function(fun=function(x) log10(x * 0.75),
                              geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
                stat_function(fun=function(x) log10(x * 0.5),
                              geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
                stat_function(fun=function(x) log10(x * 0.25),
                              geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
                stat_function(fun=function(x) log10(x* 0.1),
                              geom="line", colour="#4d4d4d", alpha=0.7, linetype="dashed") +
                annotate("text",
                         label=c("100%","75%","50%","25%", "10%"),
                         x=perc_position,
                         y=c(perc_position, perc_position*0.75, perc_position*0.50, perc_position*0.25, perc_position*0.10),
                         angle=45,
                         hjust=c( 0.9, 1.2, 1.2, 1.2, 1.2),
                         vjust=c(-0.1, 1.1, 1.1, 1.1, 1.1),
                         size=3.5);
    return(scatt);
}


main <- function() {
    cmd <- read_cmdlineargs();
    df <- merge_df(cmd[1], cmd[2]);
    breaks <- compute_breaks(df);
    perc_position <- compute_perc_position(df);
    mintotal <-compute_mintotal(df);
    theplot <- seqplot(
        df=df, 
        mintotal=mintotal, 
        perc_position=perc_position, 
        breaks=c(breaks),
        xlab=cmd[4],
        ylab=cmd[5]
        );
    ggsave(file=cmd[[3]], theplot);
}


if (!interactive()) {
    main();
}
