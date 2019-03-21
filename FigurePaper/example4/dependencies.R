## plotM plots a pie chart using ggplot
plotM <- function(x, res, motifs=motifs, motifs2=motifs2) {
    df <- data.frame(table(x))
    names(df) <- c("motif", "value")
    df <- df[order(df$value), ]
    fa <- factor(levels(df$motif), levels=c(motifs, motifs2))

    pie <- ggplot(df, aes(x="", y=value, fill=motif)) +
                    geom_bar(width = 1, stat = "identity") +
                    coord_polar("y", start=0)

    blank_theme <- theme_minimal()+
                    theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.border = element_blank(),
                    panel.grid=element_blank(),
                    axis.ticks = element_blank(),
                    plot.title=element_text(size=14, face="bold"))

    pie + scale_fill_manual(values=rainbow(length(levels(fa)))[fa]) + blank_theme +
        theme(axis.text.x=element_blank()) + labs(title=res)
}
