library(ggplot2)
library(patchwork)
library(psych)
library(GGally)
library(RColorBrewer)


head(fluorescent_signal_arcsinh)

head(fluorescent_signal_arcsinh_Znorm)

head(fluorescent_signal_full)

fluorescent_signal_arcsinh_full <- asinh(fluorescent_signal_full/150)
head(fluorescent_signal_arcsinh_full)

fluorescent_signal_arcsinh_full_ZScore <- scale(fluorescent_signal_arcsinh_full)
par(mfrow = c(6,2))

plot_list <- list()

for(i in 1:6){
  original <- as.numeric(fluorescent_signal_arcsinh_full[,i])
  normalized <- as.numeric(fluorescent_signal_arcsinh_full_ZScore[,i])
  combined <- c(original, normalized)
  combined_df <- data.frame(Value = combined ,Label = c(rep("Original", 372385), rep("Normalized", 372385)))
  xlab_ <- colnames(fluorescent_signal_arcsinh_full)[i]
  g1 <- ggplot(combined_df, aes(factor(Label), Value))  + xlab(xlab_) + geom_violin() #+ geom_density()
  plot_list[[i]] <- g1
}
wrap_plots(plot_list)

plot_list <- list()

colnames(fluorescent_signal_arcsinh_full) <- c(colnames(fluorescent_signal_arcsinh_full)[1:6], colnames(fluorescent_signal_arcsinh)[1:5], "l/d")
colnames(fluorescent_signal_arcsinh_full_ZScore) <- c(colnames(fluorescent_signal_arcsinh_full)[1:6], colnames(fluorescent_signal_arcsinh)[1:5], "l/d")


for(i in 7:12){
  original <- as.numeric(fluorescent_signal_arcsinh_full[,i])
  normalized <- as.numeric(fluorescent_signal_arcsinh_full_ZScore[,i])
  combined <- c(original, normalized)
  combined_df <- data.frame(Value = combined ,Label = c(rep("Original", 372385), rep("Normalized", 372385)))
  xlab_ <- colnames(fluorescent_signal_arcsinh_full)[i]
  g1 <- ggplot(combined_df, aes(factor(Label), Value))  + xlab(xlab_) + geom_violin() #+ geom_density()
  plot_list[[i-6]] <- g1
}

wrap_plots(plot_list)

pairs.panels(fluorescent_signal_arcsinh_full[rand_ind,1:6], lm = TRUE, stars = TRUE)

ggpairs(fluorescent_signal_arcsinh_full[rand_ind,], columns = 1:6, upper = list(continuous = wrap("cor", size = 2.5)))



cols = brewer.pal(11, "RdBu")   # goes from red to white to blue
pal = colorRampPalette(cols)
cor_colors = data.frame(correlation = seq(-1,1,0.01), 
                        correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)

panel.cor <- function(x, y, digits=2, cex.cor)
{
  par(usr = c(0, 1, 0, 1))
  u <- par('usr') 
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  r <- cor(x, y,method="spearman",use="complete.obs")
  test <- cor.test(x,y)
  bgcolor = cor_colors[2+(-r+1)*100,2]    # converts correlation into a specific color
  do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box
  
  if (test$p.value> 0.05){
    text(0.5,0.5,"Insignificant",cex=1.5)
  } else{
    text(0.5, 0.75, paste("r=",round(r,2)),cex=2.5) # prints correlatoin coefficient
    text(.5, .25, paste("p=",formatC(test$p.value, format = "e", digits = 1)),cex=2)  
    abline(h = 0.5, lty = 2) # draws a line between correlatoin coefficient and p value
  }
  
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 19, cex = 1.2, 
                        col.smooth = "blue", span = 2/3, iter = 3, ...) {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), lwd=2.5, col = col.smooth, ...)
}
panel.hist <- function(x, ...)
{
  points(x,x)
}

pairs(fluorescent_signal_arcsinh_full[rand_ind,7:12],lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,cex.labels=2)
pairs(fluorescent_signal_arcsinh_full_ZScore[rand_ind,7:12],lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist,cex.labels=2)






