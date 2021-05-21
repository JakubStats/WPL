#--------------------------------------------------------------------------------------------------------
# BiMJ_Efficiency_Hwang_et_al.R
#--------------------------------------------------------------------------------------------------------

rm(list = ls())

# Make sure all functions stored and loaded first in your R-console.

source("BiMJ_Functions_Hwang_et_al.R")

# Positive-binomial models - this should match Figure 1 in the manuscript.

Ns <- c(3, 5, 7, 10, 15, 20, 50, 75, 100)
sp <- seq(0.05, 0.95, by = 0.05)

ref <- list()

for(i in 1:length(Ns)) {
  k <- Ns[i]
  
  bin.plavar <- binpl.avar(sp, k)
  bin.pl0avar <- binpl0.avar(sp, k)
  bin.mlavar <- binmle.avar(sp, k)
  
  bin.ref <- cbind(bin.mlavar/bin.pl0avar, bin.mlavar/bin.plavar)
  
  ref <- c(ref, list(bin.ref))
}

my.title1 <- bquote("m" == .(Ns[1]))
ref1 <- ref[[1]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_1 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[2]))
ref1 <- ref[[2]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_2 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[3]))
ref1 <- ref[[3]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_3 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.6, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[4]))
ref1 <- ref[[4]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_4 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[5]))
ref1 <- ref[[5]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_5 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[6]))
ref1 <- ref[[6]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_6 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.75, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[7]))
ref1 <- ref[[7]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_7 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[8]))
ref1 <- ref[[8]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_8 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

my.title1 <- bquote("m" == .(Ns[9]))
ref1 <- ref[[9]]
colnames(ref1) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref1))
ref1 <- cbind(ref1, rep(sp, 2))
colnames(ref1) <- c("efficiency", "model", "p")
p_9 <- ggplot(ref1, aes(x = p, y = efficiency, colour = model, group = model)) + 
  geom_line(size = 1, aes(linetype = model)) + ggtitle(my.title1) + xlab("p") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

par(mfrow = c(3, 3), las=1)

multiplot(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, 
          layout = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                          ncol = 3, nrow = 3, byrow = T))

#--------------------------------------------------------------------------------------------------------

# Positive-Poisson models - this should match Figure 2 in the manuscript.

slam <- seq(0.1, 20, by = 0.1)
plavar0 <- pl.avar0(slam)
plavar <- pl.avar2(slam, 50)
mlavar <- mle.avar(slam)
ef <- cbind(mlavar, plavar0, plavar)
ref <- ef[, 1]/ef[, 2:3]

par(mfrow = c(1, 1), las = 1)

colnames(ref) <- c("PL", "WPL")
ref1 <- stack(as.data.frame(ref))
ref1 <- cbind(ref1,rep(slam, 2))
colnames(ref1) <- c("efficiency", "model", "lambda")

p_1 <- ggplot(ref1, aes(x = lambda, y = efficiency, colour = model, group = model)) +
  geom_line(size = 1, aes(linetype = model)) + xlab(expression(lambda)) + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = c("Red", "Blueviolet")) + ylim(0.8, 1.05) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed")

multiplot(p_1, layout = matrix(c(1), ncol = 1, nrow = 1, byrow = T))
