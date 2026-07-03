# cr_figure2_v2
#
# This script will create Figures 2 (version 2) for the paper 
# "Sensitivity analyses for missing data: A practical guide for planning, 
#  conducting and reporting."
# 
# Written by R Mainzer

# Load packages
library(colorspace)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Set seed
set.seed(170924)

# ------------------------------------------------------------------------------

# Generate variables
n <- 30
y <- rnorm(n)
pmy <- exp(0.1 + 1.5 * y) / (1 + exp(0.1 + 1.5 * y))
my <- ifelse(runif(n) > pmy, 1, 0)
table(my)/n
head(data.frame(y, pmy, my))

nmiss <- sum(my)

# Generate x-coordinates
x <- runif(n, min = 0.75, max = 1.25)

# ------------------------------------------------------------------------------

# Full
# x ~ 1
dat <- data.frame(name = "Full", y = y, x = x, col = 0, shape = 0, size = 1)

# Observed
# x ~ 2
yobs <- y[which(my == 0)]
xobs <- x[which(my == 0)] + 1
dat <- rbind(dat, 
             data.frame(name = "Obs", y = yobs, x = xobs, col = 1, shape = 0, size = 1),
             data.frame(name = "Obs", y = mean(yobs), x = 2.2, col = 1, shape = 1, size = 1.6))

# Missing
# x ~ 3
ymis <- y[which(my == 1)]
xmis <- x[which(my == 1)] + 2
dat <- rbind(dat, 
             data.frame(name = "Miss", y = ymis, x = xmis, col = 2, shape = 0, size = 1),
             data.frame(name = "Miss", y = mean(ymis), x = 3.2, col = 2, shape = 1, size = 1.6))

# ------------------------------------------------------------------------------

# Delta-adjusted MI
# x ~ 4, 5, 6
pm1 <- yobs
pm2 <- rnorm(nmiss)
pm3 <- pm2 - 0.765
xpm1 <- xobs + 2
xpm2 <- runif(nmiss, min = 4.75, max = 5.25)
xpm3 <- xpm2 + 1
dat <- rbind(dat, 
             data.frame(name = "delta", y = pm1, x = xpm1, col = 1, shape = 0, size = 1),
             data.frame(name = "delta", y = pm2, x = xpm2, col = 7, shape = 0, size = 1),
             data.frame(name = "delta", y = pm3, x = xpm3, col = 3, shape = 0, size = 1),
             data.frame(name = "delta", y = 0, x = 6.5, col = 3, shape = 1, size = 1.6))

# Stacked MI
# x ~ 8, 9, 10
sm1 <- yobs
sm2 <- rnorm(nmiss)
sm3 <- sm2
xsm1 <- xpm1 + 4
xsm2 <- runif(nmiss, min = 8.75, max = 9.25)
xsm3 <- xsm2 + 1

pmy_s <- exp(0.6 + 1.5 * sm2) / (1 + exp(0.6 + 1.5 * sm2))
#which(pmy_s == min(pmy_s))
#pmy_s[21] <- 0.208

mean(c(sm1/length(sm1), sm2/pmy_s))
dat <- rbind(dat, 
             data.frame(name = "stacked", y = sm1, x = xsm1, col = 1, shape = 0, size = 1),
             data.frame(name = "stacked", y = sm2, x = xsm2, col = 7, shape = 0, size = 1),
             data.frame(name = "stacked", y = sm3, x = xsm3, col = 3, shape = 0, size = 0.5/pmy_s),
             data.frame(name = "stacked", y = 0, x = 10.5, col = 3, shape = 1, size = 1.6))

# Extreme case approach
# x ~ 12, 13
ec1 <- yobs
ec2 <- rep(min(yobs), nmiss)
xec1 <- xobs + 10
xec2 <- seq(12.7, 14, length = nmiss)
dat <- rbind(dat,
             data.frame(name = "extreme", y = ec1, x = xec1, col = 1, shape = 0, size = 1),
             data.frame(name = "extreme", y = ec2, x = xec2, col = 3, shape = 0, size = 1),
             data.frame(name = "extreme", y = mean(c(ec1, ec2)), x = 13.5, col = 3, shape = 1, size = 1.6))

# ------------------------------------------------------------------------------

# Edit x axis for Full, Observed, Missing
dat$x[dat$name == "Full"] <- dat$x[dat$name == "Full"] - 3.5
dat$x[dat$name == "Obs"] <- dat$x[dat$name == "Obs"] - 2.5
dat$x[dat$name == "Miss"] <- dat$x[dat$name == "Miss"] - 1.5

# Pick colours
lighten(col = "#377EB8", amount = 0.26, method = c("absolute"), space = c("HCL"), 
        fixup = TRUE)

pal1 <- brewer.pal(n = 4, name = "Set2")
pal2 <- brewer.pal(n = 4, name = "Dark2")
pal3 <- brewer.pal(n = 6, name = "Set1")
show_col(c(pal1, pal2, pal3, "#8CC2FE"))

# Graph
fig <- ggplot(dat, aes(x = x, y = y)) + 
  geom_point(aes(color = factor(col), 
                 shape = factor(shape),
                 size = size)) +
  geom_hline(aes(yintercept = 0, col = "red"), linetype = 2) +
  geom_vline(aes(xintercept = 2.75, col = "black")) +
  annotate("curve", x = 4, y = 2.5, xend = 5.5, yend = 2.2,
           curvature = -0.5, arrow = arrow(length = unit(0.2, 'cm'))) +
  annotate("text", x = 4.6, y = 2.9, label = "Generate then shift") + 
  #annotate("curve", x = 5, y = 2.5, xend = 5.8, yend = 2.2,
  #         curvature = -0.5, arrow = arrow(length = unit(0.2, 'cm'))) +
  #annotate("text", x = 6, y = 2.8, label = "Shift") + 
  annotate("curve", x = 8, y = 2.5, xend = 9.5, yend = 2.2,
           curvature = -0.5, arrow = arrow(length = unit(0.2, 'cm'))) +
  annotate("text", x = 8.6, y = 2.9, label = "Generate then weight") + 
  #annotate("curve", x = 9, y = 2.5, xend = 9.8, yend = 2.2,
  #         curvature = -0.5, arrow = arrow(length = unit(0.2, 'cm'))) +
  #annotate("text", x = 9.5, y = 2.8, label = "Weight") + 
  scale_x_continuous(breaks = c(-2.5, -0.5, 1.5, 5, 9, 13),
                     labels = c("Full", "Observed", "Missing",
                                "Delta-adjusted MI", "Stacked MI",
                                "Extreme case")) +
  xlab(" ") +
  ylab(" ") +
  #scale_y_continuous(limits = c(-2.5, 2.2)) +
  scale_shape_manual(values = c(16, 4)) +
  scale_color_manual(values = c("#984EA3", "#377EB8", "#E41A1C", 
                                "#D95F02", "#8CC2FE", "black", "#E78AC3")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_text(hjust = 0),
        panel.grid.minor = element_blank())
fig

ggsave("results/Figure 2.jpg", plot = last_plot(), width = 9.1, height = 4.2, units = "in")
