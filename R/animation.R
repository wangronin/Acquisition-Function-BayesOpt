library(ggplot2)
# library(gganimate)
library(DiceOptim)
library(DiceKriging)
library(dplyr)
library(reshape2)
library(latex2exp)

setwd('~/Dropbox/code_base/acquisition/R')
source('./criteria.R')
source('./fitness.R')
source('./gganimate.R')
source('saving.R')
source('utilities.R')

# initialize GP
set.seed(11)
lower <- -5
upper <- 7

f <- ackley
doe <- runif(6, lower, upper) %>% data.frame(x = .)
response <- sapply(doe$x, f)

gp <- km(~1, doe, response, covtype = 'matern3_2', nugget = 1e-10,
         control = list(trace = F, pgtol = 1e-15, factr = 1e4, maxit = 1e5))

x <- c(seq(lower, upper, length.out = 600), doe$x)
res <- predict(gp, x, 'UK', checkNames = F) 
df <- data.frame(x, y = sapply(x, f), y.hat = res$mean, se = res$sd)

df.MGF <- data.frame(x)
for (t in exp(seq(-2.5, log(4), length.out = 30))) {
  df.MGF %<>% mutate(tmp = sapply(x, . %>% MGF(model = gp, t = t))) %>%
    mutate(tmp = 6 * tmp / max(tmp)) %>%
    rename_(.dots = setNames(list('tmp'), sprintf('%.5f', t)))
}

df.MGF %<>% melt(id = 'x', variable.name = 't', value.name = 'MGF') %>%
  mutate(t = as.numeric(levels(t))[t]) 

df.MGF.max <- df.MGF %>% 
  group_by(t) %>% 
  summarise(x_max = x[MGF == max(MGF)], max = max(MGF))

p <- ggplot(df.MGF) + 
  # 95% confidence interval
  geom_ribbon(data = df, aes(x = x, ymin = y.hat - 1.96 * se, 
                             ymax = y.hat + 1.96 * se),
              fill = '#26989C', show.legend = TRUE,
              alpha = 0.2) +
  # scale_fill_manual(name = NULL, values = c('95% CI' = )) + 
  # objective function
  geom_line(data = df, aes(x, y, linetype = 'objective'), 
            size = 1, alpha = 0.3) +
  # GPR model prediction
  geom_line(data = df, aes(x, y.hat, linetype = 'prediction'),
            size = 1.2, alpha = 0.3) +
  # data points
  geom_point(data = data.frame(x = doe$x, y = response), aes(x, y),
             colour = 'black', alpha = .7, shape = 20, size = 6) +
  # MGF-based acquisition function
  geom_line(aes(x, MGF, colour = t, frame = t), size = 1.8) + 
  # shaded area under the acquisition function
  geom_ribbon(aes(x, ymin = 0, ymax = MGF, fill = t, frame = t), alpha = 0.5) +
  # maximum of the acqusition function
  geom_point(data = df.MGF.max, aes(x_max, y = 0, frame = t, colour = t), 
             alpha = .8, shape = 20, size = 9) +
  # vertical line at the maximum
  geom_segment(data = df.MGF.max, aes(x_max, y = -0.5, 
                                      xend = x_max, yend = max * 1.2,
                                      frame = t, colour = t),
               linetype = 'dashed', size = 1.5) + 
  # aesthetic adjustments
  scale_colour_gradientn(name = NULL, colours = colorspace::rainbow_hcl(10), 
                         trans = "log", guide = FALSE) +
  scale_linetype_manual(name = NULL, values = c('prediction' = 'dotted', 
                                                'objective' = 'solid')) + 
  scale_fill_gradientn(name = NULL, colours = colorspace::rainbow_hcl(10), 
                       trans = "log", guide = FALSE) + 
  labs(x = TeX('$x$'), y = 'f') + 
  ggtitle("t =") + 
  theme_grey() + 
  theme(plot.title = element_text(family = "LMMath-Regular",
                                  hjust = 0.5, size = 20),
        text = element_text(family = "LMMath-Regular", size = 25),
        axis.title = element_text(family = "LMMath-Regular",
                                  face = 'italic'))

animation::ani.options(ani.width = 1600, ani.height = 1000)
gganimate(p, interval = .5, filename = "./test.mp4")
