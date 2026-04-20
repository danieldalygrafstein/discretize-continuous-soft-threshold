################################################
# Plotting simulation results

################################################
# libraries

library(tidyverse)

################################################
# read in results

results <- as_tibble(read.csv("full-sim-results.csv")[,-1])

################################################
# plot comparing type 1 error and X1/X2 cor
# split by number of categories and cat method
# using logistic method with beta = 2
# and cat algorithm of equal
# and normal X generation

fig_1_results <- results %>% 
  filter(X_gen == "normal" & y_gen == "logistic" &
           y_param == 2 & (is.na(cat_alg) | cat_alg == "equal") &
           (is.na(cat_method) | cat_method %in% c("hard", "normal_add", "triangle_add"))) %>% 
  mutate(cat_method = replace_na(cat_method, "cont")) %>% 
  group_by(sigma_x, categories, cat_method, n) %>% 
  summarise(type_1_error = mean(type_1_error))

fig_1_categories <- ggplot(fig_1_results) +
  geom_point(aes(x = sigma_x, y = type_1_error, colour = cat_method)) +
  scale_colour_manual(values = c("red4", "dodgerblue4", "gold4", "forestgreen"),
                                 labels = c("Continuous", "Hard", "Soft - normal", "Soft - triangle")) +
  geom_line(aes(x = sigma_x, y = type_1_error, colour = cat_method,
                linetype = factor(categories, levels = c("2", "3", "4", "1"),
                                  labels = c("2", "3", "4", "cont"))), size = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  labs(linetype = "Categories",
       colour = "Categorization method",
       x = bquote(sigma[.(paste0("x", ",c"))])) +
  ylab("Type-I error rate") +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(0,1)) +
  facet_wrap( ~ factor(n, levels = c("100", "1000"),
                       labels = c("n = 100", "n = 1000")), ncol = 1)

fig_1_categories

png("plots/fig_1_categories_v2.png",
    width = 8, height = 6, units = "in", res = 600)
fig_1_categories
dev.off()

################################################
# plot comparing type 1 error and X1/X2 cor
# split by sample size, association between x and y
# using quantile cat method, 2 categories,
# only hard and soft triangle no add
# X gen of lognoraml

fig_2_results <- results %>% 
  filter(X_gen == "lognormal" & y_gen == "linear" & 
           cat_alg == "quantile" &
           (cat_method %in% c("hard", "triangle_no_add", "triangle_add")) &
           categories == 2)

fig_2_sample <- ggplot(fig_2_results) +
  geom_point(aes(x = sigma_x, y = type_1_error, colour = as.factor(cat_method))) +
  geom_line(aes(x = sigma_x, y = type_1_error, colour = as.factor(cat_method),
                linetype = as.factor(y_param)), size = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  labs(linetype = "Residual variance",
       colour = "Discretization method",
       x = bquote(sigma[.(paste0("x", ",c"))])) +
  ylab("Type-I error rate") +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(0,1)) +
  facet_wrap( ~ factor(n, levels = c("100", "1000"),
                       labels = c("n = 100", "n = 1000")), ncol = 1) +
  scale_colour_manual(values = c("dodgerblue4", "forestgreen", "orange"),
                      labels = c("Hard", "Soft - triangle add groups", "Soft - triangle no add"))

fig_2_sample

png("plots/fig_2_sample_size_v2.png",
    width = 8, height = 6, units = "in", res = 600)
fig_2_sample
dev.off()


################################################
# plot comparing type 1 error and X1/X2 cor
# split by categorization method and X2 distribution
# only take hard, triangle_add, triangle_no_add
# only for logistic y with beta_2 = 2 and 3 categories
# and n=1000

fig_3_results <- results %>% 
  filter(y_gen == "linear" &
           y_param == 1.02 &
           (cat_method %in% c("hard", "triangle_add", "triangle_no_add")) &
           categories == 4 & 
           n == 1000)

fig_3_cat_alg <- ggplot(fig_3_results) +
  geom_point(aes(x = sigma_x, y = type_1_error, colour = as.factor(cat_method))) +
  geom_line(aes(x = sigma_x, y = type_1_error, colour = as.factor(cat_method),
                linetype = as.factor(X_gen)), size = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  labs(linetype = "Distribution",
       colour = "Categorization method",
       x = bquote(sigma[.(paste0("x", ",c"))])) +
  ylab("Type-I error rate") +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(0,1)) +
  facet_wrap( ~ factor(cat_alg), ncol = 1) +
  scale_colour_manual(values = c("dodgerblue4", "forestgreen", "orange"),
                      labels = c("Hard", "Soft - triangle add groups", "Soft - triangle no add"))

fig_3_cat_alg

png("plots/fig_3_cat_alg_v2.png",
    width = 8, height = 6, units = "in", res = 600)
fig_3_cat_alg
dev.off()


