# MIT License				
# 
# Copyright ©			  2021, by H.B.Wolff	
# 
# Contributors: 		H.B. Wolff1, V. Qendri1, N. Kunst1,2,3,4, F. Alarid-Escudero5, V.M.H. Coupé1	
# Affiliations: 		1.Department of Epidemiology and Data Science, Amsterdam UMC, Vrije Universiteit Amsterdam, Amsterdam Public Health, Amsterdam, Netherlands.	
#                   2. Department of Health Management and Health Economics, Faculty of Medicine, University of Oslo, Oslo, Norway.	
#                   3. Cancer Outcomes, Public Policy and Effectiveness Research (COPPER) Center, Yale University School of Medicine and Yale Cancer Center.	
#                   4. Public Health Modeling Unit, Yale University School of Public Health, New Haven, CT.	
#                   5. Division of Public Administration, Center for Research and Teaching in Economics (CIDE), Aguascalientes, AGS, Mexico, MX-AGU, Mexico.	
# Citation:			    [Methods for communicating the impact of parameter uncertainty in a multiple strategies cost-effectiveness comparison.]
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy				
# of this software and associated documentation files (the "Software"), to deal				
# in the Software without restriction, including without limitation the rights				
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell				
# copies of the Software, and to permit persons to whom the Software is				
# furnished to do so, subject to the following conditions:				
#   
#   The above copyright notice and this permission notice shall be included in all				
# copies or substantial portions of the Software.				
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR				
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,				
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE				
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER				
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,				
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE				
# SOFTWARE.



####### load libraries #######



# make sure these libraries are installed before loading
library(readxl)
library(matrixStats)
library(devtools)
library(testthat)
library(ggplot2)
library(tidyr)
library(reshape2)
library(directlabels)



####### Model Parameters Required for Plots #######



# upload you own data to 'mydata', replace path to correct folder and file
mydata <- read_excel(" ... insert path to folder here ... /uncertainty_data.xlsx")

# store the costs and effects in two dataframes
# please select the correct collumns and lines
effects<- data.matrix(mydata[,1:108])
costs<- data.matrix(mydata[,109:216])

# obtain or define variables
n_strategies <- ncol(effects)
strategies<- as.character(1:n_strategies)
n_sim <- nrow(effects)
probabilities <- seq((1/n_sim), 1, (1/n_sim))

#define your own WTP axis here (minimum, maximum, stepsize)
wtp <- seq(200, 200000, by=200)
n_wtps <- length(wtp)



####### CEAC variants #######



# relaxed criteria parameters
max_rank <- 1         # rank is smaller or equal to this value. Default: 1
acc_benefit_dif <- 0  # minimum accepteble (relevant) benefit difference. Default: 0
rel_benefit_dif <- 0  # relative accepteble (relevant) benefit difference. Default: 0

# construct empty matrices and vectors to store data
frontiers <- matrix(0, nrow = n_wtps, ncol = 5)
indicator_rs <- matrix(0, nrow = n_sim, ncol = n_strategies)
indicator_ra <- matrix(0, nrow = n_sim, ncol = n_strategies)
indicator_elc <- matrix(0, nrow = n_sim, ncol = n_strategies)
benefit <- matrix(0, nrow = n_sim, ncol = n_strategies)
ceac <- matrix(0, nrow = n_wtps, ncol = n_strategies)
inv_ceac <- matrix(0, nrow = n_wtps, ncol = n_strategies)
elc <- matrix(0, nrow = n_wtps, ncol = n_strategies)
colnames(ceac) <- strategies
colnames(inv_ceac) <- strategies
colnames(elc) <- strategies


# analysis of data
# loop through WTP thresholds
for (l in 1:n_wtps) {
  # calculate net health / monetary benefit at wtp[l]
  benefit <- (effects * wtp[l]) - costs
  #benefit <- effects - (costs / wtp[l])

  #loop through PSA
  for (sim in 1:n_sim) {

    ### store data in three matrixes at the same time: indicator_rs, indicator_rs, indicator_elc. 
    ### use comments to count either tradional or one of the relaxed indicator functions
    
    ### traditional CEAC
    Acceptable_Benefits_rs <- max(benefit[sim, ])
    indicator_rs[sim,] <- ifelse(benefit[sim,] == Acceptable_Benefits_rs, 1 , 0)

    #### inverse CEAC indicator (probability not cost effective)
    Acceptable_Benefits_ra <- min(benefit[sim, ])
    indicator_ra[sim,] <- ifelse(benefit[sim,] <= Acceptable_Benefits_ra, 1 , 0)

    #### expected loss curve
    max_Benefit <- max(benefit[sim, ])
    indicator_elc[sim, ] <- (max_Benefit - benefit[sim, ])

    #### relaxed criteria: ####

    ### CEAC using ranks
    #order_benefit <- order(benefit[sim,])
    #Acceptable_Benefits_rs <- benefit[sim, order_benefit[n_strategies - max_rank + 1]]
    #indicator_ra[sim,] <- ifelse(benefit[sim,] >= Acceptable_Benefits_rs, 1 , 0)

    ### inversed CEAC using ranks
    #order_benefit <- order(benefit[sim,])
    #Acceptable_Benefits_ra <- benefit[sim, order_benefit[max_rank]]
    #indicator_rs[sim,] <- ifelse(benefit[sim,] <= Acceptable_Benefits_ra, 1 , 0)

    ### CEAC using a fixed acceptable difference in benefit value
    #Acceptable_Benefits_rs <- max(benefit[sim, ]) - acc_benefit_dif
    #indicator_ra[sim,] <- ifelse(benefit[sim,] >= Acceptable_Benefits_rs, 1 , 0)

    ### inversed CEAC using a fixed acceptable difference in benefit value
    #Acceptable_Benefits_ra <- min(benefit[sim, ]) + acc_benefit_dif
    #indicator_ra[sim,] <- ifelse(benefit[sim,] <= Acceptable_Benefits_ra, 1 , 0)

    ### CEAC using a relative acceptable difference in benefit value
    #Acceptable_Benefits_rs <- max(benefit[sim, ]) - abs(rel_benefit_dif * max(benefit[sim, ]))
    #indicator_rs[sim,] <- ifelse(benefit[sim,] >= Acceptable_Benefits_rs, 1 , 0)

    ### inversed CEAC using a relative acceptable difference in benefit value
    #Acceptable_Benefits_ra <- min(benefit[sim, ]) + abs(rel_benefit_dif * min(benefit[sim, ]))
    #indicator_rs[sim,] <- ifelse(benefit[sim,] <= Acceptable_Benefits_ra, 1 , 0)
  }

  # calculate average of indicator functions
  ceac[l, ] <- colMeans(indicator_rs)
  inv_ceac[l, ] <- colMeans(indicator_ra)
  elc[l, ]    <- colMeans(indicator_elc)

  # calculate frontiers, alternative frontiers are also stored as options to plot.
  frontiers[l,1] <- which.min(elc[l, ]) # ELC frontier (CEAF)
  frontiers[l,2] <- which.max(ceac[l, ]) # CEAC frontier
  frontiers[l,3] <- which.min(inv_ceac[l, ]) # inverse CEAC frontier
  frontiers[l,4] <- which.max(colMins(benefit)) # maximin
  frontiers[l,5] <- which.max(colMaxs(benefit)) # maximax
}

# create melted dataframes that can be combined into one huge dataframe

#add frontiers to one of the ceac dataframes
CEAC <- data.frame(wtp, ceac, frontiers, stringsAsFactors = FALSE)
colnames(CEAC) <- c("WTP", strategies, "fstrat_elc", "fstrat_ceac", "fstrat_inv", "maximin", "maximax")
CEAC <- melt(CEAC, id.vars = c("WTP", "fstrat_elc", "fstrat_ceac", "fstrat_inv", "maximin", "maximax"), variable.name = "Strategy", value.name = "CEAC")
CEAC$Strategy <- as.character(CEAC$Strategy)
CEAC$On_Frontier_ceac <- (CEAC$fstrat_ceac == CEAC$Strategy)
CEAC$On_Frontier_elc <- (CEAC$fstrat_elc == CEAC$Strategy)
CEAC$On_Frontier_inv <- (CEAC$fstrat_inv == CEAC$Strategy)
CEAC$On_Frontier_mami <- (CEAC$maximin == CEAC$Strategy)
CEAC$On_Frontier_mama <- (CEAC$maximax == CEAC$Strategy)
CEAC$fstrat_ceac <- NULL
CEAC$fstrat_elc <- NULL
CEAC$fstrat_inv <- NULL
CEAC$maximin <- NULL
CEAC$maximax <- NULL
CEAC <- CEAC[order(CEAC$WTP), ]
rownames(CEAC) <- NULL

invCEAC <- data.frame(wtp, inv_ceac, stringsAsFactors = FALSE)
colnames(invCEAC) <- c("WTP", strategies)
invCEAC <- melt(invCEAC, id.vars = "WTP", variable.name = "Strategy", value.name = "invCEAC")
invCEAC <- invCEAC[order(invCEAC$WTP), ]
rownames(invCEAC) <- NULL

elc_df <- data.frame(wtp, elc, stringsAsFactors = FALSE)
colnames(elc_df) <- c("WTP", strategies)
elc_df <- melt(elc_df, id.vars = "WTP", variable.name = "Strategy", value.name = "ELC_DF")
elc_df <- elc_df[order(elc_df$WTP), ]
rownames(elc_df) <- NULL

# combine dataframes
CEACS_df <- CEAC[,c(1,2,4,5,6,7,8,3)]
CEACS_df[,"invCEAC"] <- invCEAC[,3]
CEACS_df[,"ELC"] <- elc_df[,3]

# calculate ICERS
dominant_strategies <- as.data.frame(as.numeric(names(table(frontiers[,1]))))[,1] # on ELC frontier
ICERS <- matrix(0, nrow = length(dominant_strategies), ncol = 4)
#calculate mean costs and effects of dominant strategies
for (i in 1:length(dominant_strategies)) {
  ICERS[i,1] <- dominant_strategies[i]
  ICERS[i,2] <- mean(costs[,dominant_strategies[i]])
  ICERS[i,3] <- mean(effects[,dominant_strategies[i]])
}

# sort on mean costs
ICERS <- ICERS[order(ICERS[,2]),]
#calculate ICERS
for (i in 2:length(dominant_strategies)) {
  ICERS[i,4] <- (ICERS[i,2] - ICERS[(i-1),2]) / (ICERS[i,3] - ICERS[(i-1),3])
}
ICERS <- data.frame(ICERS)



####### CEAC Plots #######



# plot CEAC
# van Hout BA, et al. Costs, effects and C/E-ratios alongside a clinical trial. Health Econ 1994; 3: 309-319. 1994/09/01. DOI: 10.1002/hec.4730030505.
# Fenwick E, et al. Representing uncertainty: the role of cost-effectiveness acceptability curves. Health Econ 2001; 10: 779-787. 2001/12/18. DOI: 10.1002/hec.635.
ceac_plot_RS <- ggplot(data = CEACS_df, mapping = aes_(x = CEACS_df$WTP, y = CEACS_df$CEAC, color = CEACS_df$Strategy, col = "full")) +
  geom_line() + xlab("Willingness to pay (Euro/QALY)") + ylab("Probability Cost-Effective") +
  geom_line(data=CEACS_df[CEACS_df$On_Frontier_elc==TRUE,], aes(x=WTP, y=CEAC), colour="black",size=1.5, linetype=2) +
  #geom_dl(aes(label = CEACS$Strategy), method = list(cex=0.6, dl.combine("first.points", "last.points")), color = "black") +
  #geom_rug(data = ICERS, mapping = aes_(x = ICERS[,4], y = NULL, color="black"), sides="b") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
ceac_plot_RS

ggsave("CEAC.pdf", width = 8, height = 5)

# plot ELC
# Eckermann S, et al. Health technology assessment in the cost-disutility plane. Med Decis Making 2008; 28: 172-181. 2008/03/22. DOI: 10.1177/0272989X07312474.
elc_plot <- ggplot(data = CEACS_df, mapping = aes_(x = CEACS_df$WTP, y = CEACS_df$ELC, color = invCEAC$Strategy, col = "full")) +
  geom_line() + xlab("Willingness to pay (Euro/QALY)") + ylab("Expected Loss (Euro)") +
  geom_line(data=CEACS_df[CEACS_df$On_Frontier_elc==TRUE,], aes(x=WTP, y=ELC), colour="black",size=1.5, linetype=2) +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
elc_plot

ggsave("ELC.pdf", width = 8, height = 5)

# plot elc heatmap
heatmap <- ggplot(data = CEACS_df) +
  geom_line(aes_(x = CEACS_df$WTP, y = CEACS_df$ELC, group = CEACS_df$Strategy, color = CEACS_df$CEAC)) +
  geom_line(data=CEACS_df[CEACS_df$On_Frontier_elc==TRUE,], aes(x=WTP, y=ELC), color="black",size=1.5, linetype=2) +
  scale_colour_gradient(low = "blue", high = "red", guide = "colourbar") +
  xlab("Willingness to pay (Euro/QALY)") + ylab("Expected Loss (Euro)") +
  theme_classic(base_size = 12) + theme(legend.position="right") + labs(color="Probability \ncost-effective")
heatmap

ggsave("ELC_heat.pdf", width = 8, height = 5)

# plot ceac heatmap
heatmap <- ggplot(data = CEACS_df) +
  geom_line(aes_(x = CEACS_df$WTP, y = CEACS_df$CEAC, group = CEACS_df$Strategy, color = CEACS_df$ELC)) +
  geom_line(data=CEACS_df[CEACS_df$On_Frontier_elc==TRUE,], aes(x=WTP, y=CEAC), color="black",size=1.5, linetype=2) +
  scale_colour_gradient(low = "blue", high = "red", guide = "colourbar") +
  xlab("Willingness to pay (Euro/QALY)") + ylab("Probability Cost-Effective") +
  theme_classic(base_size = 12) + theme(legend.position="right") + labs(color="Expected Loss")
heatmap

ggsave("CEAC_heat.pdf", width = 8, height = 5)



###### Net Benefit Density  #######



# Naversnik K. Output correlations in probabilistic models with multiple alternatives. Eur J Health Econ 2015; 16: 133-139. 2014/01/07. DOI: 10.1007/s10198-013-0558-0.

# calculate benefits for fixed WTP
fixed_wtp <- 50000
n_steps <- 100
smooth <- 0.5

Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)
Histogram <- matrix(0, nrow = n_steps, ncol = n_strategies)
Density <- matrix(0, nrow = n_steps, ncol = n_strategies)

# use either Net Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - costs # Net Monetary Benefits
#Benefits <- effects - (costs / fixed_wtp) # Net Health Benefits

for (strat in 1:n_strategies) {
  Benefits[ ,strat] <- sort(Benefits[ ,strat], decreasing = FALSE)
}

# calculate parameters for histogram
min_Benefit <- min(Benefits[1, ])
max_Benefit <- max(Benefits[n_sim, ])
bin_with <- (max_Benefit - min_Benefit) / n_steps

# fill histogram per strategy
for (strat in 1:n_strategies) {
  bin <- 1
  probability_per_sim <- bin_with / (sd(Benefits[, strat]) * n_sim)
  for (sim in 1:n_sim) {
    while (Benefits[sim, strat] > (min_Benefit + (bin * bin_with))) { bin <- bin + 1 }
    Histogram[bin, strat] <- (Histogram[bin, strat] + probability_per_sim)
  }
}

# smoothen histogram
for (strat in 1:n_strategies) {
  Density[1, strat] <- (Histogram[1, strat] + smooth*Histogram[2, strat]) / (2*smooth+1)
  for (step in 2:(n_steps-1)) {
    Density[step, strat] <- (smooth*Histogram[(step-1), strat] + Histogram[step, strat] + smooth*Histogram[(step+1), strat]) / (2*smooth+1)
  }
  Density[n_steps, strat] <- (smooth*Histogram[(n_steps-1), strat] + Histogram[n_steps, strat]) / (2*smooth+1)
}

# make dataframe
Density_df <- as.data.frame(Density)
colnames(Density_df) <- as.character(strategies)
Density_df["benefit"] <- seq((min_Benefit + bin_with/2), (max_Benefit - bin_with/2), bin_with)
Density_df <- melt(Density_df, id.vars = "benefit", variable.name = "Strategy", value.name = "density")

# plot
NB_density_plot <- ggplot(data = Density_df, mapping = aes_(x = Density_df$benefit, y = Density_df$density, color = Density_df$Strategy)) +
  geom_line() +
  xlab("Net Monetary Benefits") + ylab("Probability") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
NB_density_plot

ggsave("NB_density.pdf", width = 8, height = 5)




###### Incremental Benefit Density  #######



# Naversnik K. Output correlations in probabilistic models with multiple alternatives. Eur J Health Econ 2015; 16: 133-139. 2014/01/07. DOI: 10.1007/s10198-013-0558-0.

# calculate incremental benefits for fixed WTP
fixed_wtp<- 50000
n_steps <- 750
smooth <- 0.5

Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)
incBenefits <- matrix(0, nrow = n_sim, ncol = n_strategies)
Histogram <- matrix(0, nrow = n_steps, ncol = n_strategies)
Density <- matrix(0, nrow = n_steps, ncol = n_strategies)

# fill first with eitherNet Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - costs # Net Monetary Benefits
#Benefits <- effects - (costs / fixed_wtp) # Net Health Benefits

for (sim in 1:n_sim) {
  # find the maximum and second highest benefit per sim
  order_benefit <- order(Benefits[sim, ])
  simMax <- (Benefits[sim, order_benefit[n_strategies]])
  simSec <- (Benefits[sim, order_benefit[n_strategies-1]])
  
  # calculate most incremental benefits by substracting simMax from whole row
  incBenefits[sim, ] <- (Benefits[sim, ] - simMax)
  
  # replace the original maximum Benefit 'matrix location' with (simMax - simSec)
  incBenefits[sim, order_benefit[n_strategies]] <- (simMax - simSec)
}

for (strat in 1:n_strategies) {
  incBenefits[ ,strat] <- sort(incBenefits[ ,strat], decreasing = FALSE)
}

# calculate parameters for histogram
min_incBenefit <- min(incBenefits[1, ])
max_incBenefit <- max(incBenefits[n_sim, ])
bin_with <- (max_incBenefit - min_incBenefit) / n_steps


# fill histogram
for (strat in 1:n_strategies) {
  bin <- 1
  probability_per_sim <- bin_with / (sd(incBenefits[, strat]) * n_sim)
  for (sim in 1:n_sim) {
    while (incBenefits[sim ,strat] > (min_incBenefit + (bin * bin_with))) { bin <- bin + 1 }
    Histogram[bin, strat] <- (Histogram[bin, strat] + probability_per_sim)
  }
}

# smoothen histogram
for (strat in 1:n_strategies) {
  Density[1, strat] <- (Histogram[1, strat] + smooth*Histogram[2, strat]) / (2*smooth+1)
  for (step in 2:(n_steps-1)) {
    Density[step, strat] <- (smooth*Histogram[(step-1), strat] + Histogram[step, strat] + smooth*Histogram[(step+1), strat]) / (2*smooth+1)
  }
  Density[n_steps, strat] <- (smooth*Histogram[(n_steps-1), strat] + Histogram[n_steps, strat]) / (2*smooth+1)
}

# make dataframe
Density_df <- as.data.frame(Density)
colnames(Density_df) <- as.character(strategies)
Density_df["benefit"] <- seq((min_incBenefit + bin_with/2), (max_incBenefit - bin_with/2), bin_with)
Density_df <- melt(Density_df, id.vars = "benefit", variable.name = "Strategy", value.name = "density")

# plot
incB_density_plot <- ggplot(data = Density_df, mapping = aes_(x = Density_df$benefit, y = Density_df$density, color = Density_df$Strategy)) +
  geom_line() +
  geom_vline(xintercept = 0, size=1.5) +
  xlim(as.numeric(c(-3000, 2000))) +
  xlab("incremental Net Monetary Benefits") + ylab("Probability") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
incB_density_plot

ggsave("incB_density.pdf", width = 8, height = 5)



###### Stochastic Dominance  #######



# Stinnett AA, et al. Net health benefits: a new framework for the analysis of uncertainty in cost-effectiveness analysis. Med Decis Making 1998; 18: S68-80. 1998/05/05. DOI: 10.1177/0272989X98018002S09.

# calculate benefits for fixed WTP
fixed_wtp<- 50000
Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)

# use eitherNet Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - costs # Net Monetary Benefits
#Benefits <- effects - (costs / fixed_wtp) # Net Health Benefits

for (strat in 1:n_strategies) {
  Benefits[ ,strat] <- sort(Benefits[ ,strat], decreasing = TRUE)
}

# make dataframe
Benefits_df <- as.data.frame(Benefits)
colnames(Benefits_df) <- as.character(strategies)
Benefits_df["probabilities"] <- probabilities
Benefits_df <- melt(Benefits_df, id.vars = "probabilities", variable.name = "Strategy", value.name = "Benefit")

# plot
SD_plot <- ggplot(data = Benefits_df, mapping = aes_(x = Benefits_df$Benefit, y = Benefits_df$probabilities, color = Benefits_df$Strategy)) +
  geom_step() + xlab("Net Benefits") + ylab("Cumulative Probability") + scale_x_reverse() +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
SD_plot

ggsave("stoch_dom.pdf", width = 8, height = 5)



###### Incremental Benefit Plot  #######



# Bala MV, et al. Presenting results of probabilistic sensitivity analysis: the incremental benefit curve. Health Econ 2008; 17: 435-440. 2007/08/19. DOI: 10.1002/hec.1274.

# calculate incremental benefits for fixed WTP
fixed_wtp<- 50000
Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)
incBenefits <- matrix(0, nrow = n_sim, ncol = n_strategies)

# fill first with either Net Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - costs # Net Monetary Benefits
#Benefits <- effects - (costs / fixed_wtp) # Net Health Benefits

for (sim in 1:n_sim) {
  # find the maximum and second highest benefit per sim
  order_benefit <- order(Benefits[sim, ])
  simMax <- (Benefits[sim, order_benefit[n_strategies]])
  simSec <- (Benefits[sim, order_benefit[n_strategies-1]])
  
  # calculate most incremental benefits by substracting simMax from whole row
  incBenefits[sim, ] <- (Benefits[sim, ] - simMax)
  
  # replace the original maximum Benefit 'matrix location' with (simMax - simSec)
  incBenefits[sim, order_benefit[n_strategies]] <- (simMax - simSec)
}

for (strat in 1:n_strategies) {
  incBenefits[ ,strat] <- sort(incBenefits[ ,strat], decreasing = TRUE)
}

# make dataframe
incBenefits_df <- as.data.frame(incBenefits)
colnames(incBenefits_df) <- as.character(strategies)
incBenefits_df["probabilities"] <- probabilities
incBenefits_df <- melt(incBenefits_df, id.vars = "probabilities", variable.name = "Strategy", value.name = "incBenefit")

# plot
incbenefit_plot <- ggplot(data = incBenefits_df, mapping = aes_(x = incBenefits_df$incBenefit, y = incBenefits_df$probabilities, color = incBenefits_df$Strategy)) +
  geom_line() + ylab("Cumulative Probability") + xlim(as.numeric(c(-5000, 5000))) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 12) +
  theme(legend.position="none") +
  xlab("(benefit - maximum)    Incremental Benefits    (maximum - second)") + theme(axis.title.x = element_text(hjust = 0.5))
incbenefit_plot

ggsave("inc_Benefit.pdf", width = 8, height = 5)



###### Cumulative Rankogram  #######



# Epstein D. Beyond the cost-effectiveness acceptability curve: The appropriateness of rank probabilities for presenting the results of economic evaluation in multiple technology appraisal. Health Econ 2019; 28: 801-807. 2019/05/03. DOI: 10.1002/hec.3884.

# calculate benefits for fixed WTP
fixed_wtp<- 50000
Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)

# use eitherNet Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - cost # Net Monetary Benefits
#Benefits <- effectiveness - (cost / fixed_wtp) # Net Health Benefits

# store ranks per sim in Matrix
ranks_benefit <- matrix(0, nrow = n_sim, ncol = n_strategies)
for (sim in 1:n_sim) {
  # caculate the ranks per sim
  ranks_benefit[sim, ] <- order(Benefits[sim, ])
}

# calculate cumulative rank probabilites
cum_ranks_prob <- matrix(0, nrow = n_strategies, ncol = n_strategies)
for (strat in 1:n_strategies) {
  cum_prob <- 0
  for (rank in 1:n_strategies) {
    cum_prob <- (cum_prob + (sum(ranks_benefit[ ,strat]==rank)/(n_sim)))
    cum_ranks_prob[rank, strat] <- cum_prob
  }
}

# make dataframe
cum_ranks_df <- as.data.frame(cum_ranks_prob)
colnames(cum_ranks_df) <- as.character(strategies)
cum_ranks_df["ranks"] <- seq(1, n_strategies, by=1)
cum_ranks_df <- melt(cum_ranks_df, id.vars = "ranks", variable.name = "Strategy", value.name = "Probability")

# plot
rankogram_plot <- ggplot(data = cum_ranks_df, mapping = aes_(x = cum_ranks_df$ranks, y = cum_ranks_df$Probability, color = cum_ranks_df$Strategy)) +
  geom_line() + xlab("Rank greater or equal than") + ylab("Cumulative Probability") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
rankogram_plot

ggsave("cumulative_rankogram.pdf", width = 8, height = 5)



###### Return-Risk Space #######



# O'Brien BJ, et al. Building uncertainty into cost-effectiveness rankings: portfolio risk-return tradeoffs and implications for decision rules. Med Care 2000: 460-468.

# calculate incremental benefits for fixed WTP
fixed_wtp<- 50000
Benefits <- matrix(0, nrow = n_sim, ncol = n_strategies)

# fill first with eitherNet Monetary Benefits or Net Health Benefits
Benefits <- (effects * fixed_wtp) - costs # Net Monetary Benefits
#Benefits <- effects - (costs / fixed_wtp) # Net Health Benefits

# create dataframe with mu and sigma
Return_Risk <- data.frame(colMeans(Benefits), colStdevs(Benefits), as.character(strategies))
colnames(Return_Risk) <- c("mu", "sigma", "strategies")
row.names(Return_Risk) <- as.character(strategies)

Return_Risk_plot <- ggplot(data = Return_Risk, mapping = aes_(x = Return_Risk$mu, y = Return_Risk$sigma, color = Return_Risk$strategies, size = 1)) +
  geom_point() + xlab(expression(paste(mu, " Benefit (Euro)"))) + ylab(expression(paste(sigma, " Benefit (Euro)"))) +
  #xlim(as.numeric(c(444000, 449000))) + ylim(as.numeric(c(38000, 39500)))
  #+ scale_x_continuous(n.breaks = 11) +
  theme_classic(base_size = 12) + theme(legend.position="none")
Return_Risk_plot

ggsave("Return_Risk.pdf", width = 8, height = 5)



###### Isoquant CE plot #######



# Ament A, et al. The interpretation of results of economic evaluation: explicating the value of health. Health Econ 1997; 6: 625-635. 1998/02/18. DOI: 10.1002/(sici)1099-1050(199711)6:6<625::aid-hec309>3.0.co;2-o.
# define parameters
fixed_wtp <- 50000
distance <- 0.01 #QALY distance between isoquants; QALY * WTP = monetary distance

# calculate expected costs and effects
CE_df <- data.frame(as.character(strategies), colMeans(costs), colMeans(effects), stringsAsFactors = FALSE)
colnames(CE_df) <- c("strategies","costs","effects")
rownames(CE_df) <- NULL
CE_df[,"benefits"] <- CE_df$effects * fixed_wtp - CE_df$costs

# calculate variables
dom_strat <- which.max(CE_df$benefits)
intercept_from <- CE_df[dom_strat, "effects"] - CE_df[dom_strat, "costs"] / fixed_wtp # first line through max benefit
nr_lines <- ceiling((which.max(CE_df$effects) - which.min(CE_df$effects)) / distance)
intercept_to <- intercept_from - (distance * nr_lines)

#plot
CE_plane <- ggplot(data = CE_df, mapping = aes_(x = CE_df$costs, y = CE_df$effects, color = CE_df$strategies, size = 1)) +
  geom_point() +
  geom_abline(slope = (1/fixed_wtp), intercept = seq(intercept_from, intercept_to, by = -distance)) +
  ylab("effectiveness (QALY)") + xlab("costs (Euro)") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
CE_plane

ggsave("Isoquants.pdf", width = 8, height = 5)



###### Expected Benefit plot #######



# Stinnett AA, et al. Net health benefits: a new framework for the analysis of uncertainty in cost-effectiveness analysis. Med Decis Making 1998; 18: S68-80. 1998/05/05. DOI: 10.1177/0272989X98018002S09.

# construct empty matrices and vectors to store data
percentage <- 0.025 # alfa / 2
benefit <- matrix(0, nrow = n_sim, ncol = n_strategies)
benefit_mid <- matrix(0, nrow = n_wtps, ncol = n_strategies)
benefit_lower <- matrix(0, nrow = n_wtps, ncol = n_strategies)
benefit_upper <- matrix(0, nrow = n_wtps, ncol = n_strategies)
colnames(benefit) <- strategies
colnames(benefit_mid) <- strategies
colnames(benefit_lower) <- strategies
colnames(benefit_upper) <- strategies

# analysis of data loop through WTP thresholds
for (l in 1:n_wtps) {
  # calculate net health / monetary benefit at wtp[l]
  benefit <- (effects * wtp[l]) - costs
  #benefit <- effects - (costs / wtp[l])

  # perecentiles used
  for (strat in 1:n_strategies) {
    benefit_mid[l, strat] <- quantile(benefit[,strat], 0.5)
    benefit_lower[l, strat] <- quantile(benefit[,strat], percentage)
    benefit_upper[l, strat] <- quantile(benefit[,strat], (1-percentage))
  }
}

# melt into one dataframe
benefit_mid <- data.frame(wtp, benefit_mid, stringsAsFactors = FALSE)
colnames(benefit_mid) <- c("WTP", strategies)
benefit_mid <- melt(benefit_mid, id.vars = "WTP", variable.name = "Strategy", value.name = "mid")
benefit_mid <- benefit_mid[order(benefit_mid$WTP), ]
rownames(benefit_mid) <- NULL

benefit_lower <- data.frame(wtp, benefit_lower, stringsAsFactors = FALSE)
colnames(benefit_lower) <- c("WTP", strategies)
benefit_lower <- melt(benefit_lower, id.vars = "WTP", variable.name = "Strategy", value.name = "lower")
benefit_lower <- benefit_lower[order(benefit_lower$WTP), ]
rownames(benefit_lower) <- NULL

benefit_upper <- data.frame(wtp, benefit_upper, stringsAsFactors = FALSE)
colnames(benefit_upper) <- c("WTP", strategies)
benefit_upper <- melt(benefit_upper, id.vars = "WTP", variable.name = "Strategy", value.name = "upper")
benefit_upper <- benefit_upper[order(benefit_upper$WTP), ]
rownames(benefit_upper) <- NULL

# combine dataframes
Benefit_CI <- benefit_mid
Benefit_CI[,"lower"] <- benefit_lower[,3]
Benefit_CI[,"upper"] <- benefit_upper[,3]

# benefit plot
benefit_plot <- ggplot(data = Benefit_CI, mapping = aes_(x = Benefit_CI$WTP, y = Benefit_CI$mid, color = Benefit_CI$Strategy)) +
  geom_line() +
  geom_line(aes(x = WTP, y = lower), linetype=2) +
  geom_line(aes(x = WTP, y = upper), linetype=2) +
  ylab("net benefit") + xlab("willingness to pay (Euro)") +
  scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 12) +
  theme(legend.position="none")
benefit_plot

ggsave("net_Benefit_CI.pdf", width = 8, height = 5)
