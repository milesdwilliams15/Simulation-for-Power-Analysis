# Mali Treatment Simulation
# -------------------------

# Open libraries
library(dplyr)    
library(ggplot2)  

# Write village experiment simulation function
mali <- function(ind_sd=1,grp_sd=1,beta=0.25,iterations=1000,
                 villages=187,pop_min=100,pop_max=1000,sample_size=10){
        # Values denote: 
        # ind_sd = individual level standard deviation
        # grp_sd = group level standard deviation
        # beta = treatment effect
        # iterations = number of simulations
        # villages = number of villages (groups) in population
        # pop_min = minimum possible village (group) population
        # pop_max = maximum possible village (group) population
        # sample_size = number of villagers (individuals) sampled per 
        #   village (group)
  
  # Required Packages
  require(dplyr)
  require(ICC)
  require(lmtest)
  require(sandwich)
  
  m <- iterations
  b_hat1 <- c()  # Creat empty vectors to be filled with 
  b_hat2 <- c()  # treatment effect estimates and p-values.
  b_hat3 <- c()  # Vectors created for each model.
  b_hat4 <- c()
  p_val1 <- c()
  p_val2 <- c()
  p_val3 <- c()
  p_val4 <- c()
  ICC <- list()
  for(i in 1:m){
    # Generate village populations for n villages
    pop <- sample(pop_min:pop_max,size=villages,replace=T) # Populations range 
    # and are generated from a uniform
    # distribution.
    
    # Generate baseline measures per villager
    # allowing for both individual level and village
    # level variation in measures. SD = 1 at both
    # individual and village level with mean = 0.
    baseline <- list()  # Create empty list to be filled with 187 seperate vectors
    # of villager level observations per village.
    for(j in 1:villages){ 
      baseline[[j]] <- rnorm(pop[j],sd=ind_sd) + rnorm(1,sd=grp_sd) # Villager and village noise, respectively.
    }
    
    # Generate endline values following treatment,
    # again allowing for individual and village
    # level variation.
    endline <- list()
    Tr <- sample(0:1,replace=T,size=villages) # Treatment assignment indicator.
    for(j in 1:villages){           
      endline[[j]] <- baseline[[j]] + beta*Tr[j] + rnorm(pop[j],sd=ind_sd) + 
        rnorm(1, sd=grp_sd) 
      # Villager and village noise, respectively
    }
    
    # Generate data frame per village.
    # Should be a list of 187 data frames
    # per village.
    df <- list()
    for(j in 1:villages){
      df[[j]] <- data.frame(village=j,baseline=baseline[[j]],
                            endline=endline[[j]],Tr=Tr[j]) 
      df[[j]]$baseline_cat <- 0
      df[[j]]$baseline_cat[df[[j]]$baseline>quantile(df[[j]]$baseline,.25)] <- 1
      df[[j]]$baseline_cat[df[[j]]$baseline>quantile(df[[j]]$baseline,.5)] <- 3
      df[[j]]$baseline_cat[df[[j]]$baseline>quantile(df[[j]]$baseline,.75)] <- 4
      df[[j]]$endline_cat <- 0
      df[[j]]$endline_cat[df[[j]]$endline>quantile(df[[j]]$baseline,.25)] <- 1
      df[[j]]$endline_cat[df[[j]]$endline>quantile(df[[j]]$baseline,.5)] <- 3
      df[[j]]$endline_cat[df[[j]]$endline>quantile(df[[j]]$baseline,.75)] <- 4
    }
    
    # Randomly sample x baseline and endline measures
    # from each village 1,000 times and obtain treatment
    # effect estimates and p-values for 1,000 estimations
    # for each of three empirical models (baseline and endline
    # values are the categorical values):
    #   (1)  endline ~ Tr + mean(baseline) (where baseline is village level average among 
    #                                 sampled villagers)
    #   (2)  mean(endline) ~ Tr + mean(baseline) (using village level means for end and baseline)
    #   (3)  endline ~ Tr
    #   (4)  mean(endline)-mean(baseline) ~ Tr
    small_df <- list() # Empty list to be filled with sub-sampled data
    for(j in 1:villages){
      small_df[[j]] <- sample_n(df[[j]],size=sample_size)
      # Resample endline measures.
      small_df[[j]]$endline <- df[[j]]$endline[sample(1:nrow(df[[j]]),size=sample_size)]
    }
    small_df <- do.call('rbind',small_df) # Combine data frames into one data frame
    # for model estimation.
    # Estimate models:
    ## Model 1
    mod1 <- lm(endline_cat ~ Tr + baseline, small_df %>%
                  group_by(village) %>%
                  mutate(baseline=mean(baseline_cat))) 
    ## Model 2
    mod2 <- lm(endline ~ Tr + baseline, small_df %>%            
                  group_by(village) %>%
                  summarize(endline=mean(endline_cat),
                            baseline=mean(baseline_cat),
                            Tr=unique(Tr)))
    ## Model 3
    mod3 <- lm(endline_cat ~ Tr, small_df)  
    
    ## Model 4
    mod4 <- lm(change ~ Tr, small_df %>%
                 group_by(village) %>%
                 summarize(change=mean(endline_cat)-mean(baseline_cat),
                           Tr=unique(Tr)))
    
    # Use coeftest() to generate summary statistics, and use vcovCL() for 
    # generating clustered standard errors (clustered by village).
    modSum1 <- coeftest(mod1,vcovCL(mod1,cluster=small_df$village))
    modSum2 <- coeftest(mod2)
    modSum3 <- coeftest(mod3,vcovCL(mod3,cluster=small_df$village))
    modSum4 <- coeftest(mod4)
    
    # Fill vectors with each of m estimated treatment effects and
    # p-values.
    b_hat1[i] <- modSum1[2,1]
    b_hat2[i] <- modSum2[2,1]
    b_hat3[i] <- modSum3[2,1]
    b_hat4[i] <- modSum4[2,1]
    p_val1[i] <- modSum1[2,4]
    p_val2[i] <- modSum2[2,4]
    p_val3[i] <- modSum3[2,4]
    p_val4[i] <- modSum4[2,4]
    ICC[[i]] <- data.frame(ICC=c(ICCbare(x=village,
                                         y=baseline,
                                         data=small_df),
                                 ICCbare(x=village,
                                         y=endline,
                                         data=small_df)),
                           Measure=c("Baseline",
                                     "Endline"))
  }
  # End loop.
  
  # Organize results:
  
  # 1. Put estimated treatment effect and corresponding p-value
  #    replicates into a data frame.
  Output <- data.frame(b_hat=c(b_hat1,b_hat2,b_hat3,b_hat4),
                       p_val=c(p_val1,p_val2,p_val3,p_val4),
                       Model=rep(c('Model 1\nendline ~ Tr + mean(baseline)',
                                   'Model 2\nmean(endline) ~ Tr + mean(baseline)',
                                   'Model 3\nendline ~ Tr',
                                   'Model 4\nmean(endline) - mean(baseline) ~ Tr'),
                                 each=m))

  # 2. Estimate power (the proportion of times the null is
  #    rejected at the alpha = 0.05 level) for each model.
  Power <- Output %>% 
    group_by(Model) %>%
    summarize(Power=mean(p_val<=.05))
  
  # 3. Estimate the interclass correlation for both baseline and endline data.
  ICC <- do.call(rbind,ICC) %>% 
    group_by(Measure) %>%
    summarize(mean_ICC=mean(ICC))
  
  # Have the function return a list of items 1-3.
  return(list('Replicates_df'=Output,
              'Power'=Power,
              'ICC'=ICC))
}
# End Function

####################################################################

# Run experiment twice, once with ICC = .1 and again with ICC = .2
# Note: grp_sd=.5 roughly corresponds to ICC = .2 and grp_sd=.35 roughly
# corresponds to ICC = .1

## Round 1: village level variance = .35, and effect size = .25
sum_df1 <- mali(grp_sd=.35,iterations=10) 

### Treatment Effects
sum_df1$Replicates_df %>%
  ggplot(aes(b_hat)) + geom_density(size=1) +
  geom_vline(xintercept=.25) + facet_wrap(~Model) +
  labs(x='Estimated Treatment Effect',
       y='Density',
       title='1,000 Simulated Samples from Population') +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5))
### p-values
sum_df1$Replicates_df %>%
  ggplot(aes(p_val)) + geom_density(size=1) +
  geom_vline(xintercept=.05) + facet_wrap(~Model) +
  labs(x='Estimated p-value',
       y='Density',
       title='1,000 Simulated Samples from Population') +
  theme_bw()+
  theme(plot.title = element_text(hjust=.5))
### Power
sum_df1$Power
### ICC
sum_df1$ICC


## Round 2: village level variance = .5, and effect size = .25
sum_df2 <- mali(grp_sd=.5) 

### Treatment Effects
sum_df2$Replicates_df %>%
  ggplot(aes(b_hat)) + geom_density(size=1) +
  geom_vline(xintercept=.25) + facet_wrap(~Model) +
  labs(x='Estimated Treatment Effect',
       y='Density',
       title='1,000 Simulated Samples from Population') +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5))
### p-values
sum_df2$Replicates_df %>%
  ggplot(aes(p_val)) + geom_density(size=1) +
  geom_vline(xintercept=.05) + facet_wrap(~Model) +
  labs(x='Estimated p-value',
       y='Density',
       title='1,000 Simulated Samples from Population') +
  theme_bw() +
  theme(plot.title = element_text(hjust=.5))
### Power
sum_df2$Power
### ICC
sum_df2$ICC
