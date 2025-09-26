library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)  

state_mapping <- c(
  '1' = 'AL',
  '2' = 'AK',
  '4' = 'AZ',
  '5' = 'AR',
  '6' = 'CA',
  '8' = 'CO',
  '9' = 'CT',
  '10' = 'DE',
  '11' = 'DC',
  '12' = 'FL',
  '13' = 'GA',
  '15' = 'HI',
  '16' = 'ID',
  '17' = 'IL',
  '18' = 'IN',
  '19' = 'IA',
  '20' = 'KS',
  '21' = 'KY',
  '22' = 'LA',
  '23' = 'ME',
  '24' = 'MD',
  '25' = 'MA',
  '26' = 'MI',
  '27' = 'MN',
  '28' = 'MS',
  '29' = 'MO',
  '30' = 'MT',
  '31' = 'NE',
  '32' = 'NV',
  '33' = 'NH',
  '34' = 'NJ',
  '35' = 'NM',
  '36' = 'NY',
  '37' = 'NC',
  '38' = 'ND',
  '39' = 'OH',
  '40' = 'OK',
  '41' = 'OR',
  '42' = 'PA',
  '44' = 'RI',
  '45' = 'SC',
  '46' = 'SD',
  '47' = 'TN',
  '48' = 'TX',
  '49' = 'UT',
  '50' = 'VT',
  '51' = 'VA',
  '53' = 'WA',
  '54' = 'WV',
  '55' = 'WI',
  '56' = 'WY',
  '72' = 'PR'
)

## Read in results from experiments
aggr_data <- read.csv("us_opt_samples_agg_all.csv")
individual_data <- read.csv("us_opt_samples_individual_all.csv")


## Prepare df for plotting
aggr_data$MSE.Type <- ifelse(aggr_data$MSE.Type == "Pooled MSE",
                            "Equal Sampling", "Optimal Sampling") 
individual_data$state <- state_mapping[match(individual_data$state, names(state_mapping))]
individual_data$Metric <- NULL
colnames(individual_data) <- colnames(aggr_data)
plot_data <- rbind(aggr_data,individual_data)
plot_data$MSE.Type <- factor(plot_data$MSE.Type,
                             levels = c("CA", "PA", "NY", "FL", "TX", "Equal Sampling", "Optimal Sampling"))


## 650 350
ggplot(plot_data, aes(x = N_seq, y = MSE.Value, color = MSE.Type, group = MSE.Type)) +
  # Add lines
  geom_line(aes(linetype = MSE.Type), linewidth = 0.7) +
  # Add titles and labels
  labs(
    title = 'MSE Under Sampling Constraint (100 Trials)',
    x = 'Total Samples (N)',
    y = 'Average Mean Squared Error',
    color = 'Sampling Method',
    linetype = 'Sampling Method'
  ) +
  scale_color_brewer(
    palette = "Dark2" 
  ) +
  scale_linetype_manual(values = c(
    "Equal Sampling" = "solid", 
    "Optimal Sampling" = "dashed",
    "CA" = "solid", 
    "PA" = "solid", 
    "NY" = "solid", 
    "FL" = "solid", 
    "TX" = "solid"
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )
