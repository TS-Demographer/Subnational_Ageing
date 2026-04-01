# Program for Subnational Contributions to Population Ageing: Insights from Australia

# Data is sourced from ABS TableBuilder:
  # Australian Bureau of Statistics (ABS). 2011.  2011 Census of Population and Housing [Census TableBuilder]. Available: https://www.abs.gov.au/statistics/microdata-tablebuilder/tablebuilder (accessed 27 May 2025).
  # Australian Bureau of Statistics (ABS). 2021. 2021 Census of Population and Housing [Census TableBuilder]. Available: https://www.abs.gov.au/statistics/microdata-tablebuilder/tablebuilder (accessed 22 January 2026).

# Please save the data to a folder called "Data" and then set your working directory to where you have saved the folder. 
# The program will then run without any user interaction required.

# Environment setup ####
library(tidyverse)
library(ggrepel)
library(zoo)
library(gridExtra)
library(patchwork)

# Import data #
Years <- c(2011, 2021)
Variables <- c("Qual", "Inc", "Part", "CoB", "Relig", "Urban")

for (y in Years){
  for (v in c(Variables, "Pop", "Pop_age20", "Pop_age30")){
    assign(paste0(v, "_", y),
           read.csv(paste0("Data/", v, "_", y, ".csv"), 
                    sep=",", header=TRUE))
  }
}

# Define functions ####
derivative<-function(vt,vth,h){
  r<-log(vth/vt)/h
  vth2<-vt*exp(r*h/2)
  vdot=r*vth2
  return(vdot)
}
midpoint<-function(vt,vth,h){
  r<-log(vth/vt)/h
  vth2<-vt*exp(r*h/2)
  vth2[is.na(vth2)]<-0
  return(vth2)
}
relderiv <- function(vt, vth, h){
  derivv <- derivative(vt, vth, h)
  relderivv <- derivv/midpoint(vt, vth, h)
  return(relderivv)
}

# Define variables ####
Age_start <- c(30, 30, 20, 20, 20, 20)
Age_start <- c(Age_start, 0, 20, 30)
names(Age_start) <- c(Variables, "Pop", "Pop_age20", "Pop_age30")

for (v in c(Variables, "Pop", "Pop_age20", "Pop_age30")){
  data_2011 <- get(paste0(v, "_2011"))
  data_2021 <- get(paste0(v, "_2021"))
  
  ## Age-specific ####
  age <- c((Age_start[[v]]):100)
  age <- age + 0.5
  age_matrix <- matrix(replicate((ncol(data_2021)-1),age),ncol=(ncol(data_2021)-1))
  
  c_subpop_agespec_t <- matrix(0,nrow=length(age), ncol=0)
  for (i in 1:(ncol(data_2011)-1)){
    c_subpop_agespec_t <- cbind(c_subpop_agespec_t, data_2011[,i+1]/sum(data_2011[,i+1]))
  } 
  c_subpop_agespec_th <- matrix(0,nrow=length(age), ncol=0)
  for (i in 1:(ncol(data_2021)-1)){
    c_subpop_agespec_th <- cbind(c_subpop_agespec_th, data_2021[,i+1]/sum(data_2021[,i+1]))
  } 
  c_subpop_agespec_mid <- midpoint(c_subpop_agespec_t, c_subpop_agespec_th, 10)
  c_subpop_agespec_deriv <- derivative(c_subpop_agespec_t, c_subpop_agespec_th, 10)
  c_subpop_agespec_deriv[is.na(c_subpop_agespec_deriv)]<-0
  
  r_subpop_agespec <- relderiv(data_2011[,-1], data_2021[,-1], 10)
  r_subpop_agespec[is.na(r_subpop_agespec)]<-0
  
  ## Subpopulation ####
  c_subpop_t <- colSums(data_2011[ ,-1, drop = FALSE])/sum(data_2011[ ,-1])
  c_subpop_th <- colSums(data_2021[ ,-1, drop = FALSE])/sum(data_2021[ ,-1])
  c_subpop_mid <- midpoint(c_subpop_t, c_subpop_th, 10)
  c_subpop_mid_matrix <- t(matrix(replicate(length(age), c_subpop_mid), ncol=length(age)))
  
  mean_subpop_t <- colSums(age_matrix*c_subpop_agespec_t) 
  mean_subpop_th <- colSums(age_matrix*c_subpop_agespec_th)
  mean_subpop_mid <- midpoint(mean_subpop_t, mean_subpop_t, 10)
  
  c_subpop_deriv <- derivative(c_subpop_t, c_subpop_th, 10)
  mean_subpop_deriv <- derivative(mean_subpop_t, mean_subpop_th, 10)
  
  r_subpop <- relderiv(colSums(data_2011 [,-1, drop = FALSE]), colSums(data_2021[ ,-1, drop = FALSE]), 10)
  r_subpop_matrix <- t(matrix(replicate(length(age),r_subpop),ncol=length(age)))
  
  ## Total ####
  mean_t <- sum(mean_subpop_t*c_subpop_t)
  mean_th <- sum(mean_subpop_th*c_subpop_th)
  mean_mid <- midpoint(mean_t, mean_th, 10)
  
  r <- relderiv(sum(data_2011[,-1]), sum(data_2021[,-1]), 10)
  
  agespec_piece_df <- colSums(age_matrix*c_subpop_agespec_mid*c_subpop_mid_matrix*(r_subpop_agespec-r_subpop_matrix))
  subpop_piece_df <- mean_subpop_mid*c_subpop_mid*(r_subpop-r)
  
  agespec_piece_total <- sum(age_matrix*c_subpop_agespec_mid*c_subpop_mid_matrix*(r_subpop_agespec-r_subpop_matrix))
  subpop_piece_total <- sum(mean_subpop_mid*c_subpop_mid*(r_subpop-r))
  
  mean_deriv <-  agespec_piece_total + subpop_piece_total
  
  check_1_mean_deriv <- sum(mean_subpop_deriv*c_subpop_mid) + sum(mean_subpop_mid*c_subpop_deriv)
  check_2_mean_deriv <- derivative(mean_t, mean_th, 10)
  
  ## Figures ####
  subpop_df <- rbind(mean_subpop_mid,c_subpop_mid,mean_subpop_mid*c_subpop_mid)
  rownames(subpop_df)[3] <- "mean*c"
  barplot(subpop_df[1,], main = "Subpopulation means", ylim = c(0, 60))
  barplot(subpop_df[2,], main = "Subpopulations as proportion of total")
  barplot(subpop_df[3,], main = "Weighted subpopulation means")
  
  c_df <- as.data.frame(c_subpop_agespec_mid)
  c_df <- cbind(age, c_df)
  colnames(c_df) <- colnames(data_2011)
  c_df_long <- pivot_longer(c_df, cols = -c("Age"),  names_to = "Subpopulation", values_to = "c")
  
  assign(paste0(v, "_plot_c_subpop_agespec"), 
         (ggplot(data=c_df_long, aes(x=Age, y=c))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_", "Age-specific population proportions"))))
  get(paste0(v, "_plot_c_subpop_agespec"))
  
  r_df <- as.data.frame(r_subpop_agespec)
  r_df <- cbind(age, r_df)
  colnames(r_df) <- colnames(data_2011)
  r_df_long <- pivot_longer(r_df, cols = -c("Age"),  names_to = "Subpopulation", values_to = "r")
  
  assign(paste0(v, "_plot_r_subpop_agespec"), 
         (ggplot(data=r_df_long, aes(x=Age, y=r))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_", "Age-specific population growth rates"))))
  get(paste0(v, "_plot_r_subpop_agespec"))
  
  c_multiply_df <- as.data.frame(c_subpop_agespec_mid*c_subpop_mid_matrix)
  c_multiply_df <- cbind(age, c_multiply_df)
  colnames(c_multiply_df) <- colnames(data_2011)
  c_multiply_df_long <- pivot_longer(c_multiply_df, cols = -c("Age"),  names_to = "Subpopulation", values_to = "c")
  
  assign(paste0(v, "_plot_c_multiply"), 
         (ggplot(data=c_multiply_df_long, aes(x=Age, y=c))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_","Age-specific population proportions weighted by subpopulation proportions"))))
  get(paste0(v, "_plot_c_multiply"))
  
  r_difference_df <- as.data.frame(r_subpop_agespec-r_subpop_matrix)
  r_difference_df <- cbind(age, r_difference_df)
  colnames(r_difference_df) <- colnames(data_2011)
  r_difference_df_long <- pivot_longer(r_difference_df, cols = -c("Age"),  names_to = "Subpopulation", values_to = "r")
  
  assign(paste0(v, "_plot_r_difference"), 
         (ggplot(data=r_difference_df_long, aes(x=Age, y=r))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_","Age-specific population growth rates minus subpopulation growth rates"))))
  get(paste0(v, "_plot_r_difference"))
  
  c_multiply_r_difference_df <- as.data.frame((c_subpop_agespec_mid*c_subpop_mid_matrix)*(r_subpop_agespec-r_subpop_matrix))
  c_multiply_r_difference_df <- cbind(age, c_multiply_r_difference_df)
  colnames(c_multiply_r_difference_df) <- colnames(data_2011)
  c_multiply_r_difference_df_long <- pivot_longer(c_multiply_r_difference_df, cols = -c("Age"),  names_to = "Subpopulation", values_to = "r")
  
  assign(paste0(v, "_plot_c_multiply_r_difference"), 
         (ggplot(data=c_multiply_r_difference_df_long, aes(x=Age, y=r))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_","Age-specific population growth rates minus subpopulation growth rates multiplied by age-specific population proportions weighted by subpopulation proportions"))))
  get(paste0(v, "_plot_c_multiply_r_difference"))
  
  age_c_multiply_r_difference_df <- as.data.frame(age_matrix*(c_subpop_agespec_mid*c_subpop_mid_matrix)*(r_subpop_agespec-r_subpop_matrix))
  age_c_multiply_r_difference_simple_df <- as.data.frame(age_matrix*(c_subpop_agespec_mid*c_subpop_mid_matrix)*(r_subpop_agespec-r))
  age_c_multiply_r_difference_df <- cbind(age, age_c_multiply_r_difference_df)
  age_c_multiply_r_difference_simple_df <- cbind(age, age_c_multiply_r_difference_simple_df)
  
  age_c_multiply_r_difference_df_long <- pivot_longer(age_c_multiply_r_difference_df, cols = -c("age"),  names_to = "Subpopulation", values_to = "agespec_alt")
  age_c_multiply_r_difference_simple_df_long <- pivot_longer(age_c_multiply_r_difference_simple_df, cols = -c("age"),  names_to = "Subpopulation", values_to = "agespec")
  
  age_c_multiply_r_difference_df_long_test <- cbind(age_c_multiply_r_difference_simple_df_long, age_c_multiply_r_difference_df_long)
  age_c_multiply_r_difference_df_long <- cbind(age_c_multiply_r_difference_simple_df_long, age_c_multiply_r_difference_df_long$agespec_alt)
  colnames(age_c_multiply_r_difference_df_long) <- c("Age", "Subpopulation", "Agespec", "Agespec_alt")
  
  assign(paste0(v, "_plot_age_c_multiply_r_difference"), 
         (ggplot(data=age_c_multiply_r_difference_df_long, aes(x=Age, y=agespec_alt))+ 
            geom_line(aes(color=Subpopulation))+
            ggtitle(paste0(v, "_", y, "_","Weighted age-specific population growth rates minus subpopulation growth rates multiplied by age-specific population proportions weighted by subpopulation proportions"))))
  get(paste0(v, "_plot_age_c_multiply_r_difference"))
  
  # now the second part of the equation
  subpopulation_piece_df <- as.data.frame(mean_subpop_mid*c_subpop_mid*(r_subpop-r))
  subpopulation_piece_df <- cbind(colnames(data_2011)[-1], subpopulation_piece_df)
  colnames(subpopulation_piece_df) <- c("Subpopulation","values")
  
  assign(paste0(v, "_plot_subpopulation_piece_df"), 
         (ggplot(data=subpopulation_piece_df, aes(x=Subpopulation, y=values)) +
            geom_bar(stat="identity")+
            ggtitle(paste0(v, "_", y, "_","plot_subpopulation_piece_df"))))
  get(paste0(v, "_plot_subpopulation_piece_df"))
  
  subpopulation_r_piece_df <- as.data.frame((r_subpop-r))
  subpopulation_r_piece_df <- cbind(colnames(data_2011)[-1], subpopulation_r_piece_df)
  colnames(subpopulation_r_piece_df) <- c("Subpopulation","values")
  
  assign(paste0(v, "_plot_subpopulation_r_piece_df"), 
         (ggplot(data=subpopulation_r_piece_df, aes(x=Subpopulation, y=values)) +
            geom_bar(stat="identity")+
            ggtitle(paste0(v, "_", y, "_","plot_subpopulation_r_piece_df"))))
  get(paste0(v, "_plot_subpopulation_r_piece_df"))
  
  plot_subpopulation_r_piece_df <-ggplot(data=subpopulation_r_piece_df, aes(x=Subpopulation, y=values)) +
    geom_bar(stat="identity")
  
  assign(paste0(v, "_output_df"), cbind(mean_t, mean_th, mean_deriv, r, agespec_piece_total, subpop_piece_total, mean_subpop_t, mean_subpop_th, mean_subpop_deriv, mean_subpop_mid, c_subpop_mid, r_subpop, agespec_piece_df, subpop_piece_df))
  assign(paste0(v, "_agespec_output_df"), age_c_multiply_r_difference_df_long)
}  

# Manually assign the c values for the national population
prop_pop_mid <- midpoint(sum(Pop_2011[ , 2]), sum(Pop_2021[ , 2]), 10)
prop_pop_age20_mid <- midpoint(sum(Pop_age20_2011[ , 2]), sum(Pop_age20_2021[ , 2]), 10)
prop_pop_age30_mid <- midpoint(sum(Pop_age30_2011[ , 2]), sum(Pop_age30_2021[ , 2]), 10)
Pop_age20_output_df[ , 11] <- prop_pop_age20_mid/prop_pop_mid
Pop_age30_output_df[ , 11] <- prop_pop_age30_mid/prop_pop_mid  

# Figures ####
colour_scheme <- c("#f89540", "#cc4778", "#7e03a8", "#0d0887", "#22919c")

## Figure 1 ####
Variables_full <- c("Education", "Income", "Partnership", "Region of birth", "Religion", "Urbanicity")

Fig_1_df <- data.frame(matrix(0,nrow=0, ncol=5))

for (v in Variables){
  for (y in Years){
    temp <- get(paste0(v, "_", y)) %>% select(!Age)
    total <- sum(temp)
    temp <- temp/total
    age_temp <- get(paste0(v, "_", y)) %>% select(Age)
    
    temp <- cbind(age_temp, temp)
    temp <- pivot_longer(temp, 
                         cols = !Age,
                         names_to = "Group",
                         values_to = "Prop")
    temp <- temp %>% mutate(Year = y,
                            Variable = v)
    
    Fig_1_df <- rbind(Fig_1_df, temp)
  }
}

for (i in unique(Fig_1_df$Variable)){
  Fig_1_df$Variable <- replace(Fig_1_df$Variable, Fig_1_df$Variable==i, Variables_full[which(unique(Fig_1_df$Variable)==i)])
}

Figure_1_individual_plots <- lapply(unique(Fig_1_df$Variable), function(Variables_full) {
  df <- subset(Fig_1_df, Variable==Variables_full)
  df$Group <- factor(df$Group, levels = unique(df$Group))
  
  x_lab <- if (Variables_full %in% c("Region of birth", "Religion", "Urbanicity")) "Age" else NULL
  y_lab <- if (Variables_full %in% c("Education", "Region of birth")) "Proportion" else NULL
  show_strip <- Variables_full %in% c("Partnership", "Urbanicity")
  
  plot <- ggplot(data=df, aes(x=Age, y=Prop, group=Group))+
    geom_line(aes(color=Group), linewidth=1)+
    facet_grid(Year~.)+
    labs(x = x_lab, y = y_lab, color = "Subpopulation") +
    xlim(19, 101) + 
    ylim(0, 0.018) + 
    scale_colour_manual(values = colour_scheme) + 
    ggtitle(Variables_full)+
    theme_bw()+
    theme(plot.title = element_text(face="bold", size=16),
          axis.title = element_text(face="bold", size=16), 
          axis.text = element_text(colour="black", size=10), 
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
          strip.text = element_text(color = "black", face = "bold", size=16), 
          strip.background = element_rect(fill="white"), 
          legend.text = element_text(size = 10, color = "black"), 
          legend.title = element_blank(), 
          legend.position = "inside",
          legend.position.inside = c(0.79,0.86),
          legend.margin=margin(t = 0, unit='cm')
    )
  
  if (!show_strip) {
    plot <- plot + theme(strip.text = element_blank(),
                         strip.background = element_blank())
  }
  
  return(plot)
})

Figure_1 <- wrap_plots(Figure_1_individual_plots, ncol=3)
Figure_1

## Figure 2 ####
Nat_mean_deriv <- cbind(Variables, c(rep(Pop_age30_output_df[3], 2), rep(Pop_age20_output_df[3], 4)))
Nat_mean_deriv <- as.data.frame(Nat_mean_deriv)

Fig_2_df <- data.frame(matrix(0,nrow=0, ncol=5))
for (v in Variables){  
  country_result<- cbind(rep(Variables_full[which(Variables == v)], nrow(get(paste0(v, "_output_df")))), 
                         row.names(get(paste0(v, "_output_df"))), 
                         get(paste0(v, "_output_df"))[,13:14])
  country_result <- as.data.frame(country_result, stringsAsFactors = FALSE)
  row.names(country_result) <- NULL 
  country_result[,3] <- as.numeric(country_result[,3])
  country_result[,4] <- as.numeric(country_result[,4])
  country_result[,5] <- country_result[,3] + country_result[,4]
  country_result[,6] <- as.numeric(Nat_mean_deriv[which(Variables == v),2])
  colnames(country_result)<- c("Variable", "Subpopulation", "Age-specific piece", "Subpopulation piece", "Total", "National derivative")
  country_result$Subpopulation <- factor(country_result$Subpopulation, 
                                         levels = unique(country_result$Subpopulation))
  if (v == "Relig") {
    country_result$Subpopulation <- as.character(country_result$Subpopulation)
    country_result$Subpopulation[country_result$Subpopulation == "Other"] <- "Other_relig"
    country_result$Subpopulation <- factor(country_result$Subpopulation)
  }
  
  Fig_2_df <- rbind(Fig_2_df, country_result)
}

## Adjust for logarithm error ###
Fig_2_df_adj <- Fig_2_df %>% 
  group_by(Variable) %>%
  mutate(Variable_total = sum(Total))
Fig_2_df_adj <- Fig_2_df_adj %>% 
  ungroup %>% 
  mutate(
    `Age-specific piece` = Fig_2_df$`Age-specific piece`*Fig_2_df_adj$`National derivative`/Fig_2_df_adj$Variable_total,
    `Subpopulation piece`= Fig_2_df$`Subpopulation piece`*Fig_2_df_adj$`National derivative`/Fig_2_df_adj$Variable_total,
    `Total`= Fig_2_df$`Total`*Fig_2_df_adj$`National derivative`/Fig_2_df_adj$Variable_total
  ) %>%
  select(!Variable_total)

Fig_2_df <- pivot_longer(Fig_2_df_adj, 
                         cols = c("Age-specific piece", "Subpopulation piece"),  
                         names_to = "Piece", 
                         values_to = "Value")
hline_df <- Fig_2_df |>
  distinct(Variable, `National derivative`) |>
  mutate(
    colour = c("#0d0887", "#0d0887",
               "#7e03a8", "#7e03a8", "#7e03a8", "#7e03a8"),
    linetype = c(5, 5, 6, 6, 6, 6),
    linewidth = 1
  )


Figure_2 <- ggplot(data=Fig_2_df, aes(fill=Piece,x=Subpopulation, y=Value)) +
  geom_bar(stat="identity", position="stack") +
  geom_point(aes(y=Total),stat="identity",position = position_dodge(width = 0),alpha=1,size=3) +
  geom_hline(data = hline_df,aes(
    yintercept = `National derivative`, colour = colour, linetype = linetype, linewidth = linewidth))+
  scale_colour_identity()+ scale_linetype_identity()+ scale_linewidth_identity()+
  geom_hline(yintercept = 0, color = "black", linetype = 1, linewidth = 1) +
  facet_grid(~Variable, scales = "free_x") +
  labs(y="Contribution to the Change in National Mean Age",
       x="Subnational Population") +
  scale_x_discrete(labels = function(x) gsub("Other_relig", "Other", x)) +
  scale_fill_manual(values=c("#cc4778", "#f89540"), 
                    name  ="Component", 
                    breaks=c("Age-specific piece", "Subpopulation piece"),
                    labels=c("Age-specific", "Subnational"))+
  scale_y_continuous(labels = function(y) gsub("-", "\u2013 ", format(y))) + 
  guides(fill=guide_legend(override.aes=list(shape=c(NA))))+
  theme_bw()+
  theme(axis.title = element_text(face="bold", size=16), 
        axis.text = element_text(colour="black", size=14), 
        axis.text.x = element_text(angle = 50, vjust = 1.0, hjust = 1.0), 
        strip.text = element_text(color = "black", face = "bold", size=16), 
        strip.background =element_rect(fill="white"), 
        legend.position = "bottom", 
        legend.text=element_text(size = 12, color = "black"), 
        legend.title = element_text(size = 16, face = "bold"))
Figure_2

## Figure 3 ####
Fig_3_df <- data.frame(matrix(0,nrow=0, ncol=4))
for (x in Variables){
  df <- cbind(rep(Variables_full[which(Variables == x)], nrow(get(paste0(x, "_agespec_output_df")))),
              get(paste0(x, "_agespec_output_df")))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- df %>% select(!"Agespec_alt")
  colnames(df)<- c("Variable", "Age","Subpopulation", "Value")
  Fig_3_df <- rbind(Fig_3_df, df)
}

## Adjust for logarithm error ###
Fig_3_df_adj <- Fig_3_df %>% 
  group_by(Variable) %>%
  mutate(Variable_total = sum(Value)) %>%
  ungroup() %>%
  mutate(National_derivative = 
           as.numeric(c(rep(Nat_mean_deriv[1,2], 497),
                        rep(Nat_mean_deriv[3,2], 1215)
           )))
Fig_3_df <- Fig_3_df_adj %>% mutate(
  Value = Fig_3_df_adj$Value*Fig_3_df_adj$National_derivative/Fig_3_df_adj$Variable_total
)

Figure_3_individual_plots <- lapply(unique(Fig_3_df$Variable), function(Variables_full) {
  df <- subset(Fig_3_df, Variable==Variables_full)
  df$Subpopulation <- factor(df$Subpopulation, levels = unique(df$Subpopulation))
  
  x_lab <- if (Variables_full %in% c("Region of birth", "Religion", "Urbanicity")) "Age" else NULL
  y_lab <- if (Variables_full %in% c("Education", "Region of birth")) "Contribution" else NULL
  
  ggplot(data=df, aes(x=Age, y=Value, group=Subpopulation))+
    geom_line(aes(color=Subpopulation), linewidth=1)+
    geom_hline(yintercept=0, color="grey")+
    labs(x = x_lab, y = y_lab, color = "Subpopulation") +
    ylim(-0.03, 0.041)+
    scale_colour_manual(values=colour_scheme)+
    ggtitle(Variables_full)+
    theme_bw()+
    theme(plot.title = element_text(face="bold", size=16),
          axis.title = element_text(face="bold", size=16), 
          axis.text = element_text(colour="black", size=10), 
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
          strip.text = element_text(color = "black", face = "bold", size=16), 
          strip.background = element_rect(fill="white"), 
          legend.text = element_text(size = 10, color = "black"), 
          legend.title = element_blank(), 
          legend.position = "inside",
          legend.position.inside = c(0.79,0.86),
          legend.margin=margin(t = 0, unit='cm'))
})

Figure_3 <- wrap_plots(Figure_3_individual_plots, ncol=3)
Figure_3

## Table 1 ####
Table_1 <- data.frame(matrix(0,nrow=0, ncol=6))
for (x in c("Pop", "Pop_age20", "Pop_age30", Variables)){
  df <- get(paste0(x, "_output_df")) 
  df <- as.data.frame(df) %>% select(mean_subpop_t, mean_subpop_deriv, c_subpop_mid, r_subpop)
  df$mean_subpop_t <- round(df$mean_subpop_t, 1)
  df$mean_subpop_deriv <- round(df$mean_subpop_deriv, 2)
  df$c_subpop_mid <- round(df$c_subpop_mid*100, 0)
  df$r_subpop <- round(df$r_subpop*100, 1)
  df_1 <- cbind(rep(x, nrow(df)), rownames(df), df)
  Table_1 <- rbind(Table_1, df_1)
}
colnames(Table_1) <- c("Variable", "Subpopulation", "Mean age in 2011", "Annual change in mean", "Proportion of total", "Annual growth rate")
Table_1

# Paper Appendix ####
## Table A5: Component values ####
Table_A5 <- data.frame(matrix(0,nrow=0, ncol=3))
for (x in unique(Fig_2_df_adj$Variable)){
  df <- Fig_2_df_adj %>% filter(Variable == x)
  df <- as.data.frame(df) %>% select(!"National derivative")
  df$Subpopulation <- as.character(df$Subpopulation)
  df_1 <- rbind(df,
                c(paste0(x), "Sum of components", sum(df$`Age-specific piece`), sum(df$`Subpopulation piece`), sum(df$`Total`)))
  df_1 <- df_1 %>%
    mutate(across(c(`Age-specific piece`,
                    `Subpopulation piece`,
                    Total),
                  as.numeric)) 
  Table_A5 <- rbind(Table_A5, df_1)
}

Table_A5_round <- Table_A5 %>% mutate(across(c(`Age-specific piece`,
                                               `Subpopulation piece`,
                                               Total),
                                             ~ round(.x, 2)))
Table_A5_round

## Figure A1: First alternative version of Figure 3 - age-specific piece only ####  
Fig_A1_df <- data.frame(matrix(0,nrow=0, ncol=4))
for (x in Variables){  
  country_result<- cbind(rep(Variables_full[which(Variables == x)], nrow(get(paste0(x, "_agespec_output_df")))), 
                         get(paste0(x, "_agespec_output_df")))
  country_result <- as.data.frame(country_result, stringsAsFactors = FALSE)
  country_result <- country_result %>% select(!Agespec)
  colnames(country_result)<- c("Variable", "Age","Subpopulation", "Value")
  Fig_A1_df <- rbind(Fig_A1_df, country_result)
}  

## Adjust for logarithm error ###
df <- Table_A5 %>% filter(Subpopulation == "Sum of components") %>% select(!c(Subpopulation, `Subpopulation piece`, Total))
colnames(df) <- c("Variable", "True_variable_total")

Fig_A1_df_adj <- Fig_A1_df %>% 
  group_by(Variable) %>%
  mutate(Variable_total = sum(Value)) %>%
  ungroup() %>% 
  left_join(df) 

Fig_A1_df_adj <- Fig_A1_df_adj %>%
  mutate(Value = Fig_A1_df_adj$Value*Fig_A1_df_adj$True_variable_total/Fig_A1_df_adj$Variable_total
  ) %>%
  select(!c(Variable_total, True_variable_total))

Fig_A1_df <- Fig_A1_df_adj

###

Figure_A1_individual_plots <- lapply(unique(Fig_A1_df$Variable), function(Variables_full) {
  df <- subset(Fig_A1_df, Variable==Variables_full)
  df$Subpopulation <- factor(df$Subpopulation, levels = unique(df$Subpopulation))
  
  x_lab <- if (Variables_full %in% c("Region of birth", "Religion", "Urbanicity")) "Age" else NULL
  y_lab <- if (Variables_full %in% c("Education", "Region of birth")) "Contribution" else NULL
  
  ggplot(data=df, aes(x=Age, y=Value, group=Subpopulation))+
    geom_line(aes(color=Subpopulation), linewidth=1)+
    geom_hline(yintercept=0, color="grey")+
    ylim(-0.03, 0.041)+
    scale_colour_manual(values=colour_scheme)+
    labs(x = x_lab, y = y_lab, color = "Subpopulation") +
    ggtitle(Variables_full)+
    theme_bw()+
    theme(plot.title = element_text(face="bold", size=16),
          axis.title = element_text(face="bold", size=16), 
          axis.text = element_text(colour="black", size=10), 
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5), 
          strip.text = element_text(color = "black", face = "bold", size=16), 
          strip.background = element_rect(fill="white"), 
          legend.text = element_text(size = 10, color = "black"), 
          legend.title = element_blank(), 
          legend.position = "inside",
          legend.position.inside = c(0.79,0.86),
          legend.margin=margin(t = 0, unit='cm'))
})

Figure_A1 <- wrap_plots(Figure_A1_individual_plots, ncol=3)
Figure_A1