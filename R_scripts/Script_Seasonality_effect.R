#################################################################
### Fucus time -serie : testing for the effetc of seasonality ###
#################################################################

# -------------------
# Set up the analysis 
# -------------------

# Load data, packages and functions
library(vegan)
library(here)
library(tidyverse)
library(ecodist)
library(polynom)
library(purrr)

#store dataframes in alist grouped by amplicon region (16 or 18S)
raw_data=list("16S"=readRDS(here("Data/16S_time-decay_fts_df.RDS")),"18S"=readRDS(here("Data/18S_time-decay_fts_df.RDS")))

#define types of substrates
substrates=c("Water","Rock","Fucus")

#function to find the period given mantel fittec coefficient
find_period=function(x){solve(deriv(polynomial(x$coef[1:3])))}

# define object to store results in 
mantel_res=MRN_res=list()

#define the date group used here (all combinaitons amplicon * substrate)
data_group=expand.grid(substrate=substrates,amplicon=c("16S","18S")) # all type of data used here

# -------------------
# Run the Analysis
# -------------------

#run the analysis over all data groups 
for (data_g in 1:dim(data_group)[1])
{
  # retrieve type of data info
  amplicon=data_group$amplicon[data_g]
  substrate=data_group$substrate[data_g]
    
  # New time distance value: 
  data <- raw_data[[amplicon]] %>% 
    mutate(circular_date_dist=ifelse(date_dist<180,date_dist,365-date_dist))
  
  #retrieve substrate type for each sample
  substrate_type <- data %>% 
    group_by(s1) %>% 
    sample_n(1) %>% 
    select(s1=s1,Sample_type=Sample_type.1)
  
  # Create distance matrices
  bc_dist = data %>% 
    select(s1=s1,s2=s2,bc=bc) %>% 
    spread(key=s2, value=bc, fill=NA) %>% 
    column_to_rownames(var="s1") %>% 
    as.matrix
  
  time_dist = data %>% 
    select(s1=s1,s2=s2,date_dist=date_dist) %>% 
    spread(key=s2, value=date_dist, fill=NA) %>% 
    column_to_rownames(var="s1") %>% 
    as.matrix
  
  circular_time_dist = data %>% 
    select(s1=s1,s2=s2,date_dist=circular_date_dist) %>% 
    spread(key=s2, value=date_dist, fill=NA) %>% 
    column_to_rownames(var="s1") %>% 
    as.matrix
  

    #retreive samples name from this substrate
    s=substrate_type %>% 
      subset(Sample_type==data_group$substrate[data_g]) %>% 
      select(sample=s1) %>% 
      pull
    
    # prepare distanmce matric for mantel and MRM tests
    bc=as.dist(bc_dist[s,s])
    d_time=as.dist(time_dist[s,s])
    d_time2=as.dist(time_dist[s,s]*time_dist[s,s])
    d_circular_time=as.dist(circular_time_dist[s,s])
    
    # Run Mantel test
    mantel_res[[data_g]]=tibble(data=amplicon,substrate=substrate,statistics=names(mantel(bc ~ d_time)[c(1,4,5,6)]),time_difference=mantel(bc ~ d_time)[c(1,4,5,6)],circular_time_difference=mantel(bc ~ d_circular_time)[c(1,4,5,6)],)
    
    # Run MRM approach 
    MRN_res[[data_g]]=tibble(data=amplicon,substrate=substrate,term=c("Intercept","Linear","Quadratic"),coef=MRM(bc ~ d_time+d_time2)$coef[,1],p_val=MRM(bc ~ d_time+d_time2)$coef[,2])

}


# -------------------
# Output the analysis 
# -------------------

# Format result, compute period and export tables. 
Mantel_results=do.call(rbind,mantel_res) # not on paper, as a sanity check 

#Supp Table 
MRM_results=do.call(rbind,MRN_res)
write.csv(MRM_results,"Periodicity_test_results.csv")

# Period values
Period <- MRM_results %>%
            group_by(substrate,data) %>%
            group_split() %>%
            purrr::map(find_period) 

Period=tibble(amplicon=data_group$amplicon,substrate=data_group$substrate,Period=Period)
write.csv(Period,"Periodicity_values.csv")





