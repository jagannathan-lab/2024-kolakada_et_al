library(tidyverse)

x <- select(numbers_ptc_LoFbelow0_35_withtotalcount, c(minus1,aa,count)) 

normal_stop_counts <- normal_stop_counts %>%
  rename(minus1=codon,
         counts=minus1)

x<-x %>%
  rename(count_b35 =count) %>% #b35 refers to count for pLI< 0.35
  merge(numbers_ptc_LoFmid_withtotalcount[,c("minus1", "count")], by="minus1") %>%
  rename(count_mid = count) %>% # mid refers to counts between 0.35 and 0.66
  merge(numbers_ptc_LoFover0_66_withtotalcount[,c("minus1", "count")], by="minus1") %>%
  rename(count_o66 =count) %>% #o66 refers to counts with pLI > 0.66
  merge(ptc_counts, by="minus1") %>%
  rename(count_allptc =count) %>% #this refers to the total PTC counts of each amino acid
  merge(normal_stop_counts, by="minus1") %>% #this refers to the total normal stop counts for each amino acid 
  rename(count_ntc = counts)

aa_counts <- x %>% 
  select(-"minus1") %>%
  aggregate(. ~ aa, sum)

# Create the 2x2 contingency table
  contingency_table <- matrix(
    c(aa_counts$count_b35[1],
      sum(aa_counts$count_b35)-aa_counts$count_b35[1],
      aa_counts$count_o66[1],
      sum(aa_counts$count_o66)-aa_counts$count_o66[1]),
    nrow = 2
  )
  
  # Initialize an empty list to store matrices
  contingency_tables <- list()
  
  # Initialize a vector to store p-values
  p_values <- numeric()
  
  # Loop through each row in the data frame
  for (i in 1:nrow(aa_counts)) {
    # Create a contingency table for the current row
    contingency_table <- matrix(
      c(
        aa_counts$count_b35[i],
        sum(aa_counts$count_b35) - aa_counts$count_b35[i],
        aa_counts$count_o66[i],
        sum(aa_counts$count_o66) - aa_counts$count_o66[i]
      ),
      nrow = 2
    )
    
    # Store the matrix in the list
    contingency_tables[[i]] <- contingency_table
    
    # Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table)
    
    # Append the p-value to the vector
    p_values[i] <- fisher_test$p.value
  }
  
  # Add the p-values as a new column in the original data frame
  aa_counts$b35o66_pvalue <- p_values
  
  # Perform Bonferroni correction
  n_tests <- length(aa_counts$b35o66_pvalue) # Number of tests
  bonferroni_p_values <- p.adjust(aa_counts$b35o66_pvalue, method = "bonferroni") # Corrected p-values
  
  # Add the Bonferroni-corrected p-values as a new column
  aa_counts$bon_b35o66_pvalue <- bonferroni_p_values
  
  contingency_table_g <- matrix(
    c(
      aa_counts$count_b35[6],
      sum(aa_counts$count_b35) - aa_counts$count_b35[6],
      aa_counts$count_o66[6],
      sum(aa_counts$count_o66) - aa_counts$count_o66[6]
    ),
    nrow = 2
  )
  
  fisher_test_g <- fisher.test(contingency_table_g)
  
  
   
