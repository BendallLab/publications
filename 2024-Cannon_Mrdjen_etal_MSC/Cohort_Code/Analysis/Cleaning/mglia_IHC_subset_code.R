install.packages("data.table", "dplyr", "ggplot2", "viridis")
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
##### import cell table ####
  fname = "/Volumes/BryJC_Cardinal/Microglia_paper_revisions/cell_table/cell_and_objects_table_arcsinh_transformed.csv"
  mglia_IHC_RAW <- data.table::fread(fname)
  # import and add cohort data
  intermediate_cohort_annot <- data.table::fread("/Volumes/BryJC_Cardinal/Microglia_paper_revisions/IHC_info.csv")
  cohort_annot <- data.table::fread("/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/docs/cohort_info/final_merged_cohort_data.csv")
  
  mglia_IHC_annot_intermediate <- mglia_IHC_RAW %>% 
    dplyr::left_join(intermediate_cohort_annot %>% distinct(), by = "fov")
  
  cleaned_cohort_annot <- cohort_annot %>% 
    select(c("Diagnosis","Pathology","Disease_Status","patient_label")) %>% distinct() %>% slice(1:31)
  
  mglia_IHC_annot <- mglia_IHC_annot_intermediate %>% 
    dplyr::left_join(cleaned_cohort_annot, by = c("patient_label"))
  
  # channels
  channels_IHC = c("Iba1", "HLA-DR")


  # PERCENTILE normalize expression values for HLA-DR & Iba1 from 0 to 1
  mglia_IHC_annot_normalized <- mglia_IHC_annot
  normalization_vector <- apply(mglia_IHC_annot[,2:3], 2, function(x) quantile(x, 0.9999, names = F))
  mglia_IHC_annot_normalized[,2:3] <- mglia_IHC_annot_normalized[,2:3] / as.numeric(normalization_vector)
  # check whether you adjusted the range approximately from 0 to 1
  apply(mglia_IHC_annot_normalized[,2:3], 2, max)

  
  # Processed data table
  # export csv
  write.csv(mglia_IHC_annot_normalized, "/Volumes/BryJC_Cardinal/Microglia_paper_revisions/mglia_IHC.csv")

#### Plotting ####
  library(ggplot2)
  library(viridis) # For magma color palette
  
  # initial scatterplot
  mglia_IHC_annot_normalized %>%
    filter(cell_size >= 200L & cell_size <= 10650L) %>%
    filter(!(Point %in% "Caudate_TMA9-1_Pont4")) %>%
    filter(!is.na(Diagnosis)) %>%
    filter(!is.na(Pathology)) %>%
    filter(Disease_Status %in% c("AD", "Healthy")) %>%
    ggplot() +
    aes(x = `HLA-DR`, y = Iba1) +
    geom_point(shape = "circle", size = 1.5, colour = "#112446") +
    scale_x_continuous(trans = "exp") +
    scale_y_continuous(trans = "exp") +
    theme_minimal() +
    facet_grid(Disease_Status ~ Region)
    #facet_wrap(vars(Region, Disease_Status))
  
  # density plot
  mglia_IHC_annot_normalized %>%
    filter(cell_size >= 200L & cell_size <= 10650L) %>%
    filter(!(Point %in% "Caudate_TMA9-1_Pont4")) %>%
    filter(!is.na(Diagnosis)) %>%
    filter(!is.na(Pathology)) %>%
    filter(Disease_Status %in% c("AD", "Healthy")) %>%
    ggplot(aes(x = `HLA-DR`, y = `Iba1`, color = `Iba1`)) +
    geom_density_2d_filled() + 
    scale_x_continuous(trans = "exp") +
    scale_y_continuous(trans = "exp") +
    scale_color_viridis(option = "magma") +
    facet_grid(Disease_Status ~ Region) +
    theme_minimal() +
    labs(x = "HLA-DR", y = "Iba1", title = "Density Plot of HLA-DR vs Iba1") +
    theme(legend.position = "right")
  
  # for cell density w/i graph - personal
  library(ggplot2)
    # Assuming 'data' is your dataframe and it has columns 'x' and 'y' for plotting
    data <- mglia_IHC_annot_normalized %>%
      filter(cell_size >= 200L & cell_size <= 10650L) %>%
      filter(!(Point %in% "Caudate_TMA9-1_Pont4")) %>%
      filter(!is.na(Diagnosis)) %>%
      filter(!is.na(Pathology)) %>%
      filter(Disease_Status %in% c("AD", "Healthy", "Resilient"))
    # Save this data to export for Meelad. -> microglia-IHC_annot_arcsinh_transformed_99normalized_filtered data
    write.csv(data, "/Volumes/BryJC_Cardinal/AD_Resilience_Cohort_v1/Cohort_Images/ADRCohort_v1/cell_table/microglia-IHC_corrected.csv")
    
    # Calculate the geometric center
    center_HLADR <- mean(data$`HLA-DR`)
    center_Iba1 <- mean(data$Iba1)
    
    # Calculate distance from the center for each point
    data$distance <- sqrt((data$`HLA-DR` - center_HLADR)^2 + (data$Iba1 - center_Iba1)^2)
    
    # Invert distance to represent density (smaller distance = higher density)
    # Adding 1 to avoid division by zero
    data$density <- 1 / (data$distance + 1)
    
    # Plot
    plot1 <- ggplot(data, aes(x = `HLA-DR`, y = Iba1)) +
      geom_point(aes(color = ch_density_sc), size = 0.5) + 
      scale_color_viridis(option = "magma") + # Customize colors as needed
      theme_minimal() +
      labs(title = "Scatterplot Colored by Cell Density",
           
           x = "X Axis Label", # Replace with your actual x-axis label
           y = "Y Axis Label", # Replace with your actual y-axis label
           color = "Density") +
      facet_grid(Disease_Status ~ Region) +
      theme_minimal() +
      labs(x = "HLA-DR", y = "Iba1", 
           title = "Density Plot of HLA-DR vs Iba1", 
           subtitle = "Density decreases with distance from the center") +
      theme(legend.position = "right")
  
    # Density Reference - - - - - - - - - - https://slowkow.com/notes/ggplot2-color-by-density/
    # Get density of points in 2 dimensions.
    # @param x A numeric vector.
    # @param y A numeric vector.
    # @param n Create a square n by n grid to compute density.
    # @return The density within each square.
    get_density <- function(x, y, ...) {
      dens <- MASS::kde2d(x, y, ...)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    
    data <- data %>% group_by(Disease_Status, Region) %>% 
      mutate(ch_density = get_density(x = `HLA-DR`, y = Iba1, n = 50)) %>% mutate(ch_density_sc = scale(ch_density, center=FALSE)) %>% ungroup()
    
    # saving function
    picname <- paste0(getwd(),'/plots/','IHC_channel_density_scatterplots')
    ggsave(plot1, file=paste(picname,"png",sep="."), width = 14, height = 7, dpi = 320, bg = "white")
