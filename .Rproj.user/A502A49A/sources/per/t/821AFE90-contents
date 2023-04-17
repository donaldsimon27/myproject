#Load packages----------
library(tidyverse)
library(stringr)


#Initial Import dataset---------
#Imported using Jesse's pipeline
predlum0 <- readRDS("~/Documents/projectDS/rds/6_dta_symbol_remove.rds")
write.csv(predlum0, file = "~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/predlum.csv")


#Set working directory for session with Jess---------
setwd("~/Desktop/R.Projects/PredictLuminex/Scripts")


#Exported Dataset---------
#Going forward, work with this dataset
predlum <- read.csv("~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/predlum.csv")
dim(predlum)


#Impute OOR values-------------

library(tidyverse)
library(magrittr)
library(utils)

#Ncite's Impuatation script to impute OOR low and high values----------
#Assuming that your data set was exported with all the columns and consist of the "obs_conc" variable.
#Create "rds" and "csv" folders for the output data.

#Run the R functions in this script and then include your data set name
#e.g., Run this R script. Then Run "remove_symbols(dataset_project)" in the console.

#Step 1:
#Use the "remove_symbols" function to classify missing data, add metadata variables to the data set and remove symbols.

remove_symbols <- function(dataName = "dataset_project"){
  #' @importFrom dplyr mutate
  #' @importFrom magrittr %>%
  dataName <- dataName %>%
    mutate(
      meta_truemissing=ifelse(
        obs_conc == ""
        , TRUE
        , FALSE
      )
      , meta_threestar=ifelse(
        grepl("***", obs_conc,fixed=TRUE)
        , TRUE
        , FALSE
      )
      , meta_onestar=ifelse(
        grepl("^\\*\\d+", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , meta_oorgt=ifelse(
        grepl("OOR >", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , meta_oorlt=ifelse(
        grepl("OOR <", obs_conc, perl=TRUE)
        , TRUE
        , FALSE
      )
      , obs_conc_numerics=as.numeric(
        gsub(
          "\\D*(\\d+(?:\\.\\d+)*(?:E[+-]\\d+)*)"
          ,"\\1"
          ,obs_conc,perl=T
        )
      )
    )
  # Save in RDS object.
  saveRDS(dataName, "~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predlum.rds")
}

remove_symbols(dataName = predlum)
predlums1 <- readRDS("~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predlum.rds")

#Step 2:
#Impute the OOR values, oor below (<) and oor above (>) values.
#Assuming the following variables and variable names are in your data set, e.g.,
#"analyte" containing the analyte names.
#"obs_conc_numerics" created in the previouse step.

impute <- function(dataName = "dataset_project") {
  #' @importFrom dplyr arrange count group_by summarise
  #' @importFrom magrittr %>%
  #' @importFrom tidyr pivot_wider
  #' @importFrom utils write.csv
  df_rawminmaxv <- dataName %>%
    group_by(analyte) %>%
    summarise(
      Min = min(obs_conc_numerics, na.rm = TRUE),
      Max = max(obs_conc_numerics, na.rm = TRUE)
    )
  # Join experiment data and the table with the min and max analyte
  # concentrations.
  df_rawminmaxv <- left_join(dataName, df_rawminmaxv)
  # Impute OOR below (<) and OOR above (>) values.
  imp <- df_rawminmaxv %>%
    mutate(
      obs_conc_impute =
        ifelse(
          meta_oorlt == TRUE,
          Min - (Min * 0.0001),      #Min - (Min * 0.0001)
          ifelse(
            meta_oorgt == TRUE,
            Max + ((Max - Min) * 0.0001),        
            obs_conc_numerics
          )
        )
    ) %>%
    select(-Min,-Max)
  # Round Obs.Conc.Impute to two decimal places.
  imp$obs_conc_impute <- round(imp$obs_conc_impute, 2)
  # Save data to disk in RDS format.
  saveRDS(imp, "~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predstep2.rds")
  # Save data to disk in CSV format.
  write.csv(imp, file = "~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predlum.csv")
}

impute(predlums1)
final_imputed <- read_csv("~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/predlum.csv", 
                          col_names = TRUE)
write_rds(final_imputed, "~/Desktop/R.Projects/PredictLuminex/Data/rdsfold/final_imputed.rds")


#STEP 1------------
#Concatenated data points to be removed------------
predlum10 <- final_imputed %>% 
  filter(grepl(pattern = "X", type)) %>% 
  select(c("description", "analyte", "obs_conc_impute")) %>% 
  distinct() %>% 
  tidyr::pivot_wider(names_from = analyte, values_from = obs_conc_impute) %>% 
  filter(description  != "ST01064373" & 
           description != "ST01061789" & 
           description  != "ST01061782"  & 
           description  != "ST01059102" & 
           description  != "ST01046593") 


#STEP 2" Create a vector with the PDs and TPs to be removed-------------
a1 <- c("13031-W04", "13031-W08", "13031-W16", "13031-W24")

#STEP 3: Remove 13031 for correction-----------------
predlum_13031 <- predlum10 %>% 
  filter((description %in% a1)) %>% 
  view()


#STEP 4: Remove 13031 from the dataset---------------
predlum11 <- 
  predlum10 %>% 
  filter(!(description %in% a1))              #349 observations 

#Step 5.1:Mutate "15026 REC" to "15026-W24" as per luminex technician-------- 
predlum11 <- predlum11 %>% 
  mutate(description = ifelse(description == "15026 REC", "15026-W24", description)) 

#STEP 5.2: Create a vector for the observations that were mislabeled------
mislabeled <- c("130016-W24_7084169", "13009-W24", "13010-W08", "13010-W24", "13085-Dx", "13085-W04", 
                "13085-W08" , "13085-W16", "13085-W24", "13107-W16", "13108-W08")

#Step 6: Remove "mislabeled from the dataset --------------
predlum12 <- predlum11 %>% 
  filter(!(description %in% mislabeled))

#Step7: Sort out description, ie. remove detail that follow week no---------
predlum13 <- predlum12 %>% 
  mutate(description = gsub('_.*', '', description))

#Step 8: Replace NULL characters in the dataset by coding them as NA-------
#It is likely that there was not enough specimen left to run these kits 

predlum14 <- predlum13 %>% 
  as.data.frame() %>% 
  replace(.=="NULL", NA)



#Step9: Write csv------------
write_rds(predlum14, file = "~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/predMain.rds")

#Explore the data----------
getwd()
predMain <- readRDS("~/Desktop/R.Projects/PredictLuminex/Data/Exported Data/predMain.rds")
predMain <- predMain %>% mutate_at(2:51, as.numeric)

str(predMain)
class(predMain)
names(predMain)
head(predMain)
tail(predMain)
glimpse(predMain)
skimr::skim(predMain)
summary(predMain)


#Add outcome to dataset-------
data0 <- predMain |> 
  mutate(PID = description) |> 
  select(PID, description, everything()) |>
  separate(PID, c("PID", "Timepoint")) |>
  arrange(description)

data1 <- data0 |> 
  mutate(Outcome = PID) 

data2 <- data1 |> 
  mutate(Outcome = ifelse(Outcome == "11012", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11014", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11024", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11031", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11039", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11040", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11047", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11056", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11058", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11064", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11067", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11078", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "11085", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12006", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12020", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12024", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12025", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12029", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12042", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12043", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12045", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12052", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12054", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12064", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12065", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12067", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12070", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12074", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "12083", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13001", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13009", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13010", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13017", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13026", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13029", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13031", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13033", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13037", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13041", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13051", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13056", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13071", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13083", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13085", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13107", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "13108", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14010", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14025", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14026", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14029", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14031", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14046", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14051", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14053", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14057", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14069", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14091", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14099", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14110", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14112", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "14124", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15013", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15017", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15018", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15022", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15026", "PoorOutcome", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15043", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15063", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15074", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15086", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15091", "Cured", Outcome)) %>% 
  mutate(Outcome = ifelse(Outcome == "15092", "Cured", Outcome)) 


#Analyte name convention, replace hyphens with underscore----------
data3 <- data2 |> 
  rename("IL_22" = "IL-22") |> 
  rename("IL_9" = "IL-9") |> 
  rename("sIL_4R" = "sIL-4R") |> 
  rename("sIL_6R" = "sIL-6R") |> 
  rename("MMP_9" = "MMP-9") |> 
  rename("MMP_2" = "MMP-2")
write.csv(data3, "~/Desktop/R.Projects/PredictDataProcessing/Data/datalongf.csv")

#Pivot longer---
data4 <- data3 |> 
  select(description, PID, Timepoint, Outcome, everything()) |> 
  pivot_longer(cols = 5:54, names_to = "analyte", values_to = "obs_conc")

#Confirm that I can go from wider to longer easily---------
data5 <- data4 |> 
  pivot_wider(names_from = "analyte", values_from = "obs_conc")
write.csv(data5, "~/Desktop/R.Projects/PredictDataProcessing/Data/datawide.csv")


#MissForest Imputation for NA's--------
#The imputation is done on the dataset in the wide format
library(missForest)
data_final <- type.convert(data3, as.is = FALSE) %>%   #CONVERT CHARACTER VECTORS TO FACTORS, MISSForest rejects character vectors
  data.matrix() %>%                                    #mandatory transformation to a matrix for misForest
  missForest(verbose = TRUE)

data_final  |> data_final |> as_tibble()

dem <- data3 |> 
  select(1:3, 54)

data_final1 <- dem |> 
  cbind(data_final$ximp) |> 
  select(-(5:7), -58)
write.csv(data_final1, "~/Desktop/R.Projects/PredictDataProcessing/Data/datawideImp.csv")

