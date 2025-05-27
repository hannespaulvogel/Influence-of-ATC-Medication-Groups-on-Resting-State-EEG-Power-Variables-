# Influence-of-ATC-Medication-Groups-on-Resting-State-EEG-Power-Variables-

# set working directory and load environment
cd /slow/projects/eeg_meds/
conda activate envs/main

# start R
R

# attach package to current R session
library(readxl)
library(dplyr)
library(psych)

# open pv434 data
load('data/pv434/pv434/data/PV434.Rdata')

# show eegp duplicates
# - test-retest cases have been assigned twice; edat_eeg and eeg_filedate are not the same
dups = unique(eegp$pseudonym[duplicated(eegp$pseudonym)])
eegp[eegp$pseudonym %in% dups, c('pseudonym', 'edat_eeg', 'eeg_filedate')]

# remove eegp duplicates where edat_eeg and eeg_filedate are not the same
# - 13 cases with minor deviations left; keep them
eegp = eegp[!(eegp$pseudonym %in% dups & eegp$edat_eeg != as.POSIXct(eegp$eeg_filedate, format = '%d.%m.%Y')),]
sum(eegp$edat_eeg != as.POSIXct(eegp$eeg_filedate, format = '%d.%m.%Y'))
eegp[eegp$edat_eeg != as.POSIXct(eegp$eeg_filedate, format = '%d.%m.%Y'), c('pseudonym', 'edat_eeg', 'eeg_filedate')]

# merge with datajoin.add
sum(duplicated(datajoin.add$pseudonym))
sum(duplicated(eegp$pseudonym))
df = inner_join(datajoin.add, eegp, by = 'pseudonym')

### Datensatz bereinigen

# again, show df duplicates
# - test-retest cases have been assigned twice; edat_eeg and eeg_filedate are not the same
dups = unique(df$pseudonym[duplicated(df$pseudonym)])
df[df$pseudonym %in% dups, c('pseudonym', 'edat_eeg', 'eeg_filedate')]

# remove df duplicates where edat_eeg and eeg_filedate are not the same
# - 13 cases with minor deviations left; keep them
df = df[!(df$pseudonym %in% dups & df$edat_eeg != as.POSIXct(df$eeg_filedate, format = '%d.%m.%Y')),]
sum(df$edat_eeg != as.POSIXct(df$eeg_filedate, format = '%d.%m.%Y'))
df[df$edat_eeg != as.POSIXct(df$eeg_filedate, format = '%d.%m.%Y'), c('pseudonym', 'edat_eeg', 'eeg_filedate')]

# remove re-test assessments and identify remaining duplicates
df = df[is.na(df$no),]
dups = unique(df$pseudonym[duplicated(df$pseudonym)])

# keep entry with the shortest time interval between ADULT_MEDA_H_DATUM and eeg_filedate
# any duplicates left - Nope.
df[df$pseudonym %in% dups, c('pseudonym', 'ADULT_MEDA_H_DATUM', 'eeg_filedate')] 
df$interval = difftime(as.POSIXct(df$eeg_filedate, format = '%d.%m.%Y'), df$ADULT_MEDA_H_DATUM, units = "days")
for (i in 1:length(dups)) {
  mininterval = as.numeric(min(abs(df$interval[df$pseudonym %in% dups[i]]), na.rm = TRUE))
  df = df[!(df$pseudonym %in% dups[i] & (df$interval != mininterval | is.na(df$interval))),]
}
sum(duplicated(df$pseudonym))

###EEG Daten voneinander subtrahieren, Interval berechnen-> variablen ausrechnen

# not necessary for the current project: remove cases without snp data
# not necessary for the current project: df = df[!is.na(df$adult_snp.d00364_sampling_id),]

# sanity check: correlate age_eeg (LIFE) with eeg_age
# above 0.9999, one case with deviation of 3 years, one with 1 year,...
# - LIFE apparently carried out some minor database corrections, use age_eeg variable provided by LIFE
cor.test(df$eeg_age, as.numeric(df$age_eeg))
plot(df$eeg_age, as.numeric(df$age_eeg))
head(df[order(abs(df$eeg_age - as.numeric(df$age_eeg)), decreasing = TRUE), c('eeg_filename','eeg_age', 'age_eeg', 'eeg_filedate', 'edat_eeg', 'ADULT_MEDA_H_DATUM')])

# 1) calculate birthyearmonth based on age_eeg and edat_eeg
# 2) calculate exact age at eeg based on birthyearmonth and eeg_filedatetime
df$birthyearmonth = paste0(substr(as.character(df$edat_eeg-df$age_eeg*365.25),1,8), '15')
df$eeg_age_exact = as.numeric(difftime(strptime(df$eeg_filedate, format = "%d.%m.%Y"),
                                  strptime(df$birthyearmonth, format = "%Y-%m-%d"), units="days"))/365.24219052
cor(cbind(df$eeg_age, df$eeg_age_exact,as.numeric(df$age_eeg)))
plot(df$eeg_age_exact, as.numeric(df$age_eeg))

# select variables of interest (voi) ###umbenennen

voi = c('pseudonym','sex','eeg_age_exact','ADULT_MEDA_H_ATC','pow_cz_broadband','pow_cz_delta','pow_cz_theta','pow_cz_alpha','pow_cz_beta','pow_occipital_alpha','pow_occipital_peakfreq','tbr_tbr_cz')
main = df[,voi] ###im neuen dataframe, nur voi
names(main)[1:4] = c('id','sex','age','atc')


### Atc codes - tbr,tbr (teta beta ratio auf cz); alle variablen die eine N01 haben in dem Namen ATC code: # A02...#, -> beginnt dieses Medikament mit N01? - person 1-6 nutzen dieses Medikament /siehe head: #N05...#, nimmt diese Medikament

# identify cases with medication that acts on nervous system
main$N01_anesthetics = grepl('#N01', main$atc)
main$N02_analgesics = grepl('#N02', main$atc)
main$N03_antiepileptics = grepl('#N03', main$atc)
main$N04_antiparkinson = grepl('#N04', main$atc)
main$N05_psycholeptics = grepl('#N05', main$atc)
main$N06_psychoanaleptics = grepl('#N06', main$atc)
main$N07_other = grepl('#N07', main$atc)

# SSRI (https://www.whocc.no/atc_ddd_index/?code=n06ab)
# Monoamine oxidase inhibitors (https://www.whocc.no/atc_ddd_index/?code=n06af)
# Antipsychotics (https://www.whocc.no/atc_ddd_index/?code=n05a)
# Benzodiazepines (https://www.whocc.no/atc_ddd_index/?code=n05ba)
# Psychostimulants, agents used for adhd and nootropcis (https://www.whocc.no/atc_ddd_index/?code=n06b)
main$N02A_Opioids = grepl('#N02A', main$atc)
main$N05A_Antipsychotics = grepl('#N05A', main$atc)
main$N05BA_Benzodiazepines = grepl('#N05BA', main$atc)
main$N05C_Sedatives = grepl('#N05C', main$atc)
main$N06AF_MAOA = grepl('#N06AF', main$atc)
main$N06AB_SSRI = grepl('#N06AB', main$atc)
main$N06B_Stimulants = grepl('#N06B', main$atc)

# how many participants are affected?
conditions = names(main)[grep("^N", names(main))]
for (i in conditions) {
  tab = table(main[[i]])
  controls = ifelse(!is.na(tab["FALSE"]), tab["FALSE"], 0)
  cases = ifelse(!is.na(tab["TRUE"]), tab["TRUE"], 0)
  cat(sprintf("%s: %d controls vs. %d cases\n", i, controls, cases))
}

# resort data frame and save it
main = main[c('id','sex','age','atc',
  'N01_anesthetics','N02_analgesics','N03_antiepileptics','N04_antiparkinson','N05_psycholeptics','N06_psychoanaleptics','N07_other',
  'N02A_Opioids','N05A_Antipsychotics','N05BA_Benzodiazepines','N05C_Sedatives','N06AF_MAOA','N06AB_SSRI','N06B_Stimulants',
  'pow_cz_broadband','pow_cz_delta','pow_cz_theta','pow_cz_alpha','pow_cz_beta','pow_occipital_alpha','pow_occipital_peakfreq','tbr_tbr_cz')]
system('mkdir -p data/processed')
write.table(main,'data/processed/main.txt', sep = '\t', col.names = T, row.names = F, quote = F)

# ===========================================
# === create cases vs. controls variables ===
# ===========================================

###Qc check- Ausreißer herausfinden - rausrechnen

# Load the processed data
main <- read.delim('data/processed/main.txt', sep='\t')

# Define the eight dependent variables to check for outliers
dependent_vars <- c('pow_cz_broadband', 'pow_cz_delta', 'pow_cz_theta', 
                    'pow_cz_alpha', 'pow_cz_beta', 'pow_occipital_alpha', 
                    'pow_occipital_peakfreq', 'tbr_tbr_cz')

# Create a dataframe to store outlier flags
outlier_flags <- data.frame(id = main$id)

# Function to identify extreme outliers (> Q3 + 3*IQR or < Q1 - 3*IQR)
find_extreme_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 3 * iqr
  upper_bound <- q3 + 3 * iqr
  return(x < lower_bound | x > upper_bound)
}

# Apply the function to each dependent variable and store results
for(var in dependent_vars) {
  outlier_flags[[var]] <- find_extreme_outliers(main[[var]])
  # Print summary of outliers for this variable
  num_outliers <- sum(outlier_flags[[var]], na.rm = TRUE)
  cat(sprintf("Variable %s: %d extreme outliers (%.2f%%)\n", 
              var, num_outliers, 100 * num_outliers / nrow(main)))
}

# Identify participants who are outliers on any variable
outlier_flags$any_outlier <- apply(outlier_flags[, -1], 1, any)

# Print total number of outliers
total_outliers <- sum(outlier_flags$any_outlier, na.rm = TRUE)
cat(sprintf("\nTotal participants with extreme outliers: %d (%.2f%%)\n", 
            total_outliers, 100 * total_outliers / nrow(main)))

# Create clean dataset excluding extreme outliers
main_clean <- main[!outlier_flags$any_outlier, ]
cat(sprintf("Clean dataset contains %d participants\n", nrow(main_clean)))

# Optional: Save the list of outlier IDs for reference
outlier_ids <- main$id[outlier_flags$any_outlier]
writeLines(as.character(outlier_ids), "data/processed/outlier_ids.txt")

# Save the cleaned dataset
write.table(main_clean, 'data/processed/main_clean.txt', 
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)







### Dichotome Variablen erstellen

# Load the dataset
main_clean <- read.delim('data/processed/main_clean.txt', sep='\t')

# Create control group (participants takilang none of these medications)
main_clean$control_group <- !(main_clean$N03_antiepileptics | 
                             main_clean$N05_psycholeptics | 
                             main_clean$N06_psychoanaleptics)

# Create case variables for each medication
main_clean$N03_case <- NA
main_clean$N03_case[main_clean$N03_antiepileptics] <- 1
main_clean$N03_case[main_clean$control_group] <- 0

main_clean$N05_case <- NA
main_clean$N05_case[main_clean$N05_psycholeptics] <- 1
main_clean$N05_case[main_clean$control_group] <- 0

main_clean$N06_case <- NA
main_clean$N06_case[main_clean$N06_psychoanaleptics] <- 1
main_clean$N06_case[main_clean$control_group] <- 0

# Count group sizes
n_control <- sum(main_clean$control_group)
n_N03 <- sum(main_clean$N03_antiepileptics)
n_N05 <- sum(main_clean$N05_psycholeptics)
n_N06 <- sum(main_clean$N06_psychoanaleptics)

cat("Control group (takes none of N03, N05, N06):", n_control, "participants\n")
cat("N03 (antiepileptics) group:", n_N03, "participants\n")
cat("N05 (psycholeptics) group:", n_N05, "participants\n")
cat("N06 (psychoanaleptics) group:", n_N06, "participants\n\n")

# Create data subsets
N03_subset <- main_clean[!is.na(main_clean$N03_case), ]
N05_subset <- main_clean[!is.na(main_clean$N05_case), ]
N06_subset <- main_clean[!is.na(main_clean$N06_case), ]

# Simple demographic analysis for each group
demographic_analysis <- function(data, case_var, group_name) {
  cases <- data[data[[case_var]] == 1, ]
  controls <- data[data[[case_var]] == 0, ]
  
  # Age 
  age_case_mean <- mean(cases$age, na.rm=TRUE)
  age_case_sd <- sd(cases$age, na.rm=TRUE)
  age_control_mean <- mean(controls$age, na.rm=TRUE)
  age_control_sd <- sd(controls$age, na.rm=TRUE)
  
  # Sex
  sex_case_male_pct <- 100 * mean(cases$sex == 1, na.rm=TRUE)
  sex_control_male_pct <- 100 * mean(controls$sex == 1, na.rm=TRUE)
  
  cat("===", group_name, "DEMOGRAPHICS ===\n")
  cat("Sample: ", nrow(cases), " cases, ", nrow(controls), " controls\n", sep="")
  cat("Age (cases): ", round(age_case_mean, 1), " ± ", round(age_case_sd, 1), "\n", sep="")
  cat("Age (controls): ", round(age_control_mean, 1), " ± ", round(age_control_sd, 1), "\n", sep="")
  cat("Males (cases): ", round(sex_case_male_pct, 1), "%\n", sep="")
  cat("Males (controls): ", round(sex_control_male_pct, 1), "%\n\n", sep="")
}

# Run the demographic analysis for each group
cat("DEMOGRAPHIC ANALYSES\n")
cat("===================\n\n")

demographic_analysis(N03_subset, "N03_case", "N03 ANTIEPILEPTICS")
demographic_analysis(N05_subset, "N05_case", "N05 PSYCHOLEPTICS")
demographic_analysis(N06_subset, "N06_case", "N06 PSYCHOANALEPTICS")

cat("All demographic analyses completed successfully.\n")





#### safe der Case-controll daten 


# Sicherstellen, dass das Zielverzeichnis existiert
if (!dir.exists("/slow/projects/eeg_meds/results")) {
  dir.create("/slow/projects/eeg_meds/results", recursive = TRUE)
}

# 1. N03 Antiepileptics
write.table(N03_subset, '/slow/projects/eeg_meds/results/N03_case_control.txt',
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
cat("N03 Subset gespeichert: '/slow/projects/eeg_meds/results/N03_case_control.txt'\n")

# 2. N05 Psycholeptics
write.table(N05_subset, '/slow/projects/eeg_meds/results/N05_case_control.txt',
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
cat("N05 Subset gespeichert: '/slow/projects/eeg_meds/results/N05_case_control.txt'\n")

# 3. N06 Psychoanaleptics
write.table(N06_subset, '/slow/projects/eeg_meds/results/N06_case_control.txt',
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
cat("N06 Subset gespeichert: '/slow/projects/eeg_meds/results/N06_case_control.txt'\n")

cat("\nAlle Datensätze erfolgreich im results-Ordner gespeichert.\n")


# ==============================================================================================================
# === run Lin models(lm) for each medication variable (3 binary variables) x EEG-variables (8 continuous variables) ===
# ==============================================================================================================

 # This shows how you would update one row - you'll expand this into a loop
# Example:
for (group in c("N03", "N05", "N06")) {
	for (eeg_var in c("pow_cz_broadband", "pow_cz_delta", "pow_cz_theta","pow_cz_alpha", "pow_cz_beta", "pow_occipital_alpha","pow_occipital_peakfreq", "tbr_tbr_cz")) {
		result <- run_regression(group, eeg_var)
		row_idx <- which(results$medication == group & results$eeg_variable == eeg_var)
		results$t_value[row_idx] <- result$t
		results$p_value[row_idx] <- result$p
		results$df[row_idx] <- result$df
		results$partial_eta_squared[row_idx] <- result$eta 
	}
}



# Print the structure of the results dataframe
head(results)



### Tabelle per knittr Darstellen

# Prüfen, ob knitr installiert ist, und installieren falls nötig
if (!requireNamespace("knitr", quietly = TRUE)) {
  install.packages("knitr")
}

# Laden der knitr-Bibliothek
library(knitr)

# Formatieren der Tabelle mit kable
pretty_table <- kable(results, 
                     digits = 4,                   # Rundung auf 4 Nachkommastellen
                     format = "simple",            # Einfaches Tabellenformat
                     caption = "Regression Results: Medication Effects on EEG Variables")

# Anzeige der formatierten Tabelle
print(pretty_table)

# Optional: Speichern der formatierten Tabelle als Markdown-Datei
if (!dir.exists("results")) {
  dir.create("results")
}
sink("results/formatted_regression_results.md")
print(pretty_table)
sink()
cat("Formatierte Tabelle wurde in 'results/formatted_regression_results.md' gespeichert.\n")

# Optional: Speichern als HTML für noch bessere Darstellung
if (requireNamespace("knitr", quietly = TRUE) && requireNamespace("rmarkdown", quietly = TRUE)) {
  rmarkdown::render(text = c(
    "# Regression Results",
    "",
    "```{r, echo=FALSE}",
    "library(knitr)",
    "kable(results)",
    "```"
  ), output_file = "results/regression_results.html")
  cat("HTML-Version wurde in 'results/regression_results.html' gespeichert.\n")
}

