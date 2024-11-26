library(ggplot2)

## Import Data ##
#source_data <- read.csv("C:\\Users\\William Teague\\Box\\HERT_Box\\Sr90 Testing\\Sr90Y90.csv", header = FALSE, fileEncoding = "UTF-8-BOM")
source_data <- read.csv("C:\\Users\\wzt0020\\Box\\HERT_Box\\Sr90 Testing\\Sr90Y90.csv", header = FALSE, fileEncoding = "UTF-8-BOM")
#source_data <- read.csv("C:\\Users\\Will\\Box\\HERT_Box\\Sr90 Testing\\Sr90Y90.csv", header = FALSE, fileEncoding = "UTF-8-BOM")
colnames(source_data) <- c("Energy (MeV)", "Sr90 Spectrum", "Y90 Spectrum", "Combined Spectrum")

# Plot original data 
ggplot(source_data, aes(x = `Energy (MeV)`)) +
  geom_line(aes(y = `Sr90 Spectrum`, color = "Sr90 Spectrum")) +
  geom_point(aes(y = `Sr90 Spectrum`, color = "Sr90 Spectrum"), size = 2) +
  geom_line(aes(y = `Y90 Spectrum`, color = "Y90 Spectrum")) +
  geom_point(aes(y = `Y90 Spectrum`, color = "Y90 Spectrum"), size = 2) +
  geom_line(aes(y = `Combined Spectrum`, color = "Combined Spectrum")) +
  geom_point(aes(y = `Combined Spectrum`, color = "Combined Spectrum"), size = 2) +
  labs(x = "Energy (MeV)", y = "Counts/sec") +
  ggtitle("Original Sr90Y90 Energy Spectra Data") +
  theme_minimal() + guides(color = 'none')


## Normalize the spectra ##
# Select data above DART Be window threshold and that has actual particle counts
normalized_data <- source_data[source_data$`Energy (MeV)` >= 0.25 & source_data$`Combined Spectrum` != 0,]
# Find total number of particles/sec
total_counts <- sum(source_data[,4])
# Normalize flux energy spectrum
normalized_data[, 2:4] <- normalized_data[, 2:4] / total_counts

# Plot normalized data
ggplot(normalized_data, aes(x = `Energy (MeV)`)) +
  geom_line(aes(y = `Sr90 Spectrum`, color = "Sr90 Spectrum")) +
  geom_point(aes(y = `Sr90 Spectrum`, color = "Sr90 Spectrum"), size = 2) +
  geom_line(aes(y = `Y90 Spectrum`, color = "Y90 Spectrum")) +
  geom_point(aes(y = `Y90 Spectrum`, color = "Y90 Spectrum"), size = 2) +
  geom_line(aes(y = `Combined Spectrum`, color = "Combined Spectrum")) +
  geom_point(aes(y = `Combined Spectrum`, color = "Combined Spectrum"), size = 2) +
  labs(x = "Energy (MeV)", y = "Counts/sec") +
  ggtitle("Normalized Sr90Y90 Energy Spectra Data") +
  theme_minimal() + guides(color = 'none')


## Calculate expected count rates ##
flux_init <- 100/60 #counts per second
d <- 0.6  # distance between front of collimator and radiation source, cm
s <- 5  # radius of radiation source, cm

# If source radius is bigger than FOV, scale
if (s > 0.3 * (d + 3.1)) {
  flux_imod <- flux_init * 1.11^2 / s^2
}else{
  flux_imod <- flux_init
}

# expected flux with calculated conversion ratio
flux <- 2*flux_imod/(1-1/sqrt(1+s^2/(d+6.9)^2))
# multiply by flux spectrum for realistic results
expected_spectrum <- normalized_data
expected_spectrum[, 2:4] <- normalized_data[, 2:4]*flux/.02 # divide by bin width for /MeV

# Plot expected spectrum
ggplot(expected_spectrum, aes(x = `Energy (MeV)`)) +
  geom_line(aes(y = `Combined Spectrum`, color = "Combined Spectrum")) +
  geom_point(aes(y = `Combined Spectrum`, color = "Combined Spectrum"), size = 2) +
  labs(x = "Energy (MeV)", y = "Flux [#/s/sr/cm^2/MeV]") +
  scale_y_log10() +
  ggtitle("Expected Sr90Y90 Energy Spectra") +
  theme_minimal(base_size = 20) + guides(color = 'none')


## Calculate Count Rates ##
# Import geometric factor and effective energies
#geo_EC <- read.table("C:\\Users\\William Teague\\Box\\HERT_Box\\Sr90 Testing\\geofactor_EC_DARTBe.txt", header = FALSE)
geo_EC <- read.table("C:\\Users\\wzt0020\\Box\\HERT_Box\\Sr90 Testing\\geofactor_EC_DARTBe.txt", header = FALSE)
#E_eff <- read.table("C:\\Users\\William Teague\\Box\\HERT_Box\\Sr90 Testing\\effective_energies_DARTBe.txt", header = FALSE)
E_eff <- read.table("C:\\Users\\wzt0020\\Box\\HERT_Box\\Sr90 Testing\\effective_energies_DARTBe.txt", header = FALSE)
E_eff <- round(E_eff, 2)

# Create energy bins
bins <- ncol(geo_EC)
energy_edges <- seq(from = 0, to = 8, length.out = bins + 1)

# Find the indices of energy_edges that match the 'Energy (MeV)' column in expected_spectrum
min_energy <- min(expected_spectrum$`Energy (MeV)`)
max_energy <- max(expected_spectrum$`Energy (MeV)`)
matching_indices <- which(energy_edges >= min_energy & energy_edges <= max_energy)

# Select the corresponding columns
selected_geo_EC <- geo_EC[, matching_indices]
combined_spectrum <- expected_spectrum$`Combined Spectrum`

# Calculate the counts
count_rate <- apply(selected_geo_EC, 1, function(x) sum(x * combined_spectrum * 0.02))

# Filter count rates larger than 1 each 1000 seconds
expected_counts <- data.frame(Energy_Channel_Number = 1:length(count_rate[count_rate>1e-4]),
                 E_eff = E_eff[count_rate>1e-4,],
                 Count_Rate = count_rate[count_rate>1e-4])
colnames(expected_counts) <- c("Energy Channel", "Effective Energy (MeV)", "Count Rate")

# Plot count rate per energy channel
ggplot(expected_counts, aes(x = `Energy Channel`, y = `Count Rate`)) +
  geom_point(size=3) +
  scale_x_continuous(breaks = 1:max(expected_counts$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)", 
                                         breaks = 1:nrow(expected_counts), 
                                         labels = expected_counts$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0,max(expected_counts$`Count Rate`), by=ceiling(max(expected_counts$`Count Rate`)*10)/10/10)) +
  labs(x = "Energy Channel Number", y = "Count Rate") +
  theme_minimal(base_size = 20)

# Calculate the time it would take to reach 10 counts for each energy channel
expected_counts <- data.frame(expected_counts,
                        expected_time = 10/expected_counts$`Count Rate`)
colnames(expected_counts) <- c("Energy Channel", "Effective Energy (MeV)", "Count Rate", "Time to 10 Counts (s)")

# Plot how long it would take each energy channel to reach 10 counts, up to a max time
total_time <- 30*60 # total anticipated time of experiment
ggplot(expected_counts[expected_counts$`Time to 10 Counts (s)` < total_time,],
       aes(x = `Energy Channel`, y = `Time to 10 Counts (s)`/60)) +
  geom_point(size=3) +
  scale_x_continuous(breaks = 1:max(expected_counts$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)", 
                                         breaks = 1:nrow(expected_counts), 
                                         labels = expected_counts$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0,(total_time/60), by=ceiling(max(total_time/60)*1000)/1000/10)) +
  labs(x = "Energy Channel Number", y = "Expected Time to 10 counts (min)") +
  theme_minimal(base_size = 20)


## Hand Calculation ##
d_mm <- d*10  # distance between front of collimator and radiation source, mm
s_mm <- s*10  # radius of radiation source, mm

# Total flux at first collimator tooth
flux_1 <- flux_init/(2*pi*d_mm^2)/(pi*s_mm^2)
# Flux through first collimator tooth
flux_first <- flux_1 *pi^2 *d_mm^3 * (18*atan((-9+s_mm)/d_mm)+18*atan((9+s_mm)/d_mm)
                                   +d_mm*(log(d_mm^2+(-9+s_mm)^2)-log(d_mm^2+(9+s_mm)^2)))
# Total flux at last collimator tooth through first collimator tooth
flux_2 <- flux_first/(2*pi*60^2)/(pi*9^2)
# Flux through last collimator tooth
flux_last <- flux_2 *pi^2 *60^3 * (18*atan((-9+s_mm)/60)+18*atan((9+s_mm)/60)
                                +60*(log(60^2+(-9+s_mm)^2)-log(60^2+(9+s_mm)^2)))

# Create data to put though energy channel summation
instrument <- data.frame(
  # Shifted energy spectrum for after particles go through DART Be window
  Energy_MeV = expected_spectrum$`Energy (MeV)`-0.3,
  # Account for energy binning (not truly continuous)
  Flux_through_HERT = flux_last * expected_spectrum$`Combined Spectrum` * 0.02
)

# Import electron channels
#electron_channels <- read.csv("C:\\Users\\William Teague\\Box\\HERT_Box\\Energy Resolution\\electron_channels_v1.txt", header = FALSE)
electron_channels <- read.csv("C:\\Users\\wzt0020\\Box\\HERT_Box\\Energy Resolution\\electron_channels_v1.txt", header = FALSE)

# Sum bins in each energy channel
calculate_sum_in_range <- function(lower_bound, upper_bound) {
  sum(instrument$Flux_through_HERT[instrument$Energy_MeV >= lower_bound & instrument$Energy_MeV < upper_bound])
}
# Apply the function to each row of electron_channels
hc_count_rate <- apply(electron_channels, 1, function(x) calculate_sum_in_range(x[1], x[2]))

# Create data frame with energy channels, effective energies, and count rates
hc_counts <- data.frame(Energy_Channel_Number = 1:length(hc_count_rate[hc_count_rate>0]), 
                            E_eff = E_eff[hc_count_rate>0,],
                            Count_Rate = hc_count_rate[hc_count_rate>0])
colnames(hc_counts) <- c("Energy Channel", "Effective Energy (MeV)", "Count Rate")

# Plot count rate for each energy channel
ggplot(hc_counts, aes(x = `Energy Channel`, y = `Count Rate`)) +
  geom_point(size=3) +
  scale_x_continuous(breaks = 1:max(hc_counts$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)", 
                                         breaks = 1:nrow(hc_counts), 
                                         labels = hc_counts$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0,max(hc_counts$`Count Rate`), by=ceiling(max(hc_counts$`Count Rate`)*10)/10/10)) +
  labs(x = "Energy Channel Number", y = "Count Rate") +
  theme_minimal(base_size = 20)

# Calculate the time it would take to reach 10 counts for each energy channel
total_time <- 30*60 # total anticipated time of experiment
hc_counts <- data.frame(hc_counts,
                        expected_time = 10/hc_counts$`Count Rate`)
colnames(hc_counts) <- c("Energy Channel", "Effective Energy (MeV)", "Count Rate", "Time to 10 Counts (s)")

# Plot how long it would take each energy channel to reach 10 counts, up to a max time
ggplot(hc_counts[hc_counts$`Time to 10 Counts (s)` < total_time,],
                       aes(x = `Energy Channel`, y = `Time to 10 Counts (s)`/(60))) +
  geom_point(size=3) +
  scale_x_continuous(breaks = 1:max(hc_counts$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)", 
                                         breaks = 1:nrow(hc_counts), 
                                         labels = hc_counts$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0,(total_time/60), by=ceiling(max(total_time/60)*1000)/1000/10)) +
  labs(x = "Energy Channel Number", y = "Expected Time to 10 counts (min)") +
  theme_minimal(base_size = 20)

## Combine method data for plotting ##
# Combine the two data frames
expected_counts$Source <- "Geometric Factor"
hc_counts$Source <- "Hand Calculated"
combined_data <- rbind(expected_counts, hc_counts)

# Create the combined plot for count rate each energy channel, color coded by method
ggplot(combined_data, aes(x = `Energy Channel`, y = `Count Rate`, color = Source, shape = Source)) +
  geom_point(size = 6) +
  scale_x_continuous(breaks = 1:max(combined_data$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)",
                                         breaks = 1:nrow(combined_data),
                                         labels = combined_data$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0, max(combined_data$`Count Rate`), by=ceiling(max(combined_data$`Count Rate`)*10)/10/10)) +
  labs(x = "Energy Channel Number", y = "Count Rate") +
  theme_minimal(base_size = 20)

# Create the combined plot for time to reach 10 counts, color coded by method
ggplot(combined_data[combined_data$`Time to 10 Counts (s)` < total_time,], 
       aes(x = `Energy Channel`, y = `Time to 10 Counts (s)`/(60),
           color = Source, shape = Source)) +
  geom_point(size=6) +
  scale_x_continuous(breaks = 1:max(hc_counts$`Energy Channel`),
                     position = "top",
                     sec.axis = sec_axis(~ ., name = "Effective Energy (MeV)", 
                                         breaks = 1:nrow(hc_counts), 
                                         labels = hc_counts$`Effective Energy (MeV)`)) +
  scale_y_continuous(breaks = seq(0,(total_time/60), by=ceiling(max(total_time/60)*1000)/1000/10)) +
  labs(x = "Energy Channel Number", y = "Expected Time to 10 counts (min)") +
  theme_minimal(base_size = 20)

