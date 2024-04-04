# Set the file path (replace "your_file.txt" with your actual file path)
file_path <- "C:\\Users\\William Teague\\Box\\HERT_Box\\Data\\Aggregate Electron Data 1-100.txt"

# Read the data
data <- read.table(file_path, sep="", skip = 1)

# Assign header labels
column_names <- c("Einc(MeV)", "Detector1", "Detector2", "Detector3",
                  "Detector4","Detector5","Detector6",
                  "Detector7","Detector8","Detector9") 
colnames(data) <- column_names
