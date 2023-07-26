################################################################################
#
#   Risk Extrapolation project
#   
#   Download data
#
################################################################################

#----- Necessary packages
library(zip)
library(data.table)
library(readxl)
library(zen4R)

#----- Some parameters

# Dates covered
dstart <- as.Date("2011-01-01")
dend <- as.Date("2021-12-31")

# Age group limits in the Italian dataset
agebreaks <- c(1, seq(5, 100, by = 5))
agelabs <- paste(sprintf("%02i", c(0, agebreaks)), 
  c(sprintf("%02i", agebreaks - 1), "+"), sep = "")

#------------------------
# Download mortality data
#------------------------

#----- Download mortality data

# URL of data (see https://www.istat.it/it/archivio/240401)
linkd <- paste0("https://www.istat.it/storage/dati_mortalita", 
  "/decessi-comunali-giornalieri_4-21062023.zip")

# Download zip and open csv
temp <- tempfile()
download.file(linkd, temp)
ftr <- grep("csv", zip_list(temp)$filename, value = T)
deathdata <- fread(cmd = sprintf("unzip -p %s %s", temp, ftr), 
  encoding = "Latin-1", na.string = "n.d.")
unlink(temp)

# Pad province codes with zeros
deathdata[, COD_PROVCOM := sprintf("%06d", COD_PROVCOM)]

#----- Download Eurostat city lookup table

# URL of table 
# (see https://ec.europa.eu/eurostat/web/nuts/local-administrative-units)
linke <- paste0("https://ec.europa.eu/eurostat/documents/345175/501971", 
  "/EU-27-LAU-2022-NUTS-2021.xlsx")

# Open xlsx file
temp <- paste0(tempdir(), "\\lookup.xlsx")
download.file(linke, temp, mode = "wb")
lookup_euro <- read_excel(temp, sheet = "IT")
unlink(temp)

# Select variables and remove rural LAUs
lookup_euro <- subset(lookup_euro, !is.na(CITY_ID), 
  c("LAU CODE", "CITY_ID", "CITY_NAME"))

#----- Aggregate mortality by city

# Merge mortality and lookup table (selects only relevant provinces)
deathdata <- merge(deathdata, lookup_euro, 
  by.x = "COD_PROVCOM", by.y = "LAU CODE")

# Select variables (sex: total)
totalvars <- grep("T_[[:digit:]]{2}", names(deathdata), value = T)
years <- sprintf("20%s", substr(totalvars, 3, 4))

# Reshape to long
datalong <- melt.data.table(deathdata, id.vars = c("COD_PROVCOM", "REG", 
    "PROV", "CL_ETA", "GE", "CITY_ID", "CITY_NAME"),
  measure.vars = totalvars, variable.name = "year", value.name = "all")

# Rename years
datalong[, year := gsub("T_", "20", year, fixed = T)]

# Aggregate by city, date and age
datatab <- datalong[, .(all = sum(all, na.rm = T)), 
  keyby = .(CITY_ID, year, GE, CL_ETA)]

# Create date variables
datatab[, month := floor(GE / 100)]
datatab[, day := GE - month * 100]
datatab[, date := as.Date(paste(year, month, day, sep = "-"))]

# Create full list of dates
fullfactor <- as.data.table(expand.grid(CITY_CODE = unique(datatab$CITY_ID), 
  date = seq(dstart, dend, by = "day"), agegroup = 0:21))

# Merge with datatab
fulldata <- merge(fullfactor, datatab, 
  by.x = c("CITY_CODE", "date", "agegroup"), 
  by.y = c("CITY_ID", "date", "CL_ETA"), all.x = T, all.y = F)

# Fill NAs
fulldata[, all := nafill(all, fill = 0)]

# Renaming age groups, cleaning date related variables
fulldata[, ":="(agegroup = agelabs[agegroup + 1], GE = NULL, 
  year = NULL, month = NULL, day = NULL)]

# Save
fwrite(fulldata, "data/mortality.csv.gz", quote = F, compress = "gzip")

#------------------------
# Download metadata for Italy
#------------------------

#----- Metadata

# Download metadata used in The Lancet Planetary Health paper from Zenodo
download_zenodo("10.5281/zenodo.7672108", path = "data",
  files = list("metadata.csv"))

# Read metadata
metadata <- read.csv("data/metadata.csv")

# Select only Italy
metadata <- subset(metadata, CNTR_CODE == "IT", -c(CNTR_CODE, cntr_name, region, 
  nmiss, mcc_code, cityname, country, inmcc)) |>
  rename(CITY_CODE = "URAU_CODE", CITY_NAME = "URAU_NAME")

#----- Info about municipalities

# Link to the data
linkg <- paste0("https://www.istat.it/storage/codici-unita-amministrative", 
  "/Elenco-codici-statistici-e-denominazioni-delle-unita-territoriali.zip")

# Download and load into session
temp <- tempfile()
download.file(linkg, temp)
ftr <- grep("csv", zip_list(temp)$filename, value = T)
geoinfo <- fread(cmd = sprintf("unzip -p %s %s", temp, ftr), 
  encoding = "Latin-1", na.string = "n.d.")
unlink(temp)

# Link to cities, select and rename
geoinfo <- mutate(geoinfo, `LAU CODE` = sprintf("%06d", 
    .data[["Codice Comune formato numerico"]])) |>
  merge(lookup_euro) |> 
  subset(select = c(CITY_ID, `Ripartizione geografica`)) |> 
  rename(geozone = "Ripartizione geografica") |>
  unique()

# Add to metadata
metadata <- merge(metadata, geoinfo, by.x = "CITY_CODE", by.y = "CITY_ID")

# Write metadata
fwrite(metadata, "data/metadata.csv.gz", quote = F, compress = "gzip")
unlink("data/metadata.csv")

#----- Geographical data

# Download Italy map data
italymap <- gisco_get_countries(year = "2020", epsg = "4326", country = "IT")

# Save shapefile
st_write(italymap, "data/italymap.shp")
