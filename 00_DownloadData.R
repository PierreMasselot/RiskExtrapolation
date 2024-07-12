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
library(dplyr)
library(tidyr)
library(giscoR)
library(sf)

#----- Some parameters

# Dates covered
dstart <- as.Date("2011-01-01")
dend <- as.Date("2021-12-31")

# Age group limits in the Italian dataset
origagebreaks <- c(1, seq(5, 100, by = 5))

# Target age groups (for the analysis)
agebreaks <- c(0, 45, 65, 75, 85)
agelabs <- c(paste(sprintf("%02i", agebreaks[-length(agebreaks)]), 
  agebreaks[-1] - 1, sep = ""), sprintf("%i+", agebreaks[length(agebreaks)]))

#------------------------
# Download mortality data
#------------------------

#----- Download mortality data

# URL of data (see https://www.istat.it/it/archivio/240401)
linkd <- paste0("https://www.istat.it/storage/dati_mortalita", 
  "/decessi-comunali-giornalieri_4-21062023.zip")

# Download zip and open csv
temp <- tempfile(fileext = ".zip")
tdir <- gsub("\\\\[[:alnum:]]*\\.zip", "", temp)
download.file(linkd, temp)
ftr <- grep("csv", zip_list(temp)$filename, value = T)
unzip(temp, ftr, exdir = tdir)
deathdata <- fread(sprintf("%s\\%s", tdir, ftr), 
  encoding = "Latin-1", na.string = "n.d.")
unlink(temp)
unlink(sprintf("%s\\%s", tdir, ftr))

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
  measure.vars = totalvars, variable.name = "year", value.name = "deaths")

# Rename years
datalong[, year := gsub("T_", "20", year, fixed = T)]

# Aggregate by city, date and age
datatab <- datalong[, .(deaths = sum(deaths, na.rm = T)), 
  keyby = .(CITY_ID, year, GE, CL_ETA)]

# Create date variables
datatab[, month := floor(GE / 100)]
datatab[, day := GE - month * 100]
datatab[, date := as.Date(paste(year, month, day, sep = "-"))]

# Create full list of dates
fullfactor <- as.data.table(expand.grid(city_code = unique(datatab$CITY_ID), 
  date = seq(dstart, dend, by = "day"), agegroup = 0:21))

# Merge with datatab
fulldata <- merge(fullfactor, datatab, 
  by.x = c("city_code", "date", "agegroup"), 
  by.y = c("CITY_ID", "date", "CL_ETA"), all.x = T, all.y = F)

# Fill NAs
fulldata[, deaths := nafill(deaths, fill = 0)]

# Clean date related variables
fulldata[, ":="(year = NULL, month = NULL, day = NULL, GE = NULL)]

# Aggregate age-groups
fulldata[, agegroup := cut(c(0, origagebreaks)[agegroup + 1], c(agebreaks, Inf), 
  right = F, labels = agelabs)]
fulldata <- fulldata[, .(deaths = sum(deaths)), 
  by = .(city_code, date, agegroup)]

# Transform age groups as wide
fulldata <- dcast.data.table(fulldata, city_code + date ~ agegroup, 
  value.var = "deaths")
setnames(fulldata, agelabs, sprintf("deaths_%s", agelabs))

# Arrange
setorder(fulldata, city_code, date)

#----- Save

# Name
fname <- "data/mortality.csv.gz"

# Save only if not found
if (!file.exists(fname)){
  fwrite(fulldata, fname, quote = F, compress = "gzip")
}

#------------------------
# Download metadata for Italy
#------------------------

#----- Metadata

# Download data used in The Lancet Planetary Health paper from Zenodo
# Also include the metadata
# download_zenodo("10.5281/zenodo.7672108", path = tdir,
#   files = list("metadata.csv"))
download_zenodo("10.5281/zenodo.10288665", path = tdir,
  files = list("metadata.csv"))

# Read metadata
metadata <- read.csv(sprintf("%s/metadata.csv", tdir))
unlink(sprintf("%s/metadata.csv", tdir))

# Select only Italy
metadata <- subset(metadata, CNTR_CODE == "IT", -c(CNTR_CODE, cntr_name, region, 
    nmiss, mcc_code, cityname, country, inmcc, LABEL)) |>
  rename(city_code = "URAU_CODE", city_name = "URAU_NAME")

#----- Info about municipalities

# Link to the data
linkg <- paste0("https://www.istat.it/storage/codici-unita-amministrative", 
  "/Elenco-codici-statistici-e-denominazioni-delle-unita-territoriali.zip")

# Download and load into session
temp <- tempfile(fileext = ".zip")
tdir <- gsub("\\\\[[:alnum:]]*\\.zip", "", temp)
download.file(linkg, temp)
ftr <- grep("csv", zip_list(temp)$filename, value = T, useBytes = T)
unzip(temp, ftr, exdir = tdir, junkpaths = T)
geoinfo <- fread(sprintf("%s\\%s", tdir, gsub("^.*/", "", ftr, useBytes = T)), 
  encoding = "Latin-1", na.string = "n.d.")
unlink(temp)
unlink(sprintf("%s\\%s", tdir, gsub("^.*/", "", ftr, useBytes = T)))

# Link to cities, select and rename
geoinfo <- mutate(geoinfo, `LAU CODE` = sprintf("%06d", 
    .data[["Codice Comune formato numerico"]])) |>
  merge(lookup_euro) |> 
  subset(select = c(CITY_ID, `Ripartizione geografica`)) |> 
  rename(geozone = "Ripartizione geografica") |>
  unique()

# Add to metadata
metadata <- merge(metadata, geoinfo, by.x = "city_code", by.y = "CITY_ID")

# Separate spatial and age-related variables
agevars <- grep("[[:alpha:]]+_[[:digit:]]+$", names(metadata), value = T)
metadata_spatial <- metadata[, !names(metadata) %in% agevars]

# Age-related variables into long and characterise age groups
metadata_age <- metadata[, c("city_code", agevars)]
names(metadata_age) <- gsub("([^_])[[:digit:]]{2}$", "\\1", names(metadata_age))
metadata_age <- pivot_longer(metadata_age, cols = !city_code, 
    names_to = c(".value", "age"), names_sep = "_") |>
  arrange(city_code, age) |>
  mutate(age = as.numeric(age), agehigh = lead(age) - 1, .by = city_code) |>
  select(city_code, age, agehigh, prop, deathrate, lifexp)


# Write metadata
fname <- "data/metadata_spatial.csv.gz"
if (!file.exists(fname)){
  fwrite(metadata_spatial, "data/metadata_spatial.csv.gz", 
    quote = F, compress = "gzip")
  fwrite(metadata_age, "data/metadata_age.csv.gz", 
    quote = F, compress = "gzip")
}

#----- Geographical data

# Download Italy map data
italymap <- gisco_get_countries(year = "2020", epsg = "4326", country = "IT")

# Save shapefile
fname <- "data/italymap.shp"
if (!file.exists(fname)){
  st_write(italymap, fname)
}