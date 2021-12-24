gs4_deauth()
df <- read_sheet("https://docs.google.com/spreadsheets/d/1hLyjPFA2ZRpBLHnTgUnmDz7kimMZWFbz_ZGTml3-hRA/edit#gid=0")

df$DOI <- paste("<a href=",  df$DOI,">", "Link", "</a>") # this converts to clickable format

# subset data rows that are marked processed and deployed
df <- filter(df, `Deployed` == "Yes")
df <- df %>% 
  select(Authors, Year, Journal, DOI, Species, Tissue, Notes, Population, 'Cell#', 'DataID', `Article Title` )
