#downlaod the data from repo place it in data folder

if(!file.exists("data/data/decomp_with_covar.xlsx")) {
  download.file("https://figshare.com/ndownloader/files/35496569",
                "data/decomp_with_covar.xlsx")
}


