# select all and run

library(TCGAbiolinks)
library(maftools)
library(rmarkdown)



projects <- c("ACC", "BRCA", "SARC")
working_dir <- getwd()

in_dir <- paste(working_dir, "inputs", sep = "/")
maf_dir <- paste(in_dir, "maf", sep = "/")
out_dir <- paste(working_dir, "outputs", sep = "/")
report_dir <- paste(out_dir, "Report", sep = "/")


for(file in c(in_dir, maf_dir, out_dir, report_dir)) if(!file.exists(file)) dir.create(file)

for(file in projects)
{
  idir <- paste(maf_dir, file, sep = "/")
  if(!file.exists(idir)) dir.create(idir)
  
  odir <- paste(out_dir, file, sep = "/")
  if(!file.exists(odir)) dir.create(odir)
}




for(file in projects)
{
  proj_name <- paste("TCGA", file, sep = "-")
  file_name <- paste(file, "maf.RData", sep = "_")
  file_directory <- paste("C:", file, sep = "/")
  
  new_dir <- paste(maf_dir, file, sep = "/")
  setwd(new_dir)
  
  if(!file.exists(file_name))
  {
    query <- GDCquery(project = proj_name, access = "open",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation")
    
    GDCdownload(query, directory = file_directory)
    
    maf <- GDCprepare(query, save = TRUE, save.filename = file_name,
                      directory = file_directory)
  }
  setwd(working_dir)
}



for(file in projects)
{
  file_name <- paste(file, "maf.RData", sep = "_")
  new_dir <- paste(maf_dir, file, sep = "/")
  setwd(new_dir)
  load(file_name)
  data <- read.maf(data)
  pdf_dir <- paste(out_dir, file, sep = "/")
  png_dir <- paste(out_dir, file, sep = "/")
  
  
  pdf(paste(pdf_dir, paste("oncoplot", "pdf", sep = "."), sep = "/"))
  oncoplot(data, top = 10)
  dev.off()

  png(paste(png_dir, paste("oncoplot", "png", sep = "."), sep = "/"))
  oncoplot(data, top = 10)
  dev.off()


  pdf(paste(pdf_dir, paste("lollipopPlot", "pdf", sep = "."), sep = "/"))
  lollipopPlot(data, gene = "TP53", AACol = "HGVSp_Short")
  dev.off()

  png(paste(png_dir, paste("lollipopPlot", "png", sep = "."), sep = "/"))
  lollipopPlot(data, gene = "TP53", AACol = "HGVSp_Short")
  dev.off()


  pdf(paste(pdf_dir, paste("rainfallPlot", "pdf", sep = "."), sep = "/"))
  rainfallPlot(data)
  dev.off()

  png(paste(png_dir, paste("rainfallPlot", "png", sep = "."), sep = "/"))
  rainfallPlot(data)
  dev.off()


  pdf(paste(pdf_dir, paste("plotVaf", "pdf", sep = "."), sep = "/"))
  plotVaf(data)
  dev.off()

  png(paste(png_dir, paste("plotVaf", "png", sep = "."), sep = "/"))
  plotVaf(data)
  dev.off()
}



setwd(working_dir)
rmarkdown::render(input = "report.Rmd", output_format = "pdf_document",
                  output_file = "Report.pdf", output_dir = report_dir)









