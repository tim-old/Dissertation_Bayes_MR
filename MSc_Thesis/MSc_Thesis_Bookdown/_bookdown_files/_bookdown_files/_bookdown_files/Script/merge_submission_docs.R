library(pdftools)
library(here)
library(purrr)
library(stringr)
#library(doconv)
library(docxtractr)

# from https://stackoverflow.com/questions/52388945/create-pdf-in-addition-to-word-docx-using-officer
# require libreoffice https://www.libreoffice.org/download/download-libreoffice/

# Setup
set_libreoffice_path("C:/Program Files/LibreOffice/program/soffice.exe")

# List word documents
words <- list.files("Data/Thesis_Submission/",
                    pattern = "?.(docx|doc)",
                    full.names = TRUE)

# Custom function - convert word to pdf
word2pdf <- function(path){
  
  ## Extract file name
  name <- stringr::str_remove(path, "Data/Thesis_Submission/") %>% 
    stringr::str_remove(".docx")
  
  docxtractr::convert_to_pdf(path,
                             pdf_file = paste0("Data/Thesis_Submission/PDF/",
                                               name,
                                               ".pdf"))
  
}

# Convert each .doc in list to .pdf
words %>%
  map(~word2pdf(.x))

# Specify pdfs to merge
#title_page <- file.choose() # interactive alternative to select
title_page <- here("Data", "Thesis_Submission", "PDF", "B233241_MSc_DSHSC_title_page.pdf")
orig_statement <- here("Data", "Thesis_Submission", "PDF", "B233241_Declaration_of_Originality.pdf")
thesis <- here("_book", "MSc_Thesis_Bookdown.pdf")


# Merge pdfs
pdf_combine(c(title_page, orig_statement, thesis), output = here("_book", "B233241_dissertation.pdf"))
