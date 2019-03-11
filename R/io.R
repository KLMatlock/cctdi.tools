library('data.table')

#' Read Dataset
#'
#' Reads a dataset either from a list of files or from a single directory.
#'
#' Two options are available. The first option (default) reads all samples and merges them afterward while
#' the second option creates the dataset file by file which is slower but less memory intensive.
#'
#' @param fpath Path to directory of files to read from. If not set then file_list must be set.
#' @param file_list Character vector of filepaths of files to be read from.
#' @param sample_list Character vector, sample names to associate with each sample in file_list
#' @param read_col Column in the provided files to read from.
#' @param concat Whether to read all files at once and then merge (faster but more memory inefficient) or to merge file by file.
#'
#' @return Single dataframes where each col corresponds to the readings from a single sample.
#' @export
#' @importFrom data.table fread
#' @importFrom utils read.table
#'
buildSet <- function(fpath=NULL,file_list=NULL,sample_list= NULL, read_col = "TPM",concat=TRUE){

  if(!is.null(fpath)){
    file_list <- dir(fpath,full.names = T)
    if(length(file_list) == 0){stop(sprintf('Directory %s is empty.',fpath))}
    sample_list <- unlist(lapply(file_list, function(x){return(strsplit(basename(x),'.',fixed=T)[[1]][1])}))

  }else if(!is.null(file_list)){
    if(is.null(sample_list)){
      sample_list <- unlist(lapply(file_list, function(x){return(strsplit(basename(x),'.',fixed=T)[[1]][1])}))
    }else if (length(sample_list) != length(file_list)){
      stop("ERROR: passed file_list does not have corresponding sample name!")
    }
  }

  file_list <- append(file_list,NA)
  sample_list <- append(sample_list,NA)
  if(concat){
    readValid <- function(filepath, sample_name){
      #sample_name <- strsplit(basename(filepath),'.',fixed=TRUE)[[1]][1]
      df <-tryCatch({
        df <-fread(filepath, select = c('gene_id', read_col))
        colnames(df)[colnames(df)==read_col] <- sample_name
        df
      },error=function(cond){
        warning(paste(sample_name, " not a valid file. Skipping..."))
        return(NA)
      })

      return(df)
    }
    data_list <- mapply(readValid,file_list, sample_list)
    data_list <- data_list[!is.na(data_list)]

    dataset = Reduce(function(...) merge(..., by="gene_id"), data_list)
  }
  else{
    dataset <- NULL
    for (i in 1:length(file_list)){
      base_name <- basename(file_list[i])
      full_path <- file_list[i]
      sample_name <- strsplit(base_name,'.',fixed=TRUE)[[1]][1]
      dataset <- tryCatch(
        {
          if(is.null(dataset)){
            dataset <-read.table(full_path, header=TRUE, sep="\t", stringsAsFactors = F)[,c("gene_id", read_col)]

            colnames(dataset)[colnames(dataset)==read_col] <- sample_name
            dataset
          }
          else{
            temp_dataset <-read.table(full_path, header=TRUE, sep="\t", stringsAsFactors = F)[,c("gene_id", read_col)]
            colnames(temp_dataset)[colnames(temp_dataset)=="TPM"] <- sample_name
            dataset <- merge(dataset, temp_dataset, by="gene_id")
          }
        },
        error=function(cond) {
          if(!is.na(full_path)){
            warning(paste(full_path, " is not a valid file. Skipping..."))
          }
          return(dataset)
        }
      )
    }
  }
  return(dataset)
}
