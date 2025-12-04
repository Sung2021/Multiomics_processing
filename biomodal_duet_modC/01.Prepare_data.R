library(data.table)
library(dplyr)
library(methylKit)

# cxreport 읽기 함수
read_cxreport <- function(file) {
  fread(file, col.names = c("contig", "coordinate", "strand", 
                            "mc", "hmc", "c", 
                            "context", "trinucleotide", "coverage"))
}

# 5mC용 methylKit 형식 변환
convert_5mc_to_methylkit <- function(cxreport_data, sample_id) {
  cg_data <- cxreport_data %>% filter(context == "CG")
  
  mc_format <- cg_data %>%
    mutate(
      chrBase = paste0(contig, ".", coordinate),
      chr = contig,
      base = coordinate,
      coverage_mc = mc + c,
      freqC = (mc / coverage_mc) * 100,
      freqT = 100 - freqC
    ) %>%
    filter(coverage_mc > 0) %>%
    select(chrBase, chr, base, strand, coverage = coverage_mc, freqC, freqT)
  
  output_file <- paste0(sample_id, "_5mC_methylkit.txt")
  fwrite(mc_format, output_file, sep = "\t", quote = FALSE)
  
  return(output_file)
}

# 5hmC용 methylKit 형식 변환
convert_5hmc_to_methylkit <- function(cxreport_data, sample_id) {
  cg_data <- cxreport_data %>% filter(context == "CG")
  
  hmc_format <- cg_data %>%
    mutate(
      chrBase = paste0(contig, ".", coordinate),
      chr = contig,
      base = coordinate,
      coverage_hmc = hmc + c,
      freqC = (hmc / coverage_hmc) * 100,
      freqT = 100 - freqC
    ) %>%
    filter(coverage_hmc > 0) %>%
    select(chrBase, chr, base, strand, coverage = coverage_hmc, freqC, freqT)
  
  output_file <- paste0(sample_id, "_5hmC_methylkit.txt")
  fwrite(hmc_format, output_file, sep = "\t", quote = FALSE)
  
  return(output_file)
}

# 폴더 내 모든 modC 파일 처리
input_dir <- "path/to/your/modc/files"
modc_files <- list.files(input_dir, pattern = "*modC.CXreport.txt.gz$", full.names = TRUE)

# 샘플 ID 추출
sample_ids <- basename(modc_files) %>% 
  gsub("_modC.CXreport.txt.gz", "", .)

# 샘플 이름 기반 treatment 할당
# 10 = treatment2 (3개), 8 = treatment1 (4개), 9 = control (3개)
treatment_vector <- sapply(sample_ids, function(x) {
  if (grepl("10", x)) return(2)      # treatment2
  else if (grepl("8", x)) return(1)  # treatment1
  else if (grepl("9", x)) return(0)  # control
  else return(NA)
})

# 확인
cat("Sample IDs and treatments:\n")
print(data.frame(sample_id = sample_ids, treatment = treatment_vector))

# 모든 파일 변환
mc_files <- c()
hmc_files <- c()

for (i in seq_along(modc_files)) {
  cat("Processing:", sample_ids[i], "\n")
  
  data <- read_cxreport(modc_files[i])
  
  mc_file <- convert_5mc_to_methylkit(data, sample_ids[i])
  hmc_file <- convert_5hmc_to_methylkit(data, sample_ids[i])
  
  mc_files <- c(mc_files, mc_file)
  hmc_files <- c(hmc_files, hmc_file)
}

# methylKit 객체 생성
# 5mC methylRawList
mc_methylkit <- methRead(
  location = as.list(mc_files),
  sample.id = as.list(sample_ids),
  assembly = "hg38",
  treatment = treatment_vector,
  context = "CpG"
)

# 5hmC methylRawList
hmc_methylkit <- methRead(
  location = as.list(hmc_files),
  sample.id = as.list(sample_ids),
  assembly = "hg38",
  treatment = treatment_vector,
  context = "CpG"
)

# 확인
mc_methylkit
hmc_methylkit

# 그룹별 샘플 수 확인
table(treatment_vector)

