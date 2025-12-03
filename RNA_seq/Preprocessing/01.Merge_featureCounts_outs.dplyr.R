# 필요한 패키지
library(dplyr)
library(purrr)

# 1. 모든 featureCounts 파일 경로 찾기
file_paths <- list.files(
  path = ".",
  pattern = "*\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# 2. 각 파일을 읽고 필요한 컬럼만 선택
count_list <- map(file_paths, function(file) {
  sample_name <- basename(dirname(file))  # 샘플명 추출
  
  read.delim(file, skip = 1, stringsAsFactors = FALSE) %>%
    select(Geneid, Chr, Start, End, Strand, Length, 
           count = 7) %>%  # 7번째 컬럼을 count로 rename
    rename_with(~sample_name, count)  # count 컬럼을 샘플명으로 변경
})

# 3. 첫 번째 데이터프레임 시작
merged_counts <- count_list[[1]]

# 4. 나머지 데이터프레임들을 순차적으로 join
for (i in 2:length(count_list)) {
  merged_counts <- merged_counts %>%
    left_join(count_list[[i]] %>% select(Geneid, last_col()), 
              by = "Geneid")
}

# 5. 결과 저장
write.table(merged_counts, "merged_counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
