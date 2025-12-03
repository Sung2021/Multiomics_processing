# 라이브러리
library(data.table)

# 1. 모든 featureCounts 파일 경로 찾기
file_paths <- list.files(
  path = ".",  # 현재 디렉토리 (필요시 경로 변경)
  pattern = "*\\.txt$",  # featureCounts 결과 파일 패턴
  recursive = TRUE,
  full.names = TRUE
)

# 2. 첫 번째 파일 읽기 (Geneid, Chr, Start, End, Strand, Length 포함)
count_data <- fread(file_paths[1], skip = 1)  # 헤더 스킵

# 3. 나머지 파일들의 count 컬럼만 추가
for (i in 2:length(file_paths)) {
  temp <- fread(file_paths[i], skip = 1)
  count_data <- cbind(count_data, temp[, 7, with = FALSE])  # 7번째 컬럼이 count
}

# 4. 컬럼명 정리
sample_names <- basename(dirname(file_paths))  # 또는 다른 방식으로 샘플명 추출
colnames(count_data)[7:ncol(count_data)] <- sample_names

# 5. 결과 저장
write.table(count_data, "merged_counts.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
