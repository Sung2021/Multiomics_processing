draw_mat3_perfect_heatmap <- function(mat) {
  require(pheatmap)
  
  curr_min <- min(mat, na.rm = TRUE)
  
  # 1. 구간 설정 (0.7을 기준으로 나눔)
  # 0.7 이하는 0.4를 기점으로 색 변화를 억제함
  n_low <- 50
  n_high <- 100
  
  breaks_low <- seq(curr_min, 0.699, length.out = n_low)
  breaks_high <- seq(0.7, 1.0, length.out = n_high)
  all_breaks <- c(breaks_low, breaks_high)
  
  # 2. 컬러 팔레트 (핵심 튜닝)
  # colors_low: 진한 블루에서 화이트로 아주 '빨리' 변하게 설정
  # "white"를 뒤쪽에 배치하여 0.4~0.7 구간이 거의 화이트로 보이게 합니다.
  colors_low <- colorRampPalette(c("#299FFF", "white", "white", "ivory"))(length(breaks_low) - 1)
  
  # colors_high: 사용자님의 살구색(#FCD4CC) 포인트 유지
  colors_high <- colorRampPalette(c("ivory", "#FCD4CC", "red"))(length(breaks_high) - 1)
  
  all_colors <- c(colors_low, colors_high)
  
  # 3. 히트맵 그리기
  pheatmap(mat, 
           color = all_colors, 
           breaks = all_breaks,
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 6, # 샘플이 많으므로 폰트 조절
           number_format = "%.2f",
           main = "Sample Correlation: Focused Gradient (Outliers in Blue)")
}

# 실행
draw_mat3_perfect_heatmap(mat3)
