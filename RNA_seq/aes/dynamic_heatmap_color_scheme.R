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




##===========================================

# ==============================================================================
# 함수명: draw_mat3_perfect_heatmap
# 목적: 0.7 이상의 유의미한 상관관계에 해상도를 집중시키고, 아웃라이어는 블루로 격리
# 특징: 0.4~0.7 구간을 화이트로 비워 시각적 노이즈를 완벽히 제거함
# ==============================================================================
draw_mat3_perfect_heatmap <- function(mat) {
  require(pheatmap)
  
  # 1. 데이터의 전체 범위 파악
  curr_min <- min(mat, na.rm = TRUE)
  
  # 2. 색상 구간 분리 (Breaks)
  # Low 구간 (Min ~ 0.7): 아웃라이어 구간
  # High 구간 (0.7 ~ 1.0): 메인 분석 구간
  n_low <- 50
  n_high <- 100
  
  breaks_low <- seq(curr_min, 0.699, length.out = n_low)
  breaks_high <- seq(0.7, 1.0, length.out = n_high)
  all_breaks <- c(breaks_low, breaks_high)
  
  # 3. 이중 그라데이션 컬러 설계 (Colors)
  # Low: 진한 블루(#299FFF) -> 화이트(비중 높임) -> 아이보리
  # High: 아이보리 -> 살구색(#FCD4CC) -> 레드
  colors_low <- colorRampPalette(c("#299FFF", "white", "white", "ivory"))(length(breaks_low) - 1)
  colors_high <- colorRampPalette(c("ivory", "#FCD4CC", "red"))(length(breaks_high) - 1)
  
  all_colors <- c(colors_low, colors_high)
  
  # 4. 히트맵 출력 (순서 고정 및 숫자 표기)
  pheatmap(mat, 
           color = all_colors, 
           breaks = all_breaks,
           cluster_rows = FALSE,      # 샘플 순서 고정
           cluster_cols = FALSE, 
           display_numbers = TRUE,    # 상관계수 숫자 표기
           number_color = "black",
           fontsize_number = 6.5,
           number_format = "%.2f",
           main = "Custom Scaled Correlation Heatmap (Focus > 0.7)")
}

# 사용법:
# draw_mat3_perfect_heatmap(mat3)
