##----------------------------------------------------------------

library(ggplot2)
library(scales)
library(patchwork)

set.seed(42)

# ----------------------------------------------------------------
# 1. 데이터 생성 (상관관계 강제 조절)
# ----------------------------------------------------------------

# --- Omics 1: 높은 상관관계 (0.7 ~ 0.99) ---
# 모든 샘플이 공유하는 강력한 Latent Factor를 생성하여 상관관계를 높임
shared_signal1 <- rnorm(200, mean = 5, sd = 1.5)
# 각 샘플은 shared_signal과 약간의 노이즈를 섞음
omics1_raw <- t(sapply(1:10, function(i) {
  shared_signal1 + rnorm(200, sd = 0.5) # 노이즈가 작을수록 상관관계가 높아짐
}))
rownames(omics1_raw) <- paste0("Gene_S", 1:10)
colnames(omics1_raw) <- paste0("G", 1:200)


# --- Omics 2: 높은 상관관계 + 아웃라이어 2개 (0.35 ~ 0.5) ---
shared_signal2 <- rt(200, df = 3) * 5
omics2_raw <- t(sapply(1:10, function(i) {
  if (i <= 8) {
    # 일반 샘플 8개: 매우 높은 상관관계
    shared_signal2 + (rt(200, df = 3) * 0.5)
  } else {
    # 아웃라이어 샘플 2개: 신호 비중을 낮추고 노이즈를 키워 상관관계 낮춤
    (shared_signal2 * 0.3) + (rt(200, df = 3) * 4)
  }
}))
rownames(omics2_raw) <- paste0("Prot_S", 1:10)
colnames(omics2_raw) <- paste0("P", 1:200)

# ----------------------------------------------------------------
# 2. 분석 및 시각화 준비 (함수는 이전과 동일)
# ----------------------------------------------------------------

get_corr_df <- function(mat, label) {
  c_mat <- cor(t(mat), method = "pearson")
  a_mat <- abs(c_mat)
  df <- as.data.frame(as.table(a_mat))
  colnames(df) <- c("Sample1", "Sample2", "Corr")
  df$OmicsLabel <- label
  return(df)
}

df1 <- get_corr_df(omics1_raw, "Omics 1 (High Corr)")
df2 <- get_corr_df(omics2_raw, "Omics 2 (Outliers)")

common_fill_scale <- scale_fill_gradient(
  low = "white", high = "#D73027", 
  limits = c(0, 1), breaks = seq(0, 1, 0.2),
  name = "Abs Correlation", oob = squish
)

# ----------------------------------------------------------------
# 3. 개별 그래프 생성 (텍스트 색상 로직 포함)
# ----------------------------------------------------------------

p1 <- ggplot(df1, aes(Sample1, Sample2, fill = Corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Corr)), 
            size = 3.2, 
            color = ifelse(df1$Corr > 0.6, "white", "black")) +
  common_fill_scale +
  coord_equal() +
  labs(title = "Omics 1: 0.7 - 0.99 Range", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(df2, aes(Sample1, Sample2, fill = Corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Corr)), 
            size = 3.2, 
            color = ifelse(df2$Corr > 0.6, "white", "black")) +
  common_fill_scale +
  coord_equal() +
  labs(title = "Omics 2: High with 2 Outliers", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ----------------------------------------------------------------
# 4. 결과 출력
# ----------------------------------------------------------------
(p1 + p2) + plot_layout(guides = 'collect')
