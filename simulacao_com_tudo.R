## ================================
## Monte Carlo Valuation Experiment
## + CSV export + PDF figures
## + Table 1 (scenarios) + Table 2 (last5 vs all)
## ================================
set.seed(42)

## ---------- Helpers ----------
simulate_path <- function(Tt, mu, sigma, P, F0 = 100, S0 = 1) {
  S <- integer(Tt)
  F <- numeric(Tt)
  S[1] <- S0
  F[1] <- F0
  for (t in 2:Tt) {
    S[t] <- sample.int(3, 1, prob = P[S[t-1], ])
    g    <- mu[S[t]] + rnorm(1, 0, sigma)
    F[t] <- F[t-1] * (1 + g)
  }
  list(F = F, S = S)
}

# PV dos próximos H fluxos (sem perpetuidade) -> estável
value_finite_from_T <- function(F_T, g_hat, r_hat, H = 5) {
  F_proj <- F_T * cumprod(rep(1 + g_hat, H))
  disc   <- cumprod(rep(1 + r_hat, H))
  sum(F_proj / disc)
}

run_experiment <- function(
    scenario_name,
    T_obs = 20, H = 5,
    mu = c(-0.05, 0.05, 0.15),
    sigma = 0.05,
    rS = c(0.12, 0.10, 0.08),
    P = matrix(c(0.6,0.3,0.1, 0.2,0.6,0.2, 0.1,0.3,0.6), 3,3, byrow=TRUE),
    n_sim = 10000L,
    use_last5 = TRUE,
    r_hat = 0.10,
    B_boot = 2000L
) {
  Tt <- T_obs + H
  
  bias_dcf <- numeric(n_sim)
  bias_rs  <- numeric(n_sim)
  ape_dcf  <- numeric(n_sim)
  ape_rs   <- numeric(n_sim)
  Vtrue_store <- numeric(n_sim)
  
  for (s in 1:n_sim) {
    sim <- simulate_path(Tt, mu, sigma, P)
    F <- sim$F
    S <- sim$S
    
    F_obs <- F[1:T_obs]
    S_obs <- S[1:T_obs]
    F_T   <- F_obs[T_obs]
    
    # TRUE value at T_obs = PV dos próximos H fluxos realizados (sem perpetuidade)
    V_true_T <- 0
    disc <- 1
    for (h in 1:H) {
      disc <- disc * (1 + rS[S[T_obs + h]])
      V_true_T <- V_true_T + F[T_obs + h] / disc
    }
    Vtrue_store[s] <- V_true_T
    
    # DCF: média de crescimento (últimos 5 OU toda amostra observada)
    g_all <- diff(F_obs) / F_obs[-T_obs]
    if (use_last5) {
      g_last5 <- diff(F_obs[(T_obs - 5):T_obs]) / F_obs[(T_obs - 5):(T_obs - 1)]
      g_hat_dcf <- mean(g_last5)
    } else {
      g_hat_dcf <- mean(g_all)
    }
    V_dcf <- value_finite_from_T(F_T, g_hat_dcf, r_hat, H)
    
    # RS (toy): média dos mu dos estados observados (proxy simples)
    g_hat_rs <- mean(mu[S_obs])
    V_rs <- value_finite_from_T(F_T, g_hat_rs, r_hat, H)
    
    bias_dcf[s] <- V_dcf - V_true_T
    bias_rs[s]  <- V_rs  - V_true_T
    
    ape_dcf[s] <- abs(bias_dcf[s]) / abs(V_true_T)
    ape_rs[s]  <- abs(bias_rs[s])  / abs(V_true_T)
  }
  
  # Métricas
  bias_dcf_mean <- mean(bias_dcf)
  bias_rs_mean  <- mean(bias_rs)
  
  mae_dcf  <- mean(abs(bias_dcf))
  mae_rs   <- mean(abs(bias_rs))
  
  rmse_dcf <- sqrt(mean(bias_dcf^2))
  rmse_rs  <- sqrt(mean(bias_rs^2))
  
  mape_dcf <- mean(ape_dcf)
  mape_rs  <- mean(ape_rs)
  
  p10_dcf <- mean(ape_dcf > 0.10)
  p10_rs  <- mean(ape_rs  > 0.10)
  
  mae_red  <- 100 * (mae_dcf  - mae_rs)  / mae_dcf
  rmse_red <- 100 * (rmse_dcf - rmse_rs) / rmse_dcf
  mape_red <- 100 * (mape_dcf - mape_rs) / mape_dcf
  
  # Teste pareado na perda quadrática
  tt <- t.test(bias_dcf^2, bias_rs^2, paired = TRUE, alternative = "greater")
  
  # Bootstrap IC das reduções (%)
  set.seed(123)
  n <- n_sim
  rmse_red_boot <- numeric(B_boot)
  mae_red_boot  <- numeric(B_boot)
  mape_red_boot <- numeric(B_boot)
  
  for (b in 1:B_boot) {
    idx <- sample.int(n, size = n, replace = TRUE)
    
    rmse_d <- sqrt(mean(bias_dcf[idx]^2))
    rmse_r <- sqrt(mean(bias_rs[idx]^2))
    rmse_red_boot[b] <- 100 * (rmse_d - rmse_r) / rmse_d
    
    mae_d  <- mean(abs(bias_dcf[idx]))
    mae_r  <- mean(abs(bias_rs[idx]))
    mae_red_boot[b]  <- 100 * (mae_d - mae_r) / mae_d
    
    mape_d <- mean(ape_dcf[idx])
    mape_r <- mean(ape_rs[idx])
    mape_red_boot[b] <- 100 * (mape_d - mape_r) / mape_d
  }
  
  rmse_ci <- quantile(rmse_red_boot, c(0.025, 0.975))
  mae_ci  <- quantile(mae_red_boot,  c(0.025, 0.975))
  mape_ci <- quantile(mape_red_boot, c(0.025, 0.975))
  
  Vtrue_mean <- mean(abs(Vtrue_store))
  
  list(
    scenario = scenario_name,
    use_last5 = use_last5,
    # métricas
    bias_dcf = bias_dcf_mean, bias_rs = bias_rs_mean,
    mae_dcf = mae_dcf, mae_rs = mae_rs,
    rmse_dcf = rmse_dcf, rmse_rs = rmse_rs,
    mape_dcf = mape_dcf, mape_rs = mape_rs,
    p10_dcf = p10_dcf, p10_rs = p10_rs,
    mae_red = mae_red, rmse_red = rmse_red, mape_red = mape_red,
    rmse_ci_l = as.numeric(rmse_ci[1]), rmse_ci_u = as.numeric(rmse_ci[2]),
    mae_ci_l  = as.numeric(mae_ci[1]),  mae_ci_u  = as.numeric(mae_ci[2]),
    mape_ci_l = as.numeric(mape_ci[1]), mape_ci_u = as.numeric(mape_ci[2]),
    Vtrue_mean = Vtrue_mean,
    # vetores para gráficos
    bias_dcf_vec = bias_dcf, bias_rs_vec = bias_rs,
    ape_dcf_vec = ape_dcf, ape_rs_vec = ape_rs,
    # teste
    t_stat = unname(tt$statistic), p_value = tt$p.value
  )
}

to_row <- function(x) {
  data.frame(
    scenario = x$scenario,
    use_last5 = x$use_last5,
    Vtrue_mean = x$Vtrue_mean,
    bias_dcf = x$bias_dcf, bias_rs = x$bias_rs,
    mae_dcf = x$mae_dcf, mae_rs = x$mae_rs, mae_red = x$mae_red,
    rmse_dcf = x$rmse_dcf, rmse_rs = x$rmse_rs, rmse_red = x$rmse_red,
    mape_dcf = x$mape_dcf, mape_rs = x$mape_rs, mape_red = x$mape_red,
    p10_dcf = x$p10_dcf, p10_rs = x$p10_rs,
    rmse_ci_l = x$rmse_ci_l, rmse_ci_u = x$rmse_ci_u,
    mae_ci_l  = x$mae_ci_l,  mae_ci_u  = x$mae_ci_u,
    mape_ci_l = x$mape_ci_l, mape_ci_u = x$mape_ci_u,
    t_stat = x$t_stat, p_value = x$p_value,
    stringsAsFactors = FALSE
  )
}

## ---------- Scenarios ----------
P_persistent <- matrix(c(
  0.90, 0.08, 0.02,
  0.05, 0.90, 0.05,
  0.02, 0.08, 0.90
), 3,3, byrow=TRUE)

P_base <- matrix(c(
  0.6,0.3,0.1,
  0.2,0.6,0.2,
  0.1,0.3,0.6
), 3,3, byrow=TRUE)

P_switchy <- matrix(c(
  0.45, 0.35, 0.20,
  0.25, 0.45, 0.30,
  0.20, 0.35, 0.45
), 3,3, byrow=TRUE)

## ---------- Run experiments (Table 1: 3 scenarios, last5) ----------
res1 <- run_experiment("Alta persistência", P = P_persistent, sigma = 0.05, use_last5 = TRUE)
res2 <- run_experiment("Alta volatilidade", P = P_base,       sigma = 0.10, use_last5 = TRUE)
res3 <- run_experiment("Alta instabilidade (switching)", P = P_switchy, sigma = 0.05, use_last5 = TRUE)

table1_df <- rbind(to_row(res1), to_row(res2), to_row(res3))

## ---------- Table 2: last5 vs all (choose one scenario, e.g., high volatility) ----------
res2_all <- run_experiment("Alta volatilidade", P = P_base, sigma = 0.10, use_last5 = FALSE)
table2_df <- rbind(to_row(res2), to_row(res2_all))
table2_df$scenario <- paste0(table2_df$scenario, ifelse(table2_df$use_last5, " (últimos 5)", " (amostra toda)"))

## ---------- EXPORT results_df to CSV ----------
results_df <- rbind(table1_df, table2_df)
write.csv(results_df, file = "results_df.csv", row.names = FALSE)

# Also export Table 1 and Table 2 separately (optional but convenient)
write.csv(table1_df, file = "table1_scenarios.csv", row.names = FALSE)
write.csv(table2_df, file = "table2_last5_vs_all.csv", row.names = FALSE)

cat("Saved: results_df.csv, table1_scenarios.csv, table2_last5_vs_all.csv\n")

## ---------- EXPORT GRAPHS to PDF ----------
# 1) Histograms of APE for each scenario (DCF vs RS)
pdf("Figure_APE_Histogram1.pdf", width = 11, height = 7)
par(mfrow=c(3,2), mar=c(4,4,2,1))

hist(res1$ape_dcf_vec, breaks=60, main="APE - DCF (Alta persistência)", xlab="APE")
hist(res1$ape_rs_vec,  breaks=60, main="APE - RS  (Alta persistência)", xlab="APE")

hist(res2$ape_dcf_vec, breaks=60, main="APE - DCF (Alta volatilidade)", xlab="APE")
hist(res2$ape_rs_vec,  breaks=60, main="APE - RS  (Alta volatilidade)", xlab="APE")

hist(res3$ape_dcf_vec, breaks=60, main="APE - DCF (Alta instabilidade)", xlab="APE")
hist(res3$ape_rs_vec,  breaks=60, main="APE - RS  (Alta instabilidade)", xlab="APE")

par(mfrow=c(1,1))
dev.off()

# 2) Boxplot comparing APE across methods for each scenario
pdf("Figure_APE_Boxplots11.pdf", width = 11, height = 4)
par(mfrow=c(1,3), mar=c(4,4,2,1))

boxplot(res1$ape_dcf_vec, res1$ape_rs_vec, names=c("DCF","RS"),
        main="Alta persistência", ylab="APE")
boxplot(res2$ape_dcf_vec, res2$ape_rs_vec, names=c("DCF","RS"),
        main="Alta volatilidade", ylab="APE")
boxplot(res3$ape_dcf_vec, res3$ape_rs_vec, names=c("DCF","RS"),
        main="Alta instabilidade", ylab="APE")

par(mfrow=c(1,1))
dev.off()

cat("Saved: Figure_APE_Histograms1.pdf, Figure_APE_Boxplots11.pdf\n")

## ---------- Generate LaTeX Table 1 and Table 2 ----------
fmt <- function(x, d=4) sprintf(paste0("%.", d, "f"), x)
fmt2 <- function(x) sprintf("%.2f", x)

make_table1_latex <- function(df, caption="Tabela 1 -- Resultados Monte Carlo por cenário", label="tab:mc_scenarios") {
  # Keep only last5 rows expected
  df <- df[df$use_last5 == TRUE, ]
  df <- df[, c("scenario","Vtrue_mean","rmse_dcf","rmse_rs","rmse_red","mape_dcf","mape_rs","mape_red","p10_dcf","p10_rs",
               "rmse_ci_l","rmse_ci_u","mape_ci_l","mape_ci_u")]
  out <- c()
  out <- c(out, "\\begin{table}[htbp]")
  out <- c(out, "\\centering")
  out <- c(out, "\\small")
  out <- c(out, "\\begin{threeparttable}")
  out <- c(out, sprintf("\\caption{%s}", caption))
  out <- c(out, sprintf("\\label{%s}", label))
  out <- c(out, "\\begin{tabular}{lrrrrrrrr}")
  out <- c(out, "\\toprule")
  out <- c(out, "Cenário & $\\overline{|V^{true}|}$ & RMSE(DCF) & RMSE(RS) & $\\Delta$RMSE(\\%) & MAPE(DCF) & MAPE(RS) & $\\Delta$MAPE(\\%) & P(APE$>$10\\%) \\\\")
  out <- c(out, "\\midrule")
  
  for (i in 1:nrow(df)) {
    p10_pair <- paste0(fmt2(100*df$p10_dcf[i]), "/", fmt2(100*df$p10_rs[i]))  # DCF/RS in %
    line <- paste0(
      df$scenario[i], " & ",
      fmt(df$Vtrue_mean[i],2), " & ",
      fmt(df$rmse_dcf[i],2), " & ",
      fmt(df$rmse_rs[i],2), " & ",
      fmt2(df$rmse_red[i]), " & ",
      fmt(df$mape_dcf[i],4), " & ",
      fmt(df$mape_rs[i],4), " & ",
      fmt2(df$mape_red[i]), " & ",
      p10_pair, " \\\\"
    )
    out <- c(out, line)
  }
  
  out <- c(out, "\\bottomrule")
  out <- c(out, "\\end{tabular}")
  out <- c(out, "\\begin{tablenotes}")
  out <- c(out, "\\footnotesize")
  out <- c(out, "\\raggedright")
  out <- c(out, "\\item \\textit{Nota:} $\\Delta$RMSE e $\\Delta$MAPE representam redu\\c{c}\\~oes percentuais do m\\'etodo RS em rela\\c{c}\\~ao ao DCF. A coluna P(APE$>$10\\%) reporta DCF/RS (em \\%).")
  out <- c(out, "\\end{tablenotes}")
  out <- c(out, "\\end{threeparttable}")
  out <- c(out, "\\end{table}")
  paste(out, collapse="\n")
}

make_table2_latex <- function(df, caption="Tabela 2 -- Sensibilidade: janela curta vs amostra toda (cenário: alta volatilidade)", label="tab:mc_last5_vs_all") {
  df <- df[, c("scenario","rmse_dcf","rmse_rs","rmse_red","mape_dcf","mape_rs","mape_red","p10_dcf","p10_rs")]
  out <- c()
  out <- c(out, "\\begin{table}[htbp]")
  out <- c(out, "\\centering")
  out <- c(out, "\\small")
  out <- c(out, "\\begin{threeparttable}")
  out <- c(out, sprintf("\\caption{%s}", caption))
  out <- c(out, sprintf("\\label{%s}", label))
  out <- c(out, "\\begin{tabular}{lrrrrrrr}")
  out <- c(out, "\\toprule")
  out <- c(out, "Especifica\\c{c}\\~ao & RMSE(DCF) & RMSE(RS) & $\\Delta$RMSE(\\%) & MAPE(DCF) & MAPE(RS) & $\\Delta$MAPE(\\%) & P(APE$>$10\\%) \\\\")
  out <- c(out, "\\midrule")
  
  for (i in 1:nrow(df)) {
    p10_pair <- paste0(fmt2(100*df$p10_dcf[i]), "/", fmt2(100*df$p10_rs[i]))
    line <- paste0(
      df$scenario[i], " & ",
      fmt(df$rmse_dcf[i],2), " & ",
      fmt(df$rmse_rs[i],2), " & ",
      fmt2(df$rmse_red[i]), " & ",
      fmt(df$mape_dcf[i],4), " & ",
      fmt(df$mape_rs[i],4), " & ",
      fmt2(df$mape_red[i]), " & ",
      p10_pair, " \\\\"
    )
    out <- c(out, line)
  }
  
  out <- c(out, "\\bottomrule")
  out <- c(out, "\\end{tabular}")
  out <- c(out, "\\begin{tablenotes}")
  out <- c(out, "\\footnotesize")
  out <- c(out, "\\raggedright")
  out <- c(out, "\\item \\textit{Nota:} Compara estima\\c{c}\\~ao do crescimento por (i) \\`ultimos 5 per\\'iodos e (ii) toda a amostra observada. P(APE$>$10\\%) reporta DCF/RS (em \\%).")
  out <- c(out, "\\end{tablenotes}")
  out <- c(out, "\\end{threeparttable}")
  out <- c(out, "\\end{table}")
  paste(out, collapse="\n")
}

table1_tex <- make_table1_latex(table1_df)
table2_tex <- make_table2_latex(table2_df)

writeLines(table1_tex, con = "Table1_Scenarios.tex")
writeLines(table2_tex, con = "Table2_Last5_vs_All.tex")

cat("Saved: Table1_Scenarios.tex, Table2_Last5_vs_All.tex\n")

## ---------- Print LaTeX to console (optional) ----------
cat("\n\n===== TABLE 1 (LaTeX) =====\n\n")
cat(table1_tex)
cat("\n\n===== TABLE 2 (LaTeX) =====\n\n")
cat(table2_tex)




## ---------- Figura 3: CDF do APE com 3 painéis (1 PDF) ----------
pdf("figure3_cdf_3panels.pdf", width = 10, height = 4)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))  # 1 linha, 3 colunas

plot_cdf_panel <- function(res, main_suffix) {
  ape_dcf <- res$ape_dcf_vec
  ape_rs  <- res$ape_rs_vec
  
  F_dcf <- ecdf(ape_dcf)
  F_rs  <- ecdf(ape_rs)
  
  xmax <- as.numeric(quantile(c(ape_dcf, ape_rs), 0.99, na.rm = TRUE))
  
  plot(F_dcf, verticals = TRUE, do.points = FALSE,
       xlim = c(0, xmax),
       main = main_suffix,
       xlab = "APE",
       ylab = "CDF")
  
  lines(F_rs, verticals = TRUE, do.points = FALSE, lty = 2)
  legend("bottomright", legend = c("DCF", "RS"), lty = c(1, 2), bty = "n", cex = 0.9)
}

plot_cdf_panel(res1, "Alta persistência")
plot_cdf_panel(res2, "Alta vol.")
plot_cdf_panel(res3, "Alta instab.")

par(mfrow = c(1, 1))
dev.off()

cat("Figura 3 (3 painéis) salva em: figure3_cdf_3panels.pdf\n")
