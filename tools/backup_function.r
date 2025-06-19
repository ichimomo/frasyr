#' @export
#'
#'

make_hierarchical_csv <- function(data, hierarchy_cols, id_cols) {

  # 1. wide形式に変換
  wide_data <- data %>%
    arrange(factor(stat, levels=c("abc","btarget.prob", "blimit.risk", "overfishing.risk"))) %>%
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = all_of(hierarchy_cols),
      values_from = setdiff(names(data), c(hierarchy_cols, id_cols))
    ) %>%
    arrange(desc(is.na(upper_CV)), desc(beta))

  # 2. wide列名をパース
  colnames_all <- names(wide_data)
  split_names <- str_split_fixed(colnames_all[-(1:length(id_cols))], "_", length(hierarchy_cols))

  # 3. ラベル行を作成する
  label_rows <- list()

  # 1個目の追加ラベル列（ID列の直後に空白列を挿入）
  first_label_row <- c(id_cols, hierarchy_cols[1], rep(hierarchy_cols[1], nrow(split_names)))
  label_rows[[1]] <- first_label_row

  # 2個目以降のラベル行
  for (i in seq_along(hierarchy_cols)) {
    row_i <- c(id_cols, hierarchy_cols[i], split_names[, i])
    label_rows[[i + 1]] <- row_i
  }

  # 4. データ本体にも空白列を挿入
  wide_data_with_blank <- wide_data %>%
    mutate(dummy_label_col = "") %>%
    relocate(dummy_label_col, .after = all_of(id_cols[length(id_cols)]))

  # 5. 全部文字列化
  wide_data_chr <- wide_data_with_blank %>% mutate(across(everything(), as.character))

  # 6. バインドする
  final_output <- bind_rows(
    map(label_rows, ~as_tibble_row(setNames(.x, names(wide_data_with_blank)))),
    wide_data_chr
  )

  return(final_output)
}


save_excel_with_merged_headers <- function(df, file) {
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")

  # データを書き込む（文字列で）
  writeData(wb, "Summary", df, colNames = FALSE, rowNames = FALSE)

  # セル結合：ヘッダー（最初の複数行）について処理
  header_rows <- which(df[[1]] %in% c("year", "beta", "MSE", "シナリオ", "パターン"))
  n_header <- length(header_rows)

  # 列ごとに、同じ値が続いているところを結合する
  for (row in 1:n_header) {
    start_col <- 1
    while (start_col <= ncol(df)) {
      this_value <- df[[start_col]][row]
      end_col <- start_col
      while (end_col + 1 <= ncol(df) && df[[end_col + 1]][row] == this_value) {
        end_col <- end_col + 1
      }
      if (end_col > start_col) {
        mergeCells(wb, "Summary", cols = start_col:end_col, rows = row)
      }
      start_col <- end_col + 1
    }
  }

  # 保存
  saveWorkbook(wb, file, overwrite = TRUE)
}
