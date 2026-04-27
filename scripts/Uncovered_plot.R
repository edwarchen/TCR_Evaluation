library(ggplot2)
library(tidyr)
plot_sequence_grid <- function(csv_file) {
  cat("----- Reading:", csv_file, " -----\n")
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  df$Sequence <- toupper(df$Sequence)
  id_col <- ifelse("ID" %in% colnames(df), "ID", "Identifier")
  df$Full_name <- df[[id_col]]
  max_len <- max(nchar(df$Sequence))
  df$Sequence <- sprintf(paste0("%-", max_len, "s"), df$Sequence)
  seq_list <- strsplit(df$Sequence, split="")
  names(seq_list) <- df$Full_name
  
  long_df <- data.frame(
    Template = rep(names(seq_list), sapply(seq_list, length)),
    Position = unlist(lapply(seq_list, function(x) 1:length(x))),
    Base = unlist(seq_list),
    stringsAsFactors = FALSE
  )
  long_df$Base[long_df$Base == " "] <- "-"
  long_df$Template <- factor(long_df$Template, levels = rev(unique(long_df$Template)))
  
  base_colors <- c(
    "A" = "#477C54", "T" = "#FF7078", 
    "C" = "#5465AB", "G" = "#F6C971", 
    "N" = "#999999", "-" = "#F5F5F5"
  )
  
  p <- ggplot(long_df, aes(x = Position, y = Template, fill = Base)) +
    geom_tile(color = "white", size = 0.3) +
    # 关键修改：新增 geom_text 图层，将 Base 映射为文本标签
    geom_text(aes(label = Base), color = "black", size = 3.5, family = "mono", fontface = "bold") +
    scale_fill_manual(values = base_colors) +
    scale_x_continuous(position = "top", expand = c(0,0), breaks = seq(1, max_len, by=5)) +
    labs(
      title = paste("Uncovered Sequences:", basename(csv_file)),
      x = "Position (5' -> 3')", y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(family = "mono", face = "bold", size = 10),
      axis.text.x.top = element_text(size = 9),
      panel.grid = element_blank(),
      legend.position = "none", # 既然格子上已经有字母了，我们可以把底部的图例关掉！
      plot.title = element_text(face = "bold", margin = margin(b = 15))
    )
  
  # 7. 动态保存尺寸
  out_file <- sub("\\.csv$", "_heatmap.png", csv_file)
  ggsave(out_file, plot = p, 
         width = max(8, max_len * 0.15), 
         height = max(4, nrow(df) * 0.3), 
         dpi = 150, limitsize = FALSE)
  
  cat("✅ 已保存图像至:", out_file, "\n")
}

uncovered_dir <- "results/uncovered_templates/IGKV_unc"
if (!dir.exists(uncovered_dir)) {
  stop("Check folder path: ", uncovered_dir)
}

csv_files <- list.files(uncovered_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No csv files found under: ", uncovered_dir)
}

for (csv_file in csv_files) {
  plot_sequence_grid(csv_file)
}


#library(colorspace)
#colorspace::sequential_hcl(n=7, palette = 'Heat') #%>%
#  colorspace::swatchplot()

