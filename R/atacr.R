# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# My imports
#' @importFrom magrittr "%>%"


add_mean_read_count <- function(df){
  df <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(mean_read_count=mean(read_count) )
  return(df)
}

rename_control_columns <- function(df,control_name){
  df <- df %>%
    dplyr::filter( sample == control_name ) %>%
    plyr::rename(
      c("sample"="control",
        "mean_coverage"="control_bait_mean_coverage",
        "breadth"="control_bait_breadth",
        "read_count"="control_bait_read_count",
        "mean_read_count"="control_mean_read_count"
      ))
  return(df)
}

rename_sample_columns <- function(df, control_name){
  df <- df %>%
    #dplyr::filter(sample != control_name) %>%
    plyr::rename(
      c("mean_coverage"="sample_bait_mean_coverage",
        "breadth"="sample_bait_breadth",
        "read_count"="sample_bait_read_count",
        "mean_read_count"="sample_mean_read_count"
      ))
}

#' @export
normalize <- function(df,control_name){
  ## expects a dataframe with columns:
  ## sample, Bait, read_count, mean_coverage, BreadthCoverage
  ## Uses the control_name as a base set of values to normalize against
  df <- add_mean_read_count(df)
  control_df <- rename_control_columns(df, control_name)
  sample_df <- rename_sample_columns(df, control_name)

  all_data <- dplyr::left_join(sample_df, control_df, by="Bait") %>%
    dplyr::arrange(sample) %>%
    dplyr::mutate(
      read_count_ratio = sample_bait_read_count/control_bait_read_count,
      read_count_ratio = replace(read_count_ratio, read_count_ratio==Inf, NaN)
    )

  all_data <- droplevels(all_data)

  all_data <- all_data %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate( sample_mean_ratio = mean(read_count_ratio, na.rm = TRUE) )

  normalized_data <- all_data %>%
    dplyr::mutate(
      normalized_value = read_count_ratio / sample_mean_ratio,
      log_normalized_value = log2(normalized_value),
      a = log2(sample_bait_read_count * control_bait_read_count)/2
    )
  return(normalized_data )
}

#' @export
control_ratio_intensity_plot <- function(df){

  p <- ggplot2::ggplot(df) + ggplot2::aes(a, log_normalized_value) + ggplot2::geom_jitter() + ggplot2::facet_wrap( ~ sample)
  return(p)
}

#' @export
count_negatives_found <- function(df){
 x <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise( negative_control_baits = sum(negative_control) )
 return(x)
}

#' @export
find_negative_control_baits <- function(df, breadth=80, mean_depth=1.2, max_sample_percent=33, max_sample_count=5, use_max_sample_count=FALSE){
  ## find all baits with breadth < breadth, mean depth < mean_depth
  ## and those that have this property in max_sample_count or max_sample_percent of samples
  ## or fewer

  d <- df %>%
    dplyr::mutate(
      b = ifelse( sample_bait_breadth < breadth & sample_bait_mean_coverage < mean_depth,
        TRUE, FALSE)
    ) %>%
    dplyr::group_by(Bait) %>%
    dplyr::mutate( count = sum(b) )

  maximum = (max_sample_percent / 100) * dplyr::n_distinct(d$sample)
  if (use_max_sample_count){
    maximum = maximum_sample_count
  }
  d <- d %>%
    dplyr::mutate(
      negative_control = ifelse(
        b == TRUE & count <= maximum, TRUE, FALSE
      )
    ) %>%
    dplyr::mutate(count=NULL, b=NULL)

  return(d)
}



#' @export
all_versus_all_normalized_plot <- function(df){
  col_count <- dplyr::n_distinct(df$sample) + 1
  p <- df %>%
    dplyr::select(Bait, sample, log_normalized_value) %>%
    dplyr::mutate(log_normalized_value = ifelse(is.na(log_normalized_value),0, log_normalized_value)) %>%
    dplyr::mutate(log_normalized_value = ifelse(is.infinite(log_normalized_value),0, log_normalized_value)) %>%
    tidyr::spread(sample, log_normalized_value)  %>%
    GGally::ggpairs(columns=2:col_count, lower=list(continuous="points"), upper=list(continuous=GGally::wrap("cor", size = 10)))
  return(p)
}

#' @export
all_versus_all_read_count_plot <- function(df){
  col_count <- dplyr::n_distinct(df$sample) + 1
  p <- df %>%
    dplyr::select(Bait, sample, sample_bait_read_count) %>%
    dplyr::mutate(sample_bait_read_count = ifelse(is.na(sample_bait_read_count),0, sample_bait_read_count)) %>%
    dplyr::mutate(sample_bait_read_count = ifelse(is.infinite(sample_bait_read_count),0, sample_bait_read_count)) %>%
    tidyr::spread(sample, sample_bait_read_count)  %>%
    GGally::ggpairs(columns=2:col_count, lower=list(continuous="points"), upper=list(continuous=GGally::wrap("cor", size = 10)))
  return(p)
}


#' @export
negative_mean_plots <- function(df){
  p <- ggplot2::ggplot(df) + ggplot2::aes(negative_control, sample_bait_read_count)  + ggplot2::geom_violin(ggplot2::aes(colour=negative_control)) + ggplot2::geom_jitter(ggplot2::aes(colour=negative_control)) + ggplot2::facet_wrap( ~ sample, scales = "free_y")
  return(p)
}

#' @export
likelihood <- function(df){
  bait_count <- dplyr::n_distinct(df$Bait)
  a <- df %>%
    dplyr::group_by(sample) %>%
      dplyr::mutate(
          sample_negative_control_mean = mean(sample_bait_read_count[negative_control == TRUE], na.rm = TRUE)
    )

  a <- a %>% dplyr::mutate(
  p = stats::ppois(sample_bait_read_count, lambda = sample_negative_control_mean, lower = FALSE),
  corrected_p = ifelse(p * bait_count >= 1, 1, p * bait_count)
  )

 # b <- dplyr::left_join(df, a, by = "Bait")
  return(a)
}
