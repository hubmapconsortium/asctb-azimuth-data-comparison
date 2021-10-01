library(jsonlite)
library(tidyr)
library(dplyr)

extract_ct_from_json <- function(data_src) {
  column_prefix = 'AS'

  asct_unnest <- function(df) {
    udf <- unnest(df, cols = !starts_with('set'), names_repair = 'universal')
    if (exists('children', udf)) {
      return(asct_unnest(udf))
    } else {
      return(
        udf %>%
          unnest(everything(), names_repair = 'universal') %>%
          count(across(!starts_with('set')), name = 'COUNT')
        %>% select(-1)
        %>% rename_with(function(names) {
          level = 1
          n = length(names)
          for (i in 1:n) {
            if (startsWith(names[i], 'name')) {
              names[i] = paste(column_prefix, level, sep = '/')
              level = level + 1
            } else if (names[i] == 'COUNT') {
              names[i] = paste(column_prefix, level - 1, 'COUNT', sep = '/')
            }
          }
          return(names)
        }))
    }
  }

  json <- jsonlite::fromJSON(data_src)
  df <- json$tree[2,]

  asct_df <- suppressMessages(asct_unnest(df))
  return(asct_df)
}
