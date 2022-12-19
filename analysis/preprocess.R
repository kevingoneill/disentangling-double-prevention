library(tidyverse)

list.files('data/prolific/', full.names=TRUE) %>%
    read_csv() %>%
    select(prolific_id, id, vignette, factor, response, rt, measure) %>%
    filter(!is.na(measure)) %>%
    pivot_wider(names_from=measure, values_from=c(response, rt)) %>%
    rename(cause=response_vignette, rt_cause=rt_vignette,
           age=response_age,
           gender=response_gender,
           attention_check=response_attention_check,
           comments=response_comments) %>%
    write_csv('data/exp1_processed.csv')
