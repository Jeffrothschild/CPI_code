library(tidyverse)
library(patchwork)
library(timetk)
library(ggpmisc)
library(ggsignif)
library(broom)
library(tidymodels)
library(modeltime)
library(emmeans)
library(ggeffects)
library(performance)



group_data <- read_rds("cpr_data_sample_tbl.rds")
group_data %>% glimpse()


# *compliance -----------------------------------------------------------
arima_fn <- function(x){
    model_fit_auto_arima <- arima_reg() %>% 
        set_engine("auto_arima") %>% 
        fit(diet_kcal ~ date + study_day + exercise_load + weight_kg,
            data = x)
    
    model_fit_auto_arima$fit$models$model_1$coef %>% tidy() %>% pivot_wider(names_from = "names", values_from = "x") %>% 
        pull(study_day)
}


station_tbl <- group_data  %>% 
    filter(subject_weekly_training_h > 6) %>% 
    group_by(subject_id) %>% 
    nest() %>% 
    ungroup() %>%
    mutate(
        kpss = map(data, ~ tseries::kpss.test(.x$diet_kcal) ),    #stationarity  p > .05 means stationary
        kpss = map(kpss, tidy),
        arima_estimate = map_dbl(data, ~arima_fn(.x)),    
        diet.ols = map(data, ~ lm(diet_kcal ~ study_day + exercise_load + weight_kg, data = .x)),  #
        diet.breusch = map(diet.ols, ~ lmtest::bptest(.x)),  #heteroskedasticity
        diet.breusch = map(diet.breusch, broom::tidy),
        diet.ljung = map(diet.ols, ~ feasts::ljung_box(.x$residuals, lag = 7, dof = 3) ),   #Autocorrelation
        diet.ljung = map(diet.ljung, broom::tidy),
    ) %>% 
    unnest(diet.breusch) %>% select(-c(statistic, parameter, method)) %>% rename("breusch_pval" = p.value) %>% 
    unnest(diet.ljung) %>% filter(names == "lb_pvalue") %>% select(-names) %>% rename("lb_pvalue" = x) %>% 
    unnest(kpss) %>% select(-c(statistic, parameter, method)) %>% rename("kpss_p" = p.value) %>% 
    mutate(
        kpss_stationary = ifelse(kpss_p > 0.05, 1, 0),
    ) 


# to find subjects with significant trend in dietary reporting
set.seed(5012022)
ci_check_tbl <- station_tbl%>% 
    filter(kpss_stationary <1) %>%
    mutate(
        diet_resamped = map(data, ~ reg_intervals(diet_kcal ~  study_day + exercise_load +  weight_kg, data = .x)),  # 
    ) %>% 
    unnest(diet_resamped) %>% filter(term == "study_day") %>% select(-c(.alpha, .method, term)) %>% rename_with(~ str_c("boots", .), .lower:.upper) %>% 
    mutate(
        boot_ci_sig = ifelse(boots.lower * boots.upper > 1, "sig", "ns")
    ) 


subjects_with_trend <- ci_check_tbl %>% 
    select(1:2, arima_estimate, contains("boot"), breusch_pval, lb_pvalue) %>% 
    filter(
        boot_ci_sig == "sig",
        # breusch_pval > 0.05,
        # lb_pvalue > 0.05,
        boots.estimate < -5
    ) %>% 
    pull(subject_id)



low_volume_subjects <- group_data %>% filter(subject_weekly_training_h < 6) %>% 
    group_by(subject_id) %>% 
    summarise(mean_training = mean(subject_weekly_training_h, na.rm = T)) %>% 
    pull(subject_id)



group_data_tbl2 <- group_data %>% 
    filter(
        !subject_id %in% subjects_with_trend,
        !subject_id %in% low_volume_subjects,
        !is.na(diet_kcal)
    ) 

median_carbs <- group_data_tbl2 %>% 
    group_by(subject_id) %>%
    summarise(median_CHO = median(diet_carb_g_kg, na.rm = T)) 

high_carb_quantile <- quantile(median_carbs$median_CHO, 0.8, names = FALSE)
low_carb_quantile <- quantile(median_carbs$median_CHO, 0.2, names = FALSE)

group_data_tbl <- group_data_tbl2 %>% 
    mutate(subject_id = factor(subject_id)) %>% 
    group_by(subject_id) %>% 
    mutate(
        habitual_diet = case_when(
            median(diet_carb_g_kg, na.rm = T) < low_carb_quantile ~ "Low-CHO",
            median(diet_carb_g_kg, na.rm = T) > high_carb_quantile ~ "High-CHO",
            TRUE ~ "Mod-CHO"),
        habitual_diet = factor(habitual_diet, levels = c("Low-CHO", "Mod-CHO", "High-CHO"))
    ) %>% 
    ungroup()



# list tbl ----------------------------------------------------------------
pci_list_tbl <- group_data_tbl %>%
    mutate(subject_level = ifelse(subject_level == "Professional", "Elite non-professional", subject_level),
           subject_level = factor(subject_level, levels = c("Amateur", "High-level amateur", "Elite non-professional")),
           subject_id = factor(subject_id)
    ) %>% 
    group_by(subject_id) %>% 
    nest() %>% 
    mutate(
        cor = map(data, ~ cor.test(.x$diet_carb_g_kg, .x$exercise_load, , method="pearson", use = "pairwise.complete.obs")),
        cor = map(cor, broom::tidy),
    ) %>% 
    unnest(cor) %>% select(-c(method, alternative, statistic, parameter)) %>%
    rename(cor_val = "estimate",
           cor_p = "p.value"
    ) %>% 
    mutate(
        mean_carb = map_dbl(data, ~mean(.x$diet_carb_g_kg, na.rm = T)),
        sd_carb = map_dbl(data, ~sd(.x$diet_carb_g_kg, na.rm = T)),
        monotony = map_dbl(data, ~ mean_carb/sd_carb),
        carb_range = map_dbl(data, ~ max(.x$diet_carb_g_kg, na.rm = T) - min(.x$diet_carb_g_kg, na.rm = T)),
        cpi_index = carb_range * cor_val/monotony,
    )


subject_diet_level_tbl <- pci_list_tbl %>% 
    unnest(data) %>% 
    select(subject_level, habitual_diet, subject_id, subject_sex, subject_weekly_training_h, subject_primary_sport) %>% 
    summarise(habitual_diet = first(habitual_diet),
              subject_level = first(subject_level),
              sex = factor(first(subject_sex)),
              training_volume = first(subject_weekly_training_h),
              subject_primary_sport = first(subject_primary_sport)
    )



fsted_pct_tbl <- group_data_tbl %>% 
    select(subject_id, exercise_duration_min,exercise_fasted ) %>% 
    filter(exercise_duration_min >0) %>% 
    mutate(exercise_fasted = as.numeric(exercise_fasted)-1) %>% 
    group_by(subject_id) %>% 
    summarise(fasted_training = sum(exercise_fasted, na.rm = T),
              nobs = length(exercise_fasted),
              fasted_training_pct = fasted_training/nobs*100) %>% 
    select(subject_id, fasted_training_pct)



df_for_comparison <- pci_list_tbl %>% 
    left_join(subject_diet_level_tbl, by = "subject_id") %>% 
    left_join(fsted_pct_tbl, by = "subject_id") %>% 
    ungroup()


# colors ------------------------------------------------------------------
theme_JR1 <- function(axis_lines = TRUE, 
                      grid_lines = FALSE,     
                      text_size = 12,       
                      line_size = 0.2,
                      base_family= 'sans'){ 
    
    th <- ggplot2::theme_minimal(base_family = base_family, 
                                 base_size = text_size)
    th <- th + theme(panel.grid=element_blank())
    if (axis_lines) {
        th <- th + 
            theme(axis.line = element_line(size = line_size, color = "black"),
                  axis.ticks = element_line(size = line_size, color = "black"),
                  axis.text=element_text(color = "black"),
                  axis.title = element_text(size=13, colour="black", face="bold"))
    } 
    if (grid_lines) {
        th <- th + 
            theme(panel.grid.major = element_line(size = line_size))
    }
    th <- th + theme(
        axis.text.x=element_text(margin=margin(t=2)),
        axis.text.y=element_text(margin=margin(r=2)),
        axis.title.x=element_text(margin=margin(t=5)),
        axis.title.y=element_text(margin=margin(r=5)),
        plot.title=element_text(margin=margin(b=5)))
    
    return (th)
}

JR_facet_theme <- function(){
    theme_light()+
        theme(
            strip.background = element_rect(fill = "grey40", color = "grey80", size = 1),
            strip.text = element_text(colour = "white"),
            axis.title.y = element_text(face = "bold")
        )
}


ggplot2::theme_set(theme_JR1()) 
colors2 <- c("midnightblue", "#6CC458")   
colors2green <- c("#6CC458", "#6CC458")   
colors3 <- c("#6CC458", "#6CC458", "midnightblue")
colors33 <- c("#B33951", "midnightblue","#6CC458")
colors331 <- c("#B33951")
colors332 <- c("midnightblue")
colors333 <- c("#6CC458")
colors4 <- c("midnightblue", "#6CC458", "midnightblue", "#6CC458")
colors5 <- c("midnightblue", "#6CC458",  "#6CC458", "midnightblue", "grey60")
colors7 = c( "#6CC458", "grey70", "red")


# Sub-group analysis ----------------------------------------------------------

# * CPI Emmeans -----------------------------------------------------------------
df_for_mod <- df_for_comparison %>% select(cpi_index, subject_level, habitual_diet, sex, subject_id, fasted_training_pct, training_volume)

cpi_mod <- lm(cpi_index ~ subject_level + habitual_diet + sex , df_for_mod) 

check_outliers(cpi_mod)
check_model(cpi_mod)
shapiro.test(cpi_mod$residuals)
lmtest::bptest(cpi_mod)

augment(cpi_mod) %>% ggplot(aes(.fitted, .resid))+  geom_point()

#NA's due to not ebough values in code sample subset - see sex comparison for working code
emmeans(cpi_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()  
es_level_cpi <- eff_size(emmeans(cpi_mod, ~subject_level), sigma = sigma(cpi_mod), edf = cpi_mod$df.residual) %>% 
    as_tibble() %>% mutate(comp = "level") %>% 
    bind_cols(
        emmeans(cpi_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(cpi_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   
es_diet_cpi <- eff_size(emmeans(cpi_mod, ~habitual_diet), sigma = sigma(cpi_mod), edf = cpi_mod$df.residual) %>% 
    as_tibble() %>% mutate(comp = "diet") %>% 
    bind_cols(
        emmeans(cpi_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(cpi_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  
es_sex_cpi <- eff_size(emmeans(cpi_mod, ~sex), sigma = sigma(cpi_mod), edf = cpi_mod$df.residual) %>% 
    as_tibble() %>% 
    mutate(comp = "sex",
           contrast =  ifelse(contrast == "M - F", "Male - Female", contrast)) %>% 
    bind_cols(
        emmeans(cpi_mod, ~sex) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 


ES_tbl_cpi <- es_level_cpi %>% 
    bind_rows(es_diet_cpi) %>%
    bind_rows(es_sex_cpi) %>% 
    mutate(
        comp = as.factor(comp),
        contrast = str_remove_all(contrast, "\\(|\\)"),
        contrast = str_replace_all(contrast, "professional", "pro"),
        contrast = str_replace_all(contrast, "High-level amateur", "HLA"),
        contrast = str_replace_all(contrast, " - ", ":\n"),
        contrast = as.factor(contrast),
        Sign = case_when(
            p.value > 0.05 ~ "NS",
            p.value < 0.05 & effect.size > 0 ~ "+",
            p.value < 0.05 & effect.size < 0 ~ "-"
        )
    ) 

ES_plot_cpi  <- ES_tbl_cpi %>% 
    ggplot(aes(effect.size, fct_rev(fct_inorder(contrast)), xmin = lower.CL , xmax = upper.CL, color = Sign)) +
    geom_vline(xintercept = 0, lty = 3, color = "grey") +
    geom_linerange()+
    geom_point()+
    scale_color_manual(values = c("-" = "black" , 
                                  "NS" = "grey82",
                                  "+"=  "black")) + 
    labs(y = NULL, x = "Effect size")+
    theme(
        legend.position = "none")


# * training Emmeans -----------------------------------------------------------------

training_mod <- lm(training_volume ~ subject_level + habitual_diet + sex , df_for_mod)  

check_outliers(training_mod)
check_model(training_mod)
shapiro.test(training_mod$residuals)

augment(training_mod) %>% ggplot(aes(.fitted, .resid))+  geom_point()

#NA's due to not ebough values in code sample subset - see sex comparison for working code

emmeans(training_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()  
es_level_training <- eff_size(emmeans(training_mod, ~subject_level), sigma = sigma(training_mod), edf = training_mod$df.residual) %>% 
    as_tibble() %>% 
    mutate(comp = "level") %>% 
    bind_cols(
        emmeans(training_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(training_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   
es_diet_training <- eff_size(emmeans(training_mod, ~habitual_diet), sigma = sigma(training_mod), edf = training_mod$df.residual) %>% 
    as_tibble() %>% 
    mutate(comp = "diet") %>% 
    bind_cols(
        emmeans(training_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(training_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  
es_sex_training <- eff_size(emmeans(training_mod, ~sex), sigma = sigma(training_mod), edf = training_mod$df.residual) %>% 
    as_tibble() %>% 
    mutate(comp = "sex",
           contrast =  ifelse(contrast == "M - F", "Male - Female", contrast)) %>% 
    bind_cols(
        emmeans(training_mod, ~sex) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    )


ES_tbl_training <- es_level_training %>% 
    bind_rows(es_diet_training) %>%
    bind_rows(es_sex_training) %>% 
    mutate(
        comp = as.factor(comp),
        contrast = str_remove_all(contrast, "\\(|\\)"),
        contrast = str_replace_all(contrast, "professional", "pro"),
        contrast = str_replace_all(contrast, "High-level amateur", "HLA"),
        contrast = str_replace_all(contrast, " - ", ":\n"),
        contrast = as.factor(contrast),
        Sign = case_when(
            p.value > 0.05 ~ "NS",
            p.value < 0.05 & effect.size > 0 ~ "+",
            p.value < 0.05 & effect.size < 0 ~ "-"
        )
                 
    ) 

ES_plot_training  <- ES_tbl_training %>% 
    ggplot(aes(effect.size, fct_rev(fct_inorder(contrast)), xmin = lower.CL , xmax = upper.CL, color = Sign)) +
    geom_vline(xintercept = 0, lty = 3, color = "grey") +
    geom_linerange()+
    geom_point()+
    scale_color_manual(values = c("-" = "black" , 
                                  "NS" = "grey82",
                                  "+"=  "black")) + 
    labs(y = NULL, x = "Effect size")+
    theme(
        legend.position = "none")


# * fasting Emmeans -----------------------------------------------------------------

fasting_mod <- lm(fasted_training_pct ~ subject_level + habitual_diet + sex , df_for_mod)  

check_outliers(fasting_mod)
check_model(fasting_mod)
shapiro.test(fasting_mod$residuals)

augment(fasting_mod) %>% ggplot(aes(.fitted, .resid))+  geom_point()

emmeans(fasting_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()  
es_level_fasting <- eff_size(emmeans(fasting_mod, ~subject_level), sigma = sigma(fasting_mod), edf = fasting_mod$df.residual) %>% 
    as_tibble() %>% mutate(comp = "level") %>% 
    bind_cols(
        emmeans(fasting_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(fasting_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   
es_diet_fasting <- eff_size(emmeans(fasting_mod, ~habitual_diet), sigma = sigma(fasting_mod), edf = fasting_mod$df.residual) %>% 
    as_tibble() %>% mutate(comp = "diet") %>% 
    bind_cols(
        emmeans(fasting_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 

emmeans(fasting_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  
es_sex_fasting <- eff_size(emmeans(fasting_mod, ~sex), sigma = sigma(fasting_mod), edf = fasting_mod$df.residual) %>% 
    as_tibble() %>% 
    mutate(comp = "sex",
           contrast =  ifelse(contrast == "M - F", "Male - Female", contrast))  %>% 
    bind_cols(
        emmeans(fasting_mod, ~sex) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble() %>% select(p.value)
    ) 


ES_tbl_fasting <- es_level_fasting %>% 
    bind_rows(es_diet_fasting) %>%
    bind_rows(es_sex_fasting) %>% 
    mutate(
        comp = as.factor(comp),
        contrast = str_remove_all(contrast, "\\(|\\)"),
        contrast = str_replace_all(contrast, "professional", "pro"),
        contrast = str_replace_all(contrast, "High-level amateur", "HLA"),
        contrast = str_replace_all(contrast, " - ", ":\n"),
        contrast = as.factor(contrast),
        Sign = case_when(
            p.value > 0.05 ~ "NS",
            p.value < 0.05 & effect.size > 0 ~ "+",
            p.value < 0.05 & effect.size < 0 ~ "-")

    ) 

ES_plot_fasting  <- ES_tbl_fasting %>% 
    ggplot(aes(effect.size, fct_rev(fct_inorder(contrast)), xmin = lower.CL , xmax = upper.CL, color = Sign)) +
    geom_vline(xintercept = 0, lty = 3, color = "grey") +
    geom_linerange()+
    geom_point()+
    scale_color_manual(values = c("-" = "black" , # "#B33951",
                                  "NS" = "grey82",
                                  "+"=  "black")) + # "#6CC458")) +
    labs(y = NULL, x = "Effect size")+
    theme(
        legend.position = "none")



# FIGS --------------------------------------------------------------------


# *habitual diet fig -------------------------------------------------------
group_data_tbl %>% 
    ggplot(aes(fct_reorder(subject_id, diet_carb_g_kg), diet_carb_g_kg, color = habitual_diet))+
    geom_hline(yintercept = low_carb_quantile, lty = 3, color = "grey40")+    #maybe make darker?
    geom_hline(yintercept = high_carb_quantile, lty = 3, color = "grey40")+
    geom_boxplot()+
    scale_color_manual(values = colors33)+
    scale_y_continuous(breaks = c(0,3,6,9,12,15))+
    labs(x = "Individual participants", y = "Dietary carbohydrate (g/kg)", color = NULL)+
    guides(color = guide_legend(nrow = 1))+
    theme(legend.position = c(0.25, 0.95),
          legend.text = element_text(size = 12),
          axis.text.x = element_blank())


# *monotony fig -----------------------------------------------------------
max_monotony_subject <- pci_list_tbl  %>% select(subject_id, data, monotony)%>% arrange(desc(monotony)) %>% 
    mutate(subject_id = as.character(subject_id),
           subject_id = as.numeric(subject_id)) %>%
    ungroup() %>%
    slice(2) %>%
    pull("subject_id")

min_monotony_subject <- pci_list_tbl  %>% select(subject_id, data, monotony) %>% arrange(monotony) %>% 
    mutate(subject_id = as.character(subject_id),
           subject_id = as.numeric(subject_id)) %>% 
    ungroup() %>% 
    slice(1) %>% 
    pull("subject_id")


monotony_tbl<- pci_list_tbl %>% select(subject_id, monotony) %>% 
    ungroup() %>% 
    mutate(
        code = case_when(
            subject_id == max_monotony_subject ~ "High",
            subject_id == min_monotony_subject ~ "Low",
            TRUE ~ "other"),
        code = factor(code)
    ) 

max_monotony <- format(round(max(monotony_tbl$monotony, na.rm = T),1), nsmall = 1)
min_monotony <- format(round(min(monotony_tbl$monotony, na.rm = T),1), nsmall = 1)


monotony_bars <- monotony_tbl %>% 
    ggplot(aes(fct_reorder(subject_id, monotony), monotony))+
    geom_col(aes(fill = code), show.legend = F)+
    scale_fill_manual(values = colors3)+
    scale_y_continuous(expand = c(0,0))+
    labs(x = NULL, y = "Carbohydrate monotony\n(mean/SD))", title = "Carbohydrate Monotony")+
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 17))

monot_examples <- pci_list_tbl %>% 
    filter(subject_id == min_monotony_subject | subject_id == max_monotony_subject) %>% 
    mutate(
        subject_id = case_when(
            subject_id == max_monotony_subject ~ "High monotony",
            subject_id == min_monotony_subject ~ "Low monotony"),
        subject_id = factor(subject_id, levels = c("Low monotony", "High monotony"))
    ) %>% 
    unnest(data) %>% 
    ggplot(aes(study_day, diet_carb_g_kg, color = subject_id))+
    geom_point(alpha = .7, size = 2.5, show.legend = F) +
    scale_color_manual(values = colors2green)+
    scale_y_continuous(limits = c(0,10))+
    labs(x = "Study day", y = "Dietary carbohydrate\n(g/kg)", color = NULL)+
    facet_wrap(~ subject_id)+
    # facet_wrap(~ subject_id, labeller = as_labeller(monotony_facet_names))+
    theme(
        strip.text = element_text(face = "bold", size = 11),
        # strip.background = element_rect(fill = "midnightblue"),
        axis.title.y = element_text(size = 13)
    )


monotony_bars + monot_examples  +
    plot_layout(ncol = 1) +
    plot_annotation(title = "", tag_levels = 'a') &
    theme(plot.tag = element_text(size = 17, face="bold"))



# *correlation fig --------------------------------------------------------

ex_high_cor_high_range <- 513
ex_high_cor_low_range <- 508
ex_low_cor_high_range <- 514
ex_low_cor_low_range <- 506

diet_train_corr_fig <- pci_list_tbl %>% 
    mutate(
        code = case_when(
            subject_id == ex_high_cor_high_range ~ "High correlation, high range",
            subject_id == ex_high_cor_low_range ~ "High correlation, low range",
            subject_id == ex_low_cor_high_range ~ "Low correlation, high range",
            subject_id == ex_low_cor_low_range ~ "Low correlation, low range",
            TRUE ~ "other" )
    ) %>% 
    ggplot(aes(cor_val, fct_reorder(subject_id, cor_val), color = code))+
    geom_vline(xintercept=0, linetype = "dotted") +
    geom_point(aes()) +
    scale_colour_manual(values = colors5) +
    geom_errorbarh(aes(xmin=conf.low, xmax=conf.high,height = 0)) +
    theme(legend.position = "none") +
    labs(x = "Correlation between daily CHO\nintake and training load", y = "Participants") +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #    axis.line.y = element_blank(),
        legend.position = "none"
    )



diet_train_corr_facet <-  pci_list_tbl %>% 
    filter(subject_id == ex_high_cor_high_range | subject_id == ex_high_cor_low_range |  
               subject_id == ex_low_cor_high_range | subject_id == ex_low_cor_low_range) %>% 
    mutate(
        subject_id = case_when(
            subject_id == ex_high_cor_high_range ~ "High correlation, high range",
            subject_id == ex_high_cor_low_range ~ "High correlation, low range",
            subject_id == ex_low_cor_high_range ~ "Low correlation, high range",
            subject_id == ex_low_cor_low_range ~ "Low correlation, low range"),
        subject_id = factor(subject_id, levels = c("High correlation, low range", "High correlation, high range", "Low correlation, low range", "Low correlation, high range"))
    ) %>% 
    unnest(data) %>% 
    ggplot(aes(diet_carb_g_kg, exercise_load, color = subject_id))+
    geom_point(alpha = .7, size = 2.5, show.legend = F) +
    geom_smooth(se = F, method = "lm", show.legend = F)+
    scale_color_manual(values = colors4)+
    labs(x = "Dietary carbohydrate\n(g/kg)", y = "Training load (AU)", color = NULL)+
    facet_wrap(~ subject_id, scales = "free_y")+
    #JR_facet_theme()+
    theme(
        strip.text = element_text(face = "bold", size = 10, color = "white"),
        strip.background = element_rect(fill = "midnightblue"),
        axis.title = element_text(size = 13, face = "bold")
    )




diet_train_corr_fig + diet_train_corr_facet +
    plot_layout(ncol = 2,
                widths = c(1, 2)) +
    plot_annotation(title = "Diet-training correlations", tag_levels = 'a') &
    theme(plot.tag = element_text(size = 17, face="bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 19))



# *CPI fig -----------------------------------------------------------
highest_cpi <- round(max(pci_list_tbl$cpi_index),1)
lowest_cpi <- round(min(pci_list_tbl$cpi_index),1)

ex_high_CPI <- pci_list_tbl %>%  select(1,2,cpi_index) %>% 
    ungroup() %>% 
    arrange(desc(cpi_index)) %>% 
    mutate(subject_id = as.character(subject_id),
           subject_id = as.numeric(subject_id)) %>% 
    slice(1) %>% 
    pull(subject_id)


ex_low_CPI <- pci_list_tbl %>%  select(1,2,cpi_index) %>% 
    ungroup() %>% 
    arrange(abs(cpi_index)) %>% 
    mutate(subject_id = as.character(subject_id),
           subject_id = as.numeric(subject_id)) %>% 
    slice(3) %>% 
    pull(subject_id)

CPI_bars <- pci_list_tbl  %>% 
    mutate(
        code_CPI = case_when(
            subject_id == ex_high_CPI ~ "High",
            subject_id == ex_low_CPI ~ "Low",
            TRUE ~ "other"),
        code_CPI = factor(code_CPI, levels = c("High", "Low", "other"))) %>% 
    ggplot(aes(fct_reorder(subject_id, cpi_index), cpi_index))+
    geom_col(aes(fill = code_CPI), show.legend = F)+
    scale_fill_manual(values = colors3)+
    scale_y_continuous(expand = c(0,0))+
    labs(x = NULL, y = "CPI (AU)", title = "Carbohydrate Periodization Index")+
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 17))

CPI_bars

ex1_CPI <- pci_list_tbl %>%   
    filter(subject_id == ex_low_CPI | subject_id == ex_high_CPI) %>% 
    mutate(
        subject_id = case_when(
            subject_id == ex_high_CPI ~ "High CPI",
            subject_id == ex_low_CPI ~ "Low CPI"),
        subject_id = factor(subject_id, levels = c("Low CPI", "High CPI")),
        across(c(monotony, cpi_index, carb_range), ~ round(., 1)),
        cor_val = round(cor_val, 2)
    )

dat_text <- data.frame(
    label = c(paste0("Correlation: ", ex1_CPI$cor_val[1],"\n", "CHO range: ", ex1_CPI$carb_range[1], "\n", "Monotony: ", ex1_CPI$monotony[1] ), 
              paste0("Correlation: ",  ex1_CPI$cor_val[2],"\n", "CHO range: ", ex1_CPI$carb_range[2], "\n","Monotony: ", ex1_CPI$monotony[2] ) 
    ),
    subject_id   = c("Low CPI", "High CPI")
)


cpi_text <- data.frame(
    label = c(
        paste0("CPI: ",  ex1_CPI$cpi_index[1]), 
        paste0("CPI: ",  ex1_CPI$cpi_index[2])
    ),
    subject_id   = c("Low CPI", "High CPI")
)
cpi_text


CPI_examples <- pci_list_tbl %>% 
    filter(subject_id == ex_low_CPI | subject_id == ex_high_CPI) %>% 
    mutate(
        subject_id = case_when(
            subject_id == ex_high_CPI ~ "High CPI",
            subject_id == ex_low_CPI ~ "Low CPI"),
        subject_id = factor(subject_id, levels = c("Low CPI", "High CPI"))
    ) %>% 
    unnest(data) %>% 
    ggplot(aes(diet_carb_g_kg, exercise_load, color = subject_id))+
    geom_point(alpha = .7, size = 2.5, show.legend = F) +
    geom_smooth(se = F, method = "lm", show.legend = F, color = "midnightblue")+
    geom_text(data = dat_text, aes(x = 1.5, y = 2400, label = label), hjust   = -0.1,color = "black")+
    geom_text(data =  cpi_text, aes(x = 1.5, y = 1800, label = label), hjust   = -0.3,color = "black" , size = 5)+
    geom_text(data =  cpi_text, aes(x = 1.5, y = 1800, label = label), hjust   = -0.3,color = "black" , size = 5.05)+
    scale_color_manual(values = colors2green) +
    scale_x_continuous(limits = c(1.5,10))+
    scale_y_continuous(limits = c(0, 2800))+
    labs(y = "Training load (AU)", x = "Dietary carbohydrate (g/kg)", color = NULL)+
    facet_wrap(~ factor(subject_id, levels = c("Low CPI", "High CPI")))+
    theme(
        strip.text = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "midnightblue"),
        axis.title.y = element_text(size = 13),
        legend.position = "none"
    )


CPI_bars + CPI_examples  +
    plot_layout(ncol = 1) +
    plot_annotation(title = "", tag_levels = 'a') &
    theme(plot.tag = element_text(size = 17, face="bold"))



# *CPI subgroup  -----------------------------------------------------------------

min_y_cpi <- min(pci_list_tbl$cpi_index) -.2
max_y_cpi <- max(pci_list_tbl$cpi_index) +1.5



cpi_by_level <- df_for_comparison %>% 
    ggplot(aes(subject_level, cpi_index))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = subject_level, shape = subject_level), size = 1.8, alpha = .8, position = position_jitter(height = 0, width = 0.15, seed = 1))+
    scale_y_continuous(limits = c(min_y_cpi, max_y_cpi))+
    scale_color_manual(values = colors33)+
    geom_signif(comparisons = list(
        c("Amateur", "Elite non-professional")
        # c("High-level amateur", "Elite non-professional")
    ),
    map_signif_level = T, y_position = max_y_cpi*.8, annotations = c("*"), step_increase = .13, tip_length = 0.02, test = TukeyHSD) +  
    labs( x =  "Competitive level", y = "CPI (AU)", color = NULL, shape = NULL)+
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )

emmeans(cpi_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   

cpi_by_diet <- df_for_comparison %>% 
    ggplot(aes(habitual_diet, cpi_index))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = subject_level, shape = subject_level), size = 1.8, alpha = .8, position = position_jitter(height = 0, width = 0.15, seed = 2))+
    scale_y_continuous(limits = c(min_y_cpi, max_y_cpi))+
    scale_color_manual(values = colors33)+
    # geom_signif(comparisons = list(
    #   c("Mod-CHO", "Low-CHO")
    #   # c("High-CHO", "Low-CHO")
    # ),
    # map_signif_level = T, y_position = max_y_cpi*.83, annotations = c("*"), step_increase = .13, tip_length = 0.02, test = TukeyHSD) +  
    labs(x =   "Habitual diet", y = "", color = NULL, shape = NULL)+
    # guides(color = guide_legend(nrow = 1))+
    theme(
        legend.position = c(.8, .9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )

emmeans(cpi_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  

cpi_by_sex  <- df_for_comparison %>% 
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other")) %>% 
    ggplot(aes(sex, cpi_index))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = subject_level, shape = subject_level), size = 1.8, alpha = .8, position = position_jitter(height = 0, width = 0.15, seed = 1))+
    scale_y_continuous(limits = c(min_y_cpi, max_y_cpi))+
    scale_color_manual(values = colors33)+
    guides(color = guide_legend(nrow = 1))+
    # geom_signif(comparisons = list(
    #   c("Female", "Male")
    # ),
    # map_signif_level = T, y_position = max_y_cpi*.8, annotations = c("*"), step_increase = .13, tip_length = 0.02, test = TukeyHSD) + 
    labs(x = " ", y = "CPI (AU)", color = NULL, shape = NULL)+
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17))


cpi_by_level + cpi_by_diet+ cpi_by_sex + # ES_plot_cpi +
    plot_layout(ncol = 2)+
    plot_annotation(title = "", tag_levels = list(c('a', 'b', 'c', ' '))) &
    theme(plot.tag = element_text(size = 17, face="bold"))



# * diet vs load --------------------------------------------------------------

# **carb vs 3-day load ------------------------------------------------------------
plot_by_diet_fn <- function(df, dv, iv){
    
    my.formula1 <- y ~ poly(x, 1, raw = TRUE)
    my.formula2 <- y ~ poly(x, 2, raw = TRUE)
    my.formula3 <- y ~ poly(x, 3, raw=TRUE)
    
    remove_na <- df %>% filter(!is.na({{iv}}) & !is.na({{dv}})) %>% select({{dv}}, {{iv}}, subject_id, habitual_diet)
    
    
    frmla1 <- formula(paste(colnames(remove_na)[1], " ~ ", colnames(remove_na)[2], " +  (1| subject_id)", collapse = " "))
    frmla2 <- formula(paste(colnames(remove_na)[1], " ~ poly(", colnames(remove_na)[2], ",2) +  (1| subject_id)", collapse = " "))
    frmla3 <- formula(paste(colnames(remove_na)[1], " ~ poly(", colnames(remove_na)[2], ",3) +  (1| subject_id)", collapse = " "))
    
    
    diet_list_tbl <- remove_na %>%
        group_by(habitual_diet) %>% 
        nest()  %>% 
        mutate(
            load_model_1 = purrr::map(data, ~  lmerTest::lmer(frmla1, data= .x, REML = F)),
            load_model_2 = purrr::map(data, ~  lmerTest::lmer(frmla2, data= .x, REML = F)),
            load_model_3 = purrr::map(data, ~  lmerTest::lmer(frmla3, data= .x, REML = F)),
            anova = pmap(list(load_model_1, load_model_2, load_model_3), ~ lme4:::anovaLmer(..1,..2, ..3))) %>% 
        unnest(anova) %>% select(-c(npar:Df)) %>% 
        rename(p = 'Pr(>Chisq)' ) %>% 
        mutate(p = ifelse(is.na(p), 0, p),
               best_mod_load= 1:3) %>%   #number the models
        filter(p < 0.08) %>% 
        select(-p) %>% 
        slice_tail(n = 1)  %>% 
        mutate(
            RMSE_1 = map_dbl(load_model_1, sjstats::rmse),
            RMSE_2 = map_dbl(load_model_2, sjstats::rmse),
            RMSE_3 = map_dbl(load_model_3, sjstats::rmse),
            R2mar_1 = map(load_model_1, MuMIn::r.squaredGLMM),
            R2mar_2 = map(load_model_2, MuMIn::r.squaredGLMM),
            R2mar_3 = map(load_model_3, MuMIn::r.squaredGLMM),
        ) %>% 
        unnest(R2mar_1:R2mar_3) 
    
    
    LC_r2 <- round(diet_list_tbl[paste0("R2mar_", diet_list_tbl$best_mod_load[1])] %>% slice(1),2)  %>% 
        as.matrix() %>% as_tibble() %>% mutate(across(1:2, ~round(., 2))) %>%  pull(1)
    
    MC_r2 <- round(diet_list_tbl[paste0("R2mar_", diet_list_tbl$best_mod_load[2])] %>% slice(2),2)%>% 
        as.matrix() %>% as_tibble() %>% mutate(
            across(1:2, ~round(., 2)),
            across(1, ~format(., nsmall = 2))
        ) %>% 
        pull(1)
    
    HC_r2 <- round(diet_list_tbl[paste0("R2mar_", diet_list_tbl$best_mod_load[3])] %>% slice(3),2)%>% 
        as.matrix() %>% as_tibble() %>% mutate(across(1:2, ~round(., 2))) %>%  pull(1)
    
    
    diet_list_tbl %>%
        unnest(data) %>%
        ggplot(aes())+
        geom_point(data = . %>% filter(habitual_diet == "Mod-CHO"), aes({{iv}}, {{dv}}), alpha=.1, color = colors332)+
        geom_smooth(data = . %>% filter(habitual_diet == "Mod-CHO"), aes({{iv}}, {{dv}}), method = "lm", color = colors332,
                    se = F, formula = get(paste0("my.formula", diet_list_tbl$best_mod_load[2])))+
        geom_point(data = . %>% filter(habitual_diet == "Low-CHO"), aes({{iv}}, {{dv}}), alpha=.15, color = colors331)+
        geom_smooth(data = . %>% filter(habitual_diet == "Low-CHO"), aes({{iv}}, {{dv}}), method = "lm", color = colors331,
                    se = F, formula = get(paste0("my.formula", diet_list_tbl$best_mod_load[1])))+
        geom_point(data = . %>% filter(habitual_diet == "High-CHO"), aes({{iv}}, {{dv}}), alpha=.15, color = colors333)+
        geom_smooth(data = . %>% filter(habitual_diet == "High-CHO"), aes({{iv}}, {{dv}}), method = "lm", color = colors333,
                    se = F, formula = get(paste0("my.formula", diet_list_tbl$best_mod_load[3])))+
        annotate("text", x = max(remove_na[2]) *.05, y = max(remove_na[1]) *.82, hjust = 0,
                 label = paste0("Low-CHO~R^2==", LC_r2), parse=TRUE, size = 4.5,  color = colors331) +
        annotate("text", x = max(remove_na[2]) *.05, y = max(remove_na[1]) *.9, hjust = 0,
                 label = paste0("Mod-CHO~R^2==", MC_r2), parse=TRUE, size = 4.5,  color = colors332) +
        annotate("text", x = max(remove_na[2])*.05, y = max(remove_na[1])*.98, hjust = 0,
                 label = paste0("High-CHO~R^2==", HC_r2) , parse=TRUE, size = 4.5,  color = colors333)
    
    
}


carb_v_load <- plot_by_diet_fn(group_data_tbl, diet_carb_g_kg, exercise_load) +
    labs(x = "Training load (AU)", y = "Daily CHO (g/kg)"
    )

carb_v_rpe <- plot_by_diet_fn(group_data_tbl, diet_carb_g_kg, exercise_RPE_weighted) +
    labs(x = "Session RPE", y = "Daily CHO (g/kg)"
    )

carb_v_duration <- plot_by_diet_fn(group_data_tbl, diet_carb_g_kg, exercise_duration_min) +
    labs(x = "Exercise duration (min)", y = "Daily CHO (g/kg)"
    )


# **carb before duration load intensity ---------------------------------------------

pre_carb_v_load <- plot_by_diet_fn(group_data_tbl, carb_before_g_kg, exercise_load)+
    labs(y = "Pre-exercise CHO (g/kg)", x = "Training load (AU)") 

pre_carb_v_duration <- plot_by_diet_fn(group_data_tbl, carb_before_g_kg, exercise_duration_min)+
    labs(y = "Pre-exercise CHO (g/kg)", x = "Exercise duration (min)") 

pre_carb_v_rpe <- plot_by_diet_fn(group_data_tbl, carb_before_g_kg, exercise_RPE_weighted )+
    labs(y = "Pre-exercise CHO (g/kg)", x = "Session RPE") 


# ***combined  ----------------------------------------------------------

carb_v_load + carb_v_duration + carb_v_rpe + 
    pre_carb_v_load + pre_carb_v_duration + pre_carb_v_rpe +
    plot_layout(nrow = 2) +
    plot_annotation(title = "", tag_levels = 'a') &
    theme(plot.tag = element_text(size = 17, face="bold"))



# supplemental figs -------------------------------------------------------

# * simulations -----------------------------------------------------------

sim_carb_range <-  seq(0.4, 12, .2)
sim_cor <-  seq(0.1, .85, .05)
sim_monotony <-  seq(1, 6.2, .2)



sim_tbl <- expand_grid(sim_carb_range, sim_cor, sim_monotony) %>%
    mutate(
        cpi_index = sim_carb_range * sim_cor / sim_monotony,
    )

sim_a <- sim_tbl %>% 
    mutate(
        scenario = case_when(
            sim_monotony == 2 &  sim_cor == .2 ~ "Low monotony Low correlation",
            sim_monotony == 6 &  sim_cor == .2 ~ "High monotony Low correlation",
            sim_monotony == 2 &  sim_cor == .8 ~ "Low monotony High correlation",
            sim_monotony == 6 &  sim_cor == .8 ~ "High monotony High correlation",
        )
    ) %>% drop_na() %>% 
    ggplot(aes(sim_carb_range, cpi_index)) +
    geom_line(size = 1)+
    scale_y_continuous(limits = c(0, 6))+
    scale_x_continuous(breaks = c(0, 3,6,9,12)) +
    labs(x=expression(bold(Daily~CHO~range~(g~kg^-1))), y = "CPI")+
    facet_wrap(~ scenario)+
    theme(
        strip.background = element_rect(fill = "lightblue1", color = "black", size = 1),
        strip.text = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.background = element_rect(fill = "azure")
    )



sim_b <- sim_tbl %>% 
    mutate(
        scenario = case_when(
            sim_monotony == 2 &  sim_carb_range == 2 ~ "Low monotony Low range",
            sim_monotony == 6 &  sim_carb_range == 2 ~ "High monotony Low range",
            sim_monotony == 2 &  sim_carb_range == 12 ~ "Low monotony High range",
            sim_monotony == 6 &  sim_carb_range == 12 ~ "High monotony High range",
        )
    ) %>% drop_na() %>% 
    ggplot(aes(sim_cor, cpi_index)) +
    geom_line(size = 1)+
    scale_y_continuous(limits = c(0, 6))+
    labs(x = "Correlation between CHO intake and training load", y = "CPI")+
    facet_wrap(~ scenario) +
    theme(
        strip.background = element_rect(fill = "limegreen", color = "black", size = 1),
        strip.text = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.background = element_rect(fill = "honeydew1")
    )



sim_c <- sim_tbl %>% 
    mutate(
        scenario = case_when(
            sim_cor == .2 &  sim_carb_range == 2 ~ "Low correlation Low range",
            sim_cor == .8 &  sim_carb_range == 2 ~ "High correlation Low range",
            sim_cor == .2 &  sim_carb_range == 12 ~ "Low correlation High range",
            sim_cor == .8 &  sim_carb_range == 12 ~ "High correlation High range",
        )
    ) %>% drop_na() %>% 
    ggplot(aes(sim_monotony, cpi_index)) +
    geom_line(size = 1)+
    scale_y_continuous(limits = c(0, 6))+
    labs(x = "CHO monotony (mean intake / SD)", y = "CPI")+
    facet_wrap(~ scenario)+
    theme(
        strip.background = element_rect(fill = "orange", color = "black", size = 1),
        strip.text = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.background = element_rect(fill = "moccasin")
    )

sim_explain <- ggplot()+
    annotate("text", x = 1, y = 1,
             label = "CPI =  correlation * range / monotony

Fixed values are set as follows:

Low monotony: 2
High monotony: 6
Low correlation: 0.2
High correlation: 0.8
Low carb range: 2
High carb range: 12  ")+
    scale_y_continuous(limits = c(0,2))+
    theme_void()

sim_b + sim_a + sim_c + sim_explain+
    plot_layout(nrow = 2)+
    #  plot_annotation(title = "", tag_levels = 'a') &
    plot_annotation(title = "", tag_levels = list(c('a', 'b', 'c', ' '))) &
    theme(plot.tag = element_text(size = 17, face="bold"))



# *carb v study day -------------------------------------------------------

facet_tbl <- group_data_tbl %>% 
    group_by(subject_id) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(subject_id = as.character(subject_id),
           new_id = 1:nrow(.),
           new_id = paste0("ID ", new_id)
    ) 

facet_tbl %>% 
    unnest(data) %>% 
    mutate(
        exercise_fasted = factor(exercise_fasted, labels = c("Fed", "Fasted"))
    ) %>% 
    ggplot(aes(study_day, diet_carb_g_kg))+
    geom_point(aes(color = exercise_fasted))+
    scale_color_manual(values = colors2)+
    facet_wrap(~ fct_inorder(new_id) , ncol = 6)+
    labs(x = "Day of study", y = "Dietary carbohydrate (g/kg)", color = "Pre-exercise")+
    JR_facet_theme()+
    theme(legend.position = "top",
          axis.title = element_text(face = "bold"),
          legend.text =  element_text(size = 12),
          legend.title =  element_text(face = "bold", size = 12)
    )




# *carb v load ------------------------------------------------------------

facet_tbl %>% 
    unnest(data) %>% 
    ggplot(aes(diet_carb_g_kg, exercise_load))+
    geom_point(aes(), alpha = .5, color = colors333)+
    ggpubr::stat_cor(r.accuracy = 0.01, cor.coef.name = 'r', aes(label = paste(..r.label..)))+
    geom_smooth(se = F, method = "lm")+
    facet_wrap(~ fct_inorder(new_id) , ncol = 6)+
    JR_facet_theme()+
    labs(x = "Dietary carbohydrate (g/kg)", y = "Training load (AU)", color = "Pre-exercise")+
    theme(
        axis.title = element_text(face = "bold"),
    )





# * fasting -------------------------------------------------------

min_y_fasting<- min(df_for_comparison$fasted_training_pct, na.rm = T) 
max_y_fasting <- max(df_for_comparison$fasted_training_pct, na.rm = T)+20 


emmeans(fasting_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()  

fasting_by_level <- 
    df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other")) %>% 
    ggplot(aes(subject_level, fasted_training_pct))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = habitual_diet, shape = habitual_diet), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    scale_color_manual(values = colors33)+
    scale_y_continuous(limits = c(min_y_fasting, max_y_fasting*1.2), breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
    guides(color = guide_legend(nrow = 1))+
    labs( x =  "Competitive level", y = "% Training days fasted", color = NULL, shape = NULL)+
    theme(
        legend.position = c(.5, .92),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )


emmeans(fasting_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   

fasting_by_diet <- df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other")) %>% 
    ggplot(aes(habitual_diet, fasted_training_pct))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = habitual_diet, shape = habitual_diet), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    geom_signif(comparisons = list(
        # c("Mod-CHO", "Low-CHO"),
        # c("Mod-CHO", "High-CHO"),
        c("Low-CHO", "High-CHO")
    ),
    map_signif_level = T, y_position = 105, annotations = c("*"), step_increase = .1, tip_length = 0.02, test = TukeyHSD) +
    scale_y_continuous(limits = c(min_y_fasting, max_y_fasting*1.2), breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
    scale_color_manual(values = colors33)+
    labs(x =   "Habitual diet", y = "", color = NULL, shape = NULL)+
    guides(color = guide_legend(nrow = 1))+
    theme(
        legend.position = "none", #c(.5, .92),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )

emmeans(fasting_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  

fasting_by_sex <- df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other"))  %>% 
    ggplot(aes(sex, fasted_training_pct))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = habitual_diet, shape = habitual_diet), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    geom_signif(comparisons = list(
        # c("Mod-CHO", "Low-CHO"),
        # c("Mod-CHO", "High-CHO"),
        c("Male", "Female")
    ),
    map_signif_level = T, y_position = 95, annotations = c("*"), step_increase = .1, tip_length = 0.02, test = TukeyHSD) +
    scale_y_continuous(limits = c(min_y_fasting, max_y_fasting*1.2), breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
    scale_color_manual(values = colors33)+
    labs(x = " ", y = "% Training days fasted", color = NULL, shape = NULL)+
    guides(color = guide_legend(nrow = 1))+
    theme(
        legend.position = "none", #c(.5, .92),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17))


fasting_by_level + fasting_by_diet+ fasting_by_sex + #ES_plot_fasting +
    plot_layout(ncol = 2)+
    plot_annotation(title = "", tag_levels = list(c('a', 'b', 'c', ' '))) &
    theme(plot.tag = element_text(size = 17, face="bold"))






# * training volume -------------------------------------------------------

min_y_hours <- min(group_data_tbl$subject_weekly_training_h, na.rm = T) *.9
max_y_hours <- max(group_data_tbl$subject_weekly_training_h, na.rm = T) *1.1


emmeans(training_mod, ~subject_level) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()  

hours_by_level <- 
    df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other")) %>% 
    ggplot(aes(subject_level, training_volume))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = sex, shape = sex), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    scale_color_manual(values = colors33)+
    scale_y_continuous(limits = c(min_y_hours, max_y_hours)) +
    guides(color = guide_legend(nrow = 1))+
    labs( x =  "Competitive level", y = "Weekly training (h)", color = NULL, shape = NULL)+
    theme(
        legend.position = c(.5, .92),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )

emmeans(training_mod, ~habitual_diet) %>% pairs(adjust = "Holm", infer = T) %>% as_tibble()   

hours_by_diet <-  df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other")) %>% 
    ggplot(aes(habitual_diet, training_volume))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = sex, shape = sex), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    # geom_signif(comparisons = list(
    #   c("Low-CHO", "High-CHO")
    # ),
    # map_signif_level = T, y_position = 24, annotations = c("*"), step_increase = .1, tip_length = 0.02, test = TukeyHSD) +
    scale_y_continuous(limits = c(min_y_hours, max_y_hours)) +
    scale_color_manual(values = colors33)+
    labs(x =   "Habitual diet", y = "", color = NULL, shape = NULL)+
    guides(color = guide_legend(nrow = 1))+
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)
    )

emmeans(training_mod, ~sex) %>% pairs(adjust = "Holm", infer = T)  %>% as_tibble()  

hours_by_sex <- df_for_comparison %>%  
    mutate(sex = case_when(
        sex == "M" ~ "Male",
        sex == "F" ~ "Female",
        TRUE ~ "other"))  %>% 
    ggplot(aes(sex, training_volume, color = sex, shape = sex))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(color = sex, shape = sex), size = 1.8, alpha = .95, position = position_jitter(height = 0, width = 0.1, seed = 1))+
    scale_y_continuous(limits = c(min_y_hours, max_y_hours))+
    scale_color_manual(values = colors33)+
    labs(x = " ", y = "Weekly training (h)", color = NULL, shape = NULL)+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 17))



hours_by_level + hours_by_diet+ hours_by_sex + #ES_plot_training +
    plot_layout(ncol = 2)+
    plot_annotation(title = "", tag_levels = list(c('a', 'b', 'c', ' '))) &
    theme(plot.tag = element_text(size = 17, face="bold"))


# *day before after --------------------------------------------------------

yesterday_carb_v_load <- plot_by_diet_fn(group_data_tbl, lag1_diet_carb_g_kg, exercise_load)+
    labs(x = "Training load (AU)", y = "Prior day CHO (g/kg)"
    )

carb_v_load <- plot_by_diet_fn(group_data_tbl, diet_carb_g_kg, exercise_load) +
    labs(x = "Training load (AU)", y = "Training day CHO (g/kg)"
    )

tomorrow_carb_v_load <- plot_by_diet_fn(group_data_tbl, diet_carb_g_kg, lag1_exercise_load)+
    labs(x = "Training load (AU)",  y = "Next day CHO (g/kg)"
    )

yesterday_carb_v_load + carb_v_load + tomorrow_carb_v_load +
    plot_layout(nrow = 1) +
    plot_annotation(title = "", tag_levels = 'a') &
    theme(plot.tag = element_text(size = 17, face="bold"))



