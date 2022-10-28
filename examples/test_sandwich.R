old_lib_loc <- .libPaths()[1]
new_lib_loc = "/Library/Frameworks/R.framework/Versions/4.2/Resources/dev_installs"

# install.packages("PLNmodels")
# remotes::install_github("pln-team/PLNmodels", ref = "dev")

# devtools::load_all("~/OneDrive/UW/Trapnell/hooke_sandwich//")


# library(PLNmodels, lib.loc=old_lib_loc)
# library(PLNmodels) #, lib.loc=new_lib_loc)

# model <- PLN(Abundance ~ 0 + . + offset(log(O)), data = data, control = list(trace = 0))
# Theta_hat <- coef(model)
# model$get_vcov_hat("wald", Y, X)
# Theta_se_wald <- standard_error(model)
# model$get_vcov_hat("louis", Y, X)
# Theta_se_louis <- standard_error(model)
# model$get_vcov_hat("sandwich", Y, X)
# Theta_se_sandwich <- standard_error(model)


devtools::load_all("~/OneDrive/UW/Trapnell/bin/PLNmodels/")
library(monocle3)
library(splines)
library(hooke)

# ------------------------------------------------------------------------------


cds = readRDS("~/OneDrive/UW/Trapnell/zf-atlas-portal/data/full_gap_hf_ctrl_ref_mito-filt_50k_sub_anno_cds.RDS")

meso_cds = cds[, colData(cds)$major_group == "mesoderm"]
meso_cds = meso_cds[, tidyr::replace_na(between(colData(meso_cds)$timepoint, 18, 24), F)]

ccs = hooke::new_cell_count_set(meso_cds,
                                sample_group = "embryo",
                                cell_group = "cell_type_broad")


ccm = hooke::new_cell_count_model(ccs, main_model_formula_str = "~ns(timepoint, df = 2)")



best_model = ccm@best_full_model
X = ccm@full_model_family$responses
Y = ccm@full_model_family$covariates
# qr(Dn)$rank
debug(best_model$vcov_sandwich)
best_model$get_vcov_hat("sandwich",X, Y)
vcov(best_model)[1:5,1:5]
best_model$get_vcov_hat("wald",X, Y)

best_model = ccm@best_reduced_model
responses = ccm@reduced_model_family$responses
covariates = ccm@reduced_model_family$covariates
best_model$get_vcov_hat("sandwich",responses, covariates)

vcov(best_model)[1:5,1:5]


vhat_type = c("sandwich", "Wald", "Louis")


best_model = ccm@best_reduced_model
X = ccm@reduced_model_family$responses
Y = ccm@reduced_model_family$covariates
best_model$get_vcov_hat("sandwich",X, Y)


# TRY IT on SIMULATION ---------------------------------------------------------


ccs = readRDS(file = "~/OneDrive/UW/Trapnell/hooke_dev/benchmarking/sim_ccs_size1000_seed1_effect_0.25.rds")
ctrl_penalty_matrix = readRDS("~/OneDrive/UW/Trapnell/hooke_dev/benchmarking/data/ctrl_penalty_matrix.rds")
ccm = hooke::new_cell_count_model(ccs, main_model_formula_str = "~ genotype", penalty_matrix = ctrl_penalty_matrix)

undebug(new_cell_count_model)
ccm_sandwich = new_cell_count_model(ccs, main_model_formula_str = "~ genotype", penalty_matrix = ctrl_penalty_matrix, vhat_method = "wald")



best_model = ccm@best_full_model
X = ccm@full_model_family$responses
Y = ccm@full_model_family$covariates
# qr(Dn)$rank
# debug(best_model$vcov_sandwich)
best_model$get_vcov_hat("sandwich",X, Y)
vcov(best_model)[1:5,1:5]
best_model$get_vcov_hat("wald",X, Y)



# THING TO check qr(Dn)$rank

pln_data = PLNmodels::prepare_data(counts = counts(ccs) + 0,
                        covariates = colData(ccs) %>% as.data.frame,
                        offset = size_factors(ccs))

model = do.call(PLNmodels::PLNnetwork, 
                args=list("Abundance ~ 1", 
                          data = pln_data))

best_model = comb_ccm@best_full_model
X = comb_ccm@full_model_family$responses
Y = comb_ccm@full_model_family$covariates

                   