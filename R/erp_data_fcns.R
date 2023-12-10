get_electrode <- function(data, electrode_no) {
    data[electrode_no, , ]
}


idx_to_msec <- function(start_idx, end_idx, pre_stimulus) {
    time_msec <- c(start_idx, end_idx) * 2
    return(time_msec - prestimulus)
}

msec_to_idx <- function(start_msec, end_msec, pre_stimulus) {
    time_idx <- (c(start_msec, end_msec) + pre_stimulus) / 2
    return(time_idx)
}

## use this
sec_to_msec <- function(x, pre_stimulus) {
    x * 1000 - 200
}


get_unique_id <- function(data, id_len = 10) {
    all_id <- data[1, c(1:2, 5:id_len), ]
    all_id_tbl <- tbl_df(t(all_id))
    unique_id <- distinct(all_id_tbl)
    return(list(all_id_tbl = all_id_tbl, unique_id = unique_id))
}

condition_len <- function(all_id_tbl, unique_id_idx) {
    sum_id <- 0
    n_trial <- nrow(all_id_tbl)
    len_id <- length(unique_id_idx)
    for (i in 1:n_trial) {
        # print(i)
        if(sum(all_id_tbl[i, ] == unique_id_idx) == len_id) sum_id = sum_id + 1
    }
    return(sum_id)
}

condition_which <- function(all_id_tbl, unique_id_idx) {
    nr <- nrow(unique_id_idx)
    which_lst <- list()
    n_trial <- nrow(all_id_tbl)
    for(j in 1:nr) {
        which_vec <- c()
        for (i in 1:n_trial) {
            if(sum(all_id_tbl[i, ] == unique_id_idx[j, ]) == 8) which_vec <- c(which_vec, i)
        }
        which_lst[[j]] <- which_vec
    }
    return(which_lst)
}

get_VOT <- function(data, VOT_vec) {
    data_lst <- list()
    for (i in 1:length(VOT_vec)) {
        idx_VOT <- which(data[1, 5, ] == VOT_vec[i])
        data_i <- data[, , idx_VOT]
        data_lst[[i]] <- data_i
    }
    return(data_lst)
}

get_bias_cond <- function(data, bias_cond = c(0, 1)) {
    data_lst <- list()
    for (i in 1:length(bias_cond)) {
        idx_bias <- which(data[1, 6, ] == bias_cond[i])
        data_i <- data[, , idx_bias]
        data_lst[[i]] <- data_i
    }
    return(data_lst)
}

