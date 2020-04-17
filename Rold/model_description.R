#' Function for defining the elements of a model object
#'
#' @description The function returns a description, in terms of variable names, names and role, for the selected table within a Dynamic TOPMODEL object.
#'
#' @param type the names of the data frame within the Dynamic TOPMODEL object
#' @param include_states shoould states be returned
#' @param include_tmp should tempory variables used in simulation be returned
#'
#' @return A data frame descriping the properties of one of the tables found within the model
#'
model_description <- function(type=c("hillslope","channel","point_inflow","gauge"),include_states=FALSE,include_tmp=FALSE){
    type <- match.arg(type)

    if(type=="hillslope"){
        out <- data.frame(name = c("id","atb_bar","s_bar","area","delta_x","sz_dir","sf_dir", # attributes associated with catchment HSU
                                   "precip","pet", # names of input series
                                   "q_sfmax","s_rzmax","s_rz0","ln_t0","m","t_d","t_sf", # parameter names
                                   "s_sf","s_rz","s_uz","s_sz","sum_l_sz_in","l_sz","l_szmax", # states
                                   "p","e_p","e_t","l_sf","q_sf_rz","q_rz_sf","q_rz_uz","q_uz_sz","q_uz_sf","q_sz_sf","e_t","l_sf","sum_l_sz_in_t","l_sz_t","Q_minus_t","Q_plus_t","Q_minus_tDt"), ## tempory stores not needed for next timestep
                          role = c(rep("attribute",7),
                                   rep("data_series",2),
                                   rep("parameter",7),
                                   rep("state",7),
                                   rep("tmp",17)),
                          type = c("integer",rep("numeric",4),rep("list",2),
                                   rep("character",2),
                                   rep("character",7),
                                   rep("numeric",7),
                                   rep("numeric",17)
                                   ),
                          stringsAsFactors=FALSE)

    }
    if(type=="channel"){
        out <- data.frame(name= c("id","area","length","flow_dir", # states
                                  "precip","pet", # inputs
                                  "v_ch", # parameters
                                  "sum_l_sz_in", #state
                                  "p","e_p","l_sf","s_ch","sum_l_sz_in_t"),
                          role = c(rep("attribute",4),
                                   rep("data_series",2),
                                   rep("parameter",1),
                                   rep("state",1),
                                   rep("tmp",5)),
                          type = c("integer",rep("numeric",2),"list",
                                   rep("character",2),
                                   rep("character",1),
                                   rep("numeric",1),
                                   rep("numeric",5)),
                          stringsAsFactors=FALSE)
    }

    if(type=="point_inflow"){
        out <- data.frame(
            name = c("name","id","fraction"),
            type=c("character","integer","numeric"),
            role = c("data_series",rep("property",2)),
            stringsAsFactors=FALSE)
    }

    if(type=="gauge"){
        out <- data.frame(
            name = c("name","id","fraction"),
            type=c("character","integer","numeric"),
            role = c("output_label",rep("property",2)),
            stringsAsFactors=FALSE)
    }

    if(!include_states){
        out <- out[out$role!="state",,drop=FALSE]
    }
    if(!include_tmp){
        out <- out[out$role!="tmp",,drop=FALSE]
    }

    return(out)
}
