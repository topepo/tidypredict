# Parsing the model rules part of a Cubist model (i.e. without nearest-neighbor 
# adjustment). These data are contained in the `model` component of an object
# created using the Cubist package (or even the command-line version of Cubist). 

# Here is s simple example of a model file with two committees

##   id="Cubist 2.07 GPL Edition 2019-03-30"
##   prec="0" globalmean="6502.931" extrap="1" insts="0" ceiling="38989" floor="0"
##   att="outcome" mean="6502.9" sd="4707.957" min="1009" max="19999"
##   att="TotalWorkingYears" mean="11.2" sd="7.780783" min="0" max="40"
##   att="JobRole" mode="Sales_Executive"
##   att="JobLevel" mean="2" sd="1.106941" min="1" max="5"
##   redn="0.996" entries="2"
##   rules="4"
##   conds="1" cover="543" mean="2786.9" loval="1009" hival="4968" esterr="561.8"
##   type="2" att="JobLevel" cut="1" result="<="
##   coeff="2382" att="TotalWorkingYears" coeff="50"
##   conds="2" cover="124" mean="4672.2" loval="2042" hival="9724" esterr="1012.6"
##   type="3" att="JobRole" elts="Laboratory_Technician","Research_Scientist","Sales_Representative"
##   type="2" att="JobLevel" cut="1" result=">"
##   coeff="-39" att="JobLevel" coeff="2333"
##   conds="2" cover="621" mean="7136.0" loval="3886" hival="13973" esterr="988.2"
##   type="3" att="JobRole" elts="Healthcare_Representative","Human_Resources","Manufacturing_Director","Sales_Executive"
##   type="2" att="JobLevel" cut="1" result=">"
##   coeff="-1358" att="TotalWorkingYears" coeff="62" att="JobLevel" coeff="3182"
##   conds="1" cover="182" mean="16677.0" loval="11031" hival="19999" esterr="705.1"
##   type="3" att="JobRole" elts="Manager","Research_Director"
##   coeff="3355" att="TotalWorkingYears" coeff="32" att="JobLevel" coeff="3052"
##   rules="4"
##   conds="1" cover="543" mean="2786.9" loval="1009" hival="4968" esterr="569.0"
##   type="2" att="JobLevel" cut="1" result="<="
##   coeff="2672"
##   conds="2" cover="57" mean="4459.0" loval="2272" hival="5301" esterr="888.3"
##   type="2" att="TotalWorkingYears" cut="5" result="<="
##   type="2" att="JobLevel" cut="1" result=">"
##   coeff="5264" att="TotalWorkingYears" coeff="-341"
##   conds="2" cover="745" mean="6725.9" loval="2042" hival="13973" esterr="1097.2"
##   type="2" att="JobLevel" cut="1" result=">"
##   type="3" att="JobRole" elts="Healthcare_Representative","Human_Resources","Laboratory_Technician","Manufacturing_Director","Research_Scientist","Sales_Executive","Sales_Representative"
##   coeff="-1482" att="JobLevel" coeff="3620"
##   conds="1" cover="182" mean="16677.0" loval="11031" hival="19999" esterr="705.0"
##   type="3" att="JobRole" elts="Manager","Research_Director"
##   coeff="3359" att="TotalWorkingYears" coeff="32" att="JobLevel" coeff="3051"

# ------------------------------------------------------------------------------

# TODO:
# - test, test, test
# - talk to Edgar about how to structure the prediction equation
# - use of %in%? 
# - what to do about modified variable names? 

# ------------------------------------------------------------------------------

parse_model_file <- function(txt) {
  txt_rows <- stringr::str_split(txt, pattern = "\n") %>% unlist()
  
  # These are the markers for where committees start
  comm_inds <- stringr::str_which(txt_rows, "^rules=")
  num_comm <- length(comm_inds)
  # container for results for each committee
  comms <- list(length = num_comm)
  
  # Within each committee, elements are `type` (for each condition) and `coeff`
  # (for model eq). A rule starts with a `conds` element and that tells us 
  # how many lines make up the rule elements. The rule elements start right after
  # the `rules` line. Immediately after these is a single line with the information
  # on the regression equation. 
  
  for (i in seq_along(comm_inds)) {
    loc <- comm_inds[i]
    # Get the locations of the model file that encompasses the committee's rows
    if (i < num_comm) {
      uppr <- comm_inds[i + 1] - 1
    } else {
      uppr <- length(txt_rows)
    }
    num_rules <- rule_info(txt_rows[loc])
    comm_data <-
      dplyr::tibble(
        rule_num = 1:num_rules,
        rule = NA,
        reg_eq = NA
      )
    # Where are the lines that show the `conds` information?
    attr_inds <- find_cond_info(txt_rows, loc, uppr)
    cond_att <- purrr::map_dfr(attr_inds, parse_cond, txt = txt_rows)
    comm_data <- dplyr::bind_cols(comm_data, cond_att)
    
    # Loop over all of the rules and get their rule conditions
    for (j in seq_along(attr_inds)) {
      
      if (comm_data$conds[j] == 0) {
        # When the rule includes everything
        atts <- "TRUE"
      } else {
        att_loc <- attr_inds[j] + 1:comm_data$conds[j]
        atts <- purrr::map_chr(txt_rows[att_loc], make_conds)
        atts <- stringr::str_c(atts, collapse = " & ")
      }
      
      comm_data$rule[j] <- atts
    }
    
    # Get regression equations
    eq_ind <- attr_inds + comm_data$conds + 1
    comm_data$reg_eq <- purrr::map_chr(txt_rows[eq_ind], get_reg_eq)
    comm_data$committee <- i
    comms[[i]] <- comm_data
  }
  res <- 
    dplyr::bind_rows(comms) %>% 
    dplyr::mutate(
      pred_eq = stringr::str_c("(", rule, ")*(", reg_eq, ")"),
      pred_eq = rlang::parse_exprs(pred_eq)
    ) %>% 
    dplyr::select(committee, rule_num, rule, reg_eq, pred_eq, cover, mean, 
                  lower = loval, upper = hival, err = esterr)
  res
}

# ------------------------------------------------------------------------------

# Locate what line in the model file (from lines `strt` to `stp`) define the 
# conditions 
find_cond_info <- function(txt, strt = 0, stp = 0) {
  txt <- txt[(strt + 1):(stp - 1)]
  stringr::str_which(txt, "^conds=") + strt
}

# Parse the condition line 
parse_cond <- function(ind, txt) {
  entires <- stringr::str_split(txt[ind], " ") %>% unlist()
  tmp <- purrr::map(entires, ~ stringr::str_split(.x, pattern = "=") %>% unlist())
  nms <- purrr::map_chr(tmp, purrr::pluck, 1)
  vals <- purrr::map(tmp, stringr::str_remove_all, pattern = "\"")
  vals <- purrr::map_dbl(vals, ~ purrr::pluck(.x, 2) %>% as.numeric())
  names(vals) <- nms
  as.data.frame(t(vals))
}

# Determine the number of rules
rule_info <- function(txt) {
  txt <- stringr::str_remove_all(txt, "\"")
  txt <- stringr::str_remove(txt, "^rules=")
  as.integer(txt)
}


# ------------------------------------------------------------------------------

# Parse the lines that have the regression equation for the rule. For example, 
# the line might be 
#
##   coeff="1.88" att="Sepal.Width" coeff="0.65" att="Petal.Length" coeff="0.709" att="Petal.Width" coeff="-0.55"
#
# for a regression line for
#
##   outcome = 1.88 + 0.709 Petal.Length - 0.55 Petal.Width + 0.65 Sepal.Width

get_reg_eq <- function(txt) {
  entires <- stringr::str_split(txt, " ") %>% unlist()
  n <- length(entires)
  vals <- purrr::map_chr(entires, reg_terms)
  lp <- vals[1]
  vals <- vals[-1]
  if (length(vals) > 0) {
    n_elem <- length(vals)
    if (n_elem %% 2 != 0) {
      stop("number of remaining terms not even", call. = FALSE)
    }
    n_terms <- n_elem/2
    split_terms <- split(vals, rep(1:n_terms, each = 2))
    terms <- purrr::map_chr(split_terms, paste_slopes)
    terms <- stringr::str_c(terms, collapse = " + ")
    lp <- stringr::str_c(lp, " + ", terms, collapse = "")
  }
  lp
}

reg_terms <- function(txt) {
  if (stringr::str_detect(txt, "^coeff")) {
    val <- stringr::str_remove(txt[1], "coeff=\"")
    val <- stringr::str_remove(val, "\"")
  } else {
    val <- stringr::str_remove(txt[1], "att=\"")
    val <- stringr::str_remove(val, "\"") 
  }
  val
}

paste_slopes <- function(txt) {
  stringr::str_c("(", txt[2], "*", txt[1], ")", sep = " ")
}

# ------------------------------------------------------------------------------

# Parse and translate the conditions that define the rule. For example, a 
# rule like 
#
##   Year_Built > 1952
##   Bsmt_Exposure in {Av, Mn, No, No_Basement}
##   Gr_Liv_Area <= 1692
##   Kitchen_Qual = Excellent
# 
# shows up in the model file as
# 
##   conds="4" cover="35" mean="5.364152" loval="4.96379" hival="5.5214" esterr="0.039525"
##   type="3" att="Kitchen_Qual" elts="Excellent"
##   type="2" att="Gr_Liv_Area" cut="1692" result="<="
##   type="3" att="Bsmt_Exposure" elts="Av","Mn","No","No_Basement"
##   type="2" att="Year_Built" cut="1952" result=">"
#
# The type is the nature of the split (2 = numeric and 3 = qualitative). I've 
# never seen a type = 1

make_conds <- function(txt) {
  res <- purrr::map_chr(txt, single_cond)
  res <- stringr::str_c(res, collapse = " & ")
  res <- stringr::str_replace_all(res, "\"", "'")
  res
}

# parse a single condition line
single_cond <- function(txt) {
  if (stringr::str_detect(txt, "type=\"2")) {
    res <- cond_2(txt)
  } else {
    res <- cond_3(txt)
  }
  res
}

single_cond2 <- function(txt) {
  if (stringr::str_detect(txt, "type=\"2")) {
    res <- cond_2(txt)
  } else {
    res <- cond_3(txt)
  }
  res
}


# quantitative splits
cond_2 <- function(txt) {
  entires <- stringr::str_split(txt, " ") %>% unlist
  rms <- "(att=\")|(cut=\")|(result=\")"
  entires <- purrr::map_chr(entires, stringr::str_remove_all, rms)
  entires <- purrr::map_chr(entires, stringr::str_remove_all, "\"")
  stringr::str_c("(", entires[2], entires[4], entires[3], ")", sep = " ")
}

# qualitative splits
cond_3 <- function(txt) {
  entires <- stringr::str_split(txt, " ") %>% unlist
  var_name <- entires[2]
  var_name <- stringr::str_remove(var_name, "att=\"")
  var_name <- stringr::str_remove(var_name, "\"")
  elts <- entires[3]
  elts <- stringr::str_remove(elts, "elts=")
  stringr::str_c("(", var_name, " %in% c(", elts, ") )", sep = " ")
}


