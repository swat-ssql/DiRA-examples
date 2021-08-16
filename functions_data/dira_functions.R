library(stargazer)
library(foreign)
library(broom)
library(tidyverse)
# library(modelr)
library(plotly)
library(psych) #for checking direction of principal components... need to figure out my own method
library(GPArotation)


load("/Users/mfoster1/Dropbox/research/dira/DiRA-examples-master/functions_data/set_model_info.RData")

#main functions:
# dira_models <- models.dira(y ~ x1+x2+x3+x4+...+xn, data, model.names  )
# tidy.dira(dira_models)
# glance.dira(dira_models)
# stargazer.dira(formula, data, model.names)
# stargazer.dira(dira_models) # alternate format?
# plot_ly(x, y, z, data) %>%  #add markers separately
#   add_3d.directions(formula = optional, data = optional, model.names)

# plot_ly(x, y, z, data) %>%  #add markers separately
#   add_2d.directions(formula = optional, data = optional, model.names)

# Coethnic reults testing path: ~/Dropbox/fr17to_swat/2021 Spring/polsci_methods/quantProblemSet/Coethnic Results_diraTesting_2021-07-10.Rmd
# This path: ~/Dropbox/research/dira/dira_functions_2021-07-19.R
# old dira functions path: ~/Dropbox/research/congressNotFolders/createPCinfo/functionsToCreatePC_v3_createOtherExamples.R


# create models ########
#output: named list of of linear models:
#     lm_x1.predby.x2, lm_x2.predby.x1, lm_y.predby.x1, lm_y.predby.x2, lm_multiple
#     direction.x2.predby.x1.resids, direction.x1.predby.x2.resids, principal.components
#     need to add lm_y.predby.x1.orthx2,
#           lm.newdirection, lm.newdirection.orth
#calls no new functions
models.dira <- function( formula_dira , data,
                     model.names = c("lm_y.predby.x1.w.resids",
                                    "lm_y.predby.x2.w.resids",
                                    "lm_pc.w.orth"),
                     direction.start.x1 = 0, direction.end.x1 = 1,
                     direction.start.x2 = 1, direction.end.x2 = 0,
                     units.of.x1 = F){
  model.names = c("lm_multiple",  #always create directions in this order
                 "lm_x1.predby.x2", "lm_x2.predby.x1", model.names)
  model.names <- unique(model.names)

  models.dira <- list()
  models.dira[["lm_multiple"]] <- lm(formula_dira, data = data,
                                     na.action = na.exclude) #always show lm_multiple
  model.names <- model.names[-which(model.names == "lm_multiple")]

  for(model.name in model.names){
    newmodel <- lm.direction(formula = formula_dira, data = data,
                             model.name = model.name, models.dira = models.dira,
                             direction.start.x1 = 0, direction.end.x1 = 1,
                             direction.start.x2 = 1, direction.end.x2 = 0,
                             units.of.x1 )
    models.dira[[model.name]] <- newmodel
  }
  model.names <- names(models.dira) # reorder so xi~xj appears last
  short.model.names <- model.names[-which(model.names %in% c("lm_x1.predby.x2", "lm_x2.predby.x1"))]
  model.names <- c(short.model.names,c("lm_x1.predby.x2", "lm_x2.predby.x1") )
  models.dira <- models.dira[model.names]
  models.dira
}

# create regression model based on direction
lm.direction <- function(formula, data,
                           model.name, models.dira,
                         direction.start.x1 = 0, direction.end.x1 = 1,
                         direction.start.x2 = 1, direction.end.x2 = 0,
                         units.of.x1 =F ){
  x1.var <- attr(terms(formula), which = "term.labels")[1]
  x2.var <- attr(terms(formula), which = "term.labels")[2]
  y.var <- as.character(formula[[2]])[1]
  num.terms <- length(attr(terms(formula), which = "term.labels"))

  if(model.name %in% c("lm_pc.w.orth", "lm_pc","lm_new.direction")){
    rotation.matrix <- make.rotation.matrix(data, coefficients = c(x1 = x1.var,  x2 = x2.var),
                                            model.name, units.of.x1, direction.start.x1,
                                            direction.end.x1, direction.start.x2, direction.end.x2)
  }
  if(num.terms > 2){
    covariates <- attr(terms(formula), which = "term.labels")[3:length(attr(terms(formula), which = "term.labels"))]
    covariate.text <- paste(covariates, collapse = " + ")
    covariate.text <- paste0(" + ", covariate.text)
    interactionterm.text <- paste0( "+ ", x1.var,":",x2.var) #this should help all covariates
    covariate.text <- gsub(interactionterm.text, "",covariate.text , fixed = T)
  }else{
    covariate.text <- ""
  }
  if(model.name == "lm_x1.predby.x2"){
    formula <- formula(paste0(x1.var, "~", x2.var))
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if(model.name == "lm_x2.predby.x1"){
    formula <- formula(paste0(x2.var, "~", x1.var))
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if(model.name == "lm_y.predby.x1"){
    # associated with direction.x2.predby.x1
    formula.text <- paste0(y.var, "~", x1.var, covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if(model.name == "lm_y.predby.x2"){
    # associated with direction.x2.predby.x1
    formula.text <- paste0(y.var, "~", x2.var, covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if("lm_y.predby.x1.w.resids" == model.name){
    # associated with direction.x2.predby.x1 & direction.x2.predby.x1.resids
    data <- data %>% mutate(rotation.x2.resids = augment(models.dira[["lm_x2.predby.x1"]],
                                                         data = data, na.action = "na.exclude")$.resid)
    formula.text <- paste0(y.var,  "~ ", x1.var, " + rotation.x2.resids", covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data)
  }else if("lm_y.predby.x2.w.resids" == model.name){
    # associated with direction.x1.predby.x2 & direction.x1.predby.x2.resids
    data <- data %>% mutate(rotation.x1.resids = augment(models.dira[["lm_x1.predby.x2"]],
                                                         data = data, na.action = "na.exclude")$.resid)
    formula.text <- paste0(y.var,  "~ ", x2.var, " + rotation.x1.resids", covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data)
  }else if("lm_pc.w.orth" == model.name){  #associated with principal.component & principal.component.orth
    slope <- rotation.matrix[2,1]/rotation.matrix[1,1]  #b/a
    data.w.rotation <- create.rotated.variables(rotation.matrix, data, units.of.x1, slope)
    data$rotation.pc <- data.w.rotation$rotation1  #rotation.pc
    data$rotation.pc.orth  <- data.w.rotation$rotation.orth #rotation.pc.orth
    formula.text <- paste0(y.var,  "~ rotation.pc + rotation.pc.orth", covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data)
  }else if("lm_pc" == model.name){  #associated with principal.component
    slope <- rotation.matrix[2,1]/rotation.matrix[1,1]  #b/a
    data.w.rotation <- create.rotated.variables(rotation.matrix, data, units.of.x1, slope)
    data$rotation.pc <- data.w.rotation$rotation1  #rotation.pc
    formula.text <- paste0(y.var,  "~ rotation.pc", covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data)
    # tidy(newmodel)
  }else if("lm_new.direction" == model.name){# associated with new.direction & new.direction.orth
    slope <- rotation.matrix[2,1]/rotation.matrix[1,1]  #b/a
    data.w.rotation <- create.rotated.variables(rotation.matrix, data,
                                                units.of.x1 = units.of.x1, slope = slope)
    data$rotation.new.direction <- data.w.rotation$rotation1 #rotation.pc
    data$rotation.new.direction.orth  <- data.w.rotation$rotation.orth #rotation.pc.orth
    formula.text <- paste0(y.var,  "~ rotation.new.direction + rotation.new.direction.orth", covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data)
  }else if(model.name %in% c("lm_y.predby.x1.subset.x2.0",
                             "lm_y.predby.x1.subset.x2.1") ){# when subsetting data based on dummy variables
    constantval <- as.numeric(gsub("lm_y.predby.x1.subset.x2." , "", model.name))
    data <- data[ which(data[, x2.var] == constantval), ]
    formula.text <- paste0(y.var, "~", x1.var, covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if(model.name %in% c("lm_y.predby.x2.subset.x1.0",
                             "lm_y.predby.x2.subset.x1.1") ){# when subsetting data based on dummy variables
    constantval <- as.numeric(gsub("lm_y.predby.x2.subset.x1." , "", model.name))
    data <- data[ which(data[, x1.var] == constantval), ]
    formula.text <- paste0(y.var, "~", x2.var, covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else if(model.name  == "lm_y.predby.x1.subset.x1ORx2.1" ){# when subsetting data based on dummy variables
    data <- data[ which(data[, x1.var] == 1 | data[, x2.var] == 1), ]
    formula.text <- paste0(y.var, "~", x1.var, covariate.text)
    formula <- formula(formula.text)
    newmodel <- lm(formula, data = data, na.action = na.exclude)
  }else{
    stop("invalid model name", call. = FALSE)
  }
  newmodel
  # data
}

#used by lm.directions, create.x1x2seq, make.model.info.matrix
make.rotation.matrix <- function(data, coefficients, model.name, units.of.x1,
                                 direction.start.x1 = 0, direction.end.x1 = 1,
                                 direction.start.x2 = 1, direction.end.x2 = 0 ){
  if("lm_pc" == model.name | "lm_pc.w.orth" == model.name){
    pc.data <- subset(data, select = c(coefficients[["x1"]], coefficients[["x2"]]))
    rotation.matrix <- prcomp(formula(paste0("~", coefficients[["x1"]], "+",
                                             coefficients[["x2"]])),
                              data = pc.data, na.action = na.exclude)$rotation
    rotation.matrix <- adjust.rotation.matrix(rotation.matrix)
    if(units.of.x1) rotation.matrix <- rotation.matrix/rotation.matrix[1,1]
  }else if("lm_new.direction" == model.name){
    slope <- (direction.end.x2 - direction.start.x2)/(direction.end.x1 - direction.start.x1)  #new direction rotation
    rotation.matrix <- matrix(c(1, slope, -slope,1), nrow = 2)
    row.names(rotation.matrix) <- c(coefficients[["x1"]], coefficients[["x2"]])
  }else{  #filler in case an empty rotation matrix causes an error; not used anywhere
    rotation.matrix <- matrix(c(1,1,1,1), nrow = 2)
    row.names(rotation.matrix) <- c(coefficients[["x1"]], coefficients[["x2"]])
  }
  rotation.matrix
}

#make sure rotation is in correct direction (increasing in terms of x1)
# normalize matrix
#can pass it through many times with no harm bc normalized matrix will
# always renormalize to itself
adjust.rotation.matrix <- function(rotation.matrix){
  #make sure a positive a change in x1 is associated with a b change in x2
  if(rotation.matrix[1,1] < 0){
    rotation.matrix <- -1*rotation.matrix
  }
  a <- rotation.matrix[1, 1]
  b <- rotation.matrix[2, 1]
  #normalize the rotation matrix
  rotation.matrix<- rotation.matrix*(1/sqrt((a^2 + b^2) ))
  rotation.matrix
}

create.rotated.variables <- function(rotation.matrix, data, units.of.x1, slope = NULL,
                                     mean.x1 = NULL, mean.x2 = NULL){
  rotation.matrix <- adjust.rotation.matrix(rotation.matrix) #normalize and positive direction
  a <- rotation.matrix[1, 1]
  b <- rotation.matrix[2, 1]
  x1.var <- names(rotation.matrix[, 1])[1]
  x2.var <- names(rotation.matrix[, 1])[2]
  if(is.null(mean.x1) & is.null(mean.x2)){
    #calculate means when data is the original data
    #use passed in mean values when passing in predicted val data
    mean.x1 <- mean(data[, x1.var], na.rm = T)
    mean.x2 <- mean(data[, x2.var], na.rm = T)
  }

  data$rotation1 <- a*data[, x1.var] +b*data[, x2.var] - a*mean.x1 - b*mean.x2
  data$rotation.orth <- -b*data[, x1.var] +a*data[, x2.var] + b*mean.x1 - a*mean.x2
  if(units.of.x1){
    if(is.null(slope)){
      stop("must provide slope of rotation to use units.of.x1", call. = FALSE)
    }
    units.of.x1.adjustment <- 1/sqrt(1^2 + slope^2)
    data$rotation1 <- data$rotation1*units.of.x1.adjustment
    data$rotation.orth <- data$rotation.orth*units.of.x1.adjustment
  }
  data

}

# create stargazer output ####

#print results of full and both omitted models
#     calls models.dira()
stargazer.dira <- function(formula, data,
                           model.names =  c("lm_multiple", 'lm_y.predby.x1','lm_y.predby.x2' ),
                           direction.start.x1 = 0, direction.end.x1 = 1,
                           direction.start.x2 = 1, direction.end.x2 = 0,
                           units.of.x1 = F, print.direction = F,
                           print.model.descriptions = T,
                           covariate.labels = NULL,
                           dep.var.caption = NULL,
                           column.labels = NULL, ...){
  if( !("lm_multiple" %in% model.names)){
    model.names <- c("lm_multiple", model.names) #always print the multiple regression results first
  }
  dira_models <- models.dira(formula_dira = formula, data, model.names,
                                    direction.start.x1 = 0, direction.end.x1 = 1,
                                    direction.start.x2 = 1, direction.end.x2 = 0,
                                    units.of.x1 )
  coefficients <-names(dira_models[["lm_multiple"]]$model)
  names(coefficients)  <- c("y", "x1","x2")
  
  # print(paste(names(coefficients), collapse = " "))
  coefficient.info <- data.frame(var.std = names(coefficients),
                                 varnames = unname(coefficients),
                                 var.description = unname(coefficients))

  if(is.null(dep.var.caption) )  {
    dep.var.caption <-  coefficients["y"]
  }
  coefficient.info$var.description[1] <- dep.var.caption
  
  if(!is.null(covariate.labels)){
    if(length(coefficients) == 3){
      # coefficients[2:length(coefficients)] <- covariate.labels[1:length(coefficients)-1]
      coefficient.info$var.description[2:length(coefficients)] <- covariate.labels[1:length(coefficients)-1]
    }else if(length(coefficients) > 3){
      num.controls <- length(coefficients) - 3
      tmp <- c(covariate.labels[1:2],
               covariate.labels[(length(covariate.labels)-num.controls):(length(covariate.labels)-1) ])
      coefficient.info$var.description[2:length(coefficients)] <- tmp
    }
  }

  model.info <- make.model.info.matrix(data, dira_models, coefficient.info, model.names,
                                       direction.start.x1, direction.end.x1,
                                       direction.start.x2, direction.end.x2,
                                       units.of.x1)
  #associate directions with model labels and model descriptions
  model.info$models <- factor(model.info$models, levels = model.names)
  model.info <- model.info[order(model.info$models),  ]
  if(is.null(column.labels)) column.labels <- as.character(model.info$model.labels )
  stargazer(dira_models[model.names], align = T,
            dep.var.labels.include = F,
            covariate.labels = covariate.labels,
            dep.var.caption = dep.var.caption,
            column.labels = column.labels, ...)
  if(print.model.descriptions){
    for(model.name in model.names){
        cat(model.info$model.descriptions[which(model.info$models == model.name)], "\n")

    }
  }else if(units.of.x1 & sum(c("lm_pc.w.orth","lm_pc","lm_new.direction") %in% model.names) > 0) {
    scaled.models <- c("y ~ PC1 + PC2","y ~ PC1","y ~ new.direction + new.dir.orth")[c("lm_pc.w.orth","lm_pc","lm_new.direction") %in% model.names]
    numcommas <- length(scaled.models) - 2
    if(numcommas < 0) numcommas <- 0
    tmptext <- ""
    punctuationtext <- c(rep(", ", numcommas), rep(" and ", length(scaled.models)-1), "")
    for(i in 1:length(scaled.models)){
      tmptext <- paste0(tmptext, scaled.models[i], punctuationtext[i])
    }
    cat("Note: The coefficients for", tmptext, "have been scaled to be interpreted in units of x1.")
  }
  if(print.direction){
    print.models <- c("lm_x1.predby.x2","lm_x2.predby.x1")
    stargazer(dira_models[print.models], align = T, type = 'text',
              dep.var.labels.include = F,
              dep.var.caption = "",
              column.labels = c(coefficients["x1"], coefficients["x2"]) )
  }
}




#associate directions with model labels and model descriptions
#   called by stargazer.dira
#   could add it to the matrix of models, var.names, surface.directions, linecolors, linetypes if I can figure out eval() or symbol()
make.model.info.matrix <- function(data, dira_models, coefficient.info, model.names,
                                   direction.start.x1 = 0, direction.end.x1 = 1,
                                   direction.start.x2 = 1, direction.end.x2 = 0,
                                   units.of.x1 = F){
  model.info <- model.information
  model.info <- model.info %>%
    select(models, model.labels, model.descriptions, model.descriptions.not.units.of.x1) %>%
    filter(models %in% model.names)
  model.info <- model.info[!duplicated(model.info$models), ]

  for(i in 1:nrow(model.info)){
    model.name <- model.info$models[i]
    rotation.matrix <- make.rotation.matrix(data, coefficients = c(x1 = coefficient.info$varname[2],
                                                                   x2 = coefficient.info$varname[3]),
                                            model.name, units.of.x1, direction.start.x1,
                                            direction.end.x1, direction.start.x2, direction.end.x2)

    string.for.subs = model.info$model.descriptions[i]
    if( (model.name == "lm_pc" | model.name == "lm_pc.w.orth") & !units.of.x1){
      string.for.subs <- model.info$model.descriptions.not.units.of.x1[i]
    }

    model.description <- matrix.info.substitutions(string.for.subs = string.for.subs,
                                                   dira_models, coefficient.info,
                                                   model.name, rotation.matrix, units.of.x1)

    if(length(coefficients) > 3 & !grepl("no information available", model.description) ){
      model.description <- paste0(model.description, " All other covariates are assumed to be held constant.")
    }
    model.info$model.descriptions[i] <- model.description
    model.info$model.labels[i] <- matrix.info.substitutions(model.info$model.labels[i],
                                                            dira_models, coefficient.info, model.name,
                                                            rotation.matrix, units.of.x1)
  }

  model.info

}

matrix.info.substitutions <- function(string.for.subs, dira_models, coefficient.info,
                                      model.name, rotation.matrix, units.of.x1){
  slope <- rotation.matrix[2,1]/rotation.matrix[1,1]

  string.for.subs <- gsub("thismodelhere_estimate2",
                          tidy(dira_models[[model.name]])$estimate[2] %>%round(digits = 3), string.for.subs)
  string.for.subs <- gsub("rotation.matrix_a", round(rotation.matrix[1,1], digits = 3), string.for.subs)
  string.for.subs <- gsub("rotation.matrix_b", round(rotation.matrix[2,1], digits = 3), string.for.subs)
  string.for.subs <- gsub("rotation.slope", round(slope, digits = 3), string.for.subs)
  string.for.subs <- gsub("x1.var", coefficient.info$var.description[2], string.for.subs)
  string.for.subs <- gsub("x2.var", coefficient.info$var.description[3], string.for.subs)
  string.for.subs <- gsub("y.var", coefficient.info$var.description[1], string.for.subs)
  string.for.subs <- gsub("xxxx", which(model.names == model.name), string.for.subs)
  string.for.subs <- gsub("lm_multiple_estimate2",
                          tidy(dira_models[["lm_multiple"]])$estimate[2] %>%round(digits = 3),
                          string.for.subs)
  string.for.subs <- gsub("lm_multiple_estimate3",
                          tidy(dira_models[["lm_multiple"]])$estimate[3] %>%round(digits = 3),
                          string.for.subs)
  string.for.subs <- gsub("lm_x2.predby.x1_estimate2",
                          tidy(dira_models[["lm_x2.predby.x1"]])$estimate[2] %>%round(digits = 3),
                          string.for.subs)
  string.for.subs <- gsub("lm_x1.predby.x2_estimate2",
                          tidy(dira_models[["lm_x1.predby.x2"]])$estimate[2] %>%round(digits = 3),
                          string.for.subs)
  string.for.subs
}

# create models summaries with tidy() and glance() ####

tidy.dira <- function(dira_models){
  dira.data <- tibble(model.names = names(dira_models), models = dira_models)
  dira.data %>% mutate(tidied = map(models, tidy)) %>% unnest(tidied) %>% select(-models)
}

glance.dira <- function(dira_models){
  dira.data <- tibble(model.names = names(dira_models), models = dira_models)
  dira.data %>% mutate(glance = map(models, broom::glance)) %>% unnest(glance)%>% select(-models)
}
# create data for plotting #######

#returns directional data; rename to createPlottingData?
#    calls 1. models.dira, 2. base.data, 3. create.x1x2seq, 4. create.predictions
#    returns directional data
create.directions.to.plot <- function(formula, data, models = c("lm_multiple"),
                                define.constant.x1 = "mean", define.constant.x2 = "mean",
                                direction.start.x1 = 0, direction.end.x1 = 1,
                                direction.start.x2 = 1, direction.end.x2 = 0, units.of.x1 = F){
  
  model.names <- models
  
  model.info <- model.information
  model.info <- model.info %>% filter(models %in% model.names)
  dira_models <- models.dira(formula, data, model.names, units.of.x1 = F)

  if(is.null(data$direction)){ #in case data already altered in prior plotting
    combined.data <- base.data(data)
  }else{
    combined.data <- data
  }
  for(i in 1:nrow(model.info)){
    direction <- model.info$surface.directions[i]
    model.name <- model.info$models[i]
    predicteddata <- create.x1x2seq(dira_models = dira_models, data = data,
                                    direction = direction, model.name = model.name,
                                    define.constant.x1 , define.constant.x2,
                                    direction.start.x1, direction.end.x1,
                                    direction.start.x2, direction.end.x2  )
    direction.varname <- model.info$var.names[i]
    predicteddata1 <- create.predictions(data = predicteddata, dira_models = dira_models,
                                         model.name = model.name, direction = direction.varname)
    combined.data <- bind_rows(combined.data, predicteddata1)
  }
  combined.data
}


#create.directions.to.plot 2 create base data for stacked data frame with all directions
#    base data is just original data with direction = scatterplot and model = "original_data"
# calls no new functions
base.data <- function(data){
  alldata <- data %>%
    rename(rownum  = X) %>%
    mutate(direction = "scatterplot",
           model = "original_data")
  alldata
}

#create.directions.to.plot.3a: set constant value for marginals
#         called by create.directions.to.plot.3: create.x1x2seq
#     which can be called by any of the plotting commands
create.constant.value <- function(constantVal, data, coefficients, xval, model.name){
  if(constantVal == "mean"){
    constantVal <- mean(data[ ,coefficients[xval]  ], na.rm =T)
  }else if(constantVal == "median"){
    constantVal <- median(data[ ,coefficients[xval]  ], na.rm =T)
  }else{
    constantVal <- as.numeric(constantVal)
  }
  if(model.name %in% c("lm_y.predby.x1.subset.x2.0","lm_y.predby.x1.subset.x2.1")){
    constantVal <- as.numeric(gsub("lm_y.predby.x1.subset.x2." , "", model.name))
  }
  constantVal
}

has.interaction.term <- function(formula){
  terms.factor <- attr(terms(formula), which = "factors")
  interaction.term.exists <- (sum(colSums(terms.factor) > 1) ==1)
  interaction.term.exists
}

#create.directions.to.plot.3: create x1 and x2 for plotting directions of all sorts
#calls no new functions. Very long.
create.x1x2seq <- function(dira_models, data, direction = "regression.surface", model.name = "lm_multiple",
                           define.constant.x1 = "mean", define.constant.x2 = "mean",
                           direction.start.x1 = 0, direction.end.x1 = 1,
                           direction.start.x2 = 1, direction.end.x2 = 0){
  coefficients <-names(dira_models[["lm_multiple"]]$model)
  names(coefficients)  <- c("y", "x1","x2")
  direction.std <- direction
  model.info <- model.information
  model.info <- model.info %>% filter(surface.directions == direction & models == model.name) #for direction description

  if(has.interaction.term(dira_models[["lm_multiple"]]$terms)){
    model.info$direction.descriptions <- model.info$direction.descriptions.interactions
    model.info$direction.descriptions <- gsub("x1constantval", define.constant.x1,
                                              model.info$direction.descriptions)
    model.info$direction.descriptions <- gsub("x2constantval", define.constant.x2,
                                              model.info$direction.descriptions)
  } #make direction description show constant values
  if(model.name %in% c("lm_pc.w.orth", "lm_pc","lm_new.direction")){
    rotation.matrix <- make.rotation.matrix(data, coefficients,  model.name, units.of.x1 = F,
                                            direction.start.x1 = 0, direction.end.x1 = 1,
                                            direction.start.x2 = 1, direction.end.x2 = 0)
  }

  if(direction == "regression.surface"){
    x1.seq = rep(seq(min(data[ ,coefficients["x1"]  ], na.rm =T),
                     max(data[ ,coefficients["x1"]  ], na.rm =T),
                     length.out=30), each =30)
    x2.min = min(data[ ,coefficients["x2"]  ], na.rm =T)
    x2.max = max(data[ ,coefficients["x2"]  ], na.rm =T)
    x2.seq = seq(x2.min, x2.max, length.out = 30)
  }else if(direction == "marginal.x1change.x2const"){
    x1.seq = seq(min(data[ ,coefficients["x1"]  ], na.rm =T),
                 max(data[ ,coefficients["x1"]  ], na.rm =T),
                 length.out=10)
    constantNum  <- create.constant.value(constantVal = define.constant.x2, data = data,
                                          coefficients = coefficients, xval = "x2",
                                          model.name = model.name)  #rep 10
    x2.seq = rep(constantNum, each =10)
  }else if(direction =="marginal.x2change.x1const"){
    constantNum  <- create.constant.value(constantVal = define.constant.x1, data = data,
                                          coefficients = coefficients, xval = "x1",
                                          model.name = model.name)   #rep 10
    x1.seq = rep(constantNum, each =10)
    x2.seq = seq(min(data[ ,coefficients["x2"]  ], na.rm =T),
                 max(data[ ,coefficients["x2"]  ], na.rm =T),
                 length.out=10)
  }else if(direction =="direction.x2.predby.x1"){  #omit x2
    x1.seq = seq(min(data[ ,coefficients["x1"]  ], na.rm =T),
                 max(data[ ,coefficients["x1"]  ], na.rm =T),
                 length.out=10)
    tmp <- data.frame(x1.seq) %>%
      rowid_to_column("rownum")
    names(tmp) <- c("rownum", coefficients["x1"])
    x2.seq <- augment(dira_models[["lm_x2.predby.x1"]], newdata = tmp , interval = "none")$.fitted
  }else if(direction =="direction.x2.predby.x1.resids"){  #residuals direction
    x1.seq = rep(mean(data[ , coefficients["x1"]], na.rm = T), each = 10)
    x2.seq = seq(min(data[ , coefficients["x2"]], na.rm = T),
                 max(data[ , coefficients["x2"]], na.rm = T), length.out=10)
  }else if(direction =="direction.x1.predby.x2"){  #omit x1
    x2.seq = seq(min(data[ ,coefficients["x2"]  ], na.rm =T),
                 max(data[ ,coefficients["x2"]  ], na.rm =T),
                 length.out=10)
    tmp <- data.frame(x2.seq) %>%
      rowid_to_column("rownum")
    names(tmp) <- c("rownum", coefficients["x2"])
    x1.seq <- augment(dira_models[["lm_x1.predby.x2"]], newdata = tmp , interval = "none")$.fitted
  }else if(direction =="direction.x1.predby.x2.resids"){  #residuals direction
    x2.seq = rep(mean(data[ , coefficients["x2"]], na.rm = T), each = 10)
    x1.seq = seq(min(data[ , coefficients["x1"]], na.rm = T),
                 max(data[ , coefficients["x1"]], na.rm = T), length.out=10)
  }else if(direction =="new.direction"){  # user defined direction
    x2.seq = seq(direction.start.x2, direction.end.x2, length.out=10)
    x1.seq = seq(direction.start.x1, direction.end.x1, length.out=10)
    # #the next bits work *if* ensure range of x2 is reasonable for data
    # slope <- (direction.end.x2 - direction.start.x2)/(direction.end.x1 - direction.start.x1)
    # new.dir.rotation.matrix <- matrix(c(1, slope, -slope,1), nrow = 2) #normalized in create.rotated.variables
    # row.names(new.dir.rotation.matrix) <- c(coefficients["x1"], coefficients["x2"])
    # new.dir.rotation.matrix <- adjust.rotation.matrix(new.dir.rotation.matrix)
    # rotation.info <- create.rotation.slopes(rotation.matrix = new.dir.rotation.matrix, data)
    # x1.seq = seq(min(data[ , coefficients["x1"]], na.rm = T),
    #              max(data[ , coefficients["x1"]], na.rm = T), length.out=10)
    # x2.seq = x1.seq*rotation.info$slopes[1]  +rotation.info$intercepts[1]
  }else if(direction =="new.direction.orth"){  # user defined direction
    x1.seq = seq(direction.start.x1, direction.end.x1, length.out=10)
    # #the next bits need a check to see *if* range of x2 is reasonable for data
    rotation.info <- create.rotation.slopes(rotation.matrix = rotation.matrix, data)
    min_x1 <- direction.start.x1
    max_x1 <- direction.end.x1
    x2.seq = x1.seq*rotation.info$slopes[2]  +rotation.info$intercepts[2]
    if(min(x2.seq) < min(data[ ,coefficients["x2"]  ], na.rm =T) & rotation.info$slopes[2] > 0){
      min_x2 <- min(data[ ,coefficients["x2"]  ], na.rm =T)
      min_x1 <- (min_x2 - rotation.info$intercepts[2])/rotation.info$slopes[2]
    }
    if(max(x2.seq) > max(data[ ,coefficients["x2"]  ], na.rm =T) & rotation.info$slopes[2] > 0){
      max_x2 <- max(data[ ,coefficients["x2"]  ], na.rm =T)
      max_x1 <- (max_x2 - rotation.info$intercepts[2])/rotation.info$slopes[2]
    }
    if(min(x2.seq) < min(data[ ,coefficients["x2"]  ], na.rm =T) & rotation.info$slopes[2] < 0){
      min_x2 <- min(data[ ,coefficients["x2"]  ], na.rm =T)
      max_x1 <- (min_x2 - rotation.info$intercepts[2])/rotation.info$slopes[2]
    }
    if(max(x2.seq) > max(data[ ,coefficients["x2"]  ], na.rm =T) & rotation.info$slopes[2] < 0){
      max_x2 <- max(data[ ,coefficients["x2"]  ], na.rm =T)
      min_x1 <- (max_x2 - rotation.info$intercepts[2])/rotation.info$slopes[2]
    }
    x1.seq = seq(min_x1, max_x1, length.out=10)
    x2.seq = x1.seq*rotation.info$slopes[2]  +rotation.info$intercepts[2]
  }else if(direction =="principal.component"){
    pc.data <- subset(data, select = c(coefficients["x1"],coefficients["x2"]))
    x1.seq = seq(min(data[ ,coefficients["x1"]  ], na.rm =T),
                 max(data[ ,coefficients["x1"]  ], na.rm =T),
                 length.out=10)
    rotation.info <- create.rotation.slopes(rotation.matrix, data)
    x2.seq = x1.seq*rotation.info$slopes[1]  +rotation.info$intercepts[1]
    data.frame.to.add.rotated.vars <- data.frame(x1.seq, x2.seq)
    names(data.frame.to.add.rotated.vars) <- c(coefficients["x1"],coefficients["x2"])
    rotated.var <- create.rotated.variables(rotation.matrix = rotation.matrix,
                                            data = data.frame.to.add.rotated.vars,
                                            units.of.x1 = F,
                                            slope = rotation.info$slopes[1],
                                            mean.x1 = mean(pc.data[ , coefficients["x1"]], na.rm = T),
                                            mean.x2 = mean(pc.data[ , coefficients["x2"]], na.rm = T)
                                            )$rotation1
    #orthogonal component was zero for countyData, x1 = anyCollege, x2 = medianIncome1k,
    #         as it should be
    #rotated.var is never changed to units of x1 bc it doesn't affect the plot
  }else if(direction =="principal.component.orth"){  #
    x1.seq = seq(min(data[ ,coefficients["x1"]  ], na.rm =T),
                 max(data[ ,coefficients["x1"]  ], na.rm =T),
                 length.out=10)
    rotation.info <- create.rotation.slopes(rotation.matrix, data)
    x2.seq = x1.seq*rotation.info$slopes[2]  +rotation.info$intercepts[2]
    #would be good to shorten this range?
  }else{
    stop("invalid direction", call. = FALSE)
  }

  direction.name <- gsub("\\n", "\n", model.info$direction.descriptions, fixed = T)
  if(exists("rotated.var")){ #add a rotated var for omitted var models
    xvars.for.prediction <- data.frame(x1.seq, x2.seq, rotated.var,
                                       direction.name , direction.std)
  }else{
    xvars.for.prediction <- data.frame(x1.seq, x2.seq, direction.name ,
                                       direction.std)
  }
  xvars.for.prediction
}


#create slopes in x1 x2 cooridinates for rotated variables
#Use the rotation matrix
#   to calculate the slope of the first and second principal components:
#   beta_rotation = b/a
#   beta_rotation.orth = -a/b
create.rotation.slopes <- function(rotation.matrix, data){
  x1.var <- names(rotation.matrix[, 1])[1]
  x2.var <- names(rotation.matrix[, 1])[2]

  a <- rotation.matrix[x1.var,1]
  b <- rotation.matrix[x2.var,1]
  x1.mean <- mean(data[, x1.var], na.rm =T)
  x2.mean <- mean(data[, x2.var], na.rm =T)
  interceptPC1 <- x2.mean-(b/a)*x1.mean
  interceptPC2 <- (a/b)*x1.mean+x2.mean
  slopes <- c(b/a, -a/b)
  intercepts <- c(interceptPC1, interceptPC2)
  rotation.info <- data.frame(slopes, intercepts, rotationVals = c(a,b))
  rotation.info
}

#create.directions.to.plot.4:
#create predictions from x1.seq, x2.seq with confidence intervals
#requires data created by create.x1x2 function OR a data frame passed in with x1.seq, x2.seq, direction.name
#calls no new functions
create.predictions <- function(data, dira_models,
                               model.name = "lm_multiple", direction = "regression.surface"){
  coefficients <- names(dira_models[["lm_multiple"]]$model)
  model.info <- model.information
  i <- which(model.info$models == model.name & model.info$var.names == direction)
  if(model.info$include.rotated.var.for.predictions[i]) {
    rotated.var.name <- model.info$var.names[i]
    coefficients <- c(coefficients, rotated.var.name)
    names(coefficients) <- c("y", "x1","x2", "rotated.var")
  }else{
    names(coefficients)  <- c("y", "x1","x2")
    #clunky: model must be predicted from x1 x2 unless it's an omitted var model
  }

  #use lm_multiple for predicting y vals for everything
  #except the omitted variable models
  # it's correct bc the predicted vals fall
  #precisely on the regression surface in the direction
  #shown in the x1 x2 sequence
  model.name <- "lm_multiple"
  if(model.info$omitted.var.model[i]) model.name <- model.info$models[i]

  var.names <- unname(coefficients[-1])
  predicteddata <- data %>%
    rowid_to_column("rownum") %>%
    mutate(model = model.name)
  names(predicteddata) <- c("rownum", var.names,
                            "direction", "direction.std", "model")
  predicteddata <- augment(dira_models[[model.name]],
                           newdata = predicteddata , interval = "confidence")
  names(predicteddata) <- c("rownum", var.names, "direction", "direction.std",
                            "model", coefficients["y"], "lowerCI", "upperCI" )
  predicteddata
}



# create 3d plot   #######

# main plotly function
# requires either x, y, z vars (to create z ~ x + y) OR formula
# formula allows interaction terms, but otherwise any transformations must be created before plotting
# models: any model to plot
#ci: if TRUE, plots the confidence intervals associated with each direction/model. Defaults to TRUE
#plot.surface: include regression surface in plot. Defaults to TRUE.
#             can also create regression surface using add_surface()
#define.constant.x1: defaults to "mean". Can also be any value in the range of x1, or "median"
#               Only used for the marginal effects lines
#define.constant.x2: defaults to "mean". Can also be any value in the range of x2, or "median"
#               Only used for the marginal effects lines
#direction.start.x_ : for setting a new direction to plot on the regressions surface
                # Only used for new.directions. New direction moves from
#               (direction.start.x1, direction.start.x2) to (direction.end.x1, direction.end.x2)
#               Most useful for dummy variables (black & female to white & male), so defaults to
#               (0,1) to (1,0)
add_3d.directions <- function(p,  x = NULL, y = NULL, z = NULL, formula = NULL,  data = NULL,
                               model.names = "lm_multiple",
                              ci = T, plot.surface = T,
                              define.constant.x1 = "mean", define.constant.x2 = "mean",
                              direction.start.x1 = 0, direction.end.x1 = 1,
                              direction.start.x2 = 1, direction.end.x2 = 0 ){

  if(is.null(formula)) {
    formula <- create.linear.formula(p, x, y, z)
  }
  if(plot.surface){
    p <- p %>% add_surface(x, y, z, formula,  data, ci)
  }else{
    formulamessage <- paste0("Regression formula is ",
                             Reduce(paste,deparse(formula)), ".")
    message(formulamessage)
  }
  data <- data %||% plotly_data(p, id = names(p$x$visdat)[1])
  #id sets it to the first dataset used in visualization
  if(is.null(data) | length(names(data)) == 0){
    stop("Must supply data in plotly_directions", call. = FALSE)  #probably not necessary
  }

  directionalData <- create.directions.to.plot(formula, data, models = model.names,
                                         define.constant.x1 , define.constant.x2,
                                         direction.start.x1 , direction.end.x1,
                                         direction.start.x2 , direction.end.x2   )
  
  model.info <- model.information
  model.info <- model.info %>% filter( models %in% model.names &
                                        surface.directions != "regression.surface")
  # i <- 5
  for(i in 1:nrow(model.info)){
    if(model.info$omitted.var.model[i]){
      modelType <- model.info$models[i]
    }else{
      # modelname <- "linear regression\n"
      modelType <- "lm_multiple"
    }
    p <- p %>% create.subplotly(data = directionalData, modelType =modelType,
                                directionName = model.info$surface.direction[i], formula, ci = ci,
                                linetype = model.info$linetypes[i], linecolor = model.info$linecolors[i])
  }
  p
}

# calls create.directions.to.plot(), models.dira()
add_surface <- function(p, x = NULL, y = NULL, z = NULL, formula = NULL,  data = NULL, ci = T){
  data <- data %||% plotly_data(p, id = names(p$x$visdat)[1]) #id sets it to the first dataset used in visualization
  if(is.null(data) | length(names(data)) == 0){
    stop("Must supply data in add_surface", call. = FALSE)  #probably not necessary
  }
  if (is.null(formula)) {
    formula <- create.linear.formula(p, x, y, z)
  }
  formulamessage <- paste0("Regression formula is ",
                           Reduce(paste,deparse(formula)), ".")
  message(formulamessage)

  directionalData <- create.directions.to.plot(formula, data, models = c("lm_multiple"))

  dira_models <- models.dira(formula, data)
  coefficients <- names(dira_models[["lm_multiple"]]$model)
  names(coefficients)  <- c("y", "x1","x2")
  dataSubset <- subset(directionalData, direction.std == 'regression.surface')

  p <- p %>%   add_trace(data = dataSubset,
                         x = dataSubset[,coefficients["x1"]],
                         y = dataSubset[,coefficients["x2"]],
                         z = dataSubset[,coefficients["y"]],
                         type ='mesh3d',
                         intensity = dataSubset[,coefficients["y"]],
                         colorscale = list(c(0,'blue'),
                                           c(1, 'blue')),
                         opacity = .5, showscale = F,
                         legendgroup = 'regression.surface',
                         name = "Predicted regression surface")
  if(ci){
    p <- p %>%  add_trace(data = dataSubset,
                          x = dataSubset[,coefficients["x1"]],
                          y = dataSubset[,coefficients["x2"]],
                          z = dataSubset[,"upperCI"],
                          type ='mesh3d',
                          intensity = dataSubset[,"upperCI"],
                          colorscale = list(c(0,'grey'),
                                            c(1, 'grey')),
                          opacity = .5, showscale = F,
                          legendgroup = 'regression.surface', showlegend = F,
                          name = "Lower 95% CI") %>%
      add_trace(data = dataSubset,
                x = dataSubset[,coefficients["x1"]],
                y = dataSubset[,coefficients["x2"]],
                z = dataSubset[,"lowerCI"],
                type ='mesh3d',
                intensity = dataSubset[,"lowerCI"],
                colorscale = list(c(0,'grey'),  #who tf knows why color = doesn't work
                                  c(1, 'grey')),
                opacity = .5, showscale = F,
                legendgroup = 'regression.surface', showlegend = F,
                name = "Upper 95% CI")
  }
  p
}


#calls no new functions
create.subplotly <- function(p, data, directionName , modelType,
                             formula, ci,  linetype ="solid", linecolor = 'grey', ...){
  coefficients <- c(as.character(formula[[2]]),
                    as.character(formula[[3]])[2], as.character(formula[[3]])[3])
  names(coefficients)  <- c("y", "x1","x2")
  dataSubset <- subset(data, model == modelType & direction.std == directionName )

  groupingRandNum <- runif(n =1)  # to make the prediction & CIs group together
  tracename <-  dataSubset$direction[1]
  directionNames <- data.frame(zvals = c(coefficients["y"], "upperCI", "lowerCI"),
                               traceNames = c(tracename, "upper 95% CI", "lower 95% CI"),
                               linewidth = c(6,3,3), legendVal = c(T,F,F),
                               opacity = c(1,.6,.6))
  iterateLength <- 1
  if(ci) iterateLength <- nrow(directionNames)

  for(i in 1:iterateLength){  # loop through directionNames to get CIs if ci == T
    p<- p %>% add_trace( data = dataSubset,
                         x = dataSubset[,coefficients["x1"]],
                         y = dataSubset[,coefficients["x2"]],
                         z = dataSubset[,directionNames$zvals[i]], #coefficients["y"], "upperCI", "lowerCI"
                         type="scatter3d", mode='lines' , opacity = directionNames$opacity[i],
                         line = list(dash = linetype, width = directionNames$linewidth[i], color = linecolor),
                         legendgroup = groupingRandNum, showlegend = directionNames$legendVal[i],
                         name = directionNames$traceNames[i])
  }
  p
}





#creates a linear formula from passing through p and/or x, y, z
create.linear.formula <- function(p, x=NULL, y = NULL, z = NULL){
  x <-  x %||% p$x$attrs[[1]][["x"]][[2]]
  y <-  y %||% p$x$attrs[[1]][["y"]][[2]]
  z <-  z %||% p$x$attrs[[1]][["z"]][[2]]
  if(is.null(x) | is.null(y) | is.null(z)){
    stop("Must supply or pass through `x`, `y`, and `z` attributes", call. = FALSE)
  }
  if(class(x) == "formula"){
    x <- as.character(x)[2]
    y <- as.character(y)[2]
    z <- as.character(z)[2]
  }
  linearformula <- formula(paste(c(z, "~", x,"+", y), collapse = " " ))
  linearformula
}


# create 2d arranged plots ########

#visualizes rotated and/or nested variable models in 2d, stacked next to each other
#Only supports the models y~x1, y ~ x1 + x2.resids, y ~ x2, y ~ x2 + x1.resids
#p: a plotly object
#     Note: to use layout options, the layout command *must* precede this command
#model.names the model names for the rotated or nested models
#function: a function of the form y ~ x1 + x2. Only first two terms will be used
#graph.title: optional; creates a title on top of the arranged graphs
#value: returns an object of class ggarrange, which is a list of ggplot
add_2d.directions <- function( p, model.names, formula = NULL,
                               x = NULL, y = NULL, z = NULL,
                               data = NULL,
                               define.constant.x1 = "mean", define.constant.x2= "mean",
                               direction.start.x1=0 , direction.end.x1= 1,
                               direction.start.x2=1 , direction.end.x2 = 0,ci = T, categories.3 = F,
                               shareX = F, shareY = F, titleX = T, titleY= T){

  plist <- plotlist.2d(p, model.names = model.names, formula,
                              x, y, z,  data,
                              define.constant.x1 , define.constant.x2,
                              direction.start.x1, direction.end.x1,
                              direction.start.x2 , direction.end.x2,
                              ci, categories.3 )
  if(length(plist)== 1){
    returnplot <- plist[[1]]
  }else if(length(plist)== 2){
    topplot <- subplot(plist[[1]], plotly_empty(type = "scatter", mode = "markers"),
                       plist[[2]], widths = c(.46, .08, .46),
                       titleX = titleX, titleY = titleY, shareY = shareY)
    returnplot <- subplot(plotly_empty(type = "scatter", mode = "markers"),
                          topplot, plotly_empty(type = "scatter", mode = "markers"),
                          nrows = 3, heights = c(.03, .92, .05), titleX = titleX, titleY = titleY, shareY =shareY) %>%
      layout(legend = list(orientation = 'h'))
  }else if(length(plist) == 3){
    topplot <- subplot(plist[[1]], plotly_empty(type = "scatter", mode = "markers"),
                       plist[[3]], widths = c(.46, .08, .46),
                       titleX = titleX, titleY = titleY, shareY = shareY)
    midplot <- subplot(plotly_empty(type = "scatter", mode = "markers"),
                       plist[[2]], plotly_empty(type = "scatter", mode = "markers"),
                       widths = c(.25, .5, .25),
                       titleX = titleX, titleY = titleY, shareY = shareY)
    returnplot <- subplot(topplot,  plotly_empty(type = "scatter", mode = "markers"),
                          midplot, plotly_empty(type = "scatter", mode = "markers"),
                          nrows = 4, heights = c(.44, .08, .44, .04),
            titleX = titleX, titleY = titleY, shareY = shareY) %>%
      layout(legend = list(orientation = 'h'))
  }else{
    warning("model list creates more than 3 plots. To see ")
  }

    returnplot
}



plotlist.2d <- function(p, model.names, formula = NULL,
                               x = NULL, y = NULL, z = NULL,
                               data = NULL,
                               define.constant.x1 = "mean", define.constant.x2= "mean",
                               direction.start.x1=0 , direction.end.x1= 1,
                               direction.start.x2=1 , direction.end.x2 = 0,ci = T, categories.3 = F){

  plot.id = names(p$x$visdat)[1]
  data <- data %||% plotly_data(p, id = plot.id)

  if(is.null(formula)) {
    formula <- create.linear.formula(p, x, y, z)
  }
  formulamessage <- paste0("Regression formula is ",
                           Reduce(paste,deparse(formula)), ".")
  message(formulamessage)

  layout_attrs <- p$x$layoutAttrs[[p$x$cur_data]] # or plot.id?
  plot_attrs <- p$x$attrs[[p$x$cur_data]]
  if(is.null(plot_attrs$color))  plot_attrs[["colors"]] <- I("black")
  plist <- list()

  #isolate only models that can be used for 2d plot
  model.info <- model.information[which(model.information$overlayed.model != "none"), ]
  if(!categories.3){
    model.info <- model.info[which(model.info$models != "lm_multiple"), ]
  }
  if(!categories.3 & ("lm_new.direction" %in% model.names) & ("lm_y.predby.x1.subset.x1ORx2.1" %in% model.names)){
    model.info <- model.info[which(model.info$models != "lm_new.direction"), ] #don't plot new direction twice
  }

  model.names.add <- model.names
  if(sum(c("lm_y.predby.x1.w.resids", "lm_y.predby.x1") %in% model.names) > 0){
    model.names.add <- unique(c(model.names, "lm_y.predby.x1.w.resids", "lm_y.predby.x1"))
  }  #add both models if at least one exists, so the surface plot is plotted
  if(sum(c("lm_y.predby.x2.w.resids", "lm_y.predby.x2") %in% model.names) > 0){
    model.names.add <- unique(c(model.names.add, "lm_y.predby.x2.w.resids", "lm_y.predby.x2"))
  }
  model.info <- model.info[model.info$models %in% model.names.add, ]
  
  for(i in 1:nrow(model.info)){
    model.name <- model.info$models[i]
    ci.current <- ci
    plistname <- as.character(i) #only changes for lm_y.predby.xi and lm_y.predby.xi.w.resids
    # model.info.tmp <- model.info %>% filter((models == model.name))
    # model.info.tmp <- model.info.tmp[1, ] #never plot orthogonal?
    showlegend.val <- model.info$show.legend.2d[i]
    yvar <- plot_attrs$z %||% z
    yvar.label <- layout_attrs$scene$zaxis$title %||% as.character(yvar)[2]
    if(model.info$xvar.2d[i] == "x1") {
      x1var <- plot_attrs$x %||% x
      x1var.label <- layout_attrs$scene$xaxis$title %||% as.character(x1var)[2]
      x2var <- plot_attrs$y %||% y
      x2var.label <- layout_attrs$scene$yaxis$title %||% as.character(x2var)[2]
      default.constant <- as.numeric(model.info$default.x2.const[i])
    }else if(model.info$xvar.2d[i] == "x2"){
      x1var <- plot_attrs$y %||% y
      x1var.label <- layout_attrs$scene$yaxis$title %||% as.character(xvar)[2]
      x2var <- plot_attrs$x %||% x
      x2var.label <- layout_attrs$scene$xaxis$title %||% as.character(xi.predby.xj)[2]
      default.constant <- as.numeric(model.info$default.x1.const[i])
    } #set x1 and x2 axes

    default.constantx1 <- model.info$default.x1.const[i]
    default.constantx2 <- model.info$default.x2.const[i]
        # deal with special cases: subsets & regression surface plot
     if(grepl(".subset.", model.name) | categories.3) { #if we're working with subset
       if(model.name == "lm_y.predby.x1.subset.x1ORx2.1" | model.name == "lm_new.direction"){
         datatmp <- data[ which(data[ , as.character(x2var)[2]] == 1 | data[ , as.character(x1var)[2]] == 1), ]
         }else if(categories.3){
           default.constantx1 = 0
           default.constantx2 <- 0
           datatmp <- data[ which(data[ , as.character(x2var)[2]] == default.constantx1), ]
         }else{
           datatmp <- data[ which(data[ , as.character(x2var)[2]] == default.constant), ]
       }
      originalcolorlevels <- levels(as.factor((data[ , as.character(plot_attrs$color)[2]])))
      tmpcolorlevels <- as.vector(levels(as.factor((datatmp[ , as.character(plot_attrs$color)[2]]))))
      plotcolors <- plot_attrs$colors[which(originalcolorlevels %in% tmpcolorlevels)]
    }else{
      datatmp <- data
      plotcolors <- plot_attrs$colors
    } #create subset if in subset
    make.surface <- gsub("x1", model.info$xvar.2d[i],c("lm_y.predby.x1.w.resids", "lm_y.predby.x1" ))
    if(sum( (make.surface %in% model.names)) ==1){  #only deals with one omitted var model
      replacemodel <- make.surface[(make.surface %in% model.names)]
    }else{
      replacemodel <- model.name
    }  #dealing with creating regression surface plot
    if(model.name == make.surface[2]) {
      yvar <- x2var
      yvar.label <- x2var.label
      ci.current <- F
      plistname <- paste0(replacemodel, ".surface")
    }
    #baseline plot
    p <- plot_ly(data = datatmp,  x = x1var, y = yvar, text = plot_attrs$text,
                 color = plot_attrs$color, colors = plotcolors,
                 size = plot_attrs[["size"]], type = "scatter", mode = "markers",  showlegend = F) %>%
      add_markers(opacity = .6, showlegend = showlegend.val)

    overlaid.models <- c(model.name, model.info$overlayed.model[i])
    overlaid.models <- overlaid.models[overlaid.models %in% model.names]
        #exclude any models not in original list
    directionalData <- create.directions.to.plot(formula, data, models = overlaid.models,
                                                 define.constant.x1 = default.constantx1,
                                                 define.constant.x2 = default.constantx2,
                                                 direction.start.x1 , direction.end.x1, direction.start.x2 , direction.end.x2,
                                                 units.of.x1 = F)

    overlaid.models <- unique(model.information$models[which( model.information$models %in% overlaid.models)])
    #reorder so dash shows up on top
    for(overlaid.model in overlaid.models){  #add predicted lines for either or both nested and rotated models
      subplot.dir <- model.info$surface.directions[i] 
      model.info.current <- model.information %>% filter((models == overlaid.model & surface.directions == subplot.dir))
      data_1dir <- directionalData %>%
        filter(model == model.info.current$estimating.model & direction.std == subplot.dir)
      if(has.interaction.term(formula) & (overlaid.model == "lm_multiple")){
        showlegend.valtmp <- model.info$legend.model1.2d[i] #which(overlaid.model %in% overlaid.models)
        model.info.current$direction.descriptions.2d <- model.info.current$direction.descriptions.interactions.2d
        model.info.current$direction.descriptions <- model.info.current$direction.descriptions.interactions.2d
        # print(overlaid.model)
        }else{showlegend.valtmp <- showlegend.val} #adjust so interaction terms don't show lm_multiple legend twice
      tracename <- gsub("\\n","\n", model.info.current$direction.descriptions, fixed = T)
      grouping.string <- paste0( model.info.current$direction.description.2d)  # to make the prediction & CIs group together
      p <- p %>% add_trace( data = data_1dir, x = data_1dir[,as.character(x1var)[2]],
                       y = data_1dir[,as.character(yvar)[2]],
                       type = "scatter", mode='lines' ,
                       line = list(color = model.info.current$linecolors,
                                   dash = model.info.current$linetypes, width = 1.5),
                       legendgroup = grouping.string, showlegend = showlegend.valtmp,
                       name =tracename) #add predicted value line
      if(ci.current){ #add CIs for both nested and rotate models
        p <- p %>% add_ribbons(data = data_1dir, x = data_1dir[,as.character(x1var)[2]],
                               ymin = ~lowerCI,  ymax = ~upperCI, fillcolor = "grey60",
                               opacity = .1, legendgroup = grouping.string,showlegend = F, type = "scatter")
        yvar.ci <- "lowerCI"
        for(j in 1:2){ #add thin lines for upper and lower CIs
          p <- p %>% add_trace(data = data_1dir, x = data_1dir[,as.character(x1var)[2]],
                                 y = data_1dir[,yvar.ci], type = "scatter", mode='lines' ,
                               line = list(color = model.info.current$linecolors,
                                           dash = model.info.current$linetypes, width = .5,
                                           opacity = .3), legendgroup = grouping.string,
                                 showlegend = F,  name = "confidence interval")
          yvar.ci <- "upperCI"
        } #close out ci lines
       }#close out if(ci)
    } #add predicted lines for either or both nested and rotate models
    p <- p %>%
      layout(xaxis = list(title = x1var.label, showline = T),
             yaxis = list(title = yvar.label, showline = T))
    plist[[plistname]] <- p
  } #create each item in the list of plots


  plist
}


#creates the omitted variable models in 2d, stacked next to each other
#function: a function of the form y ~ x1 + x2. Only first two terms will be used
#data: the dataset with the variables specified in formula
#p: a ggplot object
#     optional; if the user passes in an existing ggplot object with
#     form ggplot(data = data, aes(x = x1var, y = x2var, z = yvar )) + geom_point
#     any colors or labels defined in that plot will be applied to
#     the output of gg_2d.directions.
#title: optional; creates a title on top of the arranged graphs
#value: returns an object of class ggarrange, which is a list of ggplot

gg_2d.directions <- function(data, formula, p = NULL, graph.title = NULL){
  x1.var <- attr(terms(formula), which = "term.labels")[1]
  x2.var <- attr(terms(formula), which = "term.labels")[2]
  y.var <- as.character(formula[[2]])[1]

  if(is.null(p)){
    p <- ggplot(data = data, aes_(x = as.name(x1.var),
                                  y = as.name(x2.var),
                                  z = as.name(y.var))) +
      geom_point()
  }
  labels <- p$labels #extract labels
  g <- ggplot_build(p)  #need g to grab color column
  data$colorcolumn <- as.vector(g$data[[1]]["colour"])[,1]

  p1 <- ggplot(data = data, aes_(x = as.name(x1.var), y =as.name(y.var))) +
    geom_point(colour = data$colorcolumn) + geom_smooth(method = lm) +
    scale_colour_identity()+
    labs(x = labels$x[1],
         y = labels$z[1])+
    theme_classic()

  p2 <- ggplot(data = data, aes_(x = as.name(x2.var), y =as.name(y.var))) +
    geom_point(colour = data$colorcolumn) + geom_smooth(method = lm) +
    labs(x = labels$y[1],
         y = labels$z[1])+
    theme_classic()

  arranged <- ggarrange(p1, p2 + rremove('ylab'), common.legend = T , legend = 'bottom')
  if(!is.null(graph.title)){
    annotate_figure(arranged,
                    top = graph.title)
  }else{
    arranged
  }
}



# specific 3d directions ####

# returns plotly object with omitted variables directions plotted
#calls create.directions.to.plot, models.dira(), add_surface(), create.subplotly()
add_directionOmitVar <- function(p, x = NULL, y = NULL, z = NULL, formula = NULL,  data = NULL,
                                 omittedVarModels = c("lm_y.predby.x1", "lm_y.predby.x1.w.resids",
                                                      "lm_y.predby.x2", "lm_y.predby.x2.w.resids"),
                                  ci = T, plot.surface = T){
  p <- p %>% add_3d.directions(x, y , z, formula,  data,
                          model.names = omittedVarModels,
                           ci, plot.surface )
  p
}

# returns plotly object with new direction plotted
add_newDirection <- function(p, x = NULL, y = NULL, z = NULL, formula = NULL,  data = NULL,
                             newDirections = c("lm_new.direction"),
                             direction.start.x1 = 0, direction.end.x1 = 1,
                             direction.start.x2 = 1, direction.end.x2 = 0,
                             ci = T, plot.surface = T){

  p <- p %>% add_3d.directions(x, y , z, formula,  data,
                                model.names = newDirections, ci = ci,
                                plot.surface = plot.surface,
                                define.constant.x1 = "mean", define.constant.x2 = "mean", #not needed here except to ensure pass through
                                direction.start.x1 , direction.end.x1,
                                direction.start.x2, direction.end.x2)
  p
}



#calls create.directions.to.plot, models.dira(), create.subplotly().
#Theoretically could move create.directions.to.plot and models.dira() into create.subplotly
add_marginals <- function(p, x = NULL, y = NULL, z = NULL, formula = NULL,  data = NULL,
                          define.constant.x1 = "mean", define.constant.x2 = "mean",
                          ci = T, plot.surface = T){

  p <- p %>% add_3d.directions(x, y , z, formula,  data,
                               model.names = c("lm_multiple"), ci = ci,
                                plot.surface = plot.surface,
                                define.constant.x1, define.constant.x2)

  p
}
