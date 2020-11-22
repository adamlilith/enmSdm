.makeLarsData <- .makeLarsData <- function (data, resp, preds, scale = TRUE, quadratic = TRUE, 
     cubic = TRUE, interaction = TRUE, interQuad = TRUE, na.rm = FALSE) {
    
	if (cubic & !quadratic) {
         quadratic <- TRUE
         warning("Forcing quadratic features because linear-by-quadratic features are requested.", 
             .immediate = TRUE)
     }
     if (!is.null(resp)) {
         if (class(resp) %in% c("integer", "numeric")) 
             resp <- names(data)[resp]
         y <- data[, resp]
     }
     else {
         y <- NULL
     }
     if (class(preds) %in% c("integer", "numeric")) 
         preds <- names(data)[preds]
     data <- data[, preds, drop = FALSE]
     if (na.rm) {
         yNa <- if (is.null(y)) {
             NULL
         }
         else {
             which(is.na(y))
         }
         predsNa <- omnibus::naRows(data)
         if (length(yNa) > 0) {
             y <- y[-yNa]
             data <- data[-yNa, ]
         }
         if (length(predsNa) > 0) {
             y <- y[-predsNa]
             data <- data[-predsNa, ]
         }
     }
     scales <- list()
     if (class(scale) == "list") {
         data <- data[, names(scale$`scaled:center`)]
         data <- as.data.frame(base::scale(data, center = scale$`scaled:center`, 
             scale = scale$`scaled:scale`))
         scales$`scaled:center` <- scale$`scaled:center`
         scales$`scaled:scale` <- scale$`scaled:scale`
     }
     else if (class(scale) == "logical") {
         if (scale) {
             data <- base::scale(data)
             scales$`scaled:center` <- attributes(data)$`scaled:center`
             scales$`scaled:scale` <- attributes(data)$`scaled:scale`
             data <- as.data.frame(data)
             if (any(scales$`scaled:scale` == 0)) {
                 zeroVar <- names(scales$`scaled:scale`[scales$`scaled:scale` == 
                   0])
                 scales$`scaled:center` <- scales$`scaled:center`[-which(names(scales$`scaled:center`) %in% 
                   zeroVar), drop = FALSE]
                 scales$`scaled:scale` <- scales$`scaled:scale`[-which(names(scales$`scaled:scale`) %in% 
                   zeroVar), drop = FALSE]
                 data <- data[, -which(names(data) %in% zeroVar), 
                   drop = FALSE]
                 preds <- preds[-which(preds %in% zeroVar)]
                 if (ncol(data) == 1 & interaction) 
                   interaction <- FALSE
                 if (ncol(data) == 1 & interQuad) 
                   interQuad <- FALSE
                 warning(paste0("Removing variables from predictor set because they have 0 variance: ", 
                   paste(zeroVar, collapse = " ")))
             }
         }
         else {
             scales$`scaled:center` <- rep(0, ncol(data))
             scales$`scaled:scale` <- rep(1, ncol(data))
             names(scales$`scaled:center`) <- names(scales$`scaled:scale`) <- names(data)
         }
     }
     inCol <- matrix(FALSE, ncol = length(preds), nrow = length(preds))
     rownames(inCol) <- 1:nrow(inCol)
     colnames(inCol) <- preds
     diag(inCol) <- TRUE
     groups <- list()
     for (i in seq_along(preds)) groups[[i]] <- preds[i]
     if (quadratic) {
         for (thisPred in preds) {
             add <- as.data.frame(data[, thisPred]^2)
             newPred <- paste0(thisPred, "_pow2")
             names(add) <- newPred
             data <- cbind(data, add)
             thisInCol <- matrix(preds %in% thisPred, ncol = length(preds))
             rownames(thisInCol) <- ncol(data)
             inCol <- rbind(inCol, thisInCol)
             vars <- c(thisPred, newPred)
             groups[[length(groups) + 1]] <- vars
         }
     }
     if (cubic) {
         for (thisPred in preds) {
             add <- as.data.frame(data[, thisPred]^3)
             newPred <- paste0(thisPred, "_pow3")
             names(add) <- newPred
             data <- cbind(data, add)
             thisInCol <- matrix(preds %in% thisPred, ncol = length(preds))
             rownames(thisInCol) <- ncol(data)
             inCol <- rbind(inCol, thisInCol)
             vars <- c(thisPred, paste0(thisPred, "_pow2"), 
                 newPred)
             groups[[length(groups) + 1]] <- vars
         }
     }
     if (interaction) {
         for (countPred1 in 1:(length(preds) - 1)) {
             for (countPred2 in (countPred1 + 1):length(preds)) {
                 thisPred <- preds[countPred1]
                 thatPred <- preds[countPred2]
                 newPred <- paste0(thisPred, "_by_", thatPred)
                 add <- as.data.frame(data[, thisPred] * data[, 
                   thatPred])
                 names(add) <- newPred
                 data <- cbind(data, add)
                 thisInCol <- matrix(preds %in% thisPred | preds %in% 
                   thatPred, ncol = length(preds))
                 rownames(thisInCol) <- ncol(data)
                 inCol <- rbind(inCol, thisInCol)
                 vars <- c(thisPred, thatPred, newPred)
                 groups[[length(groups) + 1]] <- vars
             }
         }
     }
     if (interQuad) {
         for (thisPred in preds) {
             for (thatPred in preds[!(preds %in% thisPred)]) {
                 add <- as.data.frame(data[, thisPred] * (data[, 
                   thatPred])^2)
                 names(add) <- paste0(thisPred, "_by_", 
                   thatPred, "_pow2")
                 data <- cbind(data, add)
                 thisInCol <- matrix(preds %in% thisPred | preds %in% 
                   thatPred, ncol = length(preds))
                 rownames(thisInCol) <- ncol(data)
                 inCol <- rbind(inCol, thisInCol)
                 vars <- c(thisPred, thatPred, paste0(thisPred, 
                   "_by_", thatPred), paste0(thatPred, "_by_", 
                   thisPred), paste0(thatPred, "_pow2"), 
                   paste0(thisPred, "_by_", thatPred, "_pow2"))
                 vars <- names(data)[which(names(data) %in% vars)]
                 groups[[length(groups) + 1]] <- vars
             }
         }
     }
     out <- list()
     out$resp <- resp
     out$preds <- preds
     out$y <- y
     out$data <- data
     out$scales <- scales
     out$groups <- groups
     out$inCol <- inCol
     out$features$quadratic <- quadratic
     out$features$cubic <- cubic
     out$features$interaction <- interaction
     out$features$interQuad <- interQuad
     class(out) <- c(class(out), "larsData")
     out
 }
