pivot_mapping_library <- function(mappings) {
  out <- mappings %>%
    tidyr::pivot_longer(
      cols = !`Harmonized name`,
      names_to = c(".value", "Origin file"),
      names_sep = " \\(",
      values_drop_na = TRUE
    ) %>%
    dplyr::filter(`Input name` != "-") %>%
    dplyr::mutate(`Origin file` = gsub("\\)", "", `Origin file`))
  return(out)
}

extract_missing_values <- function(appended_inputs, myinputfiles) {
  missed_mappings <- mapply(function(x, y) {
    y <- as.vector(y) %>% unlist
    names(y) <- colnames(myinputfiles)[-1]
    rows <- x %>% dplyr::filter(is.na(`Standardized name`))

    identifiers <- c(y["HMDB"], y["PubChem_CID"], y["KEGG"], y["LIPIDMAPS"], y["chebi"], y["Metabolite_Name"])
    out <- lapply(identifiers, function(z) {
      if (is.na(z)) {
        return(data.frame(rep(NA, nrow(rows))))
      } else{
        return(rows[, z])
      }
    })
    out <- do.call(cbind, out)
    colnames(out) <- c("HMDB",
                       "PubChem_CID",
                       "KEGG",
                       "LIPIDMAPS",
                       "chebi",
                       "Metabolite_Name")
    return(out)
  },
  x = appended_inputs,
  y = as.data.frame(t(myinputfiles[-1])),
  SIMPLIFY = FALSE)
  missed_mappings <- as.vector(missed_mappings)
  names(missed_mappings) <- myinputfiles$ShortFileName
  return(missed_mappings)
}

find_multimapped_metabolites <- function(mapping_library, myinputfiles) {
  inputs <- mapping_library %>%
    dplyr::select(`Harmonized name`, dplyr::starts_with("Input name"))
  out <- matrix(ncol = 3, nrow = 0)
  for (i in 2:ncol(inputs)) {
    if (any(duplicated(inputs[, i]))) {
      duplicates <- inputs[which(duplicated(inputs[, i])), c(1, i)]
      duplicates <- duplicates[!apply(duplicates == "-", 1, any), ]
      duplicates <- cbind(rep(myinputfiles$ShortFileName[i], nrow(duplicates)), duplicates)
      out <- rbind(out, duplicates)
    }
  }
  colnames(out) <- c("Input File", "Standard Name", "Input ID")
  return(out)
}

sink.reset <- function() {
  for (i in seq_len(sink.number())) {
    sink(NULL)
  }
}

write_txt_log <- function(start_time, myinputfiles) {
  sink("metLinkR_output/metLinkR_log.txt")
  cat(paste0("Run date/time:", Sys.time()), "\n")
  cat(paste0("Runtime: ", Sys.time() - start_time, " mins"), "\n")
  cat("Input file:", "\n")
  print(myinputfiles)
  cat("\n")
  cat("Session Info:", "\n")
  print(utils::sessionInfo())
  sink()
  sink.reset()
}

plot_mapping_rates <- function(mapping_rates) {
  mapping_rates <- unlist(mapping_rates)
  names(mapping_rates)[1] <- "Global"
  mapping_rates <- data.frame(mapping_rates)
  mapping_rates$dataset <- rownames(mapping_rates)
  mapping_rates$dataset <- factor(mapping_rates$dataset, levels = mapping_rates$dataset)
  colors <- c("1", rep("2", times = nrow(mapping_rates) - 1))
  p <- ggplot2::ggplot(mapping_rates,
                       ggplot2::aes(x = dataset, y = mapping_rates, fill = colors)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = c("goldenrod", "grey40")) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Dataset", y = "Mapping Rate") +
    ggplot2::guides(fill = "none")
  return(p)
}

plot_chemical_classes <- function(mapped_list_input_files,
                                  mapped_list_synonyms) {
  mapped_list_input_files <- lapply(mapped_list_input_files, function(x) {
    x %>%
      dplyr::group_by(rownum) %>%
      dplyr::filter(priority == min(priority)) %>%
      as.data.frame
  })
  refmet_classes <- lapply(mapped_list_input_files, function(x) {
    out <- x$`Super class`[which(!is.na(x$`Super class`))]
    out <- out[-which(out == "-")]
    out <- out[-which(out == "")]
    return(out)
  })
  synonym_classes <- lapply(mapped_list_synonyms, function(x) {
    if (methods::is(x, "data.frame")) {
      temp <- x %>%
        dplyr::select("Standardized name", "Super class") %>%
        dplyr::filter(`Standardized name` != "-") %>%
        unique %>%
        dplyr::pull("Super class")
    }
  })
  classes <- mapply(function(x, y) {
    return(c(x, y))
  }, x = refmet_classes, y = synonym_classes)

  plot_list <- mapply(function(x, y) {
    class_table <- as.data.frame(table(x))
    p <- ggplot2::ggplot(class_table, ggplot2::aes(x = x, y = Freq)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "ClassyFire SuperClass", y = "Count") +
      ggplot2::ggtitle(paste0(y, " Mapped Chemical Classes")) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0),
        plot.margin = margin(1, 2, 1.5, 1.2, "cm")
      )
    return(p)
  },
  x = classes,
  y = names(classes),
  SIMPLIFY = FALSE)
  return(plot_list)
}

meta_plotting <- function(chemical_plots) {
  splitted_plots <- split(chemical_plots, ceiling(seq_along(chemical_plots) /
                                                    3))
  ## if(length(splitted_plots[[length(splitted_plots)]])!=3){
  ##   n_empty <- 3 - length(splitted_plots[[length(splitted_plots)]])
  ##   splitted_plots[[length(splitted_plots)]] =
  ##     append(splitted_plots[[length(splitted_plots)]],
  ##            rep(ggplot() + theme_void(),times = n_empty))
  ## }
  lapply(splitted_plots, function(x)
    cowplot::plot_grid(plotlist = x, ncol = 1))

}

write_id_rates <- function(mapped_list_input_files,
                           mapped_list_synonyms) {
  mapped_list_input_files <- lapply(mapped_list_input_files, function(x) {
    x %>%
      dplyr::group_by(rownum) %>%
      dplyr::filter(priority == min(priority)) %>%
      as.data.frame
  })
  id_origins <- lapply(mapped_list_input_files, function(x) {
    return(x$origin[which(!is.na(x$origin))])
  })
  id_origins <- mapply(function(x, y) {
    if (methods::is(y, "data.frame")) {
      out <- c(x, rep("RaMP", times = length(unique(
        y$`Standardized name`
      )) - 1))
      return(out)
    } else{
      return(x)
    }
  }, x = id_origins, y = mapped_list_synonyms)
  table_list <- lapply(id_origins, function(x) {
    temp <- table(x)
    return(kableExtra::kable(temp, col.names = c("Identifier Type", "Count")))
  })
  return(table_list)
}


write_pdf_report <- function(mapping_rates,
                             mapped_list_input_files,
                             mapped_list_synonyms) {
  cat(
    "---
title: \"MetLinkR Report\"
author:
date: \"`r format(Sys.time(), '%d %B, %Y')`\"
output: pdf_document
---

\`\`\`{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, fig.height = 3)
\`\`\`

## Metabolite Mapping Rates
\`\`\`{r}
plot_mapping_rates(mapping_rates)
\`\`\`

## Identifier Types Used by File
\`\`\`{r, results = \'asis\'}
write_id_rates(mapped_list_input_files,mapped_list_synonyms)
\`\`\`

## Chemical Class Breakdown by File
\`\`\`{r,fig.height=12,warnings = FALSE, message = FALSE, fig.width = 12}
suppressWarnings({
meta_plotting(plot_chemical_classes(mapped_list_input_files,mapped_list_synonyms))
})
\`\`\`

",
    file = "metLinkR_output/metLinkR_report.Rmd"
  )
  rmarkdown::render("metLinkR_output/metLinkR_report.Rmd")
  if (file.exists("metLinkR_output/metLinkR_report.Rmd")) {
    #Delete file if it exists
    file.remove("metLinkR_output/metLinkR_report.Rmd")
  }
}
