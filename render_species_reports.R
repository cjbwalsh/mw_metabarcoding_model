library(quarto)

# 1. Define your loop target
species_vector <- 1:33

# 2. Build the structural YAML header for your combined document
master_lines <- c(
  "---",
  "title: ''",
  "format:",
  "  docx:",
  "    reference-doc: report_style.docx",
  "---",
  "",
  "```{r}",
  "#| echo: false",
  "#| message: false",
  "#| warning: false",
  
  "source('working_functions.R')",
  "fit <- load_model_component_bundle('ifqdtbrs_fi_qf_qd_NO_ac')$fit",
  
  paste0("current_species_list <- c(", paste(species_vector, collapse = ","), ")"),
  "```",
  ""
)

# 3. Use Quarto shortcodes to inject the template sequentially
for (i in seq_along(species_vector)) {
  master_lines <- c(
    master_lines,
    "```{r}",
    "#| echo: false",
    paste0("current_idx <- ", i),
    "```",
    "",
    "{{< include species_report_template.qmd >}}",
    ""
  )
}

# 4. Save this out as a temporary master document file
writeLines(master_lines, "auto_generated_master.qmd")

# 5. Execute the Quarto engine directly (bypasses the knitr engine loop bug)
quarto_render("auto_generated_master.qmd", output_file = "final_species_report.html")

# Optional: clean up the generated file
file.remove("auto_generated_master.qmd")
