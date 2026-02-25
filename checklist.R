# Dependencies
library(httr2)
library(dplyr)
library(jsonlite)

# Helper functions
get_property <- function(x, name) {
  if (hasName(x, name)) getElement(x, name)[[1]] else NA_character_
}

get_parent_id <- function(x) {
  for (i in rev(x$parents)) {
    if (isTRUE(taxonRank[taxonConceptID == i] %in% ranks)) break
  }

  if (length(i) < 1 || i == "MX.37600") {
    NA_character_
  } else {
    paste0("http://tun.fi/", i)
  }
}

get_invalid_names <- function(concept, type, status) {
  if (hasName(concept, type)) {
    lapply(
      concept[[type]],
      function(x, status, concept_id) {
        data.frame(
          taxonID = paste0("http://tun.fi/", x$id),
          taxonConceptID = paste0("http://tun.fi/", concept_id),
          taxonomicStatus = status,
          taxonRank = sub("MX\\.", "", x$taxonRank %||% NA),
          scientificName = x$scientificName %||% NA,
          scientificNameAuthorship = x$scientificNameAuthorship %||% NA,
          genericName = NA,
          infragenericEpithet = NA,
          specificEpithet = NA,
          infraspecificEpithet = NA,
          acceptedNameUsageID = paste0("http://tun.fi/", concept_id),
          acceptedNameUsage = NA,
          parentNameUsageID = concept$parent_id,
          parentNameUsage = NA
        )
      },
      status,
      concept$id
    )
  } else {
    list(NULL)
  }
}

flatten_concept <- function(x) {
  x$parent_id <- get_parent_id(x)

  do.call(
    rbind,
    c(
      list(
        data.frame(
          taxonID = paste0("http://tun.fi/", x$id),
          taxonConceptID = paste0("http://tun.fi/", x$id),
          taxonomicStatus = "accepted",
          taxonRank = sub("MX\\.", "", x$taxonRank %||% NA),
          scientificName = x$scientificName %||% NA,
          scientificNameAuthorship = x$scientificNameAuthorship %||% NA,
          genericName = NA,
          infragenericEpithet = NA,
          specificEpithet = NA,
          infraspecificEpithet = NA,
          acceptedNameUsageID = paste0("http://tun.fi/", x$id),
          acceptedNameUsage = NA,
          parentNameUsageID = x$parent_id,
          parentNameUsage = NA
        )
      ),
      get_invalid_names(x, "synonyms", "synonym"),
      get_invalid_names(x, "subjectiveSynonyms", "heterotypicSynonym"),
      get_invalid_names(x, "heterotypicSynonyms", "heterotypicSynonym"),
      get_invalid_names(x, "objectiveSynonyms", "homotypicSynonym"),
      get_invalid_names(x, "homotypicSynonyms", "homotypicSynonym")
    )
  )
}

has_species <- function(id) {
  children <- NameUsage[NameUsage$parentNameUsageID %in% id, ]
  if (
    any(children$taxonRank %in% c("species", "subspecies", "form", "variety"))
  ) {
    TRUE
  } else if (identical(nrow(children), 0L)) {
    FALSE
  } else {
    has_species(children$taxonID)
  }
}

ranks <- c(
  "domain",
  "kingdom",
  "subkingdom",
  "infrakingdom",
  "superphylum",
  "phylum",
  "subphylum",
  "superclass",
  "class",
  "subclass",
  "infraclass",
  "parvclass",
  "superorder",
  "order",
  "suborder",
  "infraorder",
  "parvorder",
  "superfamily",
  "family",
  "subfamily",
  "tribe",
  "subtribe",
  "genus",
  "subgenus",
  "section",
  "aggregate",
  "species",
  "subspecies",
  "form",
  "variety"
)

# Get & cache Finnish taxonomy
if (!file.exists("taxonomy.rds")) {
  req <-
    request("https://api.laji.fi") |>
    req_url_path("v0") |>
    req_url_path_append("taxa") |>
    req_url_query(
      lang = "fi",
      langFallback = "false",
      pageSize = 1000,
      onlyFinnish = "true",
      includeHidden = "true",
      selectedFields = paste(
        "id",
        "scientificNameAuthorship",
        "scientificName",
        "taxonRank",
        "parents",
        "vernacularName",
        "synonyms",
        "subjectiveSynonyms",
        "heterotypicSynonyms",
        "homotypicSynonyms",
        "objectiveSynonyms",
        "hiddenTaxon",
        sep = ","
      ),
      access_token = Sys.getenv("FINBIF_ACCESS_TOKEN"),
      page = 1
    )

  res <-
    req_perform(req) |>
    resp_body_json()

  taxonomy <- res$results

  for (i in seq_len(res$lastPage)[-1]) {
    req <- req_url_query(req, page = i)

    res <-
      req_perform(req) |>
      resp_body_json()

    taxonomy <- c(taxonomy, res$results)
  }

  saveRDS(taxonomy, "taxonomy.rds")
} else {
  taxonomy <- readRDS("taxonomy.rds")
}

# Remove hidden taxa
taxonomy <- taxonomy[!vapply(taxonomy, getElement, NA, "hiddenTaxon")]

# Extract taxon concept identfiers
taxonConceptID <- vapply(taxonomy, getElement, "", "id")

# Extract taxon ranks of taxon concepts
taxonRank <- sub("MX\\.", "", vapply(taxonomy, get_property, "", "taxonRank"))
names(taxonRank) <- taxonConceptID

# Flatten and transform the taxonomy into a single table
taxonomy_flat <- do.call(rbind, lapply(taxonomy, flatten_concept))
taxonomy_flat <- transform(
  taxonomy_flat,
  scientificName = ifelse(
    taxonRank == "aggregate",
    sub(
      paste(
        paste0(
          " ",
          c(
            "coll\\.",
            "group",
            "-group",
            "-type",
            "complex",
            "sensu lato",
            "-ryhmä",
            "sec\\.",
            "s\\.l\\.",
            "-ryhmä s\\. lat\\."
          ),
          "$"
        ),
        collapse = "|"
      ),
      " agg.",
      scientificName
    ),
    scientificName
  ),
  scientificNameAuthorship = ifelse(
    taxonRank == "aggregate",
    NA_character_,
    scientificNameAuthorship
  )
)
taxonomy_flat <- transform(
  taxonomy_flat,
  scientificName = ifelse(
    is.na(scientificNameAuthorship),
    scientificName,
    paste(scientificName, scientificNameAuthorship)
  ),
  genericName = ifelse(
    taxonRank %in%
      c(
        "genus",
        "subgenus",
        "section",
        "aggregate",
        "species",
        "subspecies",
        "form",
        "variety"
      ),
    vapply(
      strsplit(scientificName, " "),
      \(x) {
        x <- x[grepl("^[[:upper:]][[:lower:]]+$", x)]
        if (length(x) > 0) x[[1]] else NA_character_
      },
      ""
    ),
    NA
  ),
  infragenericEpithet = ifelse(
    taxonRank %in%
      c(
        "subgenus",
        "section",
        "aggregate",
        "species",
        "subspecies",
        "form",
        "variety"
      ),
    vapply(
      strsplit(scientificName, " "),
      \(x) {
        x <- gsub("\\(|\\)", "", x)
        x <- x[grepl("^[[:upper:]][[:lower:]]+$", x)]
        if (length(x) > 1) x[[2]] else NA_character_
      },
      ""
    ),
    NA
  ),
  specificEpithet = ifelse(
    taxonRank %in%
      c(
        "aggregate",
        "species",
        "subspecies",
        "form",
        "variety"
      ),
    vapply(
      strsplit(scientificName, " "),
      \(x) {
        x <- x[grepl("^[[:lower:]]([[:lower:]]|-)*$", x)]
        if (length(x) > 0) x[[1]] else NA_character_
      },
      ""
    ),
    NA
  ),
  infraspecificEpithet = ifelse(
    taxonRank %in%
      c("subspecies", "form", "variety"),
    vapply(
      strsplit(scientificName, " "),
      \(x) {
        x <- x[grepl("^[[:lower:]]([[:lower:]]|-)*$", x)]
        if (length(x) > 1) x[[2]] else NA_character_
      },
      ""
    ),
    NA
  )
)

# Extract legitimate name usages from taxonomy
NameUsage <- subset(taxonomy_flat, taxonRank %in% ranks)
NameUsage <- transform(
  NameUsage,
  taxonRank = sub("aggregate", "speciesAggregate", taxonRank)
)

# Remove GBIF taxa
NameUsage <- subset(NameUsage, !grepl("gbif", taxonID))

# Find and mark potential pro parte synonyms
proParte <-
  duplicated(NameUsage[c("scientificName", "taxonRank")]) |
  duplicated(NameUsage[c("scientificName", "taxonRank")], fromLast = TRUE)
proParte <- proParte & NameUsage$taxonomicStatus != "accepted"
NameUsage$taxonomicStatus[proParte] <- "proParteSynonym"

write.table(
  subset(NameUsage, taxonomicStatus == "proParteSynonym"),
  "duplicate-name-usage.txt",
  quote = FALSE,
  sep = "\t",
  na = "",
  row.names = FALSE
)

# Remove pro parte synonyms
NameUsage <- subset(NameUsage, taxonomicStatus != "proParteSynonym")

# Add accepted and parent name usages
rownames(NameUsage) <- NameUsage$taxonID
NameUsage$acceptedNameUsage <- NameUsage[
  NameUsage$acceptedNameUsageID,
  "scientificName"
]
NameUsage$parentNameUsage <- NameUsage[
  NameUsage$parentNameUsageID,
  "scientificName"
]

NameUsage <- subset(NameUsage, acceptedNameUsageID %in% NameUsage$taxonID)

# Remove higher taxon with no species as children
keep_taxa <- apply(
  NameUsage,
  1,
  \(x) {
    accepted <- x[["taxonomicStatus"]] == "accepted"
    higher_taxon <-
      !x[["taxonRank"]] %in% c("species", "subspecies", "form", "variety")
    if (accepted & higher_taxon) {
      has_species(x[["taxonID"]])
    } else {
      TRUE
    }
  }
)

NameUsage <- NameUsage[keep_taxa, ]
NameUsage <- subset(NameUsage, acceptedNameUsageID %in% NameUsage$taxonID)
NameUsage$acceptedNameUsage <- NULL
NameUsage$parentNameUsage <- NULL

NameUsageOld <- read.table(
  "name-usage.txt",
  TRUE,
  sep = "\t",
  quote = "",
  na.strings = ""
)

row.names(NameUsageOld) <- NameUsageOld$taxonID

new <- NameUsage[setdiff(NameUsage$taxonID, NameUsageOld$taxonID), ]

removed <- NameUsage[setdiff(NameUsageOld$taxonID, NameUsage$taxonID), ]

comparison <- compareDF::compare_df(
  NameUsage[intersect(NameUsage$taxonID, NameUsageOld$taxonID), ],
  NameUsageOld[intersect(NameUsage$taxonID, NameUsageOld$taxonID), ]
)

# Write out name usage data to a text file
write.table(
  NameUsage,
  "name-usage.txt",
  quote = FALSE,
  sep = "\t",
  na = "",
  row.names = FALSE
)

# Extract Finnish vernacular name data
vernacularName <-
  data.frame(
    taxonID = paste0("http://tun.fi/", vapply(taxonomy, getElement, "", "id")),
    vernacularName = vapply(taxonomy, get_property, "", "vernacularName"),
    language = "fi",
    countryCode = "FI"
  ) |>
  subset(
    !is.na(vernacularName) &
      vernacularName != "" &
      !grepl("\\(alalaji", vernacularName) &
      taxonID %in% NameUsage$taxonID
  )

# Write out name vernacular name to a text file
write.table(
  vernacularName,
  "vernacular-name.txt",
  quote = FALSE,
  sep = "\t",
  na = "",
  row.names = FALSE
)

# req <-
#   request("https://www.gbif.org") |>
#   req_url_path("api") |>
#   req_url_path_append("species") |>
#   req_url_path_append("search") |>
#   req_url_query(
#     dataset_key = "f95250e7-49f4-4d2e-a04e-35533dee3318",
#     issue = "PARTIALLY_PARSABLE",
#     origin = "SOURCE",
#     limit = 400,
#   )

# res <-
#   req_perform(req) |>
#   resp_body_json()

# has_issues <- NameUsage

# has_issues$issue_gbif_partially_parseable <- has_issues$taxonID %in%
#   vapply(res$results, getElement, "", "taxonID")

# req <- req_url_query(req, issue = "SCIENTIFIC_NAME_ASSEMBLED")

# res <-
#   req_perform(req) |>
#   resp_body_json()

# has_issues$issue_gbif_scientific_name_assembled <-
#   has_issues$taxonID %in% vapply(res$results, getElement, "", "taxonID")

# clb_issues <- read.csv("https://api.checklistbank.org/dataset/290762/issues")

# for (i in unique(unlist(strsplit(clb_issues$status, ";")))) {
#   has_issues[[paste0("issue_clb_", tolower(i))]] <-
#     has_issues$taxonID %in%
#     clb_issues$ID[grepl(i, clb_issues$status, fixed = TRUE)]
# }

# has_issues$issue_clb_duplicate_name <- NULL
# has_issues$issue_clb_authorship_contains_taxonomic_note <- NULL
# has_issues$issue_clb_synonym_rank_differs <- NULL
# has_issues$issue_clb_citation_unparsed <- NULL
# has_issues$issue_clb_authorship_contains_nomenclatural_note <- NULL
# has_issues$issue_clb_rank_name_suffix_conflict <- NULL

# has_issues <- filter(has_issues, if_any(starts_with("issue")))

# write.table(
#   has_issues,
#   "hasIssues.txt",
#   quote = FALSE,
#   sep = "\t",
#   na = "",
#   row.names = FALSE
# )

# duplicates <- transform(
#   NameUsage,
#   scientificName = mapply(
#     sub,
#     pattern = paste0(" ", scientificNameAuthorship),
#     x = scientificName,
#     MoreArgs = list(replacement = "", fixed = TRUE)
#   )
# )
# duplicates <- duplicates[
#   duplicated(duplicates$scientificName) |
#     duplicated(duplicates$scientificName, fromLast = TRUE),
# ]
# duplicates <- duplicates[order(duplicates$scientificName), ]
# duplicates <- transform(
#   duplicates,
#   taxonomicStatus = ifelse(
#     taxonomicStatus == "accepted",
#     taxonomicStatus,
#     "synonym"
#   )
# )

# write.table(
#   duplicates,
#   "duplicate-name-usage.txt",
#   quote = FALSE,
#   sep = "\t",
#   na = "",
#   row.names = FALSE
# )
