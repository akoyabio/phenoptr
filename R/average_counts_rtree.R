#' rtree implementation of count_within_many_impl
#' @seealso count_within_many_impl
#' @md
count_within_many_impl_rtree = function(csd, name, combos, radii,
                                        phenotype_rules) {
  # We will call count_within_detail to add a count column for
  # each 'to' phenotype. Figure out what they are.
  to_phenotypes = combos %>% purrr::map_chr(list('pair', 2)) %>% unique()
  to_phenotypes = phenotype_rules[to_phenotypes]

  # Get the categories of interest.
  categories = combos %>% purrr::map_chr('category') %>% unique()

  # Subset to what we care about, for faster distance calculation
  if (!anyNA(categories))
    csd = csd %>% dplyr::filter(`Tissue Category` %in% categories)

  # Compute individual counts for each cell by category.
  counts = purrr::map_dfr(categories, function(category) {
      if (!is.na(category))
        d_cat = csd %>% dplyr::filter(`Tissue Category`==category)
      else d_cat = csd
      dplyr::bind_cols(d_cat,
                       phenoptr:::count_within_detail(d_cat, to_phenotypes, radii))
    })


  # Compute summary stats for a single dataset, a single category and
  # (from, to) pair, and all radii
  summarize_combo = function(d, category, from, to) {
    if (!is.na(category))
      d_cat = d %>% dplyr::filter(`Tissue Category`==category)
    else {
      d_cat = d
      category = 'all'
    }

    to_count = sum(phenoptr::select_rows(d_cat, phenotype_rules[[to]]))

    # The rows of interest
    rows = d_cat %>%
      dplyr::filter(phenoptr::select_rows(., phenotype_rules[[from]]))

    # Now we can summarize per radius
    if (nrow(rows) > 0) {
      purrr::map_dfr(radii, function(radius) {
        count_col = count_col_name(to, radius)
        counts = rows[[count_col]] # Single column of interest

        tibble::tibble(category=category,
                       from=from, to=to, radius=radius,
                       from_count=nrow(rows), to_count=to_count,
                       from_with = sum(counts>0, na.rm=TRUE),
                       within_mean = mean(counts, na.rm=TRUE))
    })
    } else {
      tibble::data_frame(
        category=category,
        from=from, to=to, radius = radii,
        from_count = 0L,
        to_count = to_count,
        from_with = 0L,
        within_mean = NA
      )

    }
  }

  # Now we can do the work, computing summary stats for all combos and radii
  purrr::map_dfr(combos, function(combo) {
    summarize_combo(counts, combo$category, combo$pair[[1]], combo$pair[[2]])
  })
}


#' Compute count within for individual cells in a single field.
#'
#' Very fast version using `rtree::countWithinDistance()`.
#' @param csd Cell seg data.
#' @param radii Vector of radii to search within.
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' will use `unique_phenotypes(csd)`. Counts are from each cell to each
#' phenotype.
#' @param radii The radius or radii to search within.
#' @export
#' @md
count_within_detail = function(csd, phenotypes=NULL, radii) {
  # Check for multiple samples, this is probably an error
  if ('Sample Name' %in% names(csd) && length(unique(csd$`Sample Name`))>1)
    stop('Data appears to contain multiple samples.')

  phenotypes = validate_phenotypes(phenotypes, csd)
  field_locs = csd %>%
    dplyr::select(X=`Cell X Position`, Y=`Cell Y Position`) %>%
    as.matrix()

  result = purrr::imap(phenotypes, function(phenotype, name) {
    # Which cells are in the target phenotype?
    phenotype_cells = select_rows(csd, phenotype)

    if (sum(phenotype_cells)>0) {
      # Make an rtree of the phenotype cells
      to_cells_locs = field_locs[phenotype_cells,, drop=FALSE]
      to_cells_tree = rtree::RTree(to_cells_locs)

      # Now compute count within for each radius
      purrr::map(radii, function(radius) {
        within = rtree::countWithinDistance(to_cells_tree,
                                            field_locs, radius)

        # Subtract one for cells of type `pheno`; we don't want to count self
        within = within - phenotype_cells

        col_name = count_col_name(name, radius)
        list(within) %>% rlang::set_names(col_name)
      }) %>% purrr::flatten()
    }
    else {
      # No cells of the selected phenotype
      count_col = list(rep(0L, nrow(csd)))
      purrr::map(radii, function(radius) {
        col_name = count_col_name(name, radius)
        count_col %>% rlang::set_names(col_name)
      }) %>% purrr::flatten()
    }
  })

  tibble::as_tibble(purrr::flatten(result))
}

count_col_name = function(pheno_name, radius) {
  paste0(pheno_name, ' within ', radius)
}
