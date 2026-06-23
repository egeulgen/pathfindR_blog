---
title: "pathfindR Relese 3.0.0"
author: "Ege Ulgen"
date: 2026-06-23
categories: ["R"]
tags: ["R package", "release", "enrichment", "bioinformatics"]
editor_options: 
  markdown: 
    wrap: 72
---

## What is pathfindR?

`pathfindR` is an R package for **enrichment analysis via active
subnetworks**. Conventional enrichment methods ignore protein--protein
interaction (PPI) information, which can leave the picture incomplete.
pathfindR's main function identifies active subnetworks in a
protein-protein interaction network using a user-provided list of genes
and associated p values, then performs enrichment analyses on those
subnetworks to identify enriched terms (pathways, or more broadly, gene
sets) that may underlie the phenotype of interest.

Beyond the core workflow, pathfindR also lets you:

- cluster the enriched terms and identify a representative term in each
  cluster,
- score the enriched terms on a per-sample basis, and
- compare two pathfindR results against each other.

The methods are described in detail in Ulgen E, Ozisik O, Sezerman OU
(2019), *pathfindR: An R Package for Comprehensive Identification of
Enriched Pathways in Omics Data Through Active Subnetworks*, Front.
Genet.
([doi:10.3389/fgene.2019.00858](https://doi.org/10.3389/fgene.2019.00858)).

## What's new?

The headline of pathfindR 3.0.0 is something you mostly *won't* notice:
the Java dependency is gone. Active subnetwork search - the step that
grows candidate subnetworks over a protein--protein interaction
network - has been re-implemented entirely in R and C++ (via Rcpp),
replacing the old Java JAR that pathfindR used to shell out to.

In practice this means no JDK, and no "Java not found" issues on a fresh
machine. You install the package and it runs:

``` r
install.packages("pathfindR")        # or remotes::install_github("egeulgen/pathfindR")
library(pathfindR)
 
output <- run_pathfindR(example_pathfindR_input)
```

The greedy and simulated-annealing search methods produce results that
are numerically identical to the old JAR given the same seed and inputs,
so existing analyses remain reproducible - you're just getting them
without the extra runtime to install and configure.

## A few renames

Some exported functions were renamed to follow a more consistent,
snake_case style. If you call these directly in scripts, update the
names:

- `active_snw_search()` → `get_active_subnetworks()` (its signature
  changed too)

- `filterActiveSnws()` → `filter_active_subnetworks()`

The most common entry point, `run_pathfindR()`, is unchanged.

## Term--gene graph, split in two

The old `term_gene_graph()` did two jobs at once: build the graph and
draw it. Those are now separate, which makes it much easier to tweak the
visualization or reuse the underlying structure:

- `create_term_gene_graph()` builds the `igraph` object.
- `create_term_gene_plot()` renders a plot from that object.

``` r
g <- create_term_gene_graph(output, num_terms = 10)
create_term_gene_plot(g)
```

## Should you upgrade?

Yes! especially if Java setup has ever bitten you or your collaborators.
The science is the same, the results are the same, and the install is
now a single, self-contained step. Full notes are in the
[changelog](https://github.com/egeulgen/pathfindR/blob/master/NEWS.md).
