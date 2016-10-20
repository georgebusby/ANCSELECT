compute.manhattan_plot_positions <- function(
  chromosome,
  position,
  seperation = 1000000
) {
  O = order( chromosome, position, decreasing = F )
  chromosome = chromosome[O]
  position = position[O]
  plot_pos = position
  plot_pos[-1] = plot_pos[-1] - plot_pos[1:( length( plot_pos ) - 1 )]
  plot_pos[ which( plot_pos < -100 ) ] = seperation
  plot_pos = cumsum( plot_pos )
  return( plot_pos )
}
fix_chromosome <- function( chromosome ) {
  chromosome = as.character( chromosome )
  n = nchar( chromosome )
  w = which( n == 1 ) ;
  chromosome[w] = sprintf( "0%s", chromosome[w] )
  return( chromosome )
}


load.annotation.data <- function(
  annotations = list(
    gwas = list(
      name = "NHGRI GWAS Catalog",
      filename = "/mnt/kwiat/data/1/galton/malariagen/human/reference/NHGRI GWAS Catalog/2013-09-25_gwascatalog.txt",
      chromosome.column = "Chr_id",
      position.column = "Chr_pos",
      name.column = "Disease/Trait"
    )
    # eqtl = list(
    #   name = "GTEX eQTL",
    #   filename = "/mnt/kwiat/data/1/galton/malariagen/human/reference/GTEX/GTEX-130925-1732-19664.tab",
    #   chromosome.column = "SNP Chr.",
    #   position.column = "SNP Position",
    #   name.column = "Gene"
    # ),
    # ABS1 = list(
    #   name = "Ancient balancing selection haplotypes",
    #   filename = "/mnt/kwiat/data/1/galton/malariagen/human/reference/articles/Leffler et al Multiple Instances of Ancient Balancing Selection Shared Between Humans and Chimpanzees/TableS4_SharedHaplotypeRegions.csv",
    #   chromosome.column = "Chr (hg19)",
    #   position.column = "Position (hg19)",
    #   name.column = "rs#"
    # ),
    # ABS2 = list(
    #   name = "Ancient balancing selection coding SNPs",
    #   filename = "/mnt/kwiat/data/1/galton/malariagen/human/reference/articles/Leffler et al Multiple Instances of Ancient Balancing Selection Shared Between Humans and Chimpanzees/TableS5_SharedCodingSNPs.csv",
    #   chromosome.column = "Chr (hg19)",
    #   position.column = "Position (hg19)",
    #   name.column = "rs#"
    # )
  )
) {
  annotation.data = data.frame()
  for( name in names( annotations )) {
    annotation = annotations[[ name ]]
    if( substring( annotation$filename, nchar( annotation$filename ) - 3, nchar( annotation$filename ) ) == '.csv' ) {
      X = read.csv( annotation$filename, as.is = TRUE, header = T, check.names = F ) ;
    } else {
      X = read.delim( annotation$filename, as.is=TRUE, header=TRUE, check.names = F ) ;
    }
    annotation.data = rbind(
      annotation.data,
      data.frame(
        type = name,
        label = annotation$name,
        name = X[, annotation$name.column ],
        chromosome = fix_chromosome( X[, annotation$chromosome.column ] ),
        position = X[, annotation$position.column ]
      )
    )
  }
  return( annotation.data ) ;
}