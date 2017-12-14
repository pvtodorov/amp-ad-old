## msbb.R - Wrangling of AMP-AD Mount Sinai Brain Bank Dataset
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
suppressMessages(library( magrittr ))
library( stringr )
suppressMessages(library( synapseClient ))

## Maps Synapse IDs to their corresponding brain regions
## Uses table from wiki at https://www.synapse.org/#!Synapse:syn3157699
## For now, we withhold Affy 133 Plus2 datasets due to large batch effects
##   with the rest of the data
MSBB.regions <- function()
{
    rbind(c( "syn3191095", "BM10-FP" ),
          c( "syn3191097", "BM17-OVC" ),
          c( "syn3191099", "BM20-ITG" ),
          c( "syn3191101", "BM21-MTG" ),
          c( "syn3191103", "BM22-STG" ),
          c( "syn3191105", "BM23-PCC" ),
          c( "syn3191107", "BM32-AC" ),
          c( "syn3191109", "BM36-PHG" ),
          c( "syn3191111", "BM38-TP" ),
          c( "syn3191113", "BM4-PCG" ),
          c( "syn3191115", "BM44-IFG" ),
          c( "syn3191117", "BM46-PFC" ),
          c( "syn3191119", "BM7-SPL" ),
          c( "syn3191123", "BM8-FC" ),
##          c( "syn3191125", "BMa-AMYG" ),
          c( "syn3191127", "BMb-CD" ),
          c( "syn3191129", "BMc-HIPP" ),
##          c( "syn3191131", "BMd-NAc" ),
          c( "syn3191133", "BMe-PT" )) %>%
        set_colnames( c( "file.id", "Region" ) ) %>% as_tibble
}

## "Quiet" versions of read_tsv and read_csv
readq_tsv <- function( ... ) {suppressMessages( read_tsv( ... ) )}
readq_csv <- function( ... ) {suppressMessages( read_csv( ... ) )}

wrangle.MSBB.mRNA <- function()
{
    ## Retrieve the list of files in the expression folder
    VSyn <- synQuery( "select id, name from file where parentId==\"syn3219494\"" ) %>%
        filter( grepl( "tsv", file.name ) ) %>% inner_join( MSBB.regions(), by="file.id" ) %>%
        select( file.id, Region )

    ## Download all the files and retrieve local filenames
    cat( "Downloading data from Synapse...\n" )
    VSyn <- VSyn %>% rowwise %>%
        mutate( fnLocal = synGet( file.id, downloadLocation = local.dir )@filePath ) %>%
        ungroup

    ## Load all the expression data
    cat( "Loading individual files...\n" )
    XX <- lapply( VSyn$fnLocal, readq_tsv ) %>% setNames( VSyn$Region )

    ## Isolate the gene symbol column from other row annotations
    ## Remove transcripts that span multiple genes or don't map to any genes
    XX <- lapply( XX, function(X) {
        select( X, -ID, -GB_ACC, -ENTREZ_GENE_ID ) %>%
            rename( Gene = Gene.Symbol ) } ) %>%
        lapply( na.omit ) %>%
        lapply( function(X) {filter( X, !grepl( "///", Gene ) )} )

    ## Compute the average expression for duplicate gene entries
    cat( "Averaging multiple entries per gene...\n" )
    XX <- lapply( XX, function( X )
    { cat("."); X %>% group_by( Gene ) %>% summarize_all( funs(median) ) } )
    cat( "\n" )

    ## Temporarily store the resulting matrices
    cat( "Writing intermediate result...\n" )
    save( XX, file=file.path( local.dir, "expression.RData" ) )

    ## Reshape the matrices to be Samples by Genes
    ## Annotate each sample with corresponding brain region
    ## Collapse everything into a single data frame
    cat( "Combining everything into a single data.frame...\n" )
    ZZ <- lapply( names(XX), function(i) { XX[[i]] %>% gather(Sample, Value, -Gene) %>%
                                               mutate(Region=i) } )
    RR <- bind_rows( ZZ ) %>% spread( Gene, Value )

    ## Store the result
    write_tsv( RR, file.path( local.dir, "expression.tsv" ) )

    ## Retrieve clinical information from Synapse
    f <- synGet( "syn3205399", downloadLocation = local.dir )
    YY <- read.delim( f@filePath )
}

## Wrangles a metadata matrix from mRNA expression
wrangle.MSBB.meta <- function()
{
    stop( "Needs revision" )
    load( fn.local( "expression.RData" ) )

    ## Retrieve the brain region and patient ID from the sample names
    vl <- strsplit( colnames(X), "-" )
    vbr <- lapply( vl, function(z) {paste( z[1:2], collapse="-" )} ) %>% unlist()
    vpid <- lapply( vl, function(z) {substring(z[3], 2)} ) %>% unlist() %>% as.integer()

    ## Assign samples to their microarray platform
    pf <- c( "Affy 133AB", "Affy 133 Plus2" )
    vpf <- pf[(vbr %in% c("BMa-AMYG", "BMd-NAc")) + 1]
    
    ## Retrieve other clinical information from Synapse
    f <- synGet( "syn3205399", downloadLocation = local.dir )
    YY <- read.delim( f@filePath )

    ## Compose the metadata table
    Y <- data.frame( Sample = colnames(X), Region = vbr, Platform = vpf,
                    BrainBank = vpid ) %>% inner_join( YY )

    write.table( Y, file=fn.local("meta.tsv"), quote=FALSE, sep="\t", row.names=FALSE )
}

## Parse local directory specification
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) == 0 )
{
    cat( "NOTE: No directory specified on command line. Using default.\n" )
    local.dir <- "/data/AMP-AD/MSBB"
} else { local.dir <- argv[1] }
    

## Create directory if it doesn't exist
dir.create( local.dir, showWarnings=FALSE )
cat( "Wrangling MSBB dataset to", local.dir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapseLogin( rememberMe=TRUE )
