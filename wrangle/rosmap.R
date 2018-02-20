## Wrangling of ROS/MAP RNAseq data and matching clinical annotations
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
suppressMessages(library( synapseClient ))
library( stringr )

## Composes a mapping between ENSEMBL IDs and HUGO names
ens2hugo <- function()
{
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    tx <- ensembldb::transcripts( edb, column=c("gene_id", "gene_name") )
    data_frame( HUGO = tx$gene_name, ENSEMBL = tx$gene_id ) %>% distinct
}

## Parse local directory specification
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) == 0 )
{
    cat( "NOTE: No directory specified on command line. Using default.\n" )
    local.dir <- "/data/AMP-AD/ROSMAP"
} else { local.dir <- argv[1] }

## Create directory if it doesn't exist
dir.create( local.dir, showWarnings=FALSE )
cat( "Wrangling ROS/MAP dataset to", local.dir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapseLogin( rememberMe=TRUE )

## Read raw expression matrix
cat( "Downloading expression data...\n" )
fnX <- synGet( "syn3505720", downloadLocation = local.dir )@filePath
cat( "Loading local copy...\n" )
Xraw <- suppressMessages( read_tsv( fnX ) )

## Map ENSEMBL Gene IDs to HUGO
## Remove alternative splice forms, non-coding RNA and duplicated genes
## There are only 14 duplicates after removing non-coding transcripts
cat( "Mapping gene IDs to HUGO...\n" )
E2H <- ens2hugo()
cat( "Removing alternative splice forms, non-coding RNA and duplicate entries...\n" )
f <- function( x, pattern ) { filter( x, !grepl(pattern, HUGO) ) }
X <- Xraw %>% mutate( ENSEMBL = str_split( gene_id, "\\.", simplify=TRUE )[,1] ) %>%
    inner_join( E2H, by="ENSEMBL" ) %>% filter( !grepl("\\.", HUGO) ) %>%
    filter( !(HUGO %in% c("Y_RNA", "Metazoa_SRP", "Vault", "5S_rRNA")) ) %>%
    f( "^MIR" ) %>% f( "^RNU" ) %>% f( "^SNOR" ) %>% f( "^U[1-9]$" ) %>%
    f( "^SCARNA" ) %>% f( "^sno" ) %>% f( "^LINC" ) %>% f( "-AS[1-9]$" ) %>%
    f( "^ACA[1-9]" ) %>% filter( !duplicated( HUGO ) ) %>%
    select( -tracking_id, -gene_id, -ENSEMBL )

## Log-transform the data and combine the replicates
cat( "Additional processing...\n" )
flog <- function(v) {log2(v+1)}
fmed <- function(x) {x %>% as.matrix %>% apply( 1, median )}
XX <- X %>% mutate_at( vars(-HUGO), funs(flog) ) %>%
    mutate( `492_120515_j` = fmed(select( ., contains("492_120515") )) ) %>%
    select( -`492_120515_0`, -`492_120515_6`, -`492_120515_7` ) %>%
    gather( rnaseq_id, Value, -HUGO ) %>%
    mutate( rnaseq_id = str_sub( rnaseq_id, 0, -3 ) )

## Match sample IDs against individual identifiers
cat( "Matching sample and individual IDs...\n" )
fnZ <- synGet( "syn3382527", downloadLocation = local.dir )@filePath
XZ <- suppressMessages( read_csv(fnZ) ) %>% select( projid, rnaseq_id ) %>% na.omit %>%
    distinct %>% inner_join( XX, ., by="rnaseq_id" )

## Match expression data up against the following clinical covariates:
## ID, PMI, AOD, CDR, Braak, BrodmannArea
cat( "Matching against clinical covariates...\n" )
fnY <- synGet( "syn3191087", downloadLocation = local.dir )@filePath
Y <- suppressWarnings( suppressMessages( read_csv(fnY) ) ) %>%
    select( projid, PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc ) %>%
    mutate( BrodmannArea = "BM9/BM46" )

## Combining everything into a common data frame
cat( "Finalizing...\n" )
XY <- inner_join( Y, XZ, by="projid" ) %>% rename( ID = projid, Barcode = rnaseq_id ) %>%
    spread( HUGO, Value )

## Write out wrangled dataset to file
fnOut <- file.path( local.dir, "rosmap-wrangled.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
write_tsv( XY, fnOut )
