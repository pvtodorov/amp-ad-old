## Prototype analysis
##
## by Artem Sokolov

library( tidyverse )
library( caret )
library( randomForest )

## Creates a data slice to be used for prototype analyses
## Takes the first sample for each individual
createProtoSlice <- function()
{
    ## Load the entire dataset, removing samples with no Braak scores
    XY <- read_tsv( "/data/AMP-AD/MSBB/msbb-wrangled.tsv" ) %>% na.omit

    ## Select the first sample for each individual
    ZZ <- XY %>% group_by( individualIdentifier ) %>% slice( 1 ) %>% ungroup

    ## Clean up and save to file
    XY <- ZZ %>% select( -individualIdentifier, -CDR, -BrodmannArea, -barcode )
    save( XY, file="proto-slice.RData" )
}

## Trains a random forest model over the requested set of genes
## Returns out-of-bag error estimate across a grid of metaparameter values
## Estimates accuracy via out-of-bag error, rather than cross-validation
eval.geneset <- function( XY, vGenes )
{
    Xtrain <- XY %>% select( one_of(vGenes) ) %>% as.data.frame
    ytrain <- XY$bbscore %>% as.numeric
    cat( "Training over", length(vGenes), "genes\n" )
    mm <- train( Xtrain, ytrain, method="rf", trControl = trainControl( method="oob" ) )
    mean( mm$results$Rsquared )
}

## Given cardinality value, randomly selects a gene set and evaluates the model
eval.random <- function( XY, nGenes )
{
    setdiff( colnames(XY), "bbscore" ) %>% sample( nGenes ) %>% eval.geneset( XY, . )
}

## Evaluates nSets randomly-selected gene sets of cardinality nGenes each
eval.randSets <- function( XY, nGenes, nSets )
{
    res <- rep( NA, nSets )
    for( i in 1:nSets )
    {
        cat( "Random set", i, "  " )
        res[i] <- eval.random( XY, nGenes )
    }
    res
}

## Given a gene set, evaluates it and 100 randomly-selected background sets of same cardinality
eval.all <- function( XY, vGenes, nBk = 100 )
{
    ## Identify the set of genes in common with what's in the data
    v <- intersect( colnames(XY), vGenes )
    cat( "Of the provided", length(vGenes), "genes,", length(v), "are in the data\n" )

    ## Evaluate the true set
    vTrue <- eval.geneset( XY, v )

    ## Evaluate background sets
    vBk <- eval.randSets( XY, length(v), nBk )

    data_frame( Set = c("True", rep("Bk",nBk)), R2 = c(vTrue, vBk) )
}

main7 <- function()
{
    ## Load the data
    load( "proto-slice.RData" )

    ## Cluster 7 genes
    v <- scan( "clus7.txt", what=character() )
    eval.all( XY, v ) %>% write_csv( "res7.csv" )
}

main6 <- function()
{
    ## Load the data
    load( "proto-slice.RData" )
    
    ## Cluster 6 genes
    v <- scan( "clus6.txt", what=character() )
    eval.all( XY, v ) %>% write_csv( "res6.csv" )
}

main3 <- function()
{
    ## Load the data
    load( "proto-slice.RData" )
    
    ## Cluster 3 genes
    v <- scan( "clus3.txt", what=character() )
    eval.all( XY, v ) %>% write_csv( "res3.csv" )
}

## Plots the results of a single set
plot.set <- function( X )
{
    gg <- ggplot( filter( X, Set == "Bk" ), aes(x=R2) ) + theme_bw() +
        geom_density( size = 1, fill = "steelblue", alpha = 0.2 ) + ylab( "Density" ) +
        xlab( "Performance estimate for predicting Braak score from mRNA (R^2)" ) +
        xlim( c(0,0.3) ) +
        geom_vline( xintercept = filter( X, Set == "True" )$R2, color = "red", size=1 ) +
        theme( axis.text = element_text( face="bold", size = 12 ),
              axis.title = element_text( face="bold", size = 12 ) )
}

main.figs <- function()
{
    gg3 <- read_csv( "res3.csv" ) %>% plot.set
    gg6 <- read_csv( "res6.csv" ) %>% plot.set
    gg7 <- read_csv( "res7.csv" ) %>% plot.set

    gt <- arrangeGrob( gg3, gg6, gg7, ncol=1 )
    grid.newpage()
    grid.draw( gt )
    ggsave( "fig1.png", gt )
}
