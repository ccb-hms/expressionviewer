#' Create a shiny app to display expression values and/or alternative splicing
#'
#' This function takes a SummarizedExperiment object, SGFeatureCounts object, or
#' a list of either of these classes and produces a shiny app with a barplot
#' illustrating expression levels per gene and (in the case of SGFeatureCounts)
#' a heatmap illustrating alternative splicing per sample. The assumptions for
#' each object case are listed below.
#'
#' Summarized Experiment:
#' The sample information wishing to be displayed as the independent variable
#' is in the colData with the column name inputted into the variable "metaNames".
#' The rownames of the expression matrix are gene names.
#'
#' SGFeatureCounts Object:
#' The sample information wishing to be displayed as the independent variable
#' is in the colData with the column name inputted into the variable "metaNames".
#' The geneID() of the experiment is complete with the gene annotation style
#' desired (NCBI, HGNC Symbol, etc.).
#'
#' List of either SummarizedExperiments or SGFEatureCounts:
#' The "metaNames" are consistently named for each element of the list. Any
#' class assumptions are retained when part of a list.
#'
#' @param experiment A SummarizedExperiment object, SGFEatureCounts object, list of
#' SummarizedExperiment objects, or list of SGFEatureCounts objects.
#' @param metaNames A quoted column name for the colData column with the information
#' to describe the independent variable of the box plot ex. "sample_name"
#' @param printNames An optional name or list of names to be displayed as selection
#' choice(s) for the experiment(s) being displayed.
#' @param paper An optional paper title from which the data being displayed was
#' published.
#' @param paperlink An optional url that links to the paper with the "paper"
#' title
#' @param apptitle An optional title for the shiny app that is produced
#' @return A shiny app displaying the expression (and potential splicing) of
#' the experiment input
#' @export
ExpressionViewer <- function(experiment = counts, metaNames = "label.main",
                             printNames = " ", paper = " ",
                             paperlink = " ", apptitle = "Expression Viewer:"){
    #naming experiment list
    if(is(experiment)[1] == "list"){
        inputOps <- names(experiment)
        names(inputOps) <- printNames
        inputOps <- as.list(inputOps)}
    else if(is(experiment)[1] == "SummarizedExperiment" | is(experiment)[1] == "SGFeatureCounts"){
        inputOps <- substitute(experiment)}
    else(return(
    "input 'experiment' must be of class SummarizedExperiment, SGFeatureCounts, or a list of these classes"))

    #Checking if metaNames are present in all experiments
    if(is(experiment)[1] == "list"){
    colList <- lapply(experiment, function(x){colnames(colData(x))})
    present <- vapply(colList, function(x){metaNames %in% x}, logical(1))
    stopifnot(all(present))}

    #Get genenames present in (all) experiment(s)
    if(is(experiment)[1] == "list" & is(experiment[[1]])[1] == "SummarizedExperiment"){
        autocompletelist <- list()
        for(i in 1:length(experiment)){
            testlist[[i]] <- rownames(experiment[[i]])}
        testlist <- unlist(testlist)
        autocompletelist <- testlist[which(table(testlist) == length(experiment))]
        autocompletelist <- sort(autocompletelist)
    }
    if(is(experiment)[1] == "SummarizedExperiment"){
        autocompletelist <- rownames(experiment)
        autocompletelist <- sort(autocompletelist)
    }
    if(is(experiment)[1] == "list" & is(experiment[[1]])[1] == "SGFeatureCounts"){
        autocompletelist <- list()
        for(i in 1:length(experiment)){
            testlist[[i]] <- geneID(experiment[[i]])}
        testlist <- unlist(testlist)
        autocompletelist <- testlist[which(table(testlist) == length(experiment))]
        autocompletelist <- sort(autocompletelist)
    }
    if(is(experiment)[1] == "SGFeatureCounts"){
        autocompletelist <- geneID(experiment)
        autocompletelist <- sort(unique(autocompletelist))
    }


    #user input and page layout specification
    ui <- shiny::fluidPage(
        shiny::titlePanel(apptitle),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput(inputId = "dataset", "Dataset/Experiment:",
                            choices = inputOps, selectize = TRUE),
                br(),
                shiny::selectizeInput(inputId = "genename", label = "Gene Name:",
                               multiple = FALSE, selected = autocompletelist[1],
                               choices = autocompletelist),
                shiny::uiOutput("tab")
            ),
            shiny::mainPanel(
                shiny::plotOutput(outputId = "box", width = "100%"), position = "right",
                br(),
                shiny::plotOutput(outputId = "heat", width = "100%"), position = "right"
            )
        )
    )


    #server output/plotting specifications
    shiny::server <- function(input, output){
        getPalette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Spectral")))
        #Expression Viewer boxplot
        output$box <- shiny::renderPlot({
            if(is(experiment)[1] == "list"){
                se_index <- which(inputOps == input$dataset)
                se <- experiment[[se_index]]}
            if(is(experiment)[1] != "list"){
                se <- experiment}
            if(is(experiment[[1]])[1] == "SummarizedExperiment" | is(experiment)[1] == "SummarizedExperiment"){
                labelData <- which(names(colData(se)) %in% metaNames)
                n <- length(unique(colData(se)[,labelData]))
                cellTypeSplit <- split(x = assay(se)[input$genename,], f = colData(se)[,labelData])
                test <- names(sort(unlist(lapply(cellTypeSplit, FUN = median))))
                colData(se)[,labelData] <- factor(colData(se)[,labelData] , levels=test)
                par(mai = c(0.8,2.3,0.2,0.2))
                boxplot(formula = assay(se)[input$genename,] ~ colData(se)[,labelData],
                        data = assay(se),
                        main = input$genename,
                        xlab = "Log Expression (TPM)",
                        ylab = "",
                        cex.axis = 0.8,
                        horizontal = TRUE,
                        xlim = c(0,n+1),
                        col = getPalette(n),
                        las = 1)
                rect(xleft = min(assay(se)[input$genename,]),
                     xright = max(assay(se)[input$genename,]),
                     ybottom = seq(-0.5, n-.5),
                     ytop = seq(0.5, n+.5),
                     col= c(NA, rgb(.78,.89,1, alpha = 0.25)),
                     border = NA)
                height = 150
            }

            if(is(experiment[[1]])[1] == "SGFeatureCounts" | is(experiment)[1] == "SGFeatureCounts"){
                labelData <- which(names(colData(se)) %in% metaNames)
                n <- length(unique(colData(se)[,labelData]))
                cellTypeSplit <- split(x = assay(se)[which(geneID(sgfc) == input$genename),], f = colData(se)[,labelData])
                test <- names(sort(unlist(lapply(cellTypeSplit, FUN = median)), decreasing = TRUE))
                colData(se)[,labelData] <- factor(colData(se)[,labelData] , levels=test)
                par(mai = c(0.8,2.3,0.2,0.2))
                boxplot(formula = t(assay(se)[which(geneID(se) == input$genename),]) ~ colData(se)[,labelData],
                        data = assay(se),
                        main = input$genename,
                        xlab = "Expression (TPM)",
                        ylab = "",
                        horizontal = TRUE,
                        cex.axis = 0.8,
                        col = getPalette(n),
                        las = 1)
                rect(xleft = min(assay(se)[which(geneID(se) == input$genename),]),
                     xright = max(assay(se)[which(geneID(se) == input$genename),]),
                     ybottom = seq(-0.5, n-.5),
                     ytop = seq(0.5, n+.5),
                     col= c(NA, rgb(.78,.89,1, alpha = 0.25)),
                     border = NA)
                height = 150
                }
            })

        #SGSeq Heatmap
        if(is(experiment)[1] == "SGFeatureCounts" | is(experiment[[1]])[1] == "SGFeatureCounts"){
        output$heat <- shiny::renderPlot({
            SGSeq::plotFeatures(experiment,
                                geneID = input$genename)
            })
        }

        #Reference link
        url <- a(paper, href = paperlink)
        output$tab <- shiny::renderUI({
            tagList("Reference Link:", url)
        })
    }
    shiny::shinyApp(ui = ui, server = server)
}



