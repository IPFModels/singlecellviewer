library(shinydashboard)
library(shiny)
library(shinyBS)
library(plotly)
library(shinyjs)
library(reshape2)
library(visNetwork)
library(shinyBS)
library(bslib)
library(dashboardthemes)

options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=600*1024^2) 
ui <- dashboardPage(
  dashboardHeader(title = "Single Cell RNASeq data",titleWidth = 350,dropdownMenuOutput("userloggedin")),
  dashboardSidebar(width = 350,
                   div(style="overflow-y: scroll"),
            
                   tags$head(tags$style(HTML(".sidebar { height: 250vh; overflow-y: auto; }
                                             .shiny-notification{position: fixed;top: 33%;left: 45%;right: 30%;}
                                             " )
                   )),
                   sidebarMenu(
                     menuItem("Project List", tabName = "dashboard", icon = icon("dashboard")),
                    uiOutput("projectlist"),
                     menuItem('Project Summary', tabName = 'summ', icon = icon('hand-o-right')),
                     menuItem('Dimension Reduction Plots', tabName = 'tplot', icon = icon('hand-o-right'),
                              menuSubItem("Compare Plots", tabName = "tsneplot"),
                              menuSubItem("Interactive DR Plots", tabName = "intertsne")
                     ),
                     menuItem('Biplot', tabName = 'biplot', icon = icon('hand-o-right')),
                     menuItem('Differential Expression', tabName = 'deg', icon = icon('hand-o-right')),
                     menuItem('Seurat Heatmap', tabName = 'heatmap', icon = icon('hand-o-right')),
                     menuItem('Gene Expression Plots', tabName = 'geplot', icon = icon('hand-o-right'),
                              fluidRow(
                                column(1,h4("")),
                                column(6,menuSubItem("Gene Expression Plots", tabName = "geplots")),
                                column(2,bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                              bsTooltip(id = "q1", title = "Enter Gene name and view tSNE, violin and Ridge plots",placement = "right",trigger = "hover", options = NULL),
                              bsTooltip(id = "q2", title = "View tSNE/umap of genes expressed in a Cellgroup",placement = "right",trigger = "hover", options = NULL),
                              fluidRow(
                                column(1,h4("")),
                                column(6,menuSubItem('Gene-Cellgroup Dotplot', tabName = 'dotplot')),
                                column(2,bsButton("q3", label = "", icon = icon("question"), style = "info", size = "extra-small"))),
                              bsTooltip(id = "q3", title = "Upload genelist to view the gene expression across cellgroups as a dotplot",placement = "right",trigger = "hover", options = NULL)
                     ),
                    menuItem('Ligand Receptor Pairs', tabName = 'ligrecmenu', icon = icon('hand-o-right'),
                              menuSubItem('Ligand Receptor Bigeneplot', tabName = 'ligrec', icon = icon('hand-o-right')),
                              menuSubItem('Ligand Receptor Network', tabName = 'network', icon = icon('hand-o-right'))
                              #menuSubItem('Ligand Receptor Heatmap', tabName = 'netheatmap', icon = icon('hand-o-right'))
                     ),
                     menuItem('Help page', tabName = 'help', icon = icon('hand-o-right'))
                     
                   )#end of sidebar menu
  ),#end dashboardSidebar
  
  dashboardBody(
    shinyDashboardThemes(
      theme = "grey_light"
    ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tags$script("Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
                Shiny.onInputChange(variableName, null);
                });
                "),
    tabItems(
      tabItem(tabName = "dashboard",
              box(
                width = 12,solidHeader = T,
                title = "scRNA Data Sets",
                uiOutput("modalDTDialog"),
                DT::dataTableOutput("datasetTable")
              )
      ),
      
      ######################################################################################################################################
      
      tabItem(tabName = "summ",
              box(
                width = 12,solidHeader = TRUE,
                title = "Project Summary",
                verbatimTextOutput("prjsumm")
              )
      ),
      ######################################################################################################################################
      tabItem(tabName = "tsneplot",
              box(title = "Compare Dimension Reduction plots",solidHeader = TRUE,width=12,
                  fluidRow(
                    column(6,uiOutput("umapa")), #Dimensionality reduction method of left plot
                    column(6,uiOutput("umapb")) #Dimensionality reduction method of left plot
                  ),
                  fluidRow(
                    column(6,selectInput("categorya2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust")),
                    column(6,selectInput("categoryb2", "Select one",c('Categories' = "var",'Cluster' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"))
                  ),
                  sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  conditionalPanel(
                    condition = "input.categorya2 == 'geneexp' || input.categoryb2 == 'geneexp'",
                    fluidRow(
                      column(6,selectInput("genecolor1", "Pick color 1 for Gene Expression Plot",colors(),selected = "grey")),
                      column(6,selectInput("genecolor2", "Pick color 2 for Gene Expression Plot",colors(),selected = "blue"))
                    )
                  ),
                  fluidRow(
                    column(6,conditionalPanel(
                      condition = "input.categorya2 == 'var'",
                      uiOutput("tsnea2") # Generate list of variables for left plot
                    ),
                    conditionalPanel(
                      condition = "input.categorya2 == 'geneexp'",uiOutput("gene1aui")
                    )
                    ),
                    column(6,conditionalPanel(
                      condition = "input.categoryb2 == 'var'",
                      uiOutput("tsneb2")  # Generate list of variables for right plot
                    ),
                    conditionalPanel(
                      condition = "input.categoryb2 == 'geneexp'",uiOutput("gene2aui")
                    )
                    )),
                  fluidRow(
                    column(6,checkboxInput("subsa", label = "Check to subselect cells", value = FALSE)),
                    column(6,checkboxInput("subsb", label = "Check to subselect cells", value = FALSE))
                  ),
                  fluidRow(
                    column(6,checkboxInput("checklabel1", label = "Check for cell  group labelling", value = TRUE)),
                    column(6,checkboxInput("checklabel2", label = "Check for cell  group labelling", value = TRUE))
                  ),
                  fluidRow(
                    column(6,conditionalPanel(
                      condition = "input.subsa ==true",uiOutput("subsaui") #generate ident list for left plot
                    )),
                    column(6,conditionalPanel(
                      condition = "input.subsb ==true",uiOutput("subsbui") #generate ident list for right plot
                    )
                    ))),
              box(title = "Plots",solidHeader = TRUE,width=12,
                  plotOutput("comptsne2", height = 600),
                  downloadButton('downloadtsneplot', 'Download Dimension Reduction plot')
              )
      ),#end of degtab
      
      ###################################################################################################################################### 
      tabItem(tabName = "intertsne",
              fluidRow(
                box(
                  title = "Controls",solidHeader = TRUE,width=12,
                  fluidRow(
                    column(6,uiOutput("umapint")),
                    column(6,selectInput("intercat", "Select one",c('Categories' = "var",'Gene Expression' = "geneexp"),selected = "geneexp"))),
                  fluidRow(
                    column(6,uiOutput("setcategory")),
                    column(6,conditionalPanel(
                      condition = "input.intercat == 'var'",
                      uiOutput("intervar")
                    ),
                    conditionalPanel(
                      condition = "input.intercat == 'geneexp'",uiOutput("geneinterui")
                    )
                    )),
                  sliderInput("umap_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25)
                ),
                box(plotlyOutput("intertsne", height = 600),width=6,title = "Interactive tSNE/uMap Plot",solidHeader = TRUE),
                box(plotlyOutput("intergene", height = 600),width=6,title = "Interactive Gene Plot",solidHeader = TRUE)
                
              )
      ),#endbigeneplotTab
      ###################################################################################################################################### 
      
      tabItem(tabName = "biplot",
              fluidRow(
                box(plotOutput("bigeneplot", height = 600),width=8,title = "Bigene Plot",solidHeader = TRUE),
                
                box(
                  title = "Controls",solidHeader = TRUE,width=4,
                  uiOutput("bigenedim"),
                  uiOutput("bigene_geneaui"),
                  uiOutput("bigene_rangea"),
                  uiOutput("bigene_genebui"),
                  uiOutput("bigene_rangeb"),
                  sliderInput("bigene_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  downloadButton('downloadbiplot', 'Download Bigene plot')
                )
              ),
              box(DT::dataTableOutput('bigene_genecount'),width=12,title = "Gene Count Table",solidHeader = TRUE)
      ),#endbigeneplotTab
      
      ######################################################################################################################################
      tabItem(tabName = "deg",
              box(title = "Differentially expressed Markers",solidHeader = TRUE,width=9,
                  plotOutput("comptsne", height = 1000)
              ),
              
              fluidRow(
                box(title = "Controls",solidHeader = TRUE,width=3,
                    uiOutput("umapdeg"),
                    uiOutput("tsnea"),
                    
                    sliderInput("pointa", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                    checkboxInput("checklabel3", label = "Check for cell  group labelling", value = TRUE),
                    checkboxInput("checkviolin", label = "Check to remove points from violin plot", value = TRUE),
                    hr(),
                    uiOutput("identdef"),
                    checkboxInput("setident", label = "Check to choose a different category to compare", value = FALSE),
                    conditionalPanel(
                      condition = "input.setident ==true",uiOutput("setidentlist"),
                      uiOutput("identa"),
                      uiOutput("identb"),
                      actionButton("goButton", "Go!"),
                      sliderInput("lfc", "Log FC threshold:",min = 0.25, max = 6, value = 0.25,step=.25),
                      selectInput("test", "Select test to use",c('Wilcox' = "wilcox",'T-test' = "t", 'Poisson' = "poisson",'Negative Binomial'="negbinom"),selected = "wilcox"),
                      sliderInput("minpct", "Minimum Percent of cells:",min = 0.1, max = 10, value = 0.25)),
                    downloadButton('downloaddeg', 'Download table'),
                    downloadButton('downloadplot', 'Download Plot')
                ),
                box(DT::dataTableOutput('markergenes'),width=12,solidHeader = TRUE,title="Marker genes")
              )#End FluidRow
      ),#end of degtab
      ######################################################################################################################################
      tabItem(tabName = "heatmap",
              box(
                title = "Controls",solidHeader = TRUE,width=12,
                selectInput("shmptype", "Genes to use as input",c('Top genes across all clusters' = "topgene", 'Differentially Expressed genes' = "deggene",'Tanscription factors' = "tf",'Ligands' = "lig",'Receptors' = "rec"),selected = "deggene"),
                conditionalPanel(
                  condition = "input.shmptype =='topgene'",sliderInput("topn", "Top n genes",min = 5, max = 15, value = 2,step=1)),
                conditionalPanel(
                  condition = "input.shmptype =='deggene'",
                  fluidRow(
                    column(6,uiOutput("heatmapclust")),
                    column(6,uiOutput("heatmapgenes")))),
                uiOutput("hmpgrp"),
                #selectInput("hmpcol", "Select one",c('PurpleYellow' = "PuYl",'BlueGreen' = "BuGn", 'RedYellow' = "RdYl", 'RedBlue'="RdBu"),selected = "geneexp"),
                downloadButton('downloadheatmap', 'Download Heatmap')
              ),
              box(plotOutput("heatmap", height = 900),width=12,solidHeader = TRUE,title="Single cell heatmap of gene expression")
      ),#end of tab
      ######################################################################################################################################
      tabItem(tabName = "geplots",
              box(title = "Gene Expression Plots",solidHeader = TRUE,width=9,
                  plotOutput("geplots", height = 1000)
              ),
              
              fluidRow(
                box(title = "Controls",solidHeader = TRUE,width=3,
                    uiOutput("umapge"),
                    uiOutput("geneidui"),
                    checkboxInput("checkviolin2", label = "Check to remove points from violin plot", value = TRUE),
                    sliderInput("genenid_pointsize", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                    downloadButton('downloadplotge', 'Download Plot'))
              )#End FluidRow
      ),#end of geplot
      
      ######################################################################################################################################
      tabItem(tabName = "dotplot",
              box(title = "Controls",solidHeader = TRUE,width=12,
                  radioButtons("radiofileup", label = "File input type",choices = list("Enter gene names" = "enter", "Upload List" = "upload"),selected = "upload"),
                  fluidRow(
                    column(6,
                           conditionalPanel(
                             condition = "input.radiofileup =='upload'",fileInput('genelistfile', 'Upload Text File',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))),
                           conditionalPanel(
                             condition = "input.radiofileup =='enter'",uiOutput("enterchoice"))
                    ),
                    
                    column(6,uiOutput("setdotvar"))
                  )
              ),
              box(title = "Dot Plot",solidHeader = TRUE,width=12,
                  plotOutput("dotplot", height = 500),
                  downloadButton('downloaddotplot', 'Download Plot')
              )
      ),#end of dotplot
      
      ######################################################################################################################################
      tabItem(tabName = "ligrec",
              box(width=9, title = "Bigene Plot",solidHeader = TRUE,
                  plotOutput("bigeneplot2", height = 800)
              ),
              box(width = 3, solidHeader = TRUE,title = "Controls",
                  uiOutput("bigenedimr"),
                  uiOutput("pairby"),
                  radioButtons("clust","Select Cluster", c("All clusters"="all","Select Cluster"="clust"),selected = "all"),
                  radioButtons("gene","Select Genes", c("All genes"="allgene","Enter Genelist"="genelist"),selected = "allgene"),
                  
                  conditionalPanel(
                    condition = "input.clust == 'all' && input.gene == 'genelist'" ,
                    fileInput('genelist1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                    fileInput('genelist2', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt'))
                  ),
                  conditionalPanel(
                    condition = "input.clust == 'clust' && input.gene == 'allgene'" ,
                    uiOutput("clust1"),
                    uiOutput("clust2")
                  ),
                  conditionalPanel(
                    condition = "input.clust == 'clust' && input.gene == 'genelist'" ,
                    fileInput('genelist1.1', 'Upload Receptor Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                    fileInput('genelist2.1', 'Upload Ligand Genelist',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                    uiOutput("clust1.1"),
                    uiOutput("clust2.1")
                  ),
                  fluidRow(
                    column(6,checkboxInput("checksource", label = "Check to select by source", value = TRUE)),
                    column(6,checkboxInput("checkevi", label = "Check to select by evidence", value = TRUE)),
                    conditionalPanel(
                      condition = "input.checksource ==true",
                      column(6,uiOutput('source'))
                    ),
                    conditionalPanel(
                      condition = "input.checkevi ==true",
                      column(6,uiOutput('evidence'))
                    )
                  ),
                  uiOutput("bigene_rangea2"),
                  uiOutput("bigene_rangeb2"),
                  #sliderInput("perc_cells", "Filter by percentage of cells expressing the gene:",min = 10, max = 100, value = 50,step=5),
                  sliderInput("bigene_pointsize2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                  actionButton("lrpgo", "Change Parameters and Run"),br(),br(),
                  downloadButton('downloadlr', 'Download Results'),br(),br(),
                  downloadButton('downloadlrbigene', 'Download Plot')
              ),
              box(
                width = 12, solidHeader = TRUE,
                title = "Ligand Receptor pairs",
                DT::dataTableOutput('pairs_res')
              )#end of box
              
              
      ),#end of tabitem
      ######################################################################################################################################
      tabItem(tabName = "network",
              box(title = "Network",solidHeader = TRUE,width=8,
                  #visNetworkOutput("lrnetwork", height = 800)
                  plotOutput("lrnetwork", height = 700)
              ),#End box
              box(title = "Controls",solidHeader = TRUE,width=4,
                  uiOutput("pairbynet"),
                  #sliderInput("perc_cells2", "Filter by percentage of cells expressing the gene:",min = 10, max = 100, value = 50,step=10),
                  fluidRow(
                    column(6,checkboxInput("checksource2", label = "Check to select by source", value = TRUE)),
                    column(6,checkboxInput("checkevi2", label = "Check to select by evidence", value = TRUE)),
                    conditionalPanel(
                      condition = "input.checksource2 ==true",
                      column(6,uiOutput('source2'))
                    ),
                    conditionalPanel(
                      condition = "input.checkevi2 ==true",
                      column(6,uiOutput('evidence2'))
                    )
                  ),
                  checkboxInput("checkgrp", label = "Check to select Group(s)", value = FALSE),
                  conditionalPanel(
                    condition = "input.checkgrp ==true",
                    uiOutput('checkgrp')
                  ),
                  actionButton("lrnset", "Set filters"),
                  uiOutput("filternet"),
                  fluidRow(
                    column(6,actionButton("lrngo", "Run")),
                    column(6,downloadButton('dwldnet', 'Download Network plot')))
              ),#End box
              box(title= 'GO Analysis', solidHeader = T, width=12, 
                  uiOutput('grp1')
              ),
              box(title = "GO Terms",solidHeader = TRUE,width=12,
                  DT::dataTableOutput('gotable')
              ),
              box(title = "Ligand-Receptor pairs",solidHeader = TRUE,width=12,
                  DT::dataTableOutput('pairs_res2')
              )),#end of network
      ######################################################################################################################################
      # tabItem(tabName = "netheatmap",
      #         box(title = "Heatmap",solidHeader = TRUE,width=8,
      #             plotOutput("netheatmap")
      #         ),
      #         box(title = "Controls",solidHeader = TRUE,width=4,
      #             uiOutput("pairbyheatnet"),
      #             checkboxInput("checksourceheat", label = "Check to select by source", value = TRUE),
      #             checkboxInput("checkeviheat", label = "Check to select by evidence", value = TRUE),
      #             conditionalPanel(
      #               condition = "input.checksourceheat ==true",
      #               uiOutput('source3')
      #             ),
      #             conditionalPanel(
      #               condition = "input.checkeviheat ==true",
      #               uiOutput('evidence3')
      #             ),
      #             #sliderInput("perc_cells3", "Filter by percentage of cells expressing the gene:",min = 10, max = 100, value = 50,step=10),
      #             selectInput("hmpcolnet", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")),
      #             selectInput("clusterby", "Cluster By",c('Both'="both",'Receptor Genes' = "column",'Ligand Genes' = "row",'None' = "none")),
      #             checkboxInput("checkbox", label = "Reverse Colors", value = FALSE),
      #             actionButton("lrhgo", "Change Parameters and Run"),br(),br(),
      #             downloadButton('downloadlrheatmap', 'Download Heatmap')
      #         ),
      #         box(title = "Ligand Receptor Pairs",solidHeader = TRUE,width=12,status='primary',
      #             DT::dataTableOutput('pairs_res3')
      #         )#End box
      # ),
    
       ######################################################################################################################################
      tabItem(tabName = "help",
              h4("All datasets hosted in this website are analysed using the Seurat package developed by the Satija Lab in NYGC."), 
              h4("For more information", a("click here", href="https://satijalab.org/seurat/")),
              br(),
              h4(p(strong("1. Project Summary"))),
              h4(p(div("The Project Summary tab gives a short description about the selected dataset and basic information such as organism, total number of cells, number of cells in each cluster, total number of genes and number of dimensions used in the analysis"))),
              br(),
              h4(p(strong("2. tSNE Plots"))),
              h4(p(div("The ",em("Compare tSNE plots")," displays two plots with same options. Users can choose the Dimensionality reduction method as well as pick between categories, clusters and gene expression. There are additional options to change colors of gene expression plot (also called feature plot) as well as subselect cells based on clusters or a range of numeric values"))),
              h4(p(div("The ",em("Interactive tSNE/uMAP Plot"),"displays the same as in the ",em("Compare tSNE plots"),"tab except these plots are interactive. Users can zoom in, subselect and hover over for point information"))),
              br(),
              h4(p(strong("3. Biplot"))),
              h4(p(div("The Biplot is similar to seurat's plot with the 'overlap' option. Users specify dimensionality reduction method, 2 genes and the expression limit of each gene in terms of logUMI to view a bigene plot based on the expression of the two genes. The tab also shows a table with the number of cells in each Cell group expression the 2 genes"))),
              br(),
              h4(p(strong("4. Differential Expression"))),
              h4(p(div("The differential expression tab generates a table of marker genes specific to a cell group as well as plots showing the expression of each of those marker genes in each cell group.By default, the results displayed are each cell group again every other cell group. User can, however, find marker genes by checking the ",em("Check to choose a different category to compare")," option and selecting one cell group that can be compared against one or multiple other groups  "))),
              br(),
              h4(p(strong("5. Seurat heatmap"))),
              h4(p(div("The Seurat heatmap tap generates an expression heatmap for given cells and marker genes from the table in ",em("Differential Expression")," tab "))),
              br(),
              h4(p(strong("6. Gene Expression plots"))),
              h4(p(div("The",em("Gene Expression plots")," displays the same 3 expression, violin and ridge plots as in the ",em("Differential Expression")," tab. However, instead of having to select from a table, in this tab, users can manually specify the gene name they are interested in"))),
              h4(p(div("The ",em("Cluster-wise Gene Expression")," tab shows a general dimension reduction plot as well as a gene expression plot along with a table of list of genes present in the data ordered by the percent of cells they are present in and their average expression across these cells"))),
              h4(p(div("The ",em("Gene-Cellgroup Dotplot")," tab displays a dotplot like the one", a("here", href="https://satijalab.org/seurat/immune_alignment.html"),"where the genes are rows, cell-groups are the columns and the size of the dots show percentage of cells expressing that gene. Users please upload text files with one gene per line"))),
              br(),
              h4(p(strong("7. Ligand receptor Pairs"))), 
              h4(p(div("The ",em("Ligand Receptor Pairs")," tab generates a table of all possible ligand receptor pairs (from", a("FANTOM5", href="http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/data/PairsLigRec.txt"),") between each of the cell groups. Expression of each receptor/ligand gene per cell group is counted if it is expressed in at least 20% of the cells (by default. Users can change this and go upto 100%). This tab also displays a bigene plot for each of the ligand-receptor pair in the table"))),
              h4(p(div("The ",em("Ligand Receptor Network")," tab displays a network where each cell group is a node and the edges denote the frequency of ligand-receptor pairs between the nodes. This tab also performs GO analysis on the selected pairs and displays the results in a table"))),
              #h4(p(div("The ",em("Ligand Receptor Heatmap")," tab displays the heatmap where rows are ligand cell groups and columns are receptor cell groups and the data represents the number of ligand-receptor pairs between the groups"))),
              br(),
              h4("For unknown errors or support, email bapoorva@pennmedicine.upenn.edu") 
      )
    )#end of tabitems
  )#end of dashboard body
)#end of dashboard page

