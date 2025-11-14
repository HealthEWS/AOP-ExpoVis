library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(shinythemes)
library(heatmaply)
library(shinyjs)
library(shinyjqui)
library(tidyverse)
library(tippy)
library(vroom)
library(dplyr)
library(GOSemSim)
library(GOSim)
library(factoextra)
library(pheatmap)
library(ggplot2)
library(ggcorrplot)
library(RColorBrewer)
library(org.Hs.eg.db)
library(OrganismDbi)
library(gmp)

#### Part 1. Utility Functions ####
calc_p1 = function(N, m, nx, ny) {
  cat("N = ", N, "\tm = ",m,"\tnx = ", nx, "\tny = ", ny, "\n", sep = "")
  if(is.na(m) || is.na(nx) || is.na(ny)) {
    return(0)
  }
  result = chooseZ(N, m) * chooseZ(N - m, nx - m) * chooseZ(N - nx, ny - m) / ( chooseZ(N, nx) * chooseZ(N, ny))
  return(as.numeric(result))
}

#### Part 2. Data Loading ####
bp_phenotype_disease = readRDS("./data/bp_phenotype_disease.rds")
disease = readRDS("./data/disease.rds")
ChemicalName_df = disease %>% dplyr::select(1,2) %>% dplyr::distinct() %>% arrange(ChemicalName)
ChemicalName_display = paste(ChemicalName_df[[2]], ChemicalName_df[[1]], sep = " ")
ChemicalName_value = ChemicalName_df[[1]]
ChemicalName_value_pick = ChemicalName_value
names(ChemicalName_value_pick) = ChemicalName_display

chemical_selected <- c("2,2',4,4'-tetrabromodiphenyl ether",  
                       "2,2',4,4'-tetrabromodiphenyl ether",
                       "2,2',4,4'-tetrabromodiphenyl ether",
                       "2,2,4,4-tetrabromodiphenyl ether",
                       "PBDE-47",
                       "BDE-47",
                       "TBDP-ether",   
                       "tetrabrominated diphenyl ether 47", 
                       "BDE 47",
                       "BDE47",
                       "PBDE47",
                       "5436-43-1",
                       "PBDE 47")

CTD_disease = readRDS("./data/CTD_disease.rds")
chemicalselect_phenotype = readRDS("./data/bde_phenotype.rds")

#### Part 3. Gene and Chemical Inference Data ####
gene_phenotype = readRDS("./data/gene_phenotype.rds")
array<-table(gene_phenotype$Gene.Symbol)
gene_degree1_phenotype<-data.frame(genesymbol = as.character(names(array)),degree=as.numeric(array))
gene_degree2_disease = readRDS("./data/gene_degree2_disease.rds")

chemical_disease_curated = readRDS("./data/chemical_disease_curated.rds")
chemical_phenotype = readRDS("./data/chemical_phenotype.rds")
array1<-table(unique(chemical_phenotype$chemicalname))
chemical_degree1_phenotype<-data.frame(chemicalname = as.character(names(array1)),degree=as.numeric(array1))
array2<-table(chemical_disease_curated$ChemicalName)
chemical_degree2_disease<-data.frame(chemicalname = as.character(names(array2)),degree=as.numeric(array2))

hsGO = readRDS(file = "./data/hsGO.rds")

#### Part 4. AOP Data Loading ####
aop_ke_mie_ao<-read.csv(file="./aop_ke_mie_ao.tsv",sep = "\t",header = FALSE)
colnames(aop_ke_mie_ao)<-c("aop id", "key event id", "key event type", "key event name")

aop_ker<-read.csv(file="./aop_ke_ker.tsv",sep = "\t",header = FALSE)
colnames(aop_ker)<-c("aop id", "upstream event id", "downstream event id", "relationship id", "direct or indirect relationship", "evidence for relationship", "quantitative understanding of relationship")

aop_ec<-read.csv(file="./aop_ke_ec.tsv",sep = "\t",header = FALSE)
colnames(aop_ec)<-c("aop id", "key event id", "action", "object source", "object ontology id", "object term", "process source", "process ontology id", "process term")

#### Part 5. UI Definition ####
ui <- fluidPage( 
  theme = shinytheme("flatly"),
  titlePanel('AOP-ExpoVis'), 
  useShinydashboard(),
  use_tippy(),
  sidebarLayout(
    sidebarPanel(width = 4
                 ,multiInput(
                   inputId = "BDE47",
                   label = "Chemical Names:", 
                   selected =chemical_selected,
                   choiceNames = ChemicalName_display, 
                   choiceValues = ChemicalName_value,
                   options = list(
                     enable_search = TRUE,
                     non_selected_header = "Choose between:",
                     selected_header = "You have selected:"
                   ),
                   width = "100%"
                 ),
                 textOutput("part1_message"),
                 br(),
                 selectInput(inputId = "diseaseclass", label = "Disease Class:",multiple = FALSE, choices = NULL),
                 textOutput("part3_text1"),
                 br(),
                 selectInput(inputId = "Diseasecategory", label = "Diseasecategory",multiple = T, choices = NULL),
                 textOutput("part4_text1"),
                 br(),
                 actionButton("Start", "Refresh All...", icon = icon("rotate")),
                 p("Welcome to AOP-ExpoVis, an interactive web application designed to explore relationships between chemicals, diseases, and phenotypes. If you use AOP-ExpoVis in your research, please cite our work as follows: Author(s), Title, Journal, Year, DOI."),
    ),
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel("Chemical", 
                         DT::dataTableOutput("part1_table") %>% shinycssloaders::withSpinner(type = 5),
                         downloadButton("download_part1_table", "Download..."),
                         div(style = "width: 100%; margin-bottom: 100px;"),
                ),
                tabPanel("Disease", 
                         DT::dataTableOutput("part3_table") %>% shinycssloaders::withSpinner(type = 5),
                         downloadButton("download_part3_table", "Download..."),
                         div(style = "width: 100%; margin-bottom: 100px;"),
                ),  
                tabPanel("Phenotype", 
                         uiOutput("title_part4_table1"),
                         DT::dataTableOutput("part4_table1") %>% shinycssloaders::withSpinner(type = 5),
                         downloadButton("download_part4_table1", "Download..."),
                         hr(),
                         uiOutput("title_part4_table2"),
                         DT::dataTableOutput("part4_table2") %>% shinycssloaders::withSpinner(type = 5),
                         hr(),
                         uiOutput("title_part4_table3"),
                         DT::dataTableOutput("part4_table3") %>% shinycssloaders::withSpinner(type = 5),
                         div(style = "width: 100%; margin-bottom: 100px;"),
                ),  
                tabPanel("Heatmap", 
                         fluidRow(
                           fluidRow(
                             column(width = 4, valueBoxOutput("part5_text1", width = NULL)),   
                             column(width = 4, valueBoxOutput("part5_text2", width = NULL)),   
                           )
                         ),
                         tabsetPanel(
                           tabPanel("Inferred Phenotype Coherence", 
                                    DT::dataTableOutput("part5_table") %>% shinycssloaders::withSpinner(type = 5),
                                    downloadButton("download_part5_table", "Download..."),
                                    hr(),
                                    plotlyOutput("part5_heatfigure_plot", height = "600px") %>% jqui_resizable() %>% shinycssloaders::withSpinner(type = 5),
                                    div(style = "width: 100%; margin-bottom: 500px;"),
                           ),
                           tabPanel("Phenotype Cluster Coherence", 
                                    DT::dataTableOutput("part5_table_heatfigure") %>% shinycssloaders::withSpinner(type = 5),
                                    downloadButton("download_part5_table_heatfigure", "Download..."),
                                    hr(),
                                    fluidRow(
                                      column(width = 4,        
                                             sliderInput("thresholdvalue", HTML('Please select the exclusion ratio: <span id="tips_thresholdvalue"><i class="fa fa-info-circle" role="presentation"></i></span>'), value = 0.5, min = 0.01, max = 0.99, step = 0.01, width = "100%"), 
                                             tippy_this("tips_thresholdvalue", "<p style='font-size: 1.5em;'>The proportion chosen with the slider determines the percentage of data points to be excluded from the subsequent analysis, ranging from 0.01 to 0.99. For example, setting the exclusion ratio to 0.9 means that 90% of the data points will be excluded, and the analysis will include the remaining 10%.</p>", allowHTML = TRUE, placement = 'top'),
                                      ),
                                      column(width = 4, actionButton("refresh_heatmap", "Show/Refresh heatmap...", icon = icon("rotate")),),
                                      column(width = 4, valueBoxOutput("part5_hopkins_value", width = NULL),)
                                    ),
                                    div(style = "width: 100%; overflow-x: auto; overflow-y: auto; margin-right: 5px; margin-bottom: 5px;",
                                        plotlyOutput("part5_heatmap_plot2", height = "700px") %>% jqui_resizable()  %>% shinycssloaders::withSpinner(type = 5),
                                    ),
                           ),
                           tabPanel("AOP", 
                                    searchInput(inputId = "key_event_name_opt", label = "Please input key words of event name to search:", placeholder = "keywords of event name, separated by whitespace", btnSearch = icon("magnifying-glass"), btnReset = icon("xmark"), width = "100%"),
                                    multiInput(inputId = "aop_id", label = "Please select AOP IDs:", choices = c("NA"), width = "100%"),
                                    hr(),
                                    DT::dataTableOutput("aop_key_event_male_aop") %>% shinycssloaders::withSpinner(type = 5),
                                    downloadButton("download_aop_event_male", "Download..."),
                                    hr(),
                                    DT::dataTableOutput("aop_happy") %>% shinycssloaders::withSpinner(type = 5),
                                    downloadButton("download_aop_happy", "Download..."),
                                    div(style = "width: 100%; margin-bottom: 100px;"),
                           ),
                         )
                )
              )
    )
  )
)

#### Part 6. Server Logic ####
server <- function(input, output, session) {
  
  #### Module 1. Chemical Data Processing ####
  chemical_disease <- reactive({
    print(input$BDE47)
    if(is.null(input$BDE47)) {
      disease
    } else {
      disease %>% dplyr::filter(ChemicalName %in% input$BDE47)
    }
  })
  
  chemical1_phenotype <- reactive({
    if(is.null(input$BDE47)) {
      chemicalselect_phenotype
    } else {
      chemicalselect_phenotype %>% dplyr::filter(chemicalname %in% input$BDE47)
    }
  })
  
  output$part1_message <- renderText({
    phen_num <- length(unique(chemical1_phenotype()$phenotypeid))
    Ref_num <- length(unique(chemical1_phenotype()$pubmedids))
    paste0("Based on the selected chemical compounds, an initial match identified ", phen_num, " phenotypes involving ", Ref_num, " studies.")
  })
  
  output$part1_table <- DT::renderDataTable({
    DT::datatable(chemical1_phenotype(), selection = "single", class = "display nowrap", options = list(pageLength = 20, scrollX = TRUE))
  }) %>% bindEvent(input$refresh)
  
  output$download_part1_table <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(chemical1_phenotype(), file, row.names = F) }
  )
  
  #### Module 2. Disease Classification ####
  slimmesh_info <- reactive({
    CTD_disease[match(chemical_disease()$DiseaseID,table = CTD_disease$DiseaseID),]
  })
  
  disease_slimmesh <- reactive({
    data.frame(Disease.Name = chemical_disease()$diseasename,SlimMappings=slimmesh_info()$SlimMappings)
  })
  
  observe({
    slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE ))) %>% sort()
    freezeReactiveValue(input, "SlimMapping")
    updateSelectInput(inputId = "SlimMapping", label = "SlimMapping", choices = slimmeshcategory)
  })
  
  nervous_diease <- reactive({
    req(input$diseaseclass)
    CTD_disease[grepl(input$diseaseclass, CTD_disease$SlimMappings), ]
  })
  
  output$part2_table_slimmesh_info <- DT::renderDataTable({
    DT::datatable(slimmesh_info(), selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE)) %>%
      DT::formatStyle(columns = c("Definition"), textOverflow = 'ellipsis', overflow = 'hidden', whiteSpace = 'nowrap')
  })
  
  output$part2_table_nervous_diease <- DT::renderDataTable({
    DT::datatable(nervous_diease(), selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE))
  })
  
  output$download_part2_table1 <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(slimmesh_info(), file, row.names = F) }
  )
  
  output$download_part2_table2 <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(nervous_diease(), file, row.names = F) }
  )
  
  #### Module 3. Disease Class Management ####
  observe({
    req(slimmesh_info())
    req(disease_slimmesh())
    if (nrow(disease_slimmesh()) > 1) {
      slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE ))) %>% sort()
      disclass <- function(x){
        result<-list()
        for(i in 1:length(x)){
          slim=x[i]
          index=grepl(as.character(slim),disease_slimmesh()$SlimMappings,fixed = TRUE)
          disease= unique(subset(disease_slimmesh()$Disease.Name,index))
          result[[i]] = disease
        }
        return(result)
      }
      diseaseclass<-disclass(slimmeshcategory)
      names(diseaseclass)<-slimmeshcategory
      print(names(diseaseclass))
      freezeReactiveValue(input, "diseaseclass")
      updateSelectInput(inputId = "diseaseclass", choices = names(diseaseclass))
    }
  })
  
  diseaseclass <- reactive({
    req(slimmesh_info())
    req(disease_slimmesh())
    slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE )))
    disclass <- function(x){
      result<-list()
      for(i in 1:length(x)){
        slim=x[i]
        index=grepl(as.character(slim),disease_slimmesh()$SlimMappings,fixed = TRUE)
        disease= unique(subset(disease_slimmesh()$Disease.Name,index))
        result[[i]] = disease
      }
      return(result)
    }
    diseaseclass<-disclass(slimmeshcategory)
    names(diseaseclass)<-slimmeshcategory
    diseaseclass
  })
  
  output$part3_text1 <- renderText({
    slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE )))
    diseaseclass_select_num <- map(diseaseclass()[input$diseaseclass],length) %>% unlist() %>% sum()
    print(diseaseclass_select_num)
    paste0( "The selected chemical compounds are associated with ",length(slimmeshcategory)," categories of diseases, including ",length(diseaseclass()[[input$diseaseclass]])," instances of ",input$diseaseclass, ".")
  })
  
  output$part3_table = DT::renderDataTable({
    df = data.frame(diseaseclass()[input$diseaseclass])
    DT::datatable(df, selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE))
  }) %>% bindEvent(input$refresh)
  
  #### Module 4. Phenotype Analysis ####
  observe({
    req(input$diseaseclass)
    freezeReactiveValue(input, "Diseasecategory")
    choices = diseaseclass()[[input$diseaseclass]] 
    label = paste0(input$diseaseclass, " ", "Category:")
    updateSelectInput(inputId = "Diseasecategory", label = label, choices = choices, selected = head(choices, 3))
  })
  
  abnormal_anatomy<-c("HL-60 Cells","Cell Line,Transformed","Hep-G2","Hela Cells","K562","MCF-T Cells","HCT 116","A549","Jurkat Cells","Tumor","HT29","Caco-2 Cells","CHO","PC12","Hep G2 Cells","MCF-7 Cells")
  
  chemical1_phenotype_normal = reactive({  
    abanatomyserch<-apply(as.data.frame(abnormal_anatomy),1,grepl,chemical1_phenotype()$anatomyterms,ignore.case = TRUE)
    abanatomyindex<-apply(abanatomyserch,1,sum)
    sum(abanatomyindex)
    chemical1_phenotype_normal<-chemical1_phenotype()[!as.logical(abanatomyindex),]
    print("chemical1_phenotype_normal")
    chemical1_phenotype_normal
  })
  
  output$title_part4_table1 = renderUI({ tagList( h3("Chemical Phenotype") ) })
  
  output$part4_table1 = DT::renderDataTable({
    DT::datatable( chemical1_phenotype_normal(), selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE))
  }) %>% bindEvent(input$refresh)
  
  bp_phenotype_nervous = reactive({
    print("bp_phenotype_nervous")
    req(input$Diseasecategory)
    bp_phenotype_nervous <-bp_phenotype_disease[bp_phenotype_disease$DiseaseName %in% input$Diseasecategory,]
    bp_phenotype_nervous
  })
  
  output$title_part4_table2 = renderUI({ tagList( h3(input$diseaseclass) ) })
  
  output$part4_table2 = DT::renderDataTable({
    DT::datatable(bp_phenotype_nervous(), selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE))
  }) %>% bindEvent(input$refresh)
  
  chemical_disease_phenotype = reactive({
    print("chemical_disease_phenotype")
    chemical_disease_phenotype <-chemical1_phenotype_normal()[chemical1_phenotype_normal()$phenotypename %in% bp_phenotype_nervous()$GOName,]
    chemical_disease_phenotype
  })
  
  output$title_part4_table3 = renderUI({ tagList( h3(paste0(input$diseaseclass, " ", "Phenotype")) ) })
  
  output$part4_table3 = DT::renderDataTable({
    DT::datatable( chemical_disease_phenotype(), selection = "single", class = "display nowrap", options = list(pageLength = 10,scrollX = TRUE))
  }) %>% bindEvent(input$refresh)
  
  output$part4_text1 <- renderText({
    print(length(unique(chemical1_phenotype_normal()$phenotypeid)))
    length(unique(chemical1_phenotype_normal()$pubmedids))
    phenotype_normal = chemical1_phenotype_normal()
    phenotype_selected = bp_phenotype_nervous()
    chemicalselect_phenotype = chemical_disease_phenotype()
    paste0("In normal tissues or cells, chemicals are linked to ",length(unique(phenotype_normal$phenotypeid)), " phenotypes (",length(unique(phenotype_normal$pubmedids)), " studies). Selected diseases relate to ", length(unique(phenotype_selected$GOID)), " phenotypes, intersecting with them to yield ", length(unique(chemicalselect_phenotype$phenotypename)), " overlapping phenotypes.")
  })
  
  #### Module 5. Inference Analysis ####
  record = reactiveValues(
    N_chemical = 0, phenotype_chemical_nx = NULL, disease_chemical_ny = NULL,
    N_gene = 0, phenotype_gene_nx = NULL, disease_gene_ny = NULL,
    m_degree = NULL, m_degree_chemical = NULL,
  )
  
  observe({
    req(input$Diseasecategory)
    req(nervous_diease())
    CTD_disease_genes_curated <- gene_degree2_disease[gene_degree2_disease$DiseaseName %in% input$Diseasecategory,]
    gene_degree2_disease <- gene_degree2_disease[gene_degree2_disease$DiseaseID %in% nervous_diease()$DiseaseID,]
    array2<-table(gene_degree2_disease$GeneSymbol)
    gene_degree2_disease<-data.frame(genesymbol = as.character(names(array2)),degree=as.numeric(array2))
    colnames(gene_degree1_phenotype)<-c("genesymbol","degree1")
    colnames(gene_degree2_disease)<-c("genesymbol","degree2")
    m_degree<-merge(gene_degree1_phenotype[gene_degree1_phenotype$genesymbol %in%  CTD_disease_genes_curated$GeneSymbol,],
                    gene_degree2_disease[gene_degree2_disease$genesymbol %in% CTD_disease_genes_curated$GeneSymbol,],
                    by.x = "genesymbol",by.y = "genesymbol",all = TRUE)
    m_degree[is.na(m_degree)]<-0
    m_degree$ni <- m_degree$degree1 + m_degree$degree2
    record$m_degree = m_degree
    phenotype_degree<-table(gene_phenotype$GO.Term.ID)
    phenotype_gene_nx<-data.frame(phenotypeid=names(phenotype_degree),nx=as.numeric(phenotype_degree))
    record$phenotype_gene_nx = phenotype_gene_nx
    disease_degree<-table(CTD_disease_genes_curated$DiseaseID)
    disease_gene_ny<-data.frame(diseassid=names(disease_degree),ny=as.numeric(disease_degree))
    record$disease_gene_ny = disease_gene_ny
    N_gene = length(unique(c(as.character(gene_degree2_disease$genesymbo),as.character(gene_degree1_phenotype$genesymbol))))
    print(N_gene)
    record$N_gene = N_gene
  })
  
  output$part5_text1 <- renderValueBox({
    valueBox(value = record$N_gene, subtitle = "Gene Number", icon = icon("list-alt"))
  })
  
  CTD_disease_chems_curated = reactive({
    req(input$Diseasecategory)
    chemical_disease_curated[chemical_disease_curated$diseasename %in% input$Diseasecategory,]
  }) 
  
  observe({
    req(input$Diseasecategory)
    req(CTD_disease_chems_curated())
    m_degree_chemical<-merge(chemical_degree1_phenotype[chemical_degree1_phenotype$chemicalname %in%  CTD_disease_chems_curated()$ChemicalName,],
                             chemical_degree2_disease[chemical_degree2_disease$chemicalname %in% CTD_disease_chems_curated()$ChemicalName,],
                             by.x = "chemicalname",by.y = "chemicalname",all = TRUE)
    m_degree_chemical[is.na(m_degree_chemical)]<-0
    m_degree_chemical$ni <- m_degree_chemical$degree.x + m_degree_chemical$degree.y
    record$m_degree_chemical = m_degree_chemical
    chemical_phenotype_unique<-chemical_phenotype[!duplicated(chemical_phenotype[,c("chemicalid","phenotypeid")]),]
    phenotype_degree_chemical<-table(chemical_phenotype_unique$phenotypeid)
    phenotype_chemical_nx<-data.frame(phenotypeid=names(phenotype_degree_chemical),nx=as.numeric(phenotype_degree_chemical))
    record$phenotype_chemical_nx = phenotype_chemical_nx
    disease_degree_chemical<-chemical_disease_curated[chemical_disease_curated$diseasename %in% input$Diseasecategory,]
    disease_degree_chemical<-table(factor(disease_degree_chemical$diseasename))
    disease_chemical_ny<-data.frame(diseassid=names(disease_degree_chemical),ny=as.numeric(disease_degree_chemical))
    record$disease_chemical_ny = disease_chemical_ny
    N_chemical = length(unique(c(as.character(chemical_disease_curated$ChemicalID),as.character(chemical_phenotype$chemicalid))))
    record$N_chemical = N_chemical
  })
  
  output$part5_text2 <- renderValueBox({
    valueBox(value = record$N_chemical, subtitle = "Chemical Number", color = "orange", icon = icon("table"))
  })
  
  #### Module 6. Main Analysis Calculation ####
  part5_result = reactive({
    req(chemical_disease_phenotype())
    req(record$m_degree, record$m_degree_chemical)
    req(record$N_gene, record$N_chemical)
    req(record$phenotype_chemical_nx, record$disease_chemical_ny)
    req(record$phenotype_gene_nx, record$disease_gene_ny)
    as_disease_phenotype <- merge(chemical_disease_phenotype(), bp_phenotype_disease, by.x = "phenotypeid", by.y = "GOID" ,all.x = FALSE, all.y = TRUE)
    pheno <- unique(chemical_disease_phenotype()$phenotypeid)
    as_disease_phenotype<-bp_phenotype_disease[bp_phenotype_disease$GOID %in% pheno,]
    as_disease_phenotype[is.na(as_disease_phenotype)] <- 0
    colnames(as_disease_phenotype)
    as_disease_phenotype$chemicalname <- strsplit(as.character(as_disease_phenotype$InferenceChemicalNames),split = "|",fixed = TRUE)
    as_disease_phenotype$genename <- strsplit(as.character(as_disease_phenotype$InferenceGeneSymbols),split = "|",fixed = TRUE)
    N_chemical = record$N_chemical
    phenotype_chemical_nx = record$phenotype_chemical_nx
    disease_chemical_ny = record$disease_chemical_ny
    N_gene = record$N_gene
    phenotype_gene_nx = record$phenotype_gene_nx
    disease_gene_ny = record$disease_gene_ny
    m_degree = record$m_degree
    m_degree_chemical = record$m_degree_chemical
    as_disease_phenotype$ny_chemical <- disease_chemical_ny$ny[match(as_disease_phenotype$DiseaseName, disease_chemical_ny$diseassid)]
    as_disease_phenotype$nx_chemical <- phenotype_chemical_nx$nx[match(as_disease_phenotype$GOID, phenotype_chemical_nx$phenotypeid)]
    as_disease_phenotype$ ny_gene <- disease_gene_ny$ny[match(as_disease_phenotype$DiseaseID, disease_gene_ny$diseassid)]
    as_disease_phenotype$ nx_gene <- phenotype_gene_nx$nx[match(as_disease_phenotype$GOID, phenotype_gene_nx$phenotypeid)]
    as_disease_phenotype$p2_gene <- lapply(as_disease_phenotype$genename, function(x) match(x, m_degree$genesymbol)) %>% 
      lapply( function(x) m_degree$ni[x]*(m_degree$ni[x]-1)/(N_gene * (N_gene-1))) %>% sapply(prod)
    as_disease_phenotype$p2_chemical <- lapply(as_disease_phenotype$chemicalname, function(x) match(as.character(x),m_degree_chemical$chemicalname)) %>% 
      lapply( function(x) m_degree_chemical$ni[x]*(m_degree_chemical$ni[x]-1)/(N_chemical* (N_chemical-1))) %>% sapply(prod)
    result<- as_disease_phenotype[as_disease_phenotype$DiseaseName %in% input$Diseasecategory, ]%>% dplyr::select(c(1,2,3,4,5,7,11,12,13,14,15,16)) 
    result$p1_chemical = numeric(nrow(result))  
    for(i in 1:nrow(result)) {
      result$p1_chemical[i] = calc_p1(N_chemical, result[i, "InferenceChemicalQty"], result[i, "nx_chemical"], result[i, "ny_chemical"])
    }
    result$log10_chemical = log10(result$p1_chemical)
    result$p1_gene = numeric(nrow(result))  
    for(i in 1:nrow(result)) {
      result$p1_gene[i] = calc_p1(N_gene, result[i, "InferenceGeneQty"], result[i, "nx_gene"], result[i, "ny_gene"])
    }
    result$log10_gene = log10(result$p1_gene)
    result$log10_gene[is.na(result$log10_gene)] <- 0
    print("part5!!!!!!!!!!!!!!!!!!!!!!!!!!")
    result
  })
  
  output$part5_table = DT::renderDataTable({
    df = part5_result() %>% dplyr::select(-DiseaseID) %>% dplyr::rename(`p_(gene_comb)` = p2_gene, `p_(chemical_comb)` = p2_chemical, p_chemical = p1_chemical, p_gene = p1_gene)
    DT::datatable( df, class = "display nowrap", selection = "single", options = list(pageLength = 10,scrollX = TRUE))
  })
  
  output$download_part5_table <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(part5_result(), file, row.names = F) }
  )
  
  #### Module 7. Heatmap Visualization ####
  heatfigure = reactive({
    req(part5_result())
    scale <- part5_result()
    w_gene= 1-(exp(1)/(exp(scale$InferenceGeneQty)*2))
    w_chem= 1-(exp(1)/(exp(scale$InferenceChemicalQty)*2))
    scale$Wxya_gene<- -(w_gene*scale$log10_gene + w_gene*log10(scale$p2_gene))
    scale$Wxya_chemical <- -(w_chem*scale$log10_chemical + w_chem*log10(scale$p2_chemical))
    scale$Wxya_gene[scale$Wxya_gene<0|is.na(scale$Wxya_gene)]<-0
    scale$Wxya_chemical[scale$Wxya_chemical<0|is.na(scale$Wxya_chemical)]<-0
    scale$Wxya_chemical[scale$Wxya_chemical<0|is.infinite(scale$Wxya_chemical)]<-0
    scale$Wxya_g_c<-scale$Wxya_gene+scale$Wxya_chemical
    phenotypeid <- unique(scale$GOID)
    phenotypeTerm<-scale$GOName[match(phenotypeid,scale$GOID)]
    heatfigure<-data.frame(phenotypeid=phenotypeid,phenotypeTerm=phenotypeTerm)
    disease_df = scale|> dplyr::select(DiseaseName, DiseaseID) |> distinct() %>% arrange(DiseaseName)
    print(disease_df)
    for(i in 1:nrow(disease_df)) {
      disease_name = disease_df[i, 1]
      disease_id = disease_df[i, 2]
      heatfigure[disease_name] = scale[scale$DiseaseID==disease_id,]$Wxya_g_c[match(heatfigure$phenotypeid,scale[scale$DiseaseID==disease_id,]$GOID)]
    }
    heatfigure[is.na(heatfigure)] <- 0
    heatfigure
  })
  
  heatfigure_max = reactive({
    heatfigure = heatfigure()
    heatfigure$max<- apply(heatfigure[,3:ncol(heatfigure)],1,max,na.rm = TRUE)
    heatfigure
  })
  
  output$part5_table_heatfigure = DT::renderDataTable({
    DT::datatable( heatfigure(), class = "display nowrap", selection = "single", options = list(pageLength = 10,scrollX = TRUE))
  })
  
  output$download_part5_table_heatfigure <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(heatfigure(), file, row.names = F) }
  )
  
  heatfigure_rank_percentile_react = reactive({
    heatfigure = heatfigure()
    heatfigure_rank<-apply(heatfigure[,3:ncol(heatfigure)],2,rank,ties.method="average",na.last=FALSE)
    heatfigure_rank_percentile<-heatfigure_rank/nrow(heatfigure_rank)
    heatfigure_rank_percentile<-as.data.frame(heatfigure_rank_percentile)
    heatfigure_rank_percentile
  })
  
  output$part5_heatfigure_plot = renderPlotly({
    heatfigure = heatfigure()
    bk <- ncol(heatfigure) -2
    cat("bk =", bk, "\n")
    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
    heatfigure_rank<-apply(heatfigure[,3:ncol(heatfigure)],2,rank,ties.method="average",na.last=FALSE)
    cat("nrow(heatfigure_rank) =", nrow(heatfigure_rank), "\n")
    heatfigure_rank_percentile<-heatfigure_rank/nrow(heatfigure_rank)
    heatfigure_rank_percentile<-as.data.frame(heatfigure_rank_percentile)
    heatmaply(t(na.omit(heatfigure_rank_percentile)), colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10), k_col = min(4, nrow(heatfigure_rank_percentile)), k_row = bk, showticklabels = c(F,T))
  })
  
  heatfigure_rank_percentile_index = reactive({
    heatfigure_rank_percentile = heatfigure_rank_percentile_react()
    thresholdvalue = input$thresholdvalue
    index<-rownames(heatfigure_rank_percentile[rowSums(heatfigure_rank_percentile > thresholdvalue) == ncol(heatfigure_rank_percentile), ])
    print(index)
    index
  }) %>%bindEvent(input$refresh_heatmap)
  
  sim_df = reactive({
    heatfigure = heatfigure()
    index = heatfigure_rank_percentile_index()
    heatfigure$phenotypeTerm[as.numeric(index)]
    sim<-mgoSim(GO1 = as.character(heatfigure$phenotypeid[as.numeric(index)]), GO2 = as.character(heatfigure$phenotypeid[as.numeric(index)]), measure = "Wang",semData = hsGO ,combine=NULL)
    sim[is.na(sim)]<-0
    colnames(sim)<-heatfigure$phenotypeTerm[as.numeric(index)]
    rownames(sim)<-heatfigure$phenotypeTerm[as.numeric(index)]
    sim
  }) %>%bindEvent(input$refresh_heatmap)
  
  output$part5_heatmap_plot2 = renderPlotly({
    sim = sim_df()
    heatmaply(sim)
  })
  
  output$part5_hopkins_value <- renderValueBox({
    sim = sim_df()
    dist<-as.dist(1-sim) %>% as.matrix()
    n = ifelse(ncol(dist) >  5, 5, ncol(dist) - 1)
    cluster_test <- get_clust_tendency(dist, n, graph = FALSE)
    cat("Hopkins Statistics: ", cluster_test$hopkins_stat, "\n")
    valueBox(value = round(cluster_test$hopkins_stat, 3), subtitle = "Hopkins Statistics", color = "teal", icon = icon("table"))
  }) 
  
  #### Module 8. AOP Analysis ####
  observe({
    ao<-list()
    input_key_event_name_opt = str_squish(input$key_event_name_opt)
    if(input_key_event_name_opt == "") {
      ao$key_event_name <- aop_ke_mie_ao$`key event name` %>% unique()
    } else {
      match_cond = str_replace_all(input_key_event_name_opt, " ", "|") 
      ao$key_event_name <- aop_ke_mie_ao$`key event name`[grep(match_cond, aop_ke_mie_ao$`key event name`,fixed = FALSE,ignore.case = TRUE)] %>% unique()
    }
    key_event_male_aop = aop_ke_mie_ao %>% dplyr::filter(`key event name` %in% ao$key_event_name)
    aop_ids= key_event_male_aop$`aop id` %>%  unique() %>% sort()
    freezeReactiveValue(input, "aop_id")
    updateSelectInput(inputId = "aop_id", choices = aop_ids, selected = aop_ids)
  })
  
  key_event_male_aop_df = reactive({
    req(input$aop_id)
    cat('aop_id = ',input$aop_id, "\n")
    ao = list()
    ao$aop_id = input$aop_id
    heatfigure = heatfigure_max()
    req(heatfigure)
    key_event_male_aop <-aop_ke_mie_ao[aop_ke_mie_ao$`aop id` %in% ao$aop_id,] 
    key_event_male_aop_ec<-aop_ec[aop_ec$`aop id` %in% ao$aop_id,]    
    ias_disease_phenotype<-chemical_disease_phenotype()
    phenotype_in_aop <- aop_ec[aop_ec$`process ontology id` %in% unique(ias_disease_phenotype$phenotypeid),]
    eventbothdir<-intersect(unique(key_event_male_aop$`key event id`), unique(phenotype_in_aop$`key event id`))
    key_event_goid<-unique(aop_ec$`process ontology id`[aop_ec$`key event id` %in% key_event_male_aop $`key event id`])
    unique(key_event_goid[grepl("GO",key_event_goid)])
    intersect(phenotype_in_aop$`process ontology id`,unique(key_event_goid[grepl("GO",key_event_goid)]))
    key_event_male_aop_ec<-aop_ec[aop_ec$`aop id` %in% ao$aop_id,]
    key_event_male_aop_GO<-key_event_male_aop_ec[grepl("GO",key_event_male_aop_ec$`process ontology id`),]
    unique(key_event_male_aop_GO$`key event id`)
    unique(key_event_male_aop_ec$`key event id`)
    unique(key_event_male_aop_ec$`process ontology id`)
    getOffsprings()
    children <- GOSimEnv$children
    mapchildren<- function(eventlist){
      mapsuc<-vector()
      for (i in 1:length(eventlist)){
        eventid<-eventlist[i]
        goterm = key_event_male_aop_GO$`process ontology id`[match(eventid,key_event_male_aop_GO$`key event id`)]
        judge = ifelse(length(intersect(as.vector(ias_disease_phenotype$phenotypeid), as.vector(children[[as.character(goterm)]])) ) == 0,0, 1)
        mapsuc<-c(mapsuc,judge)
      }
      return(mapsuc)
    }
    mapchildren2<- function(eventlist){
      mapsuc<-list()
      for (i in 1:length(eventlist)){
        eventid<-eventlist[i]
        goterm = key_event_male_aop_GO$`process ontology id`[match(eventid,key_event_male_aop_GO$`key event id`)]
        judge = intersect(as.vector(ias_disease_phenotype$phenotypeid), as.vector(children[[as.character(goterm)]]))
        mapsuc[[i]]<- judge
      }
      return(mapsuc)
    }
    key_event_male_aop_GO$mapchildren<-mapchildren(key_event_male_aop_GO$`key event id`)
    eventboth<-unique(c(eventbothdir,as.character(key_event_male_aop_GO$`key event id`[key_event_male_aop_GO$mapchildren ==1])))
    aop_ke_mie_ao$`key event name`[match(as.character(key_event_male_aop_GO$`key event id`[key_event_male_aop_GO$mapchildren ==1]), aop_ke_mie_ao$`key event id`)]
    key_event_male_aop_GO$eventscore<-lapply(mapchildren2(key_event_male_aop_GO$`key event id`), match,heatfigure$phenotypeid) %>% 
      lapply(function(x) heatfigure$max[x]) %>%  
      lapply(function(x) ifelse(length(x)==0,0,max(x))) %>% 
      sapply(unique) %>% unlist
    key_event_male_aop<-aop_ke_mie_ao[aop_ke_mie_ao$`aop id` %in% ao$aop_id,]
    key_event_male_aop$event_is_in_ctd<-ifelse(key_event_male_aop$`key event id` %in% eventboth,1,0)
    index<-key_event_male_aop_ec$`process ontology id`[key_event_male_aop_ec$`key event id` %in% eventboth] %>% match(heatfigure$phenotypeid) 
    left<- distinct(heatfigure[na.omit(index),],phenotypeid,.keep_all = TRUE) 
    dim(phenotype_in_aop)
    unique(phenotype_in_aop$`process term`)
    colnames(phenotype_in_aop)
    phenotype_in_aop$eventscore<- heatfigure$max[match(phenotype_in_aop$`process ontology id`,heatfigure$phenotypeid)]
    key_event_male_aop$eventscore<-phenotype_in_aop$eventscore[match(key_event_male_aop$`key event id`,phenotype_in_aop$`key event id`)]
    if(is.null(key_event_male_aop$eventscore)) { key_event_male_aop$eventscore = 0 }
    colnames(key_event_male_aop)
    key_event_male_aop$event_children<-key_event_male_aop_GO$eventscore[match(key_event_male_aop$`key event id`,key_event_male_aop_GO$`key event id`)]
    if(is.null(key_event_male_aop$event_children)) { key_event_male_aop$event_children = 0 }
    key_event_male_aop$eventscore[is.na(key_event_male_aop$eventscore)]<-0
    key_event_male_aop$event_children[is.na(key_event_male_aop$event_children)]<-0
    key_event_male_aop$max<-apply(key_event_male_aop[,c(6,7)],1,max)
    key_event_male_aop$eventmaxscore<-key_event_male_aop$max[match(key_event_male_aop$`key event id`,key_event_male_aop$`key event id`)]
    key_event_male_aop %>% dplyr::rename(`event is in ctd` = event_is_in_ctd, `event score` = eventscore, `event children` = event_children, `event max score` = eventmaxscore)
  })
  
  output$aop_key_event_male_aop = DT::renderDataTable({
    DT::datatable( key_event_male_aop_df(), class = "display nowrap", selection = "single", options = list(pageLength = 10,scrollX = TRUE))
  })
  
  output$download_aop_event_male <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(key_event_male_aop_df(), file, row.names = F) }
  )
  
  #### Module 8. AOP Analysis (Continued) ####
  happy_df = reactive({
    key_event_male_aop = key_event_male_aop_df()
    req(key_event_male_aop)
    happy <- key_event_male_aop %>% 
      group_by(`aop id`) %>% 
      summarize(`total count` = n(), 
                proportion = sum(`event is in ctd` == 1) / `total count`, 
                `total score` = sum(case_when(`event is in ctd` == 1 ~ max, TRUE ~ 0)))
    happy
  })
  
  output$aop_happy = DT::renderDataTable({
    happy = happy_df()
    DT::datatable( happy, class = "display nowrap", selection = "single", options = list(pageLength = 10,scrollX = TRUE))
  }) 
  
  output$download_aop_happy <- downloadHandler(
    filename = function() { paste0("table.csv") },
    content = function(file) { write.csv(happy_df(), file, row.names = F) }
  )
  
}

#### Part 7. App Initialization ####
shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port=3675))
