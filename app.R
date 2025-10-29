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



calc_p1 = function(N, m, nx, ny) {
    
    cat("N = ", N, "\tm = ",m,"\tnx = ", nx, "\tny = ", ny, "\n", sep = "")
    
    if(is.na(m) || is.na(nx) || is.na(ny)) {
        return(0)
    }
    
    result = chooseZ(N, m)   * chooseZ(N - m, nx - m) * chooseZ(N - nx, ny - m) / ( chooseZ(N, nx) * chooseZ(N, ny))
    
    return(as.numeric(result))
}


#2.1 phenotype basic data
bp_phenotype_disease = readRDS("./data/bp_phenotype_disease.rds") # 4s

#2.2 BDE47_disease
disease = readRDS("./data/disease.rds")  # 18s
ChemicalName_df = disease %>% dplyr::select(1,2) %>% dplyr::distinct() %>% arrange(ChemicalName)
ChemicalName_display = paste(ChemicalName_df[[2]], ChemicalName_df[[1]], sep = " ")
ChemicalName_value = ChemicalName_df[[1]]

ChemicalName_value_pick = ChemicalName_value
names(ChemicalName_value_pick) = ChemicalName_display

BDE47_selected <- c("2,2',4,4'-tetrabromodiphenyl ether",  
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

#2.4 CTD_disease category
CTD_disease = readRDS("./data/CTD_disease.rds") # 0.1s

#2.5 as_phenotype phenotype
bde_phenotype = readRDS("./data/bde_phenotype.rds") # 2s


#5.1 relevant data in p-d inference work
#5.1.1 inference network by GENES

#gene degree  
gene_phenotype = readRDS("./data/gene_phenotype.rds")
array<-table(gene_phenotype$Gene.Symbol)
gene_degree1_phenotype<-data.frame(genesymbol = as.character(names(array)),degree=as.numeric(array))

# read subset  gene_degree2_disease
gene_degree2_disease = readRDS("./data/gene_degree2_disease.rds") # 0.08
#gene degree  
gene_phenotype = readRDS("./data/gene_phenotype.rds")
array<-table(gene_phenotype$Gene.Symbol)
gene_degree1_phenotype<-data.frame(genesymbol = as.character(names(array)),degree=as.numeric(array))

# read subset  gene_degree2_disease
gene_degree2_disease = readRDS("./data/gene_degree2_disease.rds") # 0.08


#5.1.2 inference network by CHEMICALS #####
#chmical degree 
chemical_disease_curated = readRDS("./data/chemical_disease_curated.rds")        # 0.3
# read chemical_phenotype
chemical_phenotype = readRDS("./data/chemical_phenotype.rds")        # 1.6
colnames(chemical_phenotype)

array1<-table(unique(chemical_phenotype$chemicalname))
chemical_degree1_phenotype<-data.frame(chemicalname = as.character(names(array1)),degree=as.numeric(array1))
array2<-table(chemical_disease_curated$ChemicalName)
chemical_degree2_disease<-data.frame(chemicalname = as.character(names(array2)),degree=as.numeric(array2))

#8.1 clust
hsGO = readRDS(file = "./data/hsGO.rds")



# *********PART6 INTEGRATED ANALYSIS WITH AOP Wiki

#data reading
aop_ke_mie_ao<-read.csv(file="./aop_ke_mie_ao.tsv",sep = "\t",header = FALSE)
colnames(aop_ke_mie_ao)<-c("aop id", "key event id", 
                           "key event type", "key event name")
# head(aop_ke_mie_ao)

aop_ker<-read.csv(file="./aop_ke_ker.tsv",sep = "\t",header = FALSE)
colnames(aop_ker)<-c("aop id", 
                     "upstream event id", "downstream event id", 
                     "relationship id", "direct or indirect relationship", 
                     "evidence for relationship", "quantitative understanding of relationship")
# head(aop_ker)

aop_ec<-read.csv(file="./aop_ke_ec.tsv",sep = "\t",header = FALSE)

colnames(aop_ec)<-c("aop id", "key event id", "action", "object source", "object ontology id", "object term", 
                    "process source", "process ontology id", "process term")
# head(aop_ec)

# key_event_name_opts = c("heart", "vascular", "hypertension", "Coronary")



ui <- fluidPage( 
    theme = shinytheme("flatly"), #theme = shinythemes::shinytheme("cosmo"),
    
    titlePanel('AOP-ExpoVis'), 
    useShinydashboard(),
    use_tippy(),
    sidebarLayout(
        # 侧边栏  
        sidebarPanel(width = 4,
                     
                     # part 1
                     multiInput(
                         inputId = "BDE47",
                         label = "Chemical Names:", 
                         # choices = ChemicalName,
                         selected = BDE47_selected,
                         choiceNames = ChemicalName_display, 
                         choiceValues = ChemicalName_value,
                         options = list(
                             enable_search = TRUE,
                             non_selected_header = "Choose between:",
                             selected_header = "You have selected:"
                         ),
                         width = "100%"
                     ),
                     
                     # pickerInput(
                     #     inputId = "BDE471999",
                     #     label = "Chemical Names (soluciton 2):", 
                     #     choices = ChemicalName_value_pick,
                     #     selected = BDE47_selected,
                     #     multiple = TRUE,
                     # 
                     #     options =  pickerOptions (
                     #         actionsBox = TRUE,
                     #         container = "body",
                     #     ),
                     #         
                     # 
                     #     width = "100%"
                     # ),
                     
                     textOutput("part1_message"),
                     br(),
                     
                     
                     
                     # part 3
                     selectInput(inputId = "diseaseclass",
                                 label = "Disease Class:",multiple = FALSE,
                                 choices = NULL),
                     
                     textOutput("part3_text1"),
                     br(),
                     
                     
                     
                     # part 4
                     
                     selectInput(inputId = "Urogenitalcategory",
                                 label = "Urogenitalcategory",multiple = T,
                                 choices = NULL),
                     ## 添加helptext
                     
                     textOutput("part4_text1"),
                     br(),
                     
                     actionButton("refresh", "Refresh All...", icon = icon("rotate")),
                     
                     p("Welcome to AOP-ExpoVis, an interactive web application designed to explore relationships between chemicals, diseases, and phenotypes. If you use AOP-ExpoVis in your research, please cite our work as follows: Author(s), Title, Journal, Year, DOI."),
                     
                     
                     
        ),
        
        # mainPanl
        mainPanel(width = 8,
                  
                  tabsetPanel(
                      
                      # part 1
                      tabPanel("Chemical", 
                               DT::dataTableOutput("part1_table") %>% shinycssloaders::withSpinner(type = 5),
                               
                               downloadButton("download_part1_table", "Download..."),
                               div(style = "width: 100%; margin-bottom: 100px;"),
                               
                      ),
                      
                      # part3
                      
                      tabPanel("Disease", 
                               
                               DT::dataTableOutput("part3_table") %>% shinycssloaders::withSpinner(type = 5),
                               
                               downloadButton("download_part3_table", "Download..."),
                               div(style = "width: 100%; margin-bottom: 100px;"),
                               
                      ),  
                      
                      # part 4
                      tabPanel("Phenotype", 
                               
                               # h2("xxxx"),
                               
                               # h3("bde47_phenotype_normal"),
                               uiOutput("title_part4_table1"),
                               DT::dataTableOutput("part4_table1") %>% shinycssloaders::withSpinner(type = 5),
                               downloadButton("download_part4_table1", "Download..."),
                               hr(),
                               
                               # h3("bp_phenotype_nervous"),
                               uiOutput("title_part4_table2"),
                               DT::dataTableOutput("part4_table2") %>% shinycssloaders::withSpinner(type = 5),
                               
                               hr(),
                               
                               # h3("bde_urogenital_phenotype"),
                               uiOutput("title_part4_table3"),
                               DT::dataTableOutput("part4_table3") %>% shinycssloaders::withSpinner(type = 5),
                               
                               div(style = "width: 100%; margin-bottom: 100px;"),
                               
                               
                      ),  
                      
                      
                      # part 5
                      tabPanel("Heatmap", 
                               
                               # part 5
                               fluidRow(
                                   # actionButton("refresh_heatmap", "Show/Refresh heatmap...", icon = icon("rotate")),
                                   # column(width = 3, textOutput("part5_text1"),),
                                   # column(width = 3, textOutput("part5_text2"),),
                                   
                                   fluidRow(
                                       column(width = 4,        valueBoxOutput("part5_text1", width = NULL)),   
                                       column(width = 4,        valueBoxOutput("part5_text2", width = NULL)),   
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
                                                       tippy_this("tips_thresholdvalue",
                                                                  "<p style='font-size: 1.5em;'>The proportion chosen with the slider determines the percentage of data points to be excluded from the subsequent analysis, ranging from 0.01 to 0.99. For example, setting the exclusion ratio to 0.9 means that 90% of the data points will be excluded, and the analysis will include the remaining 10%.</p>",
                                                                  allowHTML = TRUE,
                                                                  placement = 'top'
                                                                  ),
                                                ),
                                                column(width = 4,        
                                                       actionButton("refresh_heatmap", "Show/Refresh heatmap...", icon = icon("rotate")),
                                                ),
                                                column(width = 4, valueBoxOutput("part5_hopkins_value", width = NULL), # %>% shinycssloaders::withSpinner(type = 5),
                                                )
                                                
                                            ),
                                            
                                            div(style = "width: 100%; overflow-x: auto; overflow-y: auto; margin-right: 5px; margin-bottom: 5px;",

                                              plotlyOutput("part5_heatmap_plot2", height = "700px") %>% jqui_resizable()  %>% shinycssloaders::withSpinner(type = 5),
                                            
                                            ),


                                   ),

                                   tabPanel("AOP", 
                                            # selectInput(inputId = "key_event_name_opt_1",
                                            #             label = "Key Event Names:",
                                            #             multiple = T,
                                            #             choices = key_event_name_opts,
                                            #             selected = key_event_name_opts,
                                            #             width = "100%"),
                                            
                                            searchInput(
                                                inputId = "key_event_name_opt",
                                                label = "Please input key words of event name to search:", 
                                                placeholder = "keywords of event name, separated by whitespace",
                                                btnSearch = icon("magnifying-glass"), 
                                                btnReset = icon("xmark"),
                                                width = "100%"
                                            ),
                                            
                                            # selectInput
                                            multiInput(inputId = "aop_id",
                                                        label = "Please select AOP IDs:",
                                                        # multiple = TRUE,
                                                        choices = c("NA"),
                                                        width = "100%"),
                                            
                                            hr(),
                                            
                                            # DT::dataTableOutput("aop_ao") %>% shinycssloaders::withSpinner(type = 5),
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



server <- function(input, output, session) {
    
    # ###---------------------------------------------------------------------------
    # ###part1
    # ###---------------------------------------------------------------------------
    bde_disease <- reactive({
        # bde_disease <- disease[disease$ChemicalName %in% input$BDE47,]
        
        print(input$BDE47)
        
        if(is.null(input$BDE47)) {
            disease
        } else {
            disease %>%
                dplyr::filter(ChemicalName %in% input$BDE47)
        }
        
    })
    #
    bde47_phenotype <- reactive({
        # bde47_phenotype<-bde_phenotype[bde_phenotype$chemicalname
        #                                %in% input$BDE47,]
        
        if(is.null(input$BDE47)) {
            bde_phenotype
        } else {
            
            bde_phenotype %>%
                dplyr::filter(chemicalname %in% input$BDE47)
        }
        
    })
    
    output$part1_message <- renderText({
        
        phen_num <- length(unique(bde47_phenotype()$phenotypeid)) #95 根据选择的化学物初步匹配出**种表型，涉及**项研究。
        
        
        Ref_num <- length(unique(bde47_phenotype()$pubmedids)) #44 这里把bde47_phenotype数据库展示出来供用户下载
        
        paste0("Based on the selected chemical compounds, an initial match identified ", phen_num, " phenotypes involving ", Ref_num, " studies.")
    })
    # 
    output$part1_table <- DT::renderDataTable({
        DT::datatable(bde47_phenotype(), 
                      selection = "single",
                      class = "display nowrap",
                      options = list(pageLength = 20,
                                     scrollX = TRUE)
        )
    }) %>% bindEvent(input$refresh)
    
    output$download_part1_table <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(bde47_phenotype(), file, row.names = F)
        }
    )
    # 
    # 
    # ###---------------------------------------------------------------------------
    # ###part2
    # ###---------------------------------------------------------------------------
    # 
    # 
    slimmesh_info <- reactive({
        
        # bde_disease()$DiseaseID %in% CTD_disease$DiseaseID
        slimmesh_info <- CTD_disease[match(bde_disease()$DiseaseID,table = CTD_disease$DiseaseID),]
        slimmesh_info
    })
    # 
    # 
    disease_slimmesh <- reactive({
        disease_slimmesh <- data.frame(Disease.Name = bde_disease()$diseasename,SlimMappings=slimmesh_info()$SlimMappings)
    })
    # 
    # 
    
    
    observe({
        # disease_slimmesh <- data.frame(Disease.Name = bde_disease()$DiseaseName,SlimMappings=slimmesh_info()$SlimMappings)
        # strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE)
        
        slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE ))) %>% sort() #29
        
        freezeReactiveValue(input, "SlimMapping")
        updateSelectInput(
            inputId = "SlimMapping",
            label = "SlimMapping",
            choices = slimmeshcategory)
    })
    # 
    # 
    nervous_diease <- reactive({
        req(input$diseaseclass)
        
        # nervous_diease <- CTD_disease[grepl(input$SlimMapping, CTD_disease$SlimMappings), ]
        nervous_diease <- CTD_disease[grepl(input$diseaseclass, CTD_disease$SlimMappings), ]
        
        nervous_diease
    })
    # 
    # #
    # #
    output$part2_table_slimmesh_info <- DT::renderDataTable({
        DT::datatable(slimmesh_info(), 
                      selection = "single",
                      class = "display nowrap",
                      options = list(pageLength = 10,scrollX = TRUE)
        ) %>%
            DT::formatStyle(columns = c("Definition"), textOverflow = 'ellipsis', overflow = 'hidden', whiteSpace = 'nowrap')
    })
    # 
    # 
    output$part2_table_nervous_diease <- DT::renderDataTable({
        
        DT::datatable(nervous_diease(), 
                      selection = "single",
                      class = "display nowrap",
                      options = list(pageLength = 10,scrollX = TRUE)
        )
    })
    
    # 
    output$download_part2_table1 <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(slimmesh_info(), file, row.names = F)
        }
    )
    # 
    # 
    output$download_part2_table2 <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(nervous_diease(), file, row.names = F)
        }
    )
    # 
    # 
    # 
    # 
    # 
    # ###---------------------------------------------------------------------------
    # ### part3
    # ###---------------------------------------------------------------------------
    # 
    # ### 更新diseaseclass控件
    observe({
        req(slimmesh_info())
        req(disease_slimmesh())
        
        
        if (nrow(disease_slimmesh()) > 1) {
            # disease_slimmesh <- data.frame(Disease.Name = bde_disease()$DiseaseName,SlimMappings=slimmesh_info()$SlimMappings)
            # strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE)
            
            slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE ))) %>% sort() #29
            
            # print(disease_slimmesh())
            
            #3.3 grouping ias_diseases according slimmesh category
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
            updateSelectInput(
                inputId = "diseaseclass",
                choices = names(diseaseclass) )
            
        }
        
    })
    # 
    # 
    # 
    # 
    diseaseclass <- reactive({
        
        
        req(slimmesh_info())
        req(disease_slimmesh())
        # disease_slimmesh <- data.frame(Disease.Name = bde_disease()$DiseaseName,SlimMappings=slimmesh_info()$SlimMappings)
        # strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE)
        
        slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE )))  #29
        #3.3 grouping ias_diseases according slimmesh category
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
    # 
    
    output$part3_text1 <- renderText({
        
        
        slimmeshcategory<-unique(unlist(strsplit(as.character(slimmesh_info()$SlimMappings),split = "|",fixed = TRUE )))  #29
        
        # print(slimmeshcategory)
        # print(length(slimmeshcategory))
        
        diseaseclass_select_num <- map(diseaseclass()[input$diseaseclass],length) %>% unlist() %>% sum()
        print(diseaseclass_select_num)
        
        paste0( "The selected chemical compounds are associated with ",length(slimmeshcategory)," categories of diseases, including ",length(diseaseclass()[[input$diseaseclass]])," instances of ",input$diseaseclass, ".")
        
    })
    
    
    output$part3_table = DT::renderDataTable({
        
        df = data.frame(diseaseclass()[input$diseaseclass])
        
        DT::datatable(df, 
                      selection = "single",
                      class = "display nowrap",
                      options = list(pageLength = 10,scrollX = TRUE)
        )
        
    }) %>% bindEvent(input$refresh)
    
    
    # output$part3_text2 = renderText({
    #
    #   diseaseclass()[input$diseaseclass]
    #   # DT::datatable(diseaseclass()[input$diseaseclass], options = list(pageLength = 5,scrollX = TRUE))
    #
    # })
    
    
    # output$download_part3_table <- downloadHandler(
    #   # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
    #   filename = function() {
    #     paste0("part3_table.csv")
    #   },
    #   # 参数content 固定的用法 file接收filename的字符
    #   content = function(file) {
    #     write.csv(diseaseclass()[input$diseaseclass], file, row.names = F)
    #   }
    # )
    
    
    
    # ###---------------------------------------------------------------------------
    # ### part4  Phenotype Analysis
    # ###---------------------------------------------------------------------------
    # 
    # ### 输入！！
    observe({
        req(input$diseaseclass)
        
        freezeReactiveValue(input, "Urogenitalcategory")
        
        choices = diseaseclass()[[input$diseaseclass]] 
        label = paste0(input$diseaseclass, " ", "Category:")
        
        updateSelectInput(
            inputId = "Urogenitalcategory",
            label = label,
            choices = choices,
            selected = head(choices, 3))
    })
    # # Urogenitalcategory<-c("Alzheimer Disease","Parkinson Disease","Cognitive Dysfunction",'Cognition Disorders',"Memory Disorders")
    # 
    #
    abnormal_anatomy<-c("HL-60 Cells","Cell Line,Transformed",
                        "Hep-G2","Hela Cells",
                        "K562","MCF-T Cells","HCT 116","A549",
                        "Jurkat Cells","Tumor",
                        "HT29","Caco-2 Cells",
                        "CHO","PC12","Hep G2 Cells","MCF-7 Cells")              # 默认
    # 
    # 
    bde47_phenotype_normal = reactive({  
        abanatomyserch<-apply(as.data.frame(abnormal_anatomy),1,grepl,bde47_phenotype()$anatomyterms,ignore.case = TRUE)
        abanatomyindex<-apply(abanatomyserch,1,sum)
        sum(abanatomyindex)
        bde47_phenotype_normal<-bde47_phenotype()[!as.logical(abanatomyindex),]
        
        print("bde47_phenotype_normal")
        bde47_phenotype_normal
        
    })
    # 
    # 
    output$title_part4_table1 = renderUI({
        tagList(
            h3("Chemical Phenotype")
        )
        
    })
    output$part4_table1 = DT::renderDataTable({
        DT::datatable( 
            bde47_phenotype_normal(), 
            selection = "single",
            class = "display nowrap",
            options = list(pageLength = 10,scrollX = TRUE)
            
        )
        
    }) %>% bindEvent(input$refresh)
    # 
    # 
    bp_phenotype_nervous = reactive({
        print("bp_phenotype_nervous")
        
        req(input$Urogenitalcategory)
        
        bp_phenotype_nervous <-bp_phenotype_disease[
            bp_phenotype_disease$DiseaseName %in% input$Urogenitalcategory,]
        
        bp_phenotype_nervous
        
    })
    # 
    # 
    # 
    output$title_part4_table2 = renderUI({
        tagList(
            
            # h3("bp_phenotype_nervous"),
            h3(input$diseaseclass),

        )
        
    })
    output$part4_table2 = DT::renderDataTable({
        DT::datatable(
            bp_phenotype_nervous(),
            selection = "single",
            class = "display nowrap",
            options = list(pageLength = 10,scrollX = TRUE)
            
        )
        
    }) %>% bindEvent(input$refresh)
    # 
    # 
    bde_urogenital_phenotype = reactive({
        
        print("bde_urogenital_phenotype")
        
        bde_urogenital_phenotype <-bde47_phenotype_normal()[bde47_phenotype_normal()$phenotypename %in% bp_phenotype_nervous()$GOName,]
        bde_urogenital_phenotype
    })
    # 
    output$title_part4_table3 = renderUI({
        tagList(
            # h3("bde_urogenital_phenotype"),
            h3(paste0(input$diseaseclass, " ", "Phenotype")),
        )
        
    })
    output$part4_table3 = DT::renderDataTable({
        DT::datatable( 
            bde_urogenital_phenotype(), 
            selection = "single",
            class = "display nowrap",
            options = list(pageLength = 10,scrollX = TRUE)
            
        )
        
    }) %>% bindEvent(input$refresh)
    # 
    # 
    # 
    # 
    output$part4_text1 <- renderText({
        
        # print(bde47_phenotype_normal())
        print(length(unique(bde47_phenotype_normal()$phenotypeid))  )                  #92 在正常组织或细胞（非肿瘤组织或细胞）中的表型数量
        length(unique(bde47_phenotype_normal()$pubmedids))                      #38 及研究数量 输出bde47_phenotype_normal数据框
        
        phenotype_normal = bde47_phenotype_normal()
        phenotype_selected = bp_phenotype_nervous()
        bde_phenotype = bde_urogenital_phenotype()
        
        paste0("In normal (non-tumor) tissues or cells, chemicals are linked to ",length(unique(phenotype_normal$phenotypeid)),
               " phenotypes (",length(unique(phenotype_normal$pubmedids)), " studies). Selected diseases relate to ",
               length(unique(phenotype_selected$GOID)), " phenotypes, intersecting with them to yield ",
               length(unique(bde_phenotype$phenotypename)), " overlapping phenotypes.")
    })
    # 
    # 
    # ###---------------------------------------------------------------------------
    # ### part5 
    # ###---------------------------------------------------------------------------
    # 
    
    record = reactiveValues(
        N_chemical = 0,
        phenotype_chemical_nx = NULL,
        disease_chemical_ny = NULL,
        
        
        N_gene = 0,
        phenotype_gene_nx = NULL,
        disease_gene_ny = NULL,
        
        
        m_degree = NULL,
        m_degree_chemical = NULL,

    )
    
    
    
    # 
    # 
    observe({
        req(input$Urogenitalcategory)
        req(nervous_diease())
        
        # browser()
        
        CTD_disease_genes_curated <- gene_degree2_disease[gene_degree2_disease$DiseaseName %in% input$Urogenitalcategory,]
        
        # nervousid<- nervous_diease()$DiseaseID
        
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
        
        
        #nx phenotype degree inferred by genes
        phenotype_degree<-table(gene_phenotype$GO.Term.ID)
        phenotype_gene_nx<-data.frame(phenotypeid=names(phenotype_degree),nx=as.numeric(phenotype_degree))
        
        record$phenotype_gene_nx = phenotype_gene_nx
        
        #ny disease degree inferred by genes
        disease_degree<-table(CTD_disease_genes_curated$DiseaseID)
        disease_gene_ny<-data.frame(diseassid=names(disease_degree),ny=as.numeric(disease_degree))
        
        record$disease_gene_ny = disease_gene_ny
        
        
        
        #N=75247  这里得出的75247数字记录为N_gene，在后面步骤展示
        length(unique(gene_degree2_disease$genesymbol))              #1673
        length(unique(gene_degree1_phenotype$genesymbol))            #75198
        length(unique(c(as.character(gene_degree2_disease$genesymbo),as.character(gene_degree1_phenotype$genesymbol))))       # 保存这个信息  75247
        
        
        
        N_gene = length(unique(c(as.character(gene_degree2_disease$genesymbo),as.character(gene_degree1_phenotype$genesymbol))))
        print(N_gene)
        
        record$N_gene = N_gene
        
        
    })
    # 
    # 
    output$part5_text1 <- renderValueBox({
        # req(record$N_gene)
        
        # paste0( "N_gene:  ", record$N_gene)
        
        valueBox(
            value = record$N_gene,
            subtitle = "Gene Number",
            icon = icon("list-alt")
        )
    }) # %>% bindEvent(input$refresh_heatmap)
    
    
    # 
    # 
    # 
    # 
    # 
    # 
    # ### 5.1.2
    
    CTD_disease_chems_curated = reactive({
        req(input$Urogenitalcategory)
        
        CTD_disease_chems_curated <- chemical_disease_curated[chemical_disease_curated$diseasename %in% input$Urogenitalcategory,]
        CTD_disease_chems_curated
    }) 
    
    observe({
        # byx<-chemical_degree1_phenotype[chemical_degree1_phenotype$chemicalname %in%  CTD_disease_chems_curated()$ChemicalName,]
        req(input$Urogenitalcategory)
        req(CTD_disease_chems_curated())
        
        m_degree_chemical<-merge(chemical_degree1_phenotype[chemical_degree1_phenotype$chemicalname %in%  CTD_disease_chems_curated()$ChemicalName,],
                                 chemical_degree2_disease[chemical_degree2_disease$chemicalname %in% CTD_disease_chems_curated()$ChemicalName,],
                                 by.x = "chemicalname",by.y = "chemicalname",all = TRUE)
        m_degree_chemical[is.na(m_degree_chemical)]<-0
        
        m_degree_chemical$ni <- m_degree_chemical$degree.x + m_degree_chemical$degree.y
        
        record$m_degree_chemical = m_degree_chemical
        
        #nx phenotype degree inferred by chemicals
        chemical_phenotype_unique<-chemical_phenotype[!duplicated(chemical_phenotype[,c("chemicalid","phenotypeid")]),]
        phenotype_degree_chemical<-table(chemical_phenotype_unique$phenotypeid)
        phenotype_chemical_nx<-data.frame(phenotypeid=names(phenotype_degree_chemical),nx=as.numeric(phenotype_degree_chemical))
        
        record$phenotype_chemical_nx = phenotype_chemical_nx
        
        #ny disease degree inferred by chemicals
        disease_degree_chemical<-chemical_disease_curated[chemical_disease_curated$diseasename %in% input$Urogenitalcategory,]     # 这里使用了输入
        disease_degree_chemical<-table(factor(disease_degree_chemical$diseasename))
        disease_chemical_ny<-data.frame(diseassid=names(disease_degree_chemical),ny=as.numeric(disease_degree_chemical))
        
        record$disease_chemical_ny = disease_chemical_ny
        
        
        #N-chemical=14478 这里得出的14478记录为N_chemical，在后面步骤展示
        length(unique(factor(chemical_disease_curated$ChemicalID)))#10389
        length(unique(chemical_phenotype$chemicalid))#9883
        length(unique(c(as.character(chemical_disease_curated$ChemicalID),as.character(chemical_phenotype$chemicalid))))    # 保存这个信息          #14478
        
        
        N_chemical = length(unique(c(as.character(chemical_disease_curated$ChemicalID),as.character(chemical_phenotype$chemicalid))))
        
        record$N_chemical = N_chemical
        
        
    })
    # 
    output$part5_text2 <- renderValueBox({
        
        # req(record$N_chemical)
        # paste0( "N_chemical:  ",record$N_chemical)
        valueBox(
            value = record$N_chemical,
            subtitle = "Chemical Number",
            color = "orange",
            icon = icon("table")
        )
    }) # %>% bindEvent(input$refresh_heatmap)
    # 
    
    
    # 
    # ### 5.3
    # 
    part5_result = reactive({
        
        req(bde_urogenital_phenotype())
        req(record$m_degree, record$m_degree_chemical)
        req(record$N_gene, record$N_chemical)
        req(record$phenotype_chemical_nx, record$disease_chemical_ny)
        req(record$phenotype_gene_nx, record$disease_gene_ny)
        
        
        as_Urogenital_phenotype <- merge(bde_urogenital_phenotype(), bp_phenotype_disease, by.x = "phenotypeid", by.y = "GOID" ,all.x = FALSE, all.y = TRUE)
        
        pheno <- unique(bde_urogenital_phenotype()$phenotypeid)
        
        as_Urogenital_phenotype<-bp_phenotype_disease[
            bp_phenotype_disease$GOID %in% pheno,]
        
        # as_Urogenital_phenotype = as.data.frame(as_Urogenital_phenotype)
        
        as_Urogenital_phenotype[is.na(as_Urogenital_phenotype)] <- 0
        
        
        colnames(as_Urogenital_phenotype)
        
        as_Urogenital_phenotype$chemicalname <- strsplit(as.character(as_Urogenital_phenotype$InferenceChemicalNames),split = "|",fixed = TRUE)
        as_Urogenital_phenotype$genename <- strsplit(as.character(as_Urogenital_phenotype$InferenceGeneSymbols),split = "|",fixed = TRUE)
        
        #####################################################
        # restore record
        N_chemical = record$N_chemical
        phenotype_chemical_nx = record$phenotype_chemical_nx
        disease_chemical_ny = record$disease_chemical_ny
        
        
        N_gene = record$N_gene
        phenotype_gene_nx = record$phenotype_gene_nx
        disease_gene_ny = record$disease_gene_ny
        
        
        m_degree = record$m_degree
        m_degree_chemical = record$m_degree_chemical
        ####################################################
        
        
        as_Urogenital_phenotype$ny_chemical <- disease_chemical_ny$ny[
            match(as_Urogenital_phenotype$DiseaseName, disease_chemical_ny$diseassid)]
        
        as_Urogenital_phenotype$nx_chemical <- phenotype_chemical_nx$nx[
            match(as_Urogenital_phenotype$GOID, phenotype_chemical_nx$phenotypeid)]
        
        as_Urogenital_phenotype$ ny_gene <- disease_gene_ny$ny[
            match(as_Urogenital_phenotype$DiseaseID, disease_gene_ny$diseassid)]
        
        as_Urogenital_phenotype$ nx_gene <- phenotype_gene_nx$nx[
            match(as_Urogenital_phenotype$GOID, phenotype_gene_nx$phenotypeid)]
        
        
        
        as_Urogenital_phenotype$p2_gene <- lapply(as_Urogenital_phenotype$genename, function(x) match(x, m_degree$genesymbol)) %>% 
            # lapply( function(x) m_degree$ni[x]*(m_degree$ni[x]-1)/(75247*75246)) %>%
            lapply( function(x) m_degree$ni[x]*(m_degree$ni[x]-1)/(N_gene * (N_gene-1))) %>%
            sapply(prod)
        #这里的数字是N_gene及N_gene-1，不用展示
        
        as_Urogenital_phenotype$p2_chemical <- lapply(as_Urogenital_phenotype$chemicalname, function(x) match(as.character(x),m_degree_chemical$chemicalname)) %>% 
            # lapply( function(x) m_degree_chemical$ni[x]*(m_degree_chemical$ni[x]-1)/(14478*14477)) %>%
            lapply( function(x) m_degree_chemical$ni[x]*(m_degree_chemical$ni[x]-1)/(N_chemical* (N_chemical-1))) %>%
            sapply(prod)
        #这里的数字是N_chemical及N_chemical-1，不用展示
        
        
        result<- as_Urogenital_phenotype[as_Urogenital_phenotype$DiseaseName %in% input$Urogenitalcategory, ]%>%
            dplyr::select(c(1,2,3,4,5,7,11,12,13,14,15,16)) 
        
        #展示result数据框。
        
        
        result$p1_chemical = numeric(nrow(result))  
        for(i in 1:nrow(result)) {
            result$p1_chemical[i] = calc_p1(N_chemical, 
                                            result[i, "InferenceChemicalQty"],
                                            result[i, "nx_chemical"],
                                            result[i, "ny_chemical"]
            )
        }
        
        result$log10_chemical = log10(result$p1_chemical)
        
        
        result$p1_gene = numeric(nrow(result))  
        for(i in 1:nrow(result)) {
            result$p1_gene[i] = calc_p1(N_gene, 
                                        result[i, "InferenceGeneQty"],
                                        result[i, "nx_gene"],
                                        result[i, "ny_gene"]
            )
        }
        
        result$log10_gene = log10(result$p1_gene)
        
        result$log10_gene[is.na(result$log10_gene)] <- 0
        
        
        
        print("part5!!!!!!!!!!!!!!!!!!!!!!!!!!")
        # print(str(result))
        
        result
        
    })
    # 
    output$part5_table = DT::renderDataTable({
        df = part5_result() %>%
            dplyr::select(-DiseaseID) %>%
            dplyr::rename(
                `p_(gene_comb)` = p2_gene,
                `p_(chemical_comb)` = p2_chemical,
                p_chemical = p1_chemical,
                p_gene = p1_gene
            )
        
        DT::datatable( df, 
                       class = "display nowrap",
                       selection = "single",
                       options = list(pageLength = 10,scrollX = TRUE)
        )
        
    }) #  %>% bindEvent(input$refresh_heatmap)
    
    output$download_part5_table <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(part5_result(), file, row.names = F)
        }
    )
    # 
    # 5.4 combine
    # 
    heatfigure = reactive({
        req(part5_result())
        
        scale <- part5_result() # read.csv("./output/scale1.csv")
        
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
        
        heatfigure<-data.frame(phenotypeid=phenotypeid,phenotypeTerm=phenotypeTerm) # define the heatfigure df.
        
        # 选了疾病种类，对应的mesh需要相应变化
        
        disease_df = scale|> dplyr::select(DiseaseName, DiseaseID) |> distinct() %>% arrange(DiseaseName)
        
        print(disease_df)
        for(i in 1:nrow(disease_df)) {
            
            disease_name = disease_df[i, 1]
            disease_id = disease_df[i, 2]
            
            heatfigure[disease_name] = scale[scale$DiseaseID==disease_id,]$Wxya_g_c[
                match(heatfigure$phenotypeid,scale[scale$DiseaseID==disease_id,]$GOID)]
            
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
        DT::datatable( heatfigure(), 
                       class = "display nowrap",
                       selection = "single",
                       options = list(pageLength = 10,scrollX = TRUE)
        )
        
    }) # %>% bindEvent(input$refresh_heatmap)
    
    
    
    output$download_part5_table_heatfigure <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(heatfigure(), file, row.names = F)
        }
    )
    
    heatfigure_rank_percentile_react = reactive({
        heatfigure = heatfigure()

        heatfigure_rank<-apply(heatfigure[,3:ncol(heatfigure)],2,rank,ties.method="average",na.last=FALSE)
        # cor.test(heatfigure_rank[,1],heatfigure_rank[,2],method =  "kendall") ##0.7147809
        

        heatfigure_rank_percentile<-heatfigure_rank/nrow(heatfigure_rank)
        heatfigure_rank_percentile<-as.data.frame(heatfigure_rank_percentile)
        
        
        heatfigure_rank_percentile
    })
    
    
    output$part5_heatfigure_plot = renderPlotly({
        
        heatfigure = heatfigure()
        bk <- ncol(heatfigure) -2 # ncol(heatfigure[, c(3, 4, 5, 6, 7)])
        
        cat("bk =", bk, "\n")
        color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
        
        heatfigure_rank<-apply(heatfigure[,3:ncol(heatfigure)],2,rank,ties.method="average",na.last=FALSE)
        # cor.test(heatfigure_rank[,1],heatfigure_rank[,2],method =  "kendall") ##0.7147809
        
        
        cat("nrow(heatfigure_rank) =", nrow(heatfigure_rank), "\n")
        heatfigure_rank_percentile<-heatfigure_rank/nrow(heatfigure_rank)
        heatfigure_rank_percentile<-as.data.frame(heatfigure_rank_percentile)
        
        # browser()
        
        # dev.off()
        # 
        # pheatmap(t(na.omit(heatfigure_rank_percentile)),
        #          legend_breaks = seq(0.000,1,0.2),
        #          legend_labels = c("0","20","40","60","80","100"),
        #          cluster_rows = TRUE,
        #          cluster_cols = TRUE,
        #          show_rownames = TRUE,
        #          show_colnames = FALSE,
        #          cutree_rows = bk,
        #          cutree_cols = min(4, nrow(heatfigure_rank_percentile)),
        #          clustering_distance_rows="euclidean",
        #          clustering_distance_cols = "euclidean",
        #          labels_row = disease_df$DieaseName, # c("Alzheimer Disease","Parkinson Disease","Cognitive Dysfunction",'Cognition Disorders',"Memory Disorders"),
        #          fontsize_col = 5,
        #          angle_col=45,
        #          treeheight_row= 50,
        #          treeheight_col = 50,
        #          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
        #          # cellwidth = 2,
        #          # cellheight = 40,
        #          fontsize=14
        # )
        
        heatmaply(
            t(na.omit(heatfigure_rank_percentile)),
            colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
            k_col = min(4, nrow(heatfigure_rank_percentile)),
            k_row = bk,
            showticklabels = c(F,T),
            # height = 200 + bk * 40

        )
        
    }) # %>% bindEvent(input$refresh_heatmap)
    
    
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
        
        heatfigure$phenotypeTerm[as.numeric(index)] #39
        # heatfigure$max<- apply(heatfigure[,3:ncol(heatfigure)],1,max,na.rm = TRUE)
        
        # hsGO: ref the global object
        
        sim<-mgoSim(GO1 = as.character(heatfigure$phenotypeid[as.numeric(index)]),
                    GO2 = as.character(heatfigure$phenotypeid[as.numeric(index)]), 
                    measure = "Wang",semData = hsGO ,combine=NULL)
        sim[is.na(sim)]<-0
        colnames(sim)<-heatfigure$phenotypeTerm[as.numeric(index)]
        rownames(sim)<-heatfigure$phenotypeTerm[as.numeric(index)]
        
        sim
        
    }) %>%bindEvent(input$refresh_heatmap)
    
    
    output$part5_heatmap_plot2 = renderPlotly({
        sim = sim_df()
    
        heatmaply(sim,
                  # angle_col=315,
                  # cellheight =10,cellwidth = 10,
                  # fontsize_row =12,fontsize_col = 12,
            
        )
        
    }) # %>% bindEvent(input$refresh_heatmap)
    
    output$part5_hopkins_value <- renderValueBox({
        sim = sim_df()
        
        
        #8.2 statistics
        dist<-as.dist(1-sim) %>% as.matrix()
        
        #8.1 eclust package
        
        n = ifelse(ncol(dist) >  5, 5, ncol(dist) - 1)
        
        cluster_test <- get_clust_tendency(dist, n, graph = FALSE)
        cat("Hopkins Statistics: ", cluster_test$hopkins_stat, "\n")

        valueBox(
            value = round(cluster_test$hopkins_stat, 3),
            subtitle = "Hopkins Statistics",
            color = "teal",
            icon = icon("table")
        )
    }) 
    
    
    ##################################
    # AOP
    
    
    observe({
        
        ao<-list()
        #ao 
        
        input_key_event_name_opt = str_squish(input$key_event_name_opt)
        
        
        if(input_key_event_name_opt == "") {
            
            ao$key_event_name <- aop_ke_mie_ao$`key event name` %>% unique()
            
        } else {
            match_cond = str_replace_all(input_key_event_name_opt, " ", "|") 
            
            ao$key_event_name <- aop_ke_mie_ao$`key event name`[grep(match_cond,  # "heart|vascular|hypertension|Coronary",
                                                                     aop_ke_mie_ao$`key event name`,fixed = FALSE,ignore.case = TRUE)] %>% unique()
            
        }
            

        #ao screen
        # aoterm_exclude<-c("Malformation, Male reproductive tract",
        #                   "Impaired development of, Reproductive organs",
        #                   "Decreased testosterone by the fetal Leydig cells, Hypermethylation in the fetal testis",
        #                   "malformed, Male reproductive tract",
        #                   "Decreased fertility, Reduced number of oocytes ovulated",
        #                   "Induction, Male reproductive tract",
        #                   "Impaired inguinoscrotal testicular descent phase",
        #                   "Malformation, cryptorchidism - maldescended testis",
        #                   "Decrease of egg production and cummulative fecundity",
        #                   "Increased intestinal monitor peptide level",
        #                   "Decrease, Fecundity (F3)")
        # ao$key_event_name <- setdiff(ao$key_event_name,aoterm_exclude)#23 
        # ao$aop_id_1<-aop_ke_mie_ao$`aop id`[aop_ke_mie_ao$`key event name` %in% ao$key_event_name] # 61 aopid
        #key_event_male_aop
        # key_event_male_aop <-aop_ke_mie_ao[aop_ke_mie_ao$`aop id` %in% ao$aop_id_1,]
        
        # rewrite below according the 3 lines above:
        key_event_male_aop = aop_ke_mie_ao %>%
            dplyr::filter(`key event name` %in% ao$key_event_name)
        
        #aop screen
        # unique(key_event_male_aop$`aop id`[grep("ovary|ovarian|oocyte|fatal",key_event_male_aop$`key event name`)])#去除女性的生殖通路 15
        # ao$aop_id<-setdiff(ao$aop_id_1,
        #                    unique(key_event_male_aop$`aop id`[grep("ovary|ovarian|oocyte|Ovulation|fetal",key_event_male_aop$`key event name`,ignore.case = TRUE)])) %>%
        #     unique() %>% sort()
        # aop_ids = ao$aop_id
        
        aop_ids= key_event_male_aop$`aop id` %>%  unique() %>% sort()
        
        freezeReactiveValue(input, "aop_id")
        updateSelectInput(
            inputId = "aop_id",
            choices = aop_ids,
            selected = aop_ids)
    })
    
    key_event_male_aop_df = reactive({
        
        req(input$aop_id)
        
        cat('aop_id = ',input$aop_id, "\n")
        
        ao = list()
        
        ao$aop_id = input$aop_id
        heatfigure = heatfigure_max()
        
        req(heatfigure)
        
        
        key_event_male_aop <-aop_ke_mie_ao[aop_ke_mie_ao$`aop id` %in% ao$aop_id,] 
        # ao_relevant_event_ker<-aop_ker[aop_ker$`aop id` %in% ao$aop_id,] 
        key_event_male_aop_ec<-aop_ec[aop_ec$`aop id` %in% ao$aop_id,]    
        
        #6.3 key event mapping
        ias_urogenital_phenotype<-bde_urogenital_phenotype()   # ref
        phenotype_in_aop <- aop_ec[aop_ec$`process ontology id` %in% unique(ias_urogenital_phenotype$phenotypeid),]
        eventbothdir<-intersect(unique(key_event_male_aop$`key event id`),
                                unique(phenotype_in_aop$`key event id`))
        
        key_event_goid<-unique(aop_ec$`process ontology id`[aop_ec$`key event id` %in% key_event_male_aop $`key event id`])
        unique(key_event_goid[grepl("GO",key_event_goid)])#  52 GO
        intersect(phenotype_in_aop$`process ontology id`,unique(key_event_goid[grepl("GO",key_event_goid)]))# 10 matched by GO id
        
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
                judge = ifelse(length(intersect(as.vector(ias_urogenital_phenotype$phenotypeid),
                                                as.vector(children[[as.character(goterm)]])) ) == 0,0, 1)#访问列表记得双括号，此外用名称访问时，记得character化
                mapsuc<-c(mapsuc,judge)
            }
            return(mapsuc)
        }
        
        mapchildren2<- function(eventlist){
            mapsuc<-list()
            for (i in 1:length(eventlist)){
                eventid<-eventlist[i]
                goterm = key_event_male_aop_GO$`process ontology id`[match(eventid,key_event_male_aop_GO$`key event id`)]
                judge = intersect(as.vector(ias_urogenital_phenotype$phenotypeid),
                                  as.vector(children[[as.character(goterm)]]))#访问列表记得双括号，此外用名称访问时，记得character化
                mapsuc[[i]]<- judge
            }
            
            return(mapsuc)
        }
        
        # mapchildren(key_event_male_aop_GO$`key event id`)
        # mapchildren2(key_event_male_aop_GO$`key event id`)
        # sum(lengths(mapchildren2(key_event_male_aop_GO$`key event id`))!=0)
        # sum(mapchildren(key_event_male_aop_GO$`key event id`))
        key_event_male_aop_GO$mapchildren<-mapchildren(key_event_male_aop_GO$`key event id`)
        eventboth<-unique(c(eventbothdir,as.character(key_event_male_aop_GO$`key event id`[key_event_male_aop_GO$mapchildren ==1])))#28
        aop_ke_mie_ao$`key event name`[match(as.character(key_event_male_aop_GO$`key event id`[key_event_male_aop_GO$mapchildren ==1]),
                                             aop_ke_mie_ao$`key event id`)]
        
        
        key_event_male_aop_GO$eventscore<-lapply(mapchildren2(key_event_male_aop_GO$`key event id`), match,heatfigure$phenotypeid) %>% 
            lapply(function(x) heatfigure$max[x]) %>%  
            lapply(function(x) ifelse(length(x)==0,0,max(x))) %>% 
            sapply(unique) %>% unlist
        
        
        #6.4 key event info######################
        key_event_male_aop<-aop_ke_mie_ao[aop_ke_mie_ao$`aop id` %in% ao$aop_id,]
        key_event_male_aop$event_is_in_ctd<-ifelse(key_event_male_aop$`key event id` %in% eventboth,1,0)
        
        #6.5 export
        # write.table(ao_relevant_event_ker,file = "./output/ao_relevant_event_ker.txt",sep = "|",quote = FALSE,row.names = FALSE)
        
        # rski: no key_event_cytoscape table
        # write.table(key_event_cytoscape,file = "./output/key_event_cytoscape.txt",sep = "|",quote = FALSE,row.names = FALSE)
        
        
        
        
        #6.6 
        index<-key_event_male_aop_ec$`process ontology id`[key_event_male_aop_ec$`key event id` %in% eventboth] %>% match(heatfigure$phenotypeid) 
        left<- distinct(heatfigure[na.omit(index),],phenotypeid,.keep_all = TRUE) 
        
        #6.7 eventscore#####
        dim(phenotype_in_aop)
        unique(phenotype_in_aop$`process term`)#16
        colnames(phenotype_in_aop)
        phenotype_in_aop$eventscore<- heatfigure$max[match(phenotype_in_aop$`process ontology id`,heatfigure$phenotypeid)]
        key_event_male_aop$eventscore<-phenotype_in_aop$eventscore[match(key_event_male_aop$`key event id`,phenotype_in_aop$`key event id`)]
        
        
        if(is.null(key_event_male_aop$eventscore)) {
            key_event_male_aop$eventscore = 0
        }
        
        colnames(key_event_male_aop)
        key_event_male_aop$event_children<-key_event_male_aop_GO$eventscore[match(key_event_male_aop$`key event id`,key_event_male_aop_GO$`key event id`)]
        
        if(is.null(key_event_male_aop$event_children)) {
            key_event_male_aop$event_children = 0
        }
        
        key_event_male_aop$eventscore[is.na(key_event_male_aop$eventscore)]<-0
        key_event_male_aop$event_children[is.na(key_event_male_aop$event_children)]<-0
        key_event_male_aop$max<-apply(key_event_male_aop[,c(6,7)],1,max)
        key_event_male_aop$eventmaxscore<-key_event_male_aop$max[match(key_event_male_aop$`key event id`,key_event_male_aop$`key event id`)]
        
        
        key_event_male_aop %>%
            dplyr::rename(`event is in ctd` = event_is_in_ctd, 
                          `event score` = eventscore, 
                          `event children` = event_children,
                          `event max score` = eventmaxscore)
    })
    
    

    output$aop_key_event_male_aop = DT::renderDataTable({
        DT::datatable( key_event_male_aop_df(), 
                       class = "display nowrap",
                       selection = "single",
                       options = list(pageLength = 10,scrollX = TRUE)
        )
        
    }) #  %>% bindEvent(input$refresh_heatmap)
    
    output$download_aop_event_male <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(key_event_male_aop_df(), file, row.names = F)
        }
    )
    
    happy_df = reactive({
        key_event_male_aop = key_event_male_aop_df()
        
        req(key_event_male_aop)
        
        
        happy <- key_event_male_aop %>%
            group_by(`aop id`) %>%
            summarize(`total count` = n(),
                      proportion = sum(`event is in ctd` == 1) / `total count`,
                      `total score` = sum(case_when(`event is in ctd` == 1 ~ max, TRUE ~ 0)))
        
    })
    
    output$aop_happy = DT::renderDataTable({
        
        happy = happy_df()
        
        DT::datatable( happy, 
                       class = "display nowrap",
                       selection = "single",
                       options = list(pageLength = 10,scrollX = TRUE)
        )
        
        

        
    }) 
    
    output$download_aop_happy <- downloadHandler(
        # 参数filename接受一个文件名的字符串 可以是函数的返回(Reactive values and functions may be used from this function.)
        filename = function() {
            paste0("table.csv")
        },
        # 参数content 固定的用法 file接收filename的字符
        content = function(file) {
            write.csv(happy_df(), file, row.names = F)
        }
    )
    
    
    
    
}


shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port=3675))

