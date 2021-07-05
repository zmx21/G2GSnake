library(qqman)
library(filematrix)
library(grid)
library(shinydashboard)
library(shiny)
library(plotly)
library(data.table)
library(shinyWidgets)
library(latex2exp)

g2g_file_path <- dir('../../results/')[sapply(dir('../../results/'),function(x) grepl(pattern = '.allchr.txt',x = x))]
g2g_file_path <- g2g_file_path[sapply(g2g_file_path,function(x) strsplit(system(glue::glue("wc -l ../../results/{x}"),intern = T),split = ' ')[[1]][1]> 1)]

pathogen_variants <- sapply(g2g_file_path,function(x) strsplit(x=strsplit(x=x,'.allchr.txt')[[1]][1],split = '\\.')[[1]][1],USE.NAMES = F)
fm <- fm.open('../../results/G2G_Results')
pathogen_info <- data.table::fread('../../raw_data/pathogen/OUT_aa_binary_table_info.txt') %>% dplyr::filter(ID %in% pathogen_variants)

## ui.R ##
manhattan_page <- fluidPage(
  selectInput(inputId = "gene",label = 'Pathogen Gene',choices = unique(pathogen_info$GENE)),
  selectInput(inputId = "pos",label = 'Pathogen Position',choices = unique(pathogen_info$POS)),
  selectInput(inputId = "variant",label = 'Pathogen Variant',choices = pathogen_info$ID),
  mainPanel(tabsetPanel(tabPanel('Manhattan',plotlyOutput(outputId = "manhattan")),
         tabPanel('QQ',plotOutput(outputId = "qq")))))
tbl_page <- fluidPage(
  numericInput(inputId = 'p_thresh',label = 'P-Value Threshold:',value = signif(5e-8 / length(pathogen_variants),3),max = 1,min = 0,step = 1e-2),
  actionButton(inputId = 'submit_p',label = 'Submit'),
  dataTableOutput(outputId = "tbl"))

body <- dashboardBody(tabItems(tabItem(tabName = "manhattan",manhattan_page),
                               tabItem(tabName = "tbl",tbl_page)))
side_bar <- dashboardSidebar(sidebarMenu(menuItem("Results Table", tabName = "tbl", icon = icon("dashboard")),
                             menuItem("Manhattan Plot", tabName = "manhattan", icon = icon("dashboard"))))
ui <- dashboardPage(dashboardHeader(title = 'G2G Results'),side_bar,body)


server <- function(input, output) {
  observe({
    updateSelectInput(inputId = "pos", 
                      choices = pathogen_info$POS[pathogen_info$GENE == input$gene])
  })
  
  observe({
    updateSelectInput(inputId = "variant", 
                      choices = pathogen_info$ID[pathogen_info$GENE == input$gene & pathogen_info$POS == input$pos])
  })
  
  
  re <- eventReactive(input$submit_p,{
    sig_variants <- c()
    for(i in 1:ncol(fm)){
      if(min(fm[,i],na.rm = T) < input$p_thresh){
        sig_variants <- c(sig_variants,colnames(fm)[i])
      }
    }
    g2g_results <- lapply(g2g_file_path[match(sig_variants,pathogen_variants)],function(x){
      df <- data.table::fread(cmd = glue::glue("awk \'{{ if ($13 <= {input$p_thresh} || NR == 1 ) {{ print }} }}\' ../../results/{x}"),
                              select = c('CHR','POS','SNPID','Allele1','Allele2','BETA','SE','p.value','Is.SPA.converge')) %>%
        dplyr::filter(Is.SPA.converge == 1) %>% 
        dplyr::select(-Is.SPA.converge) %>% 
        dplyr::relocate(Host_SNP = SNPID,Host_Chr = CHR,Host_Pos = POS,Host_Allele1=Allele1,Host_Allele2=Allele2,BETA,SE,p.value)
      
      return(df)
    })
    names(g2g_results) <- sig_variants
    return(g2g_results)
  })
  output$tbl <- renderDataTable(data.table::rbindlist(re(),idcol = T) %>%
                                    dplyr::rename(Pathogen_Variant = `.id`) %>% dplyr::left_join(pathogen_info %>% dplyr::select(Pathogen_Variant=ID,Pathogen_Gene = GENE,Pathogen_Pos = POS)) %>%
                                    dplyr::relocate(Pathogen_Gene,Pathogen_Pos,Pathogen_Variant) %>% dplyr::arrange(p.value))
  
  output$manhattan <- renderPlotly({
    cur_result <- data.table::fread(input = glue::glue('../../results/{input$variant}.allchr.txt')) %>% dplyr::filter(p.value < 0.1)
    manhattanly::manhattanly(cur_result %>% dplyr::select(CHR=CHR,BP=POS,P=p.value,ID=SNPID),snp = 'ID',title = input$variant,genomewideline = -log10(5e-8 / length(pathogen_variants)),suggestiveline = -log10(5e-8))})
  output$qq <- renderPlot({
    cur_result <- as.vector(fm[,which(colnames(fm)==input$variant)])
    chisq1 <- qchisq(1-cur_result,1)
    lambda1 = median(chisq1)/qchisq(0.5,1)
    
    qqman::qq(cur_result,main = input$variant)
    usr = par('usr')
    text(x=1,y=usr[4] - 1,label = TeX(sprintf("$\\lambda = %.3f$", lambda1)))
    
    })
}
shinyApp(ui = ui, server = server)
