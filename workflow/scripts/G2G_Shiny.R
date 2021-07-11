library(qqman)
library(filematrix)
library(grid)
library(shinydashboard)
library(shiny)
library(plotly)
library(data.table)
library(shinyWidgets)
library(latex2exp)
library(dplyr)
library(ggtree)

g2g_file_path <- dir('../../results/')[sapply(dir('../../results/'),function(x) grepl(pattern = '.allchr.txt',x = x))]
g2g_file_path <- g2g_file_path[sapply(g2g_file_path,function(x) strsplit(system(glue::glue("wc -l ../../results/{x}"),intern = T),split = ' ')[[1]][1]> 1)]

pathogen_variants <- sapply(g2g_file_path,function(x) strsplit(x=strsplit(x=x,'.allchr.txt')[[1]][1],split = '\\.')[[1]][1],USE.NAMES = F)
fm <- fm.open('../../results/G2G_Results')
pathogen_info <- data.table::fread('../../raw_data/pathogen/OUT_aa_binary_table_filt_info.txt') %>% dplyr::filter(ID %in% pathogen_variants)
host_snps <- data.table::fread('../../results/tmp/G2G_QC.bim')


VisualizeG2G <- function(tree_file,covar_file,type = 'Pathogen',variant){
  tree <- ape::read.tree(tree_file)

  if(type == 'Pathogen'){
    dosage_vect <- data.table::fread(glue::glue('../../results/tmp/pathogen_variants/{variant}')) %>% 
      dplyr::select(IID,pathogen_dosage = any_of(variant)) %>% 
      dplyr::left_join(data.table::fread(covar_file,select = c('PID','IID')),by=c('IID'='IID')) %>% 
      dplyr::filter(!is.na(pathogen_dosage)) %>% dplyr::select(ID=PID,pathogen_dosage)
    tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,dosage_vect$ID))
    
    p1 <- ggtree(tree_filt)
    metat <- p1$data %>%
      dplyr::inner_join(dosage_vect, c('label' = 'ID'))
    p2 <- p1 +
      geom_point(data = metat,
                 aes(x = x,
                     y = y,
                     colour = factor(pathogen_dosage),
                     label = label)) + ggtitle(variant) + labs(color = '')
    return(plotly::ggplotly(p2))
    
  }else{
    host_genotype <- snpStats::read.plink(bed = '../../results/tmp/G2G_QC.bed',select.snps = variant)
    host_dosage_raw <- as(host_genotype$genotypes, Class = 'numeric')
    host_dosage <- as.vector(host_dosage_raw)  
    names(host_dosage) <- rownames(host_dosage_raw)
    host_dosage[host_dosage == 2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
    host_dosage[host_dosage == 1] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.1)
    host_dosage[host_dosage == 0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
    
    dosage_vect <- data.frame(IID = names(host_dosage),host_dosage = host_dosage) %>% 
      dplyr::left_join(data.table::fread(covar_file,select = c('PID','IID')),by=c('IID'='IID'))%>% 
      dplyr::filter(!is.na(host_dosage)) %>% dplyr::select(ID=PID,host_dosage)
    tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,dosage_vect$ID))
    
    p1 <- ggtree(tree_filt)
    metat <- p1$data %>%
      dplyr::inner_join(dosage_vect, c('label' = 'ID'))
    p2 <- p1 +
      geom_point(data = metat,
                 aes(x = x,
                     y = y,
                     colour = factor(host_dosage),
                     label = label)) + ggtitle(variant) + labs(color = '')
    return(plotly::ggplotly(p2))
    
  }
  
}


## ui.R ##
manhattan_page <- fluidPage(
  selectizeInput(inputId = "gene",label = 'Pathogen Gene',choices = unique(pathogen_info$GENE)),
  selectizeInput(inputId = "variant",label = 'Pathogen Variant',choices = pathogen_info$ID),
  mainPanel(tabsetPanel(tabPanel('Manhattan',plotlyOutput(outputId = "manhattan")),
         tabPanel('QQ',plotOutput(outputId = "qq")))))
tbl_page <- fluidPage(
  numericInput(inputId = 'p_thresh',label = 'P-Value Threshold:',value = signif(5e-8 / length(pathogen_variants),3),max = 1,min = 0,step = 1e-2),
  actionButton(inputId = 'submit_p',label = 'Submit'),
  dataTableOutput(outputId = "tbl"))

phylo_page <- fluidPage(
  radioButtons(inputId = 'var_type',label = 'Variant Type:',choices = c('Pathogen','Host')),
  selectizeInput(inputId = 'phylo_gene',label = 'Gene:',choices = NULL),
  selectizeInput(inputId = 'phylo_variant',label = 'Variant ID:',choices = NULL),
  plotlyOutput(outputId = "tree")
)

body <- dashboardBody(tabItems(tabItem(tabName = "manhattan",manhattan_page),
                               tabItem(tabName = "tbl",tbl_page),
                               tabItem(tabName = "phylo",phylo_page)))
side_bar <- dashboardSidebar(sidebarMenu(menuItem("Results Table", tabName = "tbl", icon = icon("dashboard")),
                             menuItem("Manhattan Plot", tabName = "manhattan", icon = icon("dashboard")),
                             menuItem("Phylogenetic Tree", tabName = "phylo", icon = icon("dashboard"))))
ui <- dashboardPage(dashboardHeader(title = 'G2G Results'),side_bar,body)


server <- function(input, output) {

  observe({
    updateSelectInput(inputId = "variant", 
                      choices = pathogen_info$ID[pathogen_info$GENE == input$gene])
  })
  
  observe({
    ifelse(input$var_type == 'Pathogen',choices <- unique(pathogen_info$GENE),choices <- as.character(unique(host_snps$V1)))
    ifelse(input$var_type == 'Pathogen',label <- 'Pathogen Gene',label <- 'Host Chr')
    updateSelectizeInput(inputId = "phylo_gene", 
                         choices = choices,label = label,server = TRUE)
  })
  
  observe({
    ifelse(input$var_type == 'Pathogen',choices <- pathogen_info$ID[pathogen_info$GENE == input$phylo_gene],choices <- host_snps$V2[as.character(host_snps$V1) == input$phylo_gene])
    ifelse(input$var_type == 'Pathogen',label <- 'Pathogen Variant ID',label <- 'Host SNP ID')
    
    updateSelectizeInput(inputId = "phylo_variant", 
                      choices = choices,label = label,server = TRUE)
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
  
  output$tree <- renderPlotly(VisualizeG2G('../../raw_data/pathogen/HBV_WG+og_rr.nw',
                                           '../../results/tmp/merged_covar.txt',
                                           type = input$var_type,
                                           variant = input$phylo_variant))
}
shinyApp(ui = ui, server = server)




