#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(qqman)
library(filematrix)
library(shinydashboard)
library(shiny)
library(plotly)
library(data.table)
library(readr)
library(shinyWidgets)
library(latex2exp)
library(dplyr)
library(corrplot)

#Get genes that were analyzed
results_dir <- '/home/results/'
genes <- dir(results_dir)[grepl(dir(results_dir),pattern = 'Gene_')]

#Get results
g2g_file_paths <- lapply(genes,function(x) dir(paste0(results_dir,x,'/'))[sapply(dir(paste0(results_dir,x,'/')),function(x) grepl(pattern = '.allchr.txt',x = x))])
names(g2g_file_paths) <- genes
pathogen_variants <- lapply(g2g_file_paths,function(x) sapply(x,function(y) strsplit(x=y,'.allchr.txt')[[1]][1],USE.NAMES = F))
all_g2g_file_paths <- unlist(lapply(1:length(g2g_file_paths),function(i) paste0(names(g2g_file_paths)[i],'/',g2g_file_paths[[i]])))

#Open file matrix
fm <- fm.open(glue::glue("{results_dir}G2G_Results"))

#Get variant info
pathogen_info <- do.call(rbind,lapply(genes,function(x) data.table::fread(glue::glue("{results_dir}tmp/{x}/{x}_QC.info")) %>% dplyr::filter(ID %in% unlist(pathogen_variants))))
pathogen_info$Path <- paste0(results_dir,sapply(1:nrow(pathogen_info),function(i) all_g2g_file_paths[grepl(all_g2g_file_paths,pattern = paste0(pathogen_info$ID[i],'.allchr.txt')) & grepl(all_g2g_file_paths,pattern = pathogen_info$Gene[i])]))

#Get host SNPs
host_snps <- data.table::fread(glue::glue("{results_dir}tmp/G2G_QC.bim"))

#Plot pathogen gene position against host chr position, colors as p-value
VisualizeG2G <- function(pathogen_info,pathogen_gene,host_snps,fm,p_thresh){
  cur_gene_varants <- pathogen_info$ID[pathogen_info$Gene == pathogen_gene]
  col_to_keep <- sapply(which(colnames(fm) %in% cur_gene_varants),function(i) any(fm[,i] < p_thresh,na.rm = T))
  fm_col_filt <- fm[,which(colnames(fm) %in% cur_gene_varants)][,which(col_to_keep),drop=F]
  row_to_keep <- apply(fm_col_filt,1,function(x) any(x < p_thresh,na.rm = T))
  fm_both_filt <- fm_col_filt[row_to_keep,,drop=F]
  if(nrow(fm_both_filt) == 0 | ncol(fm_both_filt) == 0){
    return(ggplot() + 
             labs(
               title = "No associations below P-value threshold",
               subtitle = "",
               caption = "none"
             ) +
             theme(panel.background = element_blank())
    )
  }
  
  colnames(fm_both_filt) <-colnames(fm)[which(colnames(fm) %in% cur_gene_varants)][col_to_keep]
  rownames(fm_both_filt) <-rownames(fm)[row_to_keep]
  
  df <- na.omit(reshape2::melt(fm_both_filt)) %>% dplyr::select(Host_SNP=Var1,Pathogen_Variant = Var2,P=value) %>% dplyr::filter(P < p_thresh)
  df_joined <- df %>% dplyr::left_join(pathogen_info %>% dplyr::select(ID,Pathogen_Gene=Gene,Pathogen_Pos = Pos,-Path),by=c('Pathogen_Variant' = 'ID')) %>%
    dplyr::left_join(host_snps %>% dplyr::select(Host_Chr = V1,rsid=V2,Host_Pos=V4),by=c('Host_SNP' = 'rsid')) %>% dplyr::mutate(Association = paste(Pathogen_Variant,'~',Host_SNP))
  p <- ggplot(df_joined,aes(x=Host_Pos,y = Pathogen_Pos,color = -log10(P),label = Association)) + facet_grid(~factor(Host_Chr,levels = unique(host_snps$V1))) + geom_point() +
    ylab(paste0('Pathogen Gene ',pathogen_gene,' Position')) + xlab('Host Chr Position') +  ylim(0,max(pathogen_info$Pos[pathogen_info$Gene==pathogen_gene])) +
    theme_bw() + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  return(ggplotly(p))
}



#Beta - Plot phylogenetic tree
# VisualizePhylo <- function(tree_file,covar_file,type = 'Pathogen',variant){
#   tree <- ape::read.tree(tree_file)
# 
#   if(type == 'Pathogen'){
#     variant_gene <- pathogen_info$Gene[pathogen_info$ID == variant]
#     dosage_vect <- data.table::fread(glue::glue("{results_dir}tmp/Gene_{variant_gene}/{variant_gene}_QC.txt")) %>% 
#       dplyr::select(IID,pathogen_dosage = any_of(variant)) %>% 
#       dplyr::left_join(covar_file %>% dplyr::select('PID','IID'),by=c('IID'='IID')) %>% 
#       dplyr::filter(!is.na(pathogen_dosage)) %>% dplyr::select(ID=PID,pathogen_dosage)
#     tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,dosage_vect$ID))
#     
#     p1 <- ggtree(tree_filt)
#     metat <- p1$data %>%
#       dplyr::inner_join(dosage_vect, c('label' = 'ID'))
#     p2 <- p1 +
#       geom_point(data = metat,
#                  aes(x = x,
#                      y = y,
#                      colour = factor(pathogen_dosage),
#                      label = label)) + ggtitle(variant) + labs(color = '')
#     return(plotly::ggplotly(p2))
#     
#   }else{
#     host_genotype <- snpStats::read.plink(bed = glue::glue("{results_dir}tmp/G2G_QC.bed"),select.snps = variant)
#     host_dosage_raw <- as(host_genotype$genotypes, Class = 'numeric')
#     host_dosage <- as.vector(host_dosage_raw)  
#     names(host_dosage) <- rownames(host_dosage_raw)
#     host_dosage[host_dosage == 2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
#     host_dosage[host_dosage == 1] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.1)
#     host_dosage[host_dosage == 0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
#     
#     dosage_vect <- data.frame(IID = names(host_dosage),host_dosage = host_dosage) %>% 
#       dplyr::left_join(covar_file %>% dplyr::select('PID','IID'),by=c('IID'='IID'))%>% 
#       dplyr::filter(!is.na(host_dosage)) %>% dplyr::select(ID=PID,host_dosage)
#     tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,dosage_vect$ID))
#     
#     p1 <- ggtree(tree_filt)
#     metat <- p1$data %>%
#       dplyr::inner_join(dosage_vect, c('label' = 'ID'))
#     p2 <- p1 +
#       geom_point(data = metat,
#                  aes(x = x,
#                      y = y,
#                      colour = factor(host_dosage),
#                      label = label)) + ggtitle(variant) + labs(color = '')
#     return(plotly::ggplotly(p2))
#     
#   }
#   
# }

#pPC vs PC plot
GeneratePCPlot <- function(pathogen_gene){
  PCs <- data.table::fread(glue::glue('{results_dir}tmp/{pathogen_gene}/merged_covar_{pathogen_gene}.txt')) %>%
    dplyr::select(contains('PC')) 
  PC_Cor <- cor(PCs)
  return(corrplot::corrplot(PC_Cor,method = 'circle',type = 'upper',title = paste0('pPC and PC Correlation - ',pathogen_gene),mar=c(0,0,1,0)))
}

## ui.R ##
results_plot_page <- fluidPage(
  selectizeInput(inputId = "result_gene",label = 'Pathogen Gene',choices = unique(pathogen_info$Gene)),
  numericInput(inputId = 'p_thresh_result',label = 'P-Value Threshold:',value = signif(5e-8,3),max = 1,min = 0,step = 1e-2),
  actionButton(inputId = 'submit_p_plot',label = 'Submit'),
  plotlyOutput(outputId = "results_plot")
)

manhattan_page <- fluidPage(
  selectizeInput(inputId = "manhattan_gene",label = 'Pathogen Gene',choices = unique(pathogen_info$Gene)),
  selectizeInput(inputId = "manhattan_variant",label = 'Pathogen Variant',choices = pathogen_info$ID),
  mainPanel(tabsetPanel(tabPanel('Manhattan',plotlyOutput(outputId = "manhattan")),
                        tabPanel('QQ',plotOutput(outputId = "qq")))))
tbl_page <- fluidPage(
  numericInput(inputId = 'p_thresh',label = 'P-Value Threshold:',value = signif(5e-8 / length(unlist(pathogen_variants)),3),max = 1,min = 0,step = 1e-2),
  actionButton(inputId = 'submit_p',label = 'Submit'),
  dataTableOutput(outputId = "tbl"))

PC_plot_page <- fluidPage(
  selectizeInput(inputId = "PC_Gene",label = 'Pathogen Gene',choices = unique(pathogen_info$Gene)),
  plotOutput(outputId = "PC_plot",height = '600px',width = '800px')
)

body <- dashboardBody(tabItems(tabItem(tabName = "results_plot",results_plot_page),
                               tabItem(tabName = "manhattan",manhattan_page),
                               tabItem(tabName = "tbl",tbl_page),
                               tabItem(tabName = "PC",PC_plot_page)))
side_bar <- dashboardSidebar(sidebarMenu(menuItem("Results Plot", tabName = "results_plot", icon = icon("dashboard")),
                                         menuItem("Results Table", tabName = "tbl", icon = icon("dashboard")),
                                         menuItem("Manhattan Plot", tabName = "manhattan", icon = icon("dashboard")),
                                         menuItem("PC vs pPC Plot", tabName = "PC", icon = icon("dashboard"))))

ui <- dashboardPage(dashboardHeader(title = 'G2G Results'),side_bar,body)


server <- function(input, output) {
  observe({
    updateSelectInput(inputId = "manhattan_variant", 
                      choices = pathogen_info$ID[pathogen_info$Gene == input$manhattan_gene])
  })
  
  observe({
    updateSelectInput(inputId = "PC_Gene",choices = unique(pathogen_info$Gene))
  })
  

  re <- eventReactive(input$submit_p,{
    sig_variants <- c()
    sig_variants_host_index <- list()
    
    for(i in 1:ncol(fm)){
      if(min(fm[,i],na.rm = T) < input$p_thresh){
        sig_variants <- c(sig_variants,colnames(fm)[i])
        cur_host_index <- which(fm[,i] < input$p_thresh)
        sig_variants_host_index <- c(sig_variants_host_index,list(cur_host_index))
        
      }
    }
    if(is.null(sig_variants)){
      df <- list(data.frame(Host_SNP= 'No associations below p-value threshold',Host_Chr = '',
                            Host_Pos = '',Host_Allele1 = '',Host_Allele2 = '',
                            BETA = '',SE = '',p.value = ''))
      names(df) <- ''
      return(df)
    }
    if(length(unlist(sig_variants_host_index)) > 500){
      df <- list(data.frame(Host_SNP= 'Cannot Display - more than 500 entires',Host_Chr = '',
                            Host_Pos = '',Host_Allele1 = '',Host_Allele2 = '',
                            BETA = '',SE = '',p.value = ''))
      names(df) <- ''
      return(df)
    }
    ReadResultsFile <- function(path,index){
        header <- colnames(data.table::fread(cmd = paste0('head -n 2 ',path)))
        if('p.value' %in% header){
            return(read.table(text = readr::read_lines(path,skip = index,n_max = 1)) %>%
                dplyr::select(CHR = V1,POS=V2,SNPID=V3,Allele1=V4,Allele2=V5,BETA = V10,SE = V11,p.value = V13,Is.SPA.converge = V15) %>%
                dplyr::filter(Is.SPA.converge == 1) %>% 
                dplyr::select(-Is.SPA.converge) %>% 
                dplyr::relocate(Host_SNP = SNPID,Host_Chr = CHR,Host_Pos = POS,
                                Host_Allele1=Allele1,Host_Allele2=Allele2,BETA,SE,p.value))
        }else if ('P' %in% header){
            return(read.table(text = readr::read_lines(path,skip = index,n_max = 1)) %>%
                dplyr::select(CHR = V1,POS=V2,SNPID=V3,REF=V4,ALT=V5,OR = V10,LOG_OR_SE = V11,p.value = V13) %>% 
                    dplyr::relocate(Host_SNP = SNPID,Host_Chr = CHR,Host_Pos = POS,
                                    Host_Allele1=REF,Host_Allele2=ALT,OR,LOG_OR_SE,p.value))
            
        }else if('LOG10P' %in% header){
          return(read.table(text = readr::read_lines(path,skip = index,n_max = 1)) %>%
                   dplyr::mutate(P=10^(-1*as.numeric(V13))) %>%
                   dplyr::select(CHR = V1,POS=V2,SNPID=V3,REF=V4,ALT=V5,BETA = V10,SE = V11,p.value = P) %>% 
                   dplyr::relocate(Host_SNP = SNPID,Host_Chr = CHR,Host_Pos = POS,
                                   Host_Allele1=REF,Host_Allele2=ALT,BETA,SE,p.value))
          
        }
    }
    
    g2g_results <- lapply(1:length(sig_variants),function(i){
      cur_path <- pathogen_info$Path[match(sig_variants[i],pathogen_info$ID)]
      cur_index <- sig_variants_host_index[[i]]
      df <- do.call(rbind,lapply(cur_index,function(x) ReadResultsFile(cur_path,x)))
      return(df)
    })
    names(g2g_results) <- sig_variants
    return(g2g_results)
  })
  output$tbl <- renderDataTable(data.table::rbindlist(re(),idcol = T) %>%
                                  dplyr::rename(Pathogen_Variant = `.id`) %>% dplyr::left_join(pathogen_info %>% dplyr::select(Pathogen_Variant=ID,Pathogen_Gene = Gene,Pathogen_Pos = Pos)) %>%
                                  dplyr::relocate(Pathogen_Gene,Pathogen_Pos,Pathogen_Variant) %>% dplyr::arrange(p.value), options = list(scrollX = TRUE))
  
  output$manhattan <- renderPlotly({
    cur_result <- data.table::fread(input = pathogen_info$Path[pathogen_info$ID == input$manhattan_variant]) 
    if('p.value' %in% colnames(cur_result)){
        cur_result <- cur_result %>% dplyr::filter(p.value < 0.1)
        manhattanly::manhattanly(cur_result %>% dplyr::select(CHR=CHR,BP=POS,P=p.value,ID=SNPID),snp = 'ID',title = input$manhattan_variant,genomewideline = -log10(5e-8 / length(unlist(pathogen_variants))),suggestiveline = -log10(5e-8))
    }else if('P' %in% colnames(cur_result)){
        cur_result <- cur_result %>% dplyr::filter(P < 0.1)
        manhattanly::manhattanly(cur_result %>% dplyr::select(CHR='#CHROM',BP=POS,P=P,ID=ID),snp = 'ID',title = input$manhattan_variant,genomewideline = -log10(5e-8 / length(unlist(pathogen_variants))),suggestiveline = -log10(5e-8))
    }else if('LOG10P' %in% colnames(cur_result)){
      cur_result <- cur_result %>% dplyr::mutate(P=10^(-1*as.numeric(LOG10P)))
      cur_result <- cur_result %>% dplyr::filter(P < 0.1)
      manhattanly::manhattanly(cur_result %>% dplyr::select(CHR=CHROM,BP=GENPOS,P=P,ID=ID),snp = 'ID',title = input$manhattan_variant,genomewideline = -log10(5e-8 / length(unlist(pathogen_variants))),suggestiveline = -log10(5e-8))
    }})
    
  output$qq <- renderPlot({
    cur_result <- as.vector(fm[,which(colnames(fm)==input$manhattan_variant)])
    chisq1 <- qchisq(1-cur_result,1)
    lambda1 = median(chisq1)/qchisq(0.5,1)
    
    qqman::qq(cur_result,main = input$manhattan_variant)
    usr = par('usr')
    text(x=1,y=usr[4] - 1,label = TeX(sprintf("$\\lambda = %.3f$", lambda1)))
    
  })
  
  output$PC_plot <- renderPlot({
    GeneratePCPlot(input$PC_Gene)
  })
  
  re_res_plot <- eventReactive(input$submit_p_plot,{
    res_plot <- VisualizeG2G(pathogen_info,input$result_gene,
                             host_snps,fm,p_thresh = input$p_thresh_result)
  })
  
  output$results_plot <- renderPlotly(re_res_plot())
  
}
shinyApp(ui = ui, server = server)
