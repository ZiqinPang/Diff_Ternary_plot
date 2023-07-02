---
title: "Diff_Ternary_plot v1.0"
author: "Ziqin Pang, Yongmei Zhou"
date: "2023/7/2"
output: html_document
editor_options: 
  chunk_output_type: console
---

options(stringsAsFactors = FALSE)
set.seed(100)
# Load package and initialization ---------------------------------------------------------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(Cairo))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(qvalue))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(tools))
suppressMessages(library(patchwork))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(grid))
suppressMessages(library(xlsx))
suppressMessages(library(UpSetR))
suppressMessages(library(configr))
suppressMessages(library(ggtern))
suppressMessages(library(writexl))
#BiocManager::install(c("ggplot2", "rlang", "Cairo"))
# Argument importation ------------------------------------------------------------------------------------------
input <- 'I22430BAC.txt'  #微生物的定量文件填写

group <- 'I22430GROUP.txt' #微生物样本的分组信息填写


#Diff='M12_vs_I20;M12_vs_I24;I20_vs_I24' #差异分析组合（3个）填写，用分号隔开
#Diff='M12_vs_I24;M12_vs_I30;I24_vs_I30' #差异分析组合（3个）填写，用分号隔开
Diff='MR12_vs_IR20;MR12_vs_IR24;IR20_vs_IR24' #差异分析组合（3个）填写，用分号隔开
#Diff='MR12_vs_IR24;MR12_vs_IR30;IR24_vs_IR30' #差异分析组合（3个）填写，用分号隔开


outdir<-'C:\\Users\\86133\\Desktop\\1' #运行脚本的文件夹路径填写F:\\桌面\\12

#ternary<-"M12,I20,I24" #三元相图的三个组名按顺序填写
#ternary<-"M12,I24,I30" #三元相图的三个组名按顺序填写
ternary<-"MR12,IR20,IR24" #三元相图的三个组名按顺序填写
#ternary<-"IR24,MR12,IR30" #三元相图的三个组名按顺序填写
setwd(outdir)

# Function defination ---------------------------------------------------------------------------------------------
loadTable <- function(inFile){
  NAVector <- c('', '#N/A', '#N/A', 'N/A', '#NA', '-1.#IND', '-1.#QNAN',
                '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA',
                'NULL', 'NaN', '<na>', 'n/a', 'na', 'null', 'nan')
  if (! file.exists(inFile)){
    stop(paste0('Input file do not exists: ', inFile))
  }
  fileFormat <- readxl::format_from_signature(inFile)
  if (fileFormat %in% c('xlsx', 'xls')){
    dataFrame <- read_excel(inFile, trim_ws = TRUE, na = NAVector)
  } else if (endsWith(inFile, '.csv')){
    tryCatch({
      dataFrame <- read_csv(inFile, trim_ws = TRUE,
                            locale = locale(encoding = "UTF-8"),
                            na = NAVector)
    }, error = function(e) {
      stop(paste0(inFile, ' should be comma (,) separated plain txt file'))
    })
  } else {
    tryCatch({
      dataFrame <- read_tsv(inFile, trim_ws = TRUE,
                            locale = locale(encoding = "UTF-8"),
                            na = NAVector)
    }, error = function(e) {
      stop(paste0(inFile, ' should be tab separated plain txt file'))
    })
  }
  if (nrow(dataFrame) == 0){
    stop('No data contained: ', inFile)
  }
  return(dataFrame)
}

loadSampleGroup <- function(inFile){
  # Load file
  sampleFrame <- loadTable(inFile)
  # Columns should contain Sample and Group
  checkColumnVector <- colnames(sampleFrame)
  if (! all(c('names', 'group') %in% checkColumnVector)){
    stop(paste0('Column names should contain [Sample, Group], ',
                'case sensitive.'))
  }
  # Sample and group should not contains NA
  if (sum(is.na(sampleFrame$names)) > 0){
    stop('Sample should not contains NA: ', inFile)
  }
  if (sum(is.na(sampleFrame$group)) > 0){
    stop('Group should not contains NA: ', inFile)
  }
  # Sample name shoule not repeat
  sampleVector <- sampleFrame$names
  repeatSampleVector <- sampleVector[duplicated(sampleVector)]
  if (length(repeatSampleVector) > 0){
    stop('Sample should not contains repeat in ', inFile, ': ',
         paste(repeatSampleVector, collapse = ', '))
  }
  # Check additional grouping
  additionalGroups <- checkColumnVector[grepl('^group[0-9+]$', checkColumnVector)]
  if (length(additionalGroups) > 0){
    # Check NA
    for (addGroup in additionalGroups){
      if (sum(is.na(sampleFrame[, addGroup])) > 0){
        stop(addGroup, ' should not contains NA: ', inFile)
      }
    }
    # Check if same group names if exist, and samples of these are the same
    ### Generate a list of standard Group and Sample
    groupSampleFrame <- sampleFrame[, c('Sample', 'Group')]
    groupSampleList <- split(groupSampleFrame, groupSampleFrame[, 'Group'])
    groupSampleList <- lapply(groupSampleList, as.list)
    ### Check other groups
    for (addGroup in additionalGroups){
      ### Generate group:sample list
      addGroupFrame <- sampleFrame[, c('Sample', addGroup)]
      addGroupList <- split(addGroupFrame, addGroupFrame[, addGroup])
      addGroupList <- lapply(addGroupList, as.list)
      #### Check if sample list of same group are the same
      for (groupName in names(addGroupList)){
        if (! (groupName %in% sampleFrame$Group)){
          groupSampleList[[groupName]] <- addGroupList[[groupName]]
        } else {
          mainGroupSampleVector <- groupSampleList[[groupName]]$Sample
          addGroupSampleVector <- addGroupList[[groupName]]$Sample
          if (length(mainGroupSampleVector) != length(addGroupSampleVector) ||
              all(sort(mainGroupSampleVector) != sort(addGroupSampleVector))){
            # Either diff length or same length but diff element
            stop(paste0('This group relate to different samples in ',
                        'different grouping: ', groupName))
          }
        }
      }
    }
  }
  # Convert to sample frame to character, and Group to factor
  sampleFrame <- as.data.frame(apply(sampleFrame, 2, as.character))
  sampleFrame <- sampleFrame %>%
    mutate(group = factor(group, levels = unique(group)))
  if (length(additionalGroups) > 0){
    for (addGroup in additionalGroups){
      uniqueGroup <- sampleFrame[, addGroup]
      sampleFrame[, addGroup] <- factor(sampleFrame[, addGroup],
                                        levels = unique(uniqueGroup))
    }
  }
  return(sampleFrame)
}

normalizeSampleFrame <- function(dataFrame, diffVector){
  checkGroupVector <- apply(dataFrame, 2, function(x){all(diffVector %in% x)})
  matchedGroupVector <- names(checkGroupVector[checkGroupVector])
  matchedGroupVector <- matchedGroupVector[matchedGroupVector != 'names']
  if (length(matchedGroupVector) == 0){
    string = 'There is no column of Group* exist these diff groups: '
    stop(paste0(string, paste0(diffVector, collapse = ', ')))
  } else {
    matchedGroup <- matchedGroupVector[1]
    dataFrame <- dataFrame[, c('names', matchedGroup)]
    colnames(dataFrame)[2] <- 'group'
  }
  return(dataFrame)
}

diffAnalysis <- function(diff,totalOTUFrame,totalGroupFrame,outdir){
  # Extract sample information of diff group
  diffGroupVector <- str_split(diff, ',')[[1]]
  groupFrame <- normalizeSampleFrame(totalGroupFrame, diffGroupVector)
  groupFrame <- groupFrame %>% select(names, group) %>% filter(group %in% diffGroupVector)
  firstColumnTitle <- colnames(totalOTUFrame)[1]
  controlGroup <- diffGroupVector[1]
  caseGroup <- diffGroupVector[2]
  # Extract data of needed samples
  totalOTUFrame <- column_to_rownames(totalOTUFrame, var = firstColumnTitle)
  dataFrame <- totalOTUFrame[, groupFrame$names, drop = FALSE]
  dataFrame  = round(dataFrame)
  dataFrame = dataFrame[rowSums(cpm(dataFrame) > 1) >= 2,]
  dataFrame = dataFrame[,c(groupFrame[groupFrame$group==controlGroup,]$names,groupFrame[groupFrame$group==caseGroup,]$names)]
  conditions = data.frame(conditions=factor(c(rep(controlGroup, dim(groupFrame[groupFrame$group==controlGroup,])[1]), rep(caseGroup, dim(groupFrame[groupFrame$group==caseGroup,])[1]))))
  rownames(conditions) = colnames(dataFrame)
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = dataFrame, colData = conditions, design = ~ conditions)
  dds = DESeq(ddsFullCountTable)
  contrast=c("conditions",controlGroup,caseGroup)
  res = results(dds, contrast)
  baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == controlGroup])
  baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == caseGroup])
  res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
  res = cbind(sampleA=controlGroup, sampleB=caseGroup, as.data.frame(res))
  res$padj[is.na(res$padj)]  <- 1
  res$pvalue[is.na(res$pvalue)]  <- 1
  res = as.data.frame(res[order(res$pvalue),])
  res2 = res[,c(5,3,4,6,9,10)]
  res2$log2FoldChange=res2$log2FoldChange*-1
  colnames(res2)=c("baseMean",controlGroup,caseGroup,"log2FoldChange","pvalue","padj")
  res2$regulated=ifelse((as.numeric(res2$pvalue) < 0.05 & as.numeric(res2$log2FoldChange) > 0),"up",
                        ifelse((as.numeric(res2$pvalue) < 0.05 & as.numeric(res2$log2FoldChange) < 0),"down","nosig"))
  resdiff=res2[res2$regulated!="nosig",]
  write.table(rbind(OTU=colnames(res2),res2), file=paste(outdir, paste0(controlGroup,"_vs_",caseGroup,".all.xls"),sep="\\"), quote = F, sep = '\t', row.names=T, col.names = F)
  write.table(rbind(OTU=colnames(resdiff),resdiff), file=paste(outdir, paste0(controlGroup,"_vs_",caseGroup,".diff.xls"),sep="\\"), quote = F, sep = '\t', row.names=T, col.names = F)
}

str_comb <- function(name,x,upset_files){
  wb <- xlsx::createWorkbook()
  n <- length(name)
  num=0
  col=1
  for (i in 1:n) {
    num=num+choose(n,i)
  }
  com_data <- unique(Data$ID)
  id <- c("id")
  for (j in 1:n) {
    comb_res <- combn(name,j)
    m <- ncol(comb_res)
    for (l in 1:m) {
      outfile <- paste(comb_res[,l],collapse = '-')
      res <- Reduce(intersect,x[c(comb_res[,l])])
      id <- c(id,outfile)
      com_data <- cbind(com_data, all_id %in% res)
      col=col+1
      if(col==num)break
    }
    break
  }
  colnames(com_data) <- id
  sheet <- createSheet(wb,sheetName = "intersection")
  addDataFrame(com_data,sheet,row.names = F,col.names = T)
  saveWorkbook(wb,file = upset_files)
}

ternaryPlot <- function(ternaryFrame, prefix, cols){
  p <- ggtern(data = ternaryFrame, 
              aes_string(x = colnames(ternaryFrame)[2], y = colnames(ternaryFrame)[3], z = colnames(ternaryFrame)[4])) +
    geom_point(aes(size = Abundance, color = group), alpha = 0.8) +
    labs(x=paste(colnames(ternaryFrame)[2],' (Group1)',sep = "   \n"),y=paste(colnames(ternaryFrame)[3],' (Group2)',sep = "\n"),z=paste(colnames(ternaryFrame)[4],'(Group3)',sep = "\n"))+
    scale_colour_manual(values = cols) +
    # theme_rgbw(base_size = 12 ) +
    theme_bw(base_size = 12 ) +
    theme(panel.background = element_rect(fill='white', colour='black'))+
    # labs(title = "Ternary plot") +
    labs(title = "") +
    theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'))+
    theme(axis.line = element_line(colour = "black"))+
    theme(legend.text=element_text(size=9,face = "bold"))+
    theme(axis.title = element_text(color='black',size=10,face = "bold"))
  pdfFigure <- paste0(prefix, ".ternary_plot.pdf")
  ggsave(pdfFigure, p, height = 7, width = 7)
  pngFigure <- str_replace(pdfFigure, ".pdf$", ".png")
  ggsave(pngFigure, p, height = 8, width = 8)
}

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5","#F0E442","#CCCCCC")#all

# ---------------------------------------------------------------------------------------------------------------------------------------------
# Difference analysis -------------------------------------------------------------------------------------------------------------------------
# Load original OTU and group（N0B0=A;N0B1=B;N0B2=C）
OTUFrame <- loadTable(input)
GroupFrame <- loadSampleGroup(group)

dif_num <- str_split(Diff, ';')[[1]]
for(diff in 1:length(dif_num)){
  diffgroup=gsub('_vs_',',',dif_num[diff])
  diffAnalysis(diffgroup,OTUFrame,GroupFrame,outdir)
}


# Data statistics------------------------------------------------------------------------------------------------------------------------------
AB <- read.table(paste(outdir, paste0(dif_num[1],".all.xls"),sep="\\"),header = TRUE,row.names = 1,sep = "\t")
AC <- read.table(paste(outdir, paste0(dif_num[2],".all.xls"),sep="\\"),header = TRUE,row.names = 1,sep = "\t")
BC <- read.table(paste(outdir, paste0(dif_num[3],".all.xls"),sep="\\"),header = TRUE,row.names = 1,sep = "\t")

AB_AE=AB[AB$regulated=="down",]
AB_BE=AB[AB$regulated=="up",]
AC_AE=AC[AC$regulated=="down",]
AC_CE=AC[AC$regulated=="up",]
BC_BE=BC[BC$regulated=="down",]
BC_CE=BC[BC$regulated=="up",]

Group1_Enriched=intersect(row.names(AB_AE),row.names(AC_AE)) %>% unique()
Group2_Enriched=intersect(row.names(AB_BE),row.names(BC_BE)) %>% unique()
Group3_Enriched=intersect(row.names(AC_CE),row.names(BC_CE)) %>% unique()
Group1_Depleted=intersect(row.names(AB_BE),row.names(AC_CE)) %>% unique()
Group2_Depleted=intersect(row.names(AB_AE),row.names(BC_CE)) %>% unique()
Group3_Depleted=intersect(row.names(AC_AE),row.names(BC_BE)) %>% unique()
ABAE=setdiff(row.names(AB_AE),row.names(AC_AE))%>% unique()#1VS2_1E
ABBE=setdiff(row.names(AB_BE),row.names(BC_BE))%>% unique()#1VS2_2E
ACAE=setdiff(row.names(AC_AE),row.names(AB_AE))%>% unique()#1VS3_1E
ACCE=setdiff(row.names(AC_CE),row.names(BC_CE))%>% unique()#1VS3_3E
BCBE=setdiff(row.names(BC_BE),row.names(AB_BE))%>% unique()#2VS3_2E
BCCE=setdiff(row.names(BC_CE),row.names(AC_CE))%>% unique()#2VS3_3E
alldiff <- Reduce(union, list(Group1_Enriched,Group2_Enriched,Group3_Enriched,Group1_Depleted,Group2_Depleted,Group3_Depleted,ABAE,ABBE,ACAE,ACCE,BCBE,BCCE))%>% unique()
allnosig <- Reduce(union, list(row.names(AB[AB$regulated=="nosig",]), row.names(AC[AC$regulated=="nosig",]), row.names(BC[BC$regulated=="nosig",])))%>% unique()
nosig= setdiff(allnosig,alldiff)%>% unique()#nosig

Data=rbind(data.frame("ID"=Group1_Enriched,Set_name=rep("Group1_Enriched",length(Group1_Enriched))),
           data.frame("ID"=Group2_Enriched,Set_name=rep("Group2_Enriched",length(Group2_Enriched))),
           data.frame("ID"=Group3_Enriched,Set_name=rep("Group3_Enriched",length(Group3_Enriched))),
           data.frame("ID"=Group1_Depleted,Set_name=rep("Group1_Depleted",length(Group1_Depleted))),
           data.frame("ID"=Group2_Depleted,Set_name=rep("Group2_Depleted",length(Group2_Depleted))),
           data.frame("ID"=Group3_Depleted,Set_name=rep("Group3_Depleted",length(Group3_Depleted))),
           data.frame("ID"=ABAE,Set_name=rep("ABAE",length(ABAE))),
           data.frame("ID"=ABBE,Set_name=rep("ABBE",length(ABBE))),
           data.frame("ID"=ACAE,Set_name=rep("ACAE",length(ACAE))),
           data.frame("ID"=ACCE,Set_name=rep("ACCE",length(ACCE))),
           data.frame("ID"=BCBE,Set_name=rep("BCBE",length(BCBE))),
           data.frame("ID"=BCCE,Set_name=rep("BCCE",length(BCCE))))

# upsetR ----------------------------------------------------------------------------------------------------------------------------------------------------------------
sets <- unique(Data$Set_name)
listInput <- list() 
for (i in 1:length(sets)) {listInput[[i]] <- Data$ID[which(Data$Set_name == sets[i])]}
names(listInput) <- sets
upsetR_pdf=paste(outdir, "upsetR.pdf", sep = "\\")
upsetR_png=paste(outdir, "upsetR.png", sep = "\\")
upset_files=paste(outdir, "upset_intersection.xlsx",sep="\\")

# upsetR  Statistics information
all_id <- unique(Data$ID)
str_comb(names(listInput), listInput, upset_files)

#upsetR draw
pdf(upsetR_pdf, height = 6, width =  8, onefile = F) 
upset(fromList(listInput), order.by = "freq", nsets = length(sets),
      mainbar.y.label = 'Intersection Size', sets.x.label = "Set Size", matrix.color = "grey23",
      main.bar.color = "grey23", sets.bar.color = "grey23")
dev.off()
png(upsetR_png, height =  6, width =  8, res = 300, units = "in") 
upset(fromList(listInput), order.by = "freq", nsets = length(sets),
      mainbar.y.label = 'Intersection Size', sets.x.label = "Set Size", matrix.color = "grey23",
      main.bar.color = "grey23", sets.bar.color = "grey23")
dev.off()

# ead upsetR data,statistics on the data,对各分组信息进行统计，及重复名字问题长名短名一起选短名， 短名短名一起，选择短名集合数量少的留下，数量多的重复短名删除----------
dataFrame <- read_xlsx(upset_files)
dataFrame <- column_to_rownames(dataFrame, var = "id")

onlyone=NULL #获得仅在一类中存在的微生物
more=NULL  #获得多类中存在的微生物，方便后面的重新分类
for(i in 1:dim(dataFrame)[1]){
  # i=1
  File_i<- dataFrame[i,]
  if(as.numeric(t(table(File_i[1,]==TRUE))[,"TRUE"])==1){
    onlyone=rbind(onlyone,File_i) 
  }else{
    more=rbind(more,File_i) 
  }
}
# 统计每类中仅在一类中存在的微生物的个数，方便后面微生物重新归类是对判断
onlyone_num=data.frame(Group1_Enriched=ifelse((length(t(table(onlyone$Group1_Enriched)))==2),as.numeric(t(table(onlyone$Group1_Enriched==TRUE))[,"TRUE"]),0),
                       Group2_Enriched=ifelse((length(t(table(onlyone$Group2_Enriched)))==2),as.numeric(t(table(onlyone$Group2_Enriched==TRUE))[,"TRUE"]),0),
                       Group3_Enriched=ifelse((length(t(table(onlyone$Group3_Enriched)))==2),as.numeric(t(table(onlyone$Group3_Enriched==TRUE))[,"TRUE"]),0),
                       Group1_Depleted=ifelse((length(t(table(onlyone$Group1_Depleted)))==2),as.numeric(t(table(onlyone$Group1_Depleted==TRUE))[,"TRUE"]),0),
                       Group2_Depleted=ifelse((length(t(table(onlyone$Group2_Depleted)))==2),as.numeric(t(table(onlyone$Group2_Depleted==TRUE))[,"TRUE"]),0),
                       Group3_Depleted=ifelse((length(t(table(onlyone$Group3_Depleted)))==2),as.numeric(t(table(onlyone$Group3_Depleted==TRUE))[,"TRUE"]),0),
                       ABAE=ifelse((length(t(table(onlyone$ABAE)))==2),as.numeric(t(table(onlyone$ABAE==TRUE))[,"TRUE"]),0),
                       ABBE=ifelse((length(t(table(onlyone$ABBE)))==2),as.numeric(t(table(onlyone$ABBE==TRUE))[,"TRUE"]),0),
                       ACAE=ifelse((length(t(table(onlyone$ACAE)))==2),as.numeric(t(table(onlyone$ACAE==TRUE))[,"TRUE"]),0),
                       ACCE=ifelse((length(t(table(onlyone$ACCE)))==2),as.numeric(t(table(onlyone$ACCE==TRUE))[,"TRUE"]),0),
                       BCBE=ifelse((length(t(table(onlyone$BCBE)))==2),as.numeric(t(table(onlyone$BCBE==TRUE))[,"TRUE"]),0),
                       BCCE=ifelse((length(t(table(onlyone$BCCE)))==2),as.numeric(t(table(onlyone$BCCE==TRUE))[,"TRUE"]),0)) 

## 多类中存在的微生物重新归类
MORE_DATA=NULL
head(more)
for(a in 1:dim(more)[1]){
  short_i<- more[a,c(1:6)]
  long_i<- more[a,c(7:12)]
  if(grepl("TRUE",paste(colnames(t(table(short_i==TRUE))), collapse = '_'))){
    
    test=data.frame("ID"=row.names(t(short_i[1,]==TRUE)),t(short_i[1,]==TRUE))
    colnames(test)=c("ID","GEGE")
    
    if(dim(test[test$GEGE==TRUE,])[1]>1){
      
      min=colnames(short_i[,short_i[1,]==TRUE])[1]
      for( b in 2:(length(colnames(short_i[,short_i[1,]==TRUE]))) ){
        #b = "Group2_Depleted"
        b_tmp=colnames(short_i[,short_i[1,]==TRUE])[b]
        if(as.numeric(onlyone_num[,min]) > as.numeric(onlyone_num[,b_tmp])){
          min=b_tmp
        }
      }
      MORE_DATA_i=data.frame("ID"=row.names(short_i),group=min)
      MORE_DATA=rbind(MORE_DATA,MORE_DATA_i)
    }else{
      MORE_DATA_i=data.frame("ID"=row.names(short_i),group=test[test$GEGE==TRUE,]$ID)
      MORE_DATA=rbind(MORE_DATA,MORE_DATA_i)
    }
  }else{
    
    test=data.frame("ID"=row.names(t(long_i[1,]==TRUE)),t(long_i[1,]==TRUE))
    colnames(test)=c("ID","GEGE")
    
    if(dim(test[test$GEGE==TRUE,])[1]>1){      
      min=colnames(long_i[,long_i[1,]==TRUE])[1]
      for( b in 2:(length(colnames(long_i[,long_i[1,]==TRUE]))) ){
        #b = "Group2_Depleted"
        b_tmp=colnames(long_i[,long_i[1,]==TRUE])[b]
        if(as.numeric(onlyone_num[,min]) > as.numeric(onlyone_num[,b_tmp])){
          min=b_tmp
        }
      }
      MORE_DATA=data.frame("ID"=row.names(long_i),group=min) 
      MORE_DATA=rbind(MORE_DATA,MORE_DATA_i)
    }else{
      MORE_DATA_i=data.frame("ID"=row.names(long_i),group=test[test$GEGE==TRUE,]$ID)
      MORE_DATA=rbind(MORE_DATA,MORE_DATA_i)
    }
  }
}

test2=data.frame("ID"=row.names(onlyone),onlyone)
one_DATA=NULL
for(l in 2:dim(test2)[2]){
  # l= 3
  test_tmp=test2[,c(1,l)]
  one_DATA_i=data.frame("ID"=unique(test_tmp[test_tmp[,2]=="TRUE",]$ID),group=rep(colnames(test_tmp)[2],length(unique(test_tmp[test_tmp[,2]=="TRUE",]$ID))))
  one_DATA=rbind(one_DATA,one_DATA_i)
}

nosig_DATA=data.frame("ID"=nosig,group=rep("nosig",length(nosig)))
plot_data=rbind(one_DATA,MORE_DATA,nosig_DATA)

row.names(plot_data)=plot_data$ID
id=intersect(row.names(AB),intersect(row.names(AC),row.names(BC)))%>%unique()
plot_data2 <- cbind(AB[id,],AB[id,],AC[id,],BC[id,],plot_data[id,])
abundance <- plot_data2[c(str_split(ternary, ',')[[1]],"group")]
abundance <- abundance[order(abundance$group),]
# 对微生物在组中的丰度均值以及分类数据进行保存
write.table(rbind(OTU=colnames(abundance),abundance), file=paste(outdir, "ternaryPlot.xls",sep="\\"), quote = F, sep = '\t', row.names=T, col.names = F)

# 绘制三元相图
abundance <- data.frame("Taxonomy"=row.names(abundance),abundance)
abundance$Abundance <- apply(abundance[,2:4],1,mean)
abundance <- abundance[order(abundance[,"Abundance"],decreasing = T),]
#修改分类名
abundance$group=ifelse(abundance$group== "ABAE","1VS2_1E(2D)",ifelse(abundance$group=="ABBE",'1VS2_2E(1D)',
                                                                     ifelse(abundance$group=="ACAE","1VS3_1E(3D)",ifelse(abundance$group=="ACCE","1VS3_3E(1D)",
                                                                                                                         ifelse(abundance$group=="BCBE","2VS3_2E(3D)",ifelse(abundance$group=="BCCE","2VS3_3E(2D)",abundance$group))))))
level=c("Group1_Enriched","Group2_Enriched","Group3_Enriched","Group1_Depleted","Group2_Depleted","Group3_Depleted","1VS2_1E(2D)","1VS2_2E(1D)","1VS3_1E(3D)","1VS3_3E(1D)","2VS3_2E(3D)","2VS3_3E(2D)","nosig")
abundance <- transform(abundance,group=factor(group,levels =level)) 
# Ternary plot
prefix <- (paste0(outdir, '\\',gsub(',','_vs_',ternary)))
ternaryPlot(abundance, prefix, Palette)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

