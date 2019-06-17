## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, cache = FALSE)


## ----warning=FALSE, message=FALSE----------------------------------------
#es necessita tenir instalat (https://github.com/rich-iannone/DiagrammeR/issues/133)
#install.packages("webshot")
#webshot::install_phantomjs()


library(DiagrammeR)
m <- mermaid("
        gantt
        dateFormat  YYYY-MM-DD
        title Diagrama de Gantt
        
        section PEC 1
        Plan de trabajo             :active,          first_1,    2019-03-05, 2019-03-17
        Entrega PEC1    :crit,    first_36, 2019-03-18, 24h
        
        section PEC 2
        Revisión Bibliográfica métodos integración       :active,        first_2,    2019-03-19, 2019-03-31
        Revisión Bibliográfica métodos preprocesado      :active,    import_1,   2019-03-25, 2019-03-30
        Elección del conjunto de datos :active,  import_3,   after import_1, 7d
        Preparación de los datos: active, import_5, after import_3, 9d
        Entrega PEC 2            :crit,          import_4,   2019-04-23, 2019-04-24
        
        section PEC 3
        Prueba de diferentes métodos de preprocesado  :active,        extras_1,   after import_4,  2d
        Prueba de diferentes métodos de integración   :  active             extras_2,   after extras_1, 12d
        Creación del modelo predictivo            : active              extras_3,  2019-05-11, 2019-05-15
        Evaluación del modelo predictivo            : active              extras_4,   2019-05-16, 1d
        Valoración de las variables en el modelo    : active              extras_5,   2019-05-18, 1d
        Prueba de creación del modelo con otros datos         : active              extras_6,   2019-05-19, 1d
        Entrega de la PEC 3: crit,  extras_7, 2019-05-20, 24h

        section PEC 4
        Redacción de la memoria                 :active,        extras_10,   2019-04-01,  2019-06-05
        Entrega de la PEC 4               : crit,              extras_12,   after extras_10, 24h

        section PEC 5
        Preparación de la presentación                 :active,        extras_15,   after extras_12,  7d
        Entrega de la PEC 5a                :crit,              extras_22,   after extras_15, 24h
        Defensa del TFM            :    crit,           extras_32,   after extras_22, 9d
        ")

m$x$config = list(ganttConfig = list(
  axisFormatter = list(list(
    "%b %d" 
    ,htmlwidgets::JS(
      'function(d){ return d.getDay() == 1 }' 
    )
  ))
))

m


## ----  echo=TRUE---------------------------------------------------------
mainDir <-getwd()
workingDir <- mainDir
celDir <- file.path(workingDir, "celfiles")


## ---- echo=TRUE, message=FALSE-------------------------------------------
library(xtable)
library(Biobase)
library(oligo)
library(ggplot2)
library(ggrepel)
source("https://raw.githubusercontent.com/uebvhir/UEB_PCA/master/UEB_plotPCA3.R")
library(genefilter)
library(clariomdhumantranscriptcluster.db)
library(limma)


## ------------------------------------------------------------------------
targets <-read.csv2(file.path(workingDir,"targets.csv"), header = TRUE, 
                    sep = ";", row.names = 1) 
 

 
pd <-read.AnnotatedDataFrame(file.path(workingDir,"targets.csv"),
                             header = TRUE, row.name = "FileName", sep=";")


## ---- results='asis'-----------------------------------------------------
x.big<-xtable(targets, caption="Archivo targets que muestra la asociación entre las muestras y sus covariables",
               label="targettable1", include.rownames = TRUE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## celFiles <- list.celfiles(celDir, full.names = TRUE)
## rawData <- read.celfiles(celFiles, phenoData = pd)
## colnames(rawData) <-rownames(pData(rawData)) <- pd@data$ShortName


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## colores <-  pd@data$colores
## grupos <- as.factor(pd@data$Group)
## colorPCA <- c("green", "yellow")
## 
## boxplot(rawData, which="all", cex.axis=0.6, col = colores,  las = 2,
##           main="Boxplot for arrays intensity: Raw Data")
## 
## plotPCA3(exprs(rawData), labels = colnames(rawData), factor = grupos,
##          title="Raw data", scale = FALSE, colores = colorPCA, size = 2.5)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## eset_rma <- rma(rawData)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## boxplot(eset_rma,main="Boxplot of Normalized data", cex.axis=0.5, col=colores, las=2)
## plotPCA3(exprs(eset_rma), labels = colnames(eset_rma), factor = grupos,
##          title="Normalized data", scale = FALSE, colores = colorPCA, size = 2.5)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## annotation(eset_rma) <- "clariomdhumantranscriptcluster.db"
## filtered <- nsFilter(eset_rma, require.entrez = TRUE,
##          var.func=IQR, remove.dupEntrez = TRUE, require.GOBP = FALSE,
##          require.GOCC = FALSE, require.GOMF = FALSE,
##          var.filter = TRUE, var.cutoff = 0.66,
##          filterByQuantile = TRUE, feature.exclude = "^AFFX")
## dim(filtered$eset) #8239 60
## eset_filtered <-filtered$eset


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## grupo <- as.factor(targets$Group)
## design2 <- model.matrix( ~ 0 + grupo)
## colnames(design2)<-c("CASE","CONTROL")
## rownames(design2)<-targets$ShortName
## print(design2)
## dim(design2) #60 2
## 
## contrastsMatrix2 <- makeContrasts(CASEvsCONTROL = CASE - CONTROL,
##                                  levels = design2)
## fit2<-lmFit(eset_filtered, design2)
## fit.main2<-contrasts.fit(fit2, contrastsMatrix2)
## fit.main2<-eBayes(fit.main2)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## topTab_CASEvsCONTROL <- topTable (fit.main2, number = nrow(fit.main2),
##                                   coef="CASEvsCONTROL", adjust="fdr")


## ---- echo=FALSE---------------------------------------------------------
top <-read.csv2(file.path(workingDir,"toptab.csv"), header = TRUE, 
                    sep = ";", row.names = 1) 


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(top[1:50, 1:8], caption = "50 primeros genes más diferencialmente expresados",
               label = "topt", include.rownames = TRUE,comment=FALSE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny",comment=FALSE)


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## #se lee archivo de anotaciones que se ha bajado de la web de affymetrix
## anotacions <- read.csv(file.path(dataDir, "miRNA-4_0-st-v1.annotations.mod.csv"),
##                        sep=",",header=TRUE)
## 
## #se redefinen el nombre de las columnas
## library(dplyr)
## anotacions2 <- add_rownames(anotacions,"Probe.Set.ID2")
## colnames(anotacions2) <- c("Probe.Set.ID","Probe.Set.Name","Accession",
##                            "Transcript.ID.Array.Design.","Sequence.Type",
##                            "Species.Scientific.Name","Alignments","Sequence.Length",
##                            "Sequence","Genome.Context","Clustered.miRNAs.within.10kb",
##                            "Target.Genes","GeneChip.Array","Annotation.Date","Year",
##                            "Sequence.Source")
## 
## #se seleccionan las de la especie humana
## Hanotacions <- anotacions2[which(anotacions2$Species.Scientific.Name == "Homo sapiens"),]


## ---- eval=FALSE,echo=TRUE-----------------------------------------------
## grupo <- as.factor(targets1$Group)
## design <- model.matrix( ~ 0 + grupo)
## colnames(design)<-c("CASE","CONTROL")
## rownames(design)<-targets1$ShortName
## 
## contrastsMatrix <- makeContrasts(CASEvsCONTROL = CASE - CONTROL,
##                                  levels = design)
## 
## fit <- lmFit(end.data1, design)
## fit.main <- contrasts.fit(fit, contrastsMatrix)
## fit.main <- eBayes(fit.main)
## 
## topTab_CASEvsCONTROL <- topTable (fit.main, number = nrow(fit.main),
##                                   coef="CASEvsCONTROL", adjust="fdr")


## ---- echo=FALSE---------------------------------------------------------
top2 <-read.csv2(file.path(workingDir,"toptab2.csv"), header = TRUE, 
                    sep = ";") 


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(top2[1:50, 1:7], caption = "50 primeros miRNA más diferencialmente expresados",
               label = "topt2")
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE, include.rownames = FALSE)


## ---- eval=FALSE, echo=TRUE----------------------------------------------
## dades <- read.table(file.path("dades.csv"), sep = ";", head = T,row.names = 1)
## library(compareGroups)
## #githubinstall::githubinstall("mmotaF")
## library(mmotaF)
## #githubinstall::githubinstall("anaStatsUEB")
## library(anaStatsUEB)


## ---- eval=FALSE, echo=TRUE----------------------------------------------
## res <- compareGroups(~. , data = dades[,-c(1,2,9:11,13,36,56,83,84)], max.xlev = 50,
##                      method = 2)
## restab <- createTable(res)
## export2csv(restab, "descriptive.csv", which.table="descr", sep=";", nmax = TRUE)
## descriptive <- read.table(file.path("descriptive.csv"), sep = ";", head = T)


## ---- results='asis'-----------------------------------------------------
descriptive <- read.table(file.path("descriptive.csv"), sep = ";", head = T)
x.big <- xtable(descriptive[1:50,], caption = "Análisis descriptivo de las primeras 50 variables",
               label = "desc", include.rownames = FALSE, comment=FALSE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- eval=FALSE, echo=TRUE----------------------------------------------
## pdf("descriptivePlot.pdf")
## desc_plot(dades[,-c(1,2,9:11,13,36,56,83,84)], rowcol = c(1,3),show.lg = TRUE,
##           cex.lab = 0.01)
## dev.off()


## ---- eval=FALSE, echo=TRUE----------------------------------------------
## res <- compareGroups(type~. , data = dades[,-c(1,2,9:11,13,36,56,83,84)],
##                      max.xlev = 10, method = 2 )
## restab1 <- createTable(res)
## export2csv(restab1, "comp.csv", sep=";", nmax = TRUE)


## ---- results='asis'-----------------------------------------------------
comp <- read.table(file.path("comp.csv"), sep = ";", head = T)
x.big <- xtable(comp[1:50,], caption = "Análisis estadístico de las primeras 50 variables",
               label = "desc", include.rownames = FALSE, comment=FALSE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- results='asis'-----------------------------------------------------
descrip <- read.table(file.path("celulas.csv"), sep = ";", head = T)
x.big <- xtable(descrip, caption = "Análisis descriptivo de las poblaciones celulares",
               label = "desc", include.rownames = FALSE, comment=FALSE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- results='asis'-----------------------------------------------------
compcel <- read.table(file.path("compcel.csv"), sep = ";", head = T)
x.big <- xtable(compcel, caption = "Análisis estadístico de las poblaciones celulares",
               label = "desc", include.rownames = FALSE, comment=FALSE )
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- eval=TRUE, echo=TRUE, results='hide', message=FALSE----------------
library(mixOmics)


## ---- echo=TRUE----------------------------------------------------------
genes <- read.csv2(file.path("ClariomD_TopTab.csv"), sep = ";", header = TRUE, dec = ".")


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(genes[1:20, 1:8], caption = "toptable del análisis de genes",
               label = "desc")
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE, include.rownames = FALSE)


## ---- echo=TRUE----------------------------------------------------------
mirna <- read.csv2(file.path("mirna_expr.csv"), sep = ";", dec = ",", header = TRUE)


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(mirna[1:20, 1:8], caption = "toptable del análisis de miRNA",
               label = "desc")
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE, include.rownames = FALSE)


## ---- echo=TRUE----------------------------------------------------------
clinical <- read.csv2(file.path("ClinicalData_1_clean.csv"), header = TRUE, 
                      dec = ",", sep = ";")

clinical <- clinical[order(clinical$Row.names),] 
rownames(clinical) <- clinical[, 1]
clinical <- clinical[, -c(1, 11)]

#se eliminan las filas que no están en las ómicas
clinical.sel <- clinical[-c(2, 3, 5),]
#se eliminan aquellas columnas no numéricas
clinical.sel <- clinical.sel[, -c(1, 5:9)]


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(clinical.sel[1:15, 1:10], caption = "Variables clínicas seleccionadas",
               label = "desc", include.rownames = FALSE)
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ------------------------------------------------------------------------
colnames(clinical.sel)


## ------------------------------------------------------------------------
cell <- read.csv2(file.path("PoblacionesCelulares2.csv"), sep = ";", dec = ",", header = TRUE)
cell <- cell[order(cell$X),]
rownames(cell) <- cell[, 1]

cell <- cell [, -c(1, 2)]
cell.sel <-  cell[-c(2, 3, 5),]


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(cell.sel[1:15, 1:4], caption = "Poblaciones celulares seleccionadas",
               label = "desc", include.rownames = FALSE)
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- echo=TRUE----------------------------------------------------------
grupos <- read.csv2(file.path("grupos.csv"), sep = ";", header = TRUE)

Y <- grupos$Group
summary(Y) 


## ---- results='asis'-----------------------------------------------------
x.big <- xtable(grupos, caption = "Correspondencia entre los pacientes y la condición experimental a la que pertenecen",
               label = "desc", include.rownames = FALSE)
print(x.big, tabular.environment='longtable', floating=FALSE, size = "tiny", comment=FALSE)


## ---- echo=TRUE----------------------------------------------------------
gen.sel <- subset(genes, P.Value < 0.01)
#se quita la muestra P05 que en los miRNA no está
gen.sel <- gen.sel[, c(2, 11, 12, 14:70)] 
rownames(gen.sel) <- gen.sel$Gene.Symbol
gen.sel <- gen.sel[, -1]


## ---- echo=TRUE----------------------------------------------------------
mirna.sel <- subset(mirna, P.Value < 0.01)
mirna.sel <- mirna.sel[, c(2, 11:69)]
rownames(mirna.sel) <- mirna.sel$Probe.Set.Name
mirna.sel <- mirna.sel[, -1]


## ---- echo=TRUE----------------------------------------------------------
colnames(gen.sel) %in% colnames(mirna.sel)
colnames(gen.sel) %in% rownames(clinical.sel)
rownames(clinical.sel) %in% rownames(cell.sel)


## ---- echo=TRUE----------------------------------------------------------
gen.selt <- t(gen.sel)
mirna.selt <- t(mirna.sel)


## ---- echo=TRUE----------------------------------------------------------
set.seed(123)
ind <- sample(nrow(grupos), 0.67*dim(grupos)[1])

train.gen.selt <- gen.selt[ind, ]
test.gen.selt <- gen.selt[-ind, ]
train.mirna.selt <- mirna.selt[ind, ]
test.mirna.selt <- mirna.selt[-ind, ]
train.clinical.sel <- clinical.sel[ind, ]
test.clinical.sel <- clinical.sel[-ind, ]
train.cell.sel <- cell.sel[ind, ]
test.cell.sel <- cell.sel[-ind, ]
train.grupos <- grupos[ind, ]
test.grupos <- grupos[-ind, ]


## ---- echo=TRUE----------------------------------------------------------
rownames(train.gen.selt) %in% rownames(train.mirna.selt)
rownames(train.cell.sel) %in% rownames(train.mirna.selt)
rownames(train.gen.selt) %in% rownames(train.clinical.sel)
train.grupos$X %in% rownames(train.cell.sel)


## ---- echo=TRUE----------------------------------------------------------
X <- list(mRNA = train.gen.selt, 
          miRNA = train.mirna.selt,
          cell = train.cell.sel,
          clinical = train.clinical.sel)

Y <- train.grupos$Group


## ---- echo=TRUE----------------------------------------------------------
list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), cell = c(5, 5), clinical = c(6 ,6))


## ---- echo=TRUE----------------------------------------------------------
MyResult.diablo <- block.splsda(X, Y, keepX = list.keepX, ncomp = 2, 
                                scale = TRUE, mode = "regression")


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
plotIndiv(MyResult.diablo,
            ind.names = TRUE,
            legend = TRUE, cex=c(3, 3),
            title = 'PLOT OF SAMPLES')



## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
plotVar(MyResult.diablo, 
          var.names = c(TRUE, TRUE, TRUE, TRUE),
          pch = c(5, 5, 5, 5))


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
plotDiablo(MyResult.diablo, ncomp = 1)


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
circosPlot(MyResult.diablo, cutoff = 0.65)


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
cimDiablo(MyResult.diablo,
            color.blocks = c('darkorchid', 'brown1', 'lightgreen', "lightblue"),
            comp = 1,
            margin = c(10,15),
            legend.position = "right")


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
plotLoadings(MyResult.diablo, comp = 1, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")


## ---- echo=TRUE, fig.align='center', fig.height=7, fig.width=8-----------
network(MyResult.diablo, blocks = c(1, 2, 3, 4),
        color.node = c('darkorchid', 'brown1', 'lightgreen', "lightblue"),
        cutoff = 0.95,
        save = 'pdf', 
        name.save = file.path('CorNetwork95'))


## ---- echo=TRUE----------------------------------------------------------
X.test <- list(mRNA = test.gen.selt, 
               miRNA = test.mirna.selt,
               cell = test.cell.sel,
               clinical = test.clinical.sel)
Mypredict.diablo <- predict(MyResult.diablo, newdata = X.test)


## ---- echo=TRUE----------------------------------------------------------
confusion.mat <- get.confusion_matrix(
    truth = test.grupos$Group, 
  predicted = Mypredict.diablo$MajorityVote$centroids.dist[,2])

confusion.mat

