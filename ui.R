
check.package=function(package){
	if (!require(package,character.only = TRUE)){
		install.packages(package)
		library(package,character.only = TRUE)
	}
}

if (!require("DT")){
	install.packages('DT', repos = 'http://cran.rstudio.com')
	library("DT")
}
check.package("gtools")
check.package("limma")
check.package("glmnet")
check.package("parallel")
check.package("survival")
check.package("FactoMineR")
check.package("matrixStats") 

check.package("affy")
check.package("genefilter")
check.package("annotate")
check.package("shinyFiles")


#if (!require("MASS")) install.packages('MASS')

shinyUI(navbarPage("Robust Gene Selection for qPCR Validation", theme="bootstrap.css",
	tabPanel("Main",
		sidebarPanel(

			# READ DATA

			div(
				fileInput('file', h5('Predictors File')),
				style="height:80px"		
			),

			div(
				fileInput('adjust.factors', h5('Outcome + Adjust Factors File')),
				style="height:70px"		
			),

			div(hr(),style="height:5px"),
			#radioButtons("aux.param", label = h5("aux.param"),inline=T,
				#choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "TRUE"),

			selectInput(inputId = "response",label = h5("Type of Response"),
				choices = list("Dichotomous (Logistic Regression)" ="binomial", 
						"Survival (Cox Regression)" = "survival",
						"Continuous (Linear Model)" = "gaussian",
						"Counts (Poisson Regression)"="poisson", 
						"Nominal (Multinomial Regression)"="multinomial"),
				selected = "binomial"),


			actionButton("goButton", "Start Resampling",icon=icon("forward", lib = "font-awesome")),


			div(hr(),style="height:5px"),

			selectInput(inputId = "whichpanel",label = h5("Options"),
				choices = list("Data Options" ="1", "Resampling Options" = "2", "Candidate Options" = "3", "Lasso Options"="4", "Plot Options"="5"),
				selected = "1"),


			# data options: weights, scale, log

			wellPanel(

				conditionalPanel(
					condition = "input.whichpanel == 1",

					h4("Data Options"),
					radioButtons("log", label = h5("Logarithmic Transformation (Base 2)"),inline=T,
						choices = list("YES" ="log2", "NO" = "FALSE"),selected = "FALSE"),

					radioButtons("stan", label = h5("Standardize Features (Z-score)"),inline=T,
						choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "FALSE"),
					#conditionalPanel(
						#condition = "input.stan == TRUE",
						helpText("Note: standardizing is performed at each iteration independently"),
					#),


					fileInput('weights', h5('Sample Weights'))

				),

				conditionalPanel(
					condition = "input.whichpanel == 2",

					h4("Resampling Options"),	
					sliderInput("M",h5("Number of Resamplings"), min = 10,max = 2000,value = 100,step=10),
					sliderInput("OOB.percent",h5("% OOB samples in each Resampling"),min = 0,max = 95,value = 10,step=1),
					sliderInput("var.percent",h5("% of Features in each Resampling"),min = 5,max = 100,value = 30,step=1)

				),

				conditionalPanel(
					condition = "input.whichpanel == 3",

					h4("Candidate Options"),

					# pregunta filosòfica: és correcte que l'usuari ho pugui triar? És complicat de repondre, potser millor deixar l'opció, ja que 
						# en algun cas es podria decidir agafar gens amb un cert FC (sense ajustar) per detectar a la qPCR
					radioButtons(inputId = "filteradjusted",label = h5("Factor Adjusting?"), # seran els resiuds del model lineal
						choices = list("Yes" =TRUE, "No" = FALSE),
						selected = FALSE,inline=TRUE),

					conditionalPanel(
						condition = "input.filteradjusted == 'FALSE'",
						sliderInput("mfc",h5("Minimum Fold Change Required"), min = 1,max = 10,value = 1,step=0.1)
					),

					sliderInput("keep.filter",h5("Number of Lasso Candidate Features"),min = 20,max = 1000, value = 200,step=10),

					conditionalPanel(
						condition = "input.filteradjusted == 'FALSE'",
						radioButtons("type.ranking.no.adj", label = h5("Ranking Method for Lasso Candidate Features"),
							choices = list("Limma's P-Value" ="limma", "Wilcoxon (time consuming)" = "non.parametric", "Fold Change" = "fc"), 
							selected = "limma")
					),

					conditionalPanel(
						condition = "input.filteradjusted == 'TRUE'",
						radioButtons("type.ranking.adj", label = h5("Ranking Method for Lasso Candidate Features"),
							choices = list("Limma's P-Value" ="limma"), 
							selected = "limma")
					)
				),


				conditionalPanel(
					condition = "input.whichpanel == 4",

					h4("Lasso Options"),	
					sliderInput("alpha",h5("Elasticnet Mixing Parameter"), min = 0,max = 1,value = 1,step=0.01), # 0=ridge penalty, 1: lasso

					radioButtons("lambda", label = h5("Cross-Validated Lambda"),inline=FALSE,
						choices = list("1SE Min CV error" ="lambda.1se", "Min CV error" = "lambda.min"),selected = "lambda.1se"),

					numericInput("dfmax", h5("dfmax"),value=NA),
					numericInput("pmax", h5("pmax (No sé el que fa)"),value=NA) # hi ha algun problema, no entenc aquest parametre

				),

				conditionalPanel(
					condition = "input.whichpanel == 5",

					h4("Plot Options"),
					uiOutput("ui"),

					radioButtons("alternative.plot", label = h5("Type Plot"),inline=T,
						choices = list("CV Lasso Errors" ="lasso.errors", 
							"Biomarkers Length" = "length.biomarkers"),
						selected = "lasso.errors"),

					radioButtons("ordre", label = h5("Feature Order"),inline=T,
						choices = list("% Times Selected" ="perctimes", "Uniqueness" = "uniqueness",
							"LOOCV(BA) *"="uni.error","Fold-Change *"="FC"),
						selected = "perctimes"),

					helpText("* using features selected at least one time"),

					selectInput(inputId = "bivariate.combination",label = h5("Bivariate Plot"),
						choices = list("% Selected vs LOOCV(BA)" ="1", "% Selected vs FC" = "2",
								"Uniqueness vs LOOCV(BA)" ="3", "Uniqueness vs FC" = "4"),
						selected = "1")

				)


			),


		width=3),
		column(9,plotOutput("main_plot"))
	),


	tabPanel("Common Biomarkers",
		sidebarPanel(
			sliderInput("n",h5("Biomarker Length"),min = 0,max = 50,value = c(0,0),step=1),
			textInput("llista",h5("Iterations without:"),value=""),
			textInput("forced.features",h5("Search biomarkers with any of:"),value=""),
		width=3),
		column(8,DT::dataTableOutput("biomarkers"))
	),

	tabPanel("Data Visualization",
		sidebarPanel(
			#textInput("geneX", h5("gene axis X"),value=""),
			#textInput("geneY", h5("gene axis Y"),value=""),
			#textInput("genes", h5("Genes"),value=""),

			radioButtons("gener_or_biomarker", label = h5("Represent Feature or Biomarker"),inline=T,
				choices = list("Gene" ="Gene", "Common Biomarker" = "Common Biomarker"),selected = "Gene"),

			uiOutput("ui.gene.selection"),	
			uiOutput("ui.biomarker.selection"),			
			hr(),

			radioButtons("scatter_label", label = h5("Sample ID"),inline=T,
				choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "FALSE"),
			hr(),

			uiOutput("ui.plot"),
			uiOutput("ui.plot2"),

			#conditionalPanel(
			#	condition = "input['type.plot'] == 'boxplot'",

			#	radioButtons(inputId = "hide.points",label = h5("Sample Points"),
			#		choices = list("Hide" ="hide", "Show" = "show"),selected = "show",inline=TRUE)

			#),

			radioButtons(inputId = "graph.adjusted",label = h5("Factor Adjusting?"), # seran els resiuds del model lineal
				choices = list("Yes" =TRUE, "No" = FALSE),
				selected = FALSE,inline=TRUE),
			hr(),
			selectInput("legend_loc", label = h5("Legend Position"),
				choices = list("bottomleft"="bottomleft","bottom"="bottom", "bottomright" ="bottomright","left"="left","center"="center",
				"right"="right","topleft" = "topleft","top"="top","topright" = "topright","none"="none"), selected = "topright"),
			selectInput("sel_pal", label = h5("Palette"),
				choices = list("1"="1","2"="2", "3" ="3","4"="4","5"="5","6"="6"), selected = "1"),
		width=3),
		column(9,plotOutput("gene_plot"))
	),

	tabPanel("DGE",
		sidebarPanel(

			# possible method? SAM, limma, Rank Product, Wilcoxon...

			# test

			selectInput("DGEtest", label = h5("Test"),
				choices = list("limma" ="limma","T-Test (Equal Var)" = "ttest.eq","T-Test (Unequal Var)" = "ttest.uneq", 
				"Wilcox" = "wilcox"), 
				selected = "limma"),
	
			radioButtons(inputId = "limmaadjusted",label = h5("Factor Adjusting?"),
				choices = list("Yes" ="Yes", "No" = "No"),
				selected = "No",inline=TRUE),


			h4("Filtering Options"),

			#sliderInput("num_genes",h5("Number of Genes"),min = 10,max = 1000,value = 100,step=10), # number of genes
			numericInput("min_p_val", h5("Adjusted P-Value"),value=0.05,step=0.01), # p-value minim
			#sliderInput("min_fc",h5("Minimum Limma's |Log Fold Change|"), min = 0,max = 10,value = 0,step=0.2),# limma's fold change minim
			#sliderInput("min_fc",h5("Fold Change"), min = 1,max = 20,value = 1,step=0.2),# limma's fold change minim

			conditionalPanel(
				condition = "input.limmaadjusted == 'No'",	
				sliderInput("min_fc",h5("Fold Change"), min = 1,max = 20,value = 1,step=0.2)

			),

			conditionalPanel(
				condition = "input.limmaadjusted == 'No'",	
				selectInput("sort.by2", label = h5("Sort By"),
					choices = list("name" ="name", "Fold Change" = "FC", "Abs(log2FC(Ratio))" = "log2FC(Ratio)" , 
					"P-Value" = "P.Value", "AveExpr" = "AveExpr"), 
					selected = "P.Value")

			),

			conditionalPanel(
				condition = "input.limmaadjusted == 'Yes'",	
				selectInput("sort.by1", label = h5("Sort By"),
					choices = list("name" ="name","logFC(Diff)" = "logFC(Diff)","Abs(logFC(Diff)]" = "Abs(logFC(Diff))", 
					"P-Value" = "P.Value", "AveExpr" = "AveExpr"), 
					selected = "P.Value")

			),

			#selectInput("sort.by", label = h5("Sort By"),
				#choices = list("name" ="name", "Fold Change" = "FC", "P-Value" = "P.Value", "AveExpr" = "AveExpr"), 
				#selected = "P.Value"),

		width=3),
		column(8,DT::dataTableOutput("limma"))
	),

	tabPanel("Data Preprocessing",
		sidebarPanel(

			# READ DATA

			div(
				downloadButton('downloadData.P', h5('Start Preprocessing'),class="orange"),
				style="height:80px"		
			),

			#helpText("- Data Preprocessing might take some time"),
			#helpText("- Data will be downloaded to the downalod directory??"),

			div(
				fileInput('cel.file', h5('CEL Files'),multiple=TRUE),
				style="height:70px"		
			),

	
			#shinyDirButton('directory', h5('Nice Folder Selection'),'Select Folder with class.txt and CEL files'),

			div(hr(),style="height:5px"),

			# normalization method
			selectInput("norm.method", label = h5("Normalization Method"),
				choices = list(" "=" ","RMA" ="rma", "MAS5.0 (log2)" = "mas5"), 
				selected = " "),

			# probeset or gene at output matrix
			radioButtons("type.array", label = h5("Output Matrix Feature Name"),
				choices = list("Probeset"="Probeset", "Gene" = "Gene"), 
				selected = "Probeset",inline=TRUE),

			# Probesets with gene annotated
			radioButtons("entrez.id", label = h5("Required Entrez Gene ID"),inline=T, 
					choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "TRUE"),

			# exclude probesets QC "AFFX"
			radioButtons("feat.exclude", label = h5("Remove Affymetrix QC Probesets"),inline=T,
					choices = list("YES" ="^AFFX", "NO" = "something.impossible"),selected = "^AFFX"),

			# one probeset for gene
			radioButtons("dup.entrez", label = h5("Statistic-based One Probeset for Gene"),inline=T,
					choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "TRUE"),				

			# select % of probests by Statistic
			radioButtons("performVARfiltering", label = h5("Statistic-based Feature Filtering"),inline=T,
					choices = list("YES" ="TRUE", "NO" = "FALSE"),selected = "TRUE"),
			
			# Filter Criterion Statistic
			conditionalPanel(
				condition = "input.performVARfiltering == 'TRUE' || input['dup.entrez'] == 'TRUE'",
				radioButtons("func.filter", label = h5("Statistic"),
					choices = list("IQR" ="IQR", "Standard Deviation" = "sd"), 
					selected = "IQR",inline=TRUE)
			),

			# % of filter out probesets
			conditionalPanel(
				condition = "input.performVARfiltering == 'TRUE'",
				sliderInput("filter.percent",h5("% of Filtered Features by Statistic"),min = 1,max = 99,value = 25,step=1)
			),


		width=3),
		column(9,verbatimTextOutput("summary"))

	)

))



