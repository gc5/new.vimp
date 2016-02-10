

options(shiny.maxRequestSize=-1,warn=-1) 

source("data/funcions.R")

# funció adicional per fer plot

shinyServer(
	function(input, output,session){

		# Per si vull fer lo de seleccionar folder, crec que facilitaria lo de llegir data (evitant copiar-ho a temp file)
		#volumes <- getVolumes() #c('R Installation'=R.home())
		#shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))

		dataInput <- reactive({

			if (!input$response %in% c("binomial")){
				stop("Sorry! Work in progress... Only dichotomic response is supported")
			}

			if(is.null(input$file$datapath)) stop("Upload Predictors File")
			if(is.null(input$adjust.factors$datapath)) stop("Upload Outcome File")

			# data = predictors
			em=read.table(input$file$datapath,h=T,sep="\t",na.strings=c("","NA"),fill=T,row.names=1,
				check.names=FALSE) # repassar
			em=na.omit(em) # que fer amb els missings?
			
			# outcome + adjust
			
			adjust.factors=read.table(input$adjust.factors$datapath,h=T,sep="\t",fill=T,row.names=1)

			if (!identical(sort(rownames(adjust.factors)),sort(colnames(em)))){
				message=paste("\n\nMissing samples in Outcome File:\n",
					paste(colnames(em)[!colnames(em) %in% rownames(adjust.factors)],collapse="\n"),
					"\n\nMissing samples in Data File:\n",
					paste(rownames(adjust.factors)[!rownames(adjust.factors) %in% colnames(em)],collapse="\n"),sep="")
				stop(message)
			}

			adjust.factors=adjust.factors[colnames(em),,drop=F]

			# fer a algun lloc make.names

			if (input$response == "binomial"){
				g=factor(make.names(adjust.factors[,1]))

				if (nlevels(g) != 2){
					stop("First column in Outcome File have more than two groups")
				}

				if (min(table(g))<4){stop("Required a minimum of 4 observations at each class")}

				adjust.factors=adjust.factors[,-1,drop=F]

			}

			if (input$response == "survival"){

				if (ncol(adjust.factors)<=1) stop("Outcome File should have at least two columns (time, status)")

				time=adjust.factors[,1]
				status=adjust.factors[,2]

				if (!class(time) %in% c("numeric","integer")){
					stop("First column in Outcome File (time) is not numeric")
				}

				if (!all(status %in% c(0,1))){
					stop("Second column in Outcome File (status) can contain only {0,1}")
				}

				if (any(time<=0)){
					stop("Time to event must be positive (>0)")
				}

				g=Surv(time,status)

				adjust.factors=adjust.factors[,-c(1,2),drop=F]

			}

			if (input$log=="log2"){
				em=log(em,2)
			}

			#rownames(em)=make.names(rownames(em),unique=T)
			rownames(em)=make.unique(rownames(em),sep=" &")

			if(ncol(adjust.factors)==0){
				adjust.factors=data.frame(matrix(ncol=0,nrow=ncol(em)))
				rownames(adjust.factors)=colnames(em)
			} 

			# FI comprovacions de data

			return(list(em=em,g=g,adjust.factors=adjust.factors))

		})
		
		dataPlotar <- reactive({

			withProgress(message="Loading Data",detail="...",value=0.1,{			

				data=dataInput()

			})

			withProgress(message="Performing Resamplings",detail="",value=0,{			

				mc.cores=detectCores()-1
				#mc.cores=1

				if (mc.cores<1) mc.cores=1

				# llegir weights
				if(is.null(input$weights$datapath)){
					weights.vec=rep(1,ncol(data$em))
					names(weights.vec)=colnames(data$em)
				} else {
					weights.vec=t(read.table(input$weights$datapath,h=F,sep="\t",fill=T,row.names=1))[1,]

					if (!identical(sort(names(weights.vec)),sort(colnames(data$em)))){
						message=paste("\n\nMissing samples in Weights File:\n",
							paste(colnames(data$em)[!colnames(data$em) %in% names(weights.vec)],collapse="\n"),
							"\n\nMissing samples in Data File:\n",
							paste(names(weights.vec)[!names(weights.vec) %in% colnames(data$em)],collapse="\n"),sep="")
						stop(message)
					}
				}

				# 
				type.ranking=c(input$type.ranking.no.adj,input$type.ranking.adj)[as.logical(input$filteradjusted)+1]

				# algun check amb els adjust.factors
				if ( as.logical(input$filteradjusted) & (ncol(data$adjust.factors)==0) ){
						stop("Missing Adjusting Factors File")					
				}				

				resultats=bioselect.paralel.lasso(XX=t(data$em),Y=data$g,M=input$M,keep.filter=input$keep.filter,
					method.filter=type.ranking,adjust.filter=input$filteradjusted,OOB.percent=input$OOB.percent/100,var.percent=input$var.percent/100,
					stan=input$stan,weights.vec=weights.vec,adjust.factors=data$adjust.factors,mfc=input$mfc,alpha=input$alpha,
					dfmax=input$dfmax,pmax=input$pmax,lambda=input$lambda,mc=mc.cores) # scale(t(data$em))

			})

			if (ncol(resultats$bios)<=1){stop("All iterations have 0 selected variables. Perform more iterations or assumue your outcome is unpredictable :(")}

			# ficar aqui el calcul del power de forma individual, pels 40 primers gens directament

			withProgress(message="Last calculations",detail="",value=0,{			

				calc=calculs(quant=min(ncol(resultats$bios[,-1]),40),sortida.lasso=resultats,em=data$em,g=data$g,weights.vec=weights.vec,
					adjust.factors=data$adjust.factors)

			})

			calc$em=data$em
			calc$g=data$g
			calc$adjust.factors=data$adjust.factors

			return(calc)
		
		})

		update.data <- eventReactive(input$goButton, {

			calc=dataPlotar()
		})

		output$main_plot <- renderPlot({
	
			value=input$quant
			if (is.null(value)){
				value=min(25,length(update.data()$perc.sel))
			}

			data=update.data()

			withProgress(message="Creating Plot",detail="",value=1,{
				grafic(data,value,alternative.plot=input$alternative.plot,ordre=input$ordre,type.bivariate=input$bivariate.combination)
			})
	
		},res=100,width=1236,height=1236*10/12) #  ,height=500,res=80 # ,width=1248*0.9,height=1132.8*0.9

		output$ui <- renderUI({

			what=length(update.data()$perc.sel)
			sliderInput("quant", "",min = 1,max = min(40,what), value = min(25,what), step = 1)

		})



		biomark.calc <- reactive({

			calc<-update.data()

			coefs=!is.na(calc$coefs) # TRUE=l'ha fet servir, FALSE= no l'ha fet servir
			vars=calc$vars

			#### IMPORTANT A SOLUCIONAR AQUI!!!!

			# si triem una variable, quines iteracions hem de treure?
			# 1. simplement no ensenyar els biomarkers que han sortit, però no sembla molar del tot...
			# 2. treure totes les iteracions on la variable era candidata, així suposes que aquella falla segur...


			# treiem les iteracions que han inclos alguna de les variables de la llista (fet segons 1.)
			#if (input$llista!=""){
			#	llista.orig=do.call(c,strsplit(input$llista,split=","))
			#	llista=make.names(llista.orig)
			#	if (!all(llista %in% colnames(coefs))){
			#		missatge=paste("{'",paste(llista.orig[!(llista %in% colnames(coefs))],collapse="','"),"'} never used",sep="")
			#		stop(missatge)
			#	}
			#	w=apply(subset(coefs,select=llista),1,sum)==0 # els que són 0 és que no han fet servir cap de les variables
			#	coefs=coefs[w,]
			#}

			# fet segons 2.
			if (input$llista!=""){
				llista.orig=do.call(c,strsplit(input$llista,split=","))
				llista=make.names(llista.orig)
				if (!all(llista %in% colnames(vars))){
					missatge=paste("{'",paste(llista.orig[!(llista %in% colnames(vars))],collapse="','"),"'} not in dataset",sep="")
					stop(missatge)
				}
				w=apply(subset(vars,select=llista),1,sum)==0 # els que són 0 és que no han fet servir cap de les variables
				coefs=coefs[w,]
				vars=vars[w,]
			}

			biomarkers=sapply(1:nrow(coefs),compt.biomarkers,x=coefs)
			biomarkers=table(biomarkers)

			# trec biomarkers que han sortit una vegada només (repassar si provoca problemes fer-ho, suposo que no)
				# excepte si tots han sortit només una vegada
			if (!all(biomarkers==1)) biomarkers=biomarkers[biomarkers>1]

			# numero total de vegades que el biomarker en questio podia sortir

			# trec iteracions que han inclós les variables de input$llista?
			biomark.total=sapply(names(biomarkers),compt.biomarker.possible,vars=vars)
			biomark.total[is.na(biomark.total)]=nrow(coefs)

			# Nota: el biomark.total s'ha de veure afectat per les iteracions que es treuen que incloeien alguna variable?
			# Estudiar que passa en aquest cas
			# No afecta ja que el total no es calcula sobre la matriu coefs. S'hauria de fer? Potser si, ja que treure la variable es podria interpretar com
			# treure-la al 100%

			biomarkers=data.frame(biomarker=names(biomarkers),times=as.numeric(biomarkers),
				perc=round(100*as.numeric(biomarkers)/as.numeric(biomark.total),2),length=sapply(strsplit(names(biomarkers),split=","),length))

			biomarkers=biomarkers[order(biomarkers$times,decreasing=T),]

			if (class(calc$g)=="factor") biomarkers=data.frame(biomarkers,T_2=round(sapply(as.character(biomarkers[,1]),calc.T2,g=calc$g,em=calc$em),2),
								check.names=F)

			return(list(biomarkers=biomarkers,vars=vars))
		})

		biomark.tune <- reactive({

			biomarkers=biomark.calc()

			vars<-biomarkers$vars
			biomarkers=biomarkers$biomarkers

			lim.n1=min(input$n)
			lim.n2=max(input$n)

			#if (input$n>0){ biomarkers=biomarkers[biomarkers$length<=input$n,]}
			if ((lim.n1>0)|(lim.n2>0)){ biomarkers=biomarkers[(biomarkers$length>=lim.n1) & (biomarkers$length<=lim.n2),]}

			colnames(biomarkers)=c("Biomarker","# Times","Percentage Possible (%)","Length Biomark","T^2")

			if (input$forced.features!=""){

				forced.orig=do.call(c,strsplit(input$forced.features,split=","))
				forced=make.names(forced.orig)
				if (!all(forced %in% colnames(vars))){
					missatge=paste("{'",paste(forced.orig[!(forced %in% colnames(vars))],collapse="','"),"'} not in dataset",sep="")
					stop(missatge)
				}

				biomarkers=biomarkers[apply(sapply(forced,grepl,x=biomarkers$Biomarker),1,any),]
			}

			return(biomarkers)

		})

		output$biomarkers <- DT::renderDataTable({
		
			if(input$goButton>0){

				biomarkers=biomark.tune()

				biomarkers
			} else {
				stop("Waiting for resamplings")
			}

		},selection="none",rownames=FALSE, options = list(lengthMenu = c(10, 25, 50), pageLength = 25))

		output$biomarkersAAAAAAA <- DT::renderDataTable({

			if(input$goButton>0){

				biomarkers=biomark.calc()

				vars<-biomarkers$vars
				biomarkers=biomarkers$biomarkers

				lim.n1=min(input$n)
				lim.n2=max(input$n)

				#if (input$n>0){ biomarkers=biomarkers[biomarkers$length<=input$n,]}
				if ((lim.n1>0)|(lim.n2>0)){ biomarkers=biomarkers[(biomarkers$length>=lim.n1) & (biomarkers$length<=lim.n2),]}

				colnames(biomarkers)=c("Biomarker","# Times","Percentage Possible (%)","Length Biomark","T^2")

				if (input$forced.features!=""){

					forced.orig=do.call(c,strsplit(input$forced.features,split=","))
					forced=make.names(forced.orig)
					if (!all(forced %in% colnames(vars))){
						missatge=paste("{'",paste(forced.orig[!(forced %in% colnames(vars))],collapse="','"),"'} not in dataset",sep="")
						stop(missatge)
					}

					biomarkers=biomarkers[apply(sapply(forced,grepl,x=biomarkers$Biomarker),1,any),]
				}

				biomarkers
			} else {
				stop("Waiting for resamplings")
			}

		},selection="none",rownames=FALSE, options = list(lengthMenu = c(10, 25, 50), pageLength = 25))


		output$ui.gene.selection <- renderUI({

			if(is.null(input$file$datapath) | is.null(input$adjust.factors$datapath)){
				#div(id='myDiv', class='simpleDiv','Waiting for Predictor and Outcome Files',
				#	style="color:blue;font-size:16px;background-color: white")
				selectizeInput('no.used1', h5("Genes"), choices = "Waiting for Predictor and Outcome Files", multiple = FALSE)

			} else {
				features.list=rownames(dataInput()$em)
				selectizeInput('features', h5("Genes"), choices = features.list, multiple = TRUE)
			}

		})

		output$ui.biomarker.selection <- renderUI({

			if(input$goButton==0){
				#div(id='myDiv', class='simpleDiv','Waiting for resamplings',
				#	style="color:blue;font-size:16px;background-color: white")
				#div(selectizeInput('no.used2', h5("Common Biomarkers"), choices = "Waiting for resamplings", multiple = FALSE),
				#	style="background-color:blue")
				selectInput('no.used2', h5("Common Biomarkers"), choices = "Waiting for resamplings",selected="Waiting for resamplings")
				#bsButton("no.used2", label = "Waiting for resamplings", icon = icon("ban"),disabled=TRUE)

			} else {

				biomark.list=as.character(biomark.tune()$Biomarker)
				biomark.list=biomark.list[!biomark.list %in% ""]
				selectizeInput('common.biomarkers', h5("Common Biomarkers"), choices = biomark.list, multiple = FALSE)
			}

		})

		# gene plot quan hi ha factor adjusting, dos possbilitats:
			# 1- residus del model lineal vs grup (lo que hi ha fet)
			# 2- model lineal amb tot, prediccions vs grup? no ho ensenyaria bé. residus calculats amb els coeficients d'aquí

		# deixo 1, però s'ha de tenir en compte que el efecte dels adjustig factors en Data Visualization i en DGE és diferent!
		# si hi ha informació compartida entre x i adjusting factors, serà força diferent, ja que al gene.plot es treu tota la info de una primer i després es mira lo que
		# queda en l'altra (lm(residus(adjust.factors) ~ grup)). En canvi al DGE es modelitza conjuntament, que podria tenir més sentit

		# es pot fer només pujar la data, ja que no fa el update.calc()
		output$gene_plot <- renderPlot({

			data=dataInput()

			grey="grey90" # grey85, #E0E0E0

			#gene.orig=strsplit(input$genes,split=",")[[1]]
			if (input$gener_or_biomarker=="Gene"){
				if (length(input$features)==0){stop("")}
				gene.orig=input$features
				
			} else {
				if (length(input$common.biomarkers)==0){stop("")}
				gene.orig=strsplit(input$common.biomarkers,split=",")[[1]]
			}

	

			#if ((input$geneX=="")&(input$geneY=="")){
			#	stop("")
			#}

			#if (length(gene.orig)==0){
			#	stop("")
			#}

			if ( as.logical(input$graph.adjusted) & (ncol(data$adjust.factors)==0) ){
				stop("Missing Adjusting Factors File")
			}


			#if (input$geneX=="" | input$geneY==""){

			if ((length(gene.orig)==1)&(class(data$g)=="factor")){

				# univariate plot

				#gene.orig=c(input$geneX,input$geneY)
				#gene.orig=gene.orig[!gene.orig %in% ""]
				#gene=make.names(gene.orig)
				gene=gene.orig

				if (!(gene %in% rownames(data$em))){
					stop(paste("'",gene.orig,"'", "not in dataset",sep=""))
				}

				x=as.numeric(data$em[gene,])
				if (input$graph.adjusted){
					x=resid(lm(x~.,data=data.frame(x,data$adjust.factors)))
				}

				if (input$type.plot=="Index"){
					par(mar=c(6,6,3,2))
					plot(x,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
					rect(-1000,-1000,1000,1000,col=grey)
			
					xat=seq(1,length(x),by=floor(length(x)/7))
					yat=seq(min(x),max(x),length.out=7)
					abline(h=yat,col="white",lwd=1.5)

					clrs=data$g
					if (input$sel_pal=="1") palette=c("#56A8A4","#EBA94D")
					if (input$sel_pal=="2") palette=c("#63074E","#BFC421")
					if (input$sel_pal=="3") palette=c("#97A861","#9D6187") # c("#47BBD1","#840505")
					if (input$sel_pal=="4") palette=c("#60A3EC","#FCC29C")
					if (input$sel_pal=="5") palette=c("#BC7C82","#A9FF05")
					if (input$sel_pal=="6") palette=c("#3B69F7","#B4B701")

					levels(clrs)=palette
					clrs=as.character(clrs)

					if (input$legend_loc!="none") legend(input$legend_loc,levels(data$g),col=palette,cex=1.5,text.font=2,pch=19)

					points(1:length(x),x,col="black",pch=21,bg=clrs,cex=1.5)

					if (input$scatter_label){
						text(1:length(x),x,colnames(data$em),font=2,cex=1,pos=3)
					}

					axis(2,at=yat,round(yat,2),las=2,cex.axis=1.5,font=2)
					axis(1,at=xat,round(xat,2),las=1,cex.axis=1.5,font=2)
				
					title(xlab="Sample Index",ylab=gene.orig,cex.lab=1.5,font.lab=2,line=4.5)
				}

				if (input$type.plot=="boxplot"){
					par(mar=c(6,6,3,2))
					boxplot(x~data$g,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
					rect(-1000,-1000,1000,1000,col=grey)
			
					xat=seq(1,length(x),by=floor(length(x)/7))
					yat=seq(min(x),max(x),length.out=7)
					abline(h=yat,col="white",lwd=1.5)

					clrs=data$g
					if (input$sel_pal=="1") palette=c("#56A8A4","#EBA94D")
					if (input$sel_pal=="2") palette=c("#63074E","#BFC421")
					if (input$sel_pal=="3") palette=c("#97A861","#9D6187") # c("#47BBD1","#840505")
					if (input$sel_pal=="4") palette=c("#60A3EC","#FCC29C")
					if (input$sel_pal=="5") palette=c("#BC7C82","#A9FF05")
					if (input$sel_pal=="6") palette=c("#3B69F7","#B4B701")

					levels(clrs)=palette
					clrs=as.character(clrs)

					b=boxplot(x~data$g,border=palette,add=T,range=c(0,1.5)[(input$hide.points=="hide")+1],
						col="white",pars=list(boxwex=0.5),yaxt="n",
						names=rep("",nlevels(data$g)),lwd=2,pch=21,outbg="white")

					if(input$hide.points=="show"){
						x.point=as.numeric(data$g)+runif(length(data$g),-0.2,0.2)
						points(x.point,x,col="black",pch=21,bg=clrs,cex=1.5)

						if (input$scatter_label){
							text(x.point,x,colnames(data$em),font=2,cex=1,pos=3)
						}
					}

					axis(2,at=yat,round(yat,2),las=2,cex.axis=1.5,font=2)
					axis(1,at=1:nlevels(data$g),levels(data$g),las=1,cex.axis=1.5,font=2)
				
					title(ylab=gene.orig,cex.lab=1.5,font.lab=2,line=4.5)
				}

			}

			#} else {

			if ((length(gene.orig)==2)&(class(data$g)=="factor")){
			
				# bivariate plot

				#geneX.orig=input$geneX
				#geneY.orig=input$geneY
				
				geneX.orig=gene.orig[1]
				geneY.orig=gene.orig[2]

				#geneX=make.names(geneX.orig)
				#geneY=make.names(geneY.orig)

				geneX=geneX.orig
				geneY=geneY.orig
		
				if (!(geneX %in% rownames(data$em))){
					stop(paste("'",geneX.orig,"'", "not in dataset",sep=""))
				}

				if (!(geneY %in% rownames(data$em))){
					stop(paste("'",geneY.orig,"'", "not in dataset",sep=""))
				}

				x=as.numeric(data$em[geneX,])
				y=as.numeric(data$em[geneY,])

				if (input$graph.adjusted){
					x=resid(lm(x~.,data=data.frame(x,data$adjust.factors)))
					y=resid(lm(y~.,data=data.frame(y,data$adjust.factors)))
				}
			
				par(mar=c(6,6,3,2))
				plot(x,y,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
				rect(-1000,-1000,1000,1000,col=grey)
			
				xat=seq(min(x),max(x),length.out=7)
				yat=seq(min(y),max(y),length.out=7)
				abline(h=yat,v=xat,col="white",lwd=1.5)

				clrs=data$g
				if (input$sel_pal=="1") palette=c("#56A8A4","#EBA94D")
				if (input$sel_pal=="2") palette=c("#63074E","#BFC421")
				if (input$sel_pal=="3") palette=c("#97A861","#9D6187") # c("#47BBD1","#840505")
				if (input$sel_pal=="4") palette=c("#60A3EC","#FCC29C")
				if (input$sel_pal=="5") palette=c("#BC7C82","#A9FF05")
				if (input$sel_pal=="6") palette=c("#3B69F7","#B4B701")

				levels(clrs)=palette
				clrs=as.character(clrs)

				if (input$legend_loc!="none") legend(input$legend_loc,levels(data$g),col=palette,cex=1.5,text.font=2,pch=19)

				points(x,y,col="black",pch=21,bg=clrs,cex=1.5)

				if (input$scatter_label){
					text(x,y,colnames(data$em),font=2,cex=1,pos=3)
				}


				axis(2,at=yat,round(yat,2),las=2,cex.axis=1.5,font=2)
				axis(1,at=xat,round(xat,2),las=1,cex.axis=1.5,font=2)
				
				title(xlab=geneX.orig,ylab=geneY.orig,cex.lab=1.5,font.lab=2,line=4.5)


			}

			if ((length(gene.orig)>2)&(class(data$g)=="factor")){
			
				# PCA plot

				#gene.orig=make.names(gene.orig)

				if (!all(gene.orig %in% rownames(data$em))){
					stop(paste("'",paste(gene.orig[!gene.orig %in% rownames(data$em)],collapse="','"),"'", "not in dataset",sep=""))
				}

				pca.matrix=t(data$em[gene.orig,])

				if (input$graph.adjusted){
					for (i in 1:ncol(pca.matrix)){
						pca.matrix[,i]=resid(lm(pca.matrix[,i]~.,data=data$adjust.factors))
					}
				}

				pca=PCA(pca.matrix,graph=F)

				x=pca$ind$coord[,1]
				y=pca$ind$coord[,2]
			
				par(mar=c(6,6,3,2))
				plot(x,y,type="n",xlab="",ylab="",xaxt="n",yaxt="n",asp=1)
				rect(-1000,-1000,1000,1000,col=grey)
			
				#xat=seq(min(x),max(x),length.out=7)
				#yat=seq(min(y),max(y),length.out=7)

				range=diff(c(min(c(x,y)),max(c(x,y))))
				by=min(c(0.2,0.5,1,2,5,10,Inf)[c(0.2,0.5,1,2,5,10,Inf) >(range/7)])
				xat=seq(-1000,1000,by=by)
				yat=seq(-1000,1000,by=by)

				abline(h=yat,v=xat,col="white",lwd=1.5)

				abline(h=0,v=0,lty=2,col="black",lwd=1.2)

				clrs=data$g
				if (input$sel_pal=="1") palette=c("#56A8A4","#EBA94D")
				if (input$sel_pal=="2") palette=c("#63074E","#BFC421")
				if (input$sel_pal=="3") palette=c("#97A861","#9D6187") # c("#47BBD1","#840505")
				if (input$sel_pal=="4") palette=c("#60A3EC","#FCC29C")
				if (input$sel_pal=="5") palette=c("#BC7C82","#A9FF05")
				if (input$sel_pal=="6") palette=c("#3B69F7","#B4B701")

				levels(clrs)=palette
				clrs=as.character(clrs)

				if (input$legend_loc!="none") legend(input$legend_loc,levels(data$g),col=palette,cex=1.5,text.font=2,pch=19)

				points(x,y,col="black",pch=21,bg=clrs,cex=1.5)

				if (input$scatter_label){
					text(x,y,colnames(data$em),font=2,cex=1,pos=3)
				}


				axis(2,at=yat,round(yat,2),las=2,cex.axis=1.5,font=2)
				axis(1,at=xat,round(xat,2),las=1,cex.axis=1.5,font=2)
				
				title(main="PCA",
					xlab=paste("Dim 1 (",round(pca$eig[1,2],2),"%)",sep=""),
					ylab=paste("Dim 2 (",round(pca$eig[2,2],2),"%)",sep=""),
					cex.lab=1.5,font.lab=2,line=4.5)

			}

			if (class(data$g)=="Surv"){
				# fer alguna kaplan-meier? com resumir el model? cox score partit per quartils i fent kaplan?

				df=subset(data.frame(merged.data),select=gene.orig)
				fit=coxph(data$g~.,data=df)
				qs=cut(predict(fit),breaks=4)

				levels(qs)=c("Q1","Q2","Q3","Q4")
				qs=factor(qs)

				par(mar=c(6,6,3,2))
				plot(survfit(data$g~qs),col=1:nlevels(qs),lwd=2,xlab="",ylab="",xaxt="n",yaxt="n")
				rect(-1000,-1000,1000,1000,col=grey)

				xat=seq(0,max(data$g[,1]),length.out=7)
				yat=seq(0,1,by=0.1)
				abline(h=yat,v=xat,col="white",lwd=1.5)

				lines(survfit(data$g~qs),col=1:nlevels(qs),lwd=2)

				axis(2,at=yat,paste(yat*100,"%",sep=""),las=2,cex.axis=1.5,font=2)
				axis(1,at=xat,round(xat,2),las=1,cex.axis=1.5,font=2)

				if (input$legend_loc!="none") legend(input$legend_loc,paste(levels(qs)," (n=",table(qs),")",sep=""),
									col=1:nlevels(qs),lty=1,lwd=2,cex=1.5,text.font=2)

				title(main="Kaplan-Meier curves of the Cox linear predictor quartiles",
					xlab="Time",ylab="Survival Probability",cex.lab=1.5,font.lab=2,line=4.5)


			}

		},res=100,width=800,height=800) #  ,height=500,res=80

		output$ui.plot <- renderUI({

			#gene.orig=strsplit(input$genes,split=",")[[1]]

			if (input$gener_or_biomarker=="Gene"){
				gene.orig=input$features
			} else {
				if (length(input$common.biomarkers)==0){
					gene.orig=NULL
				} else {
					gene.orig=strsplit(input$common.biomarkers,split=",")[[1]]
				}
			}

			what=length(gene.orig)
			if (what==1){
				radioButtons(inputId = "type.plot",label = h5("Type of Plot"),
					choices = list("Index Plot" ="Index", "Boxplot" = "boxplot"),selected = "boxplot",inline=TRUE)

			}

		})

		# es podria fer amb conditionalPanel normalment, però al dependre del ui.plot de damunt també ho he de fer així i poder posar el what
		output$ui.plot2 <- renderUI({

			#gene.orig=strsplit(input$genes,split=",")[[1]]

			if (input$gener_or_biomarker=="Gene"){
				gene.orig=input$features
			} else {
				if (length(input$common.biomarkers)==0){
					gene.orig=NULL
				} else {
					gene.orig=strsplit(input$common.biomarkers,split=",")[[1]]
				}
			}

			what=length(gene.orig)
			if (what==1){
				conditionalPanel(
					condition = "input['type.plot'] == 'boxplot'",

					radioButtons(inputId = "hide.points",label = h5("Sample Points"),
						choices = list("Hide" ="hide", "Show" = "show"),selected = "show",inline=TRUE)

				)
			}
		})


		# es pot fer només pujar la data, ja que no fa el update.calc()
		limma.calc <- reactive({

			withProgress(message="Calculating",detail="",value=0.5,{	

				data=dataInput()

				y=data$g
				x=data$em

				if (input$DGEtest=="limma"){
					if(input$limmaadjusted == "Yes"){

						if (ncol(data$adjust.factors)==0){
							stop("Missing Adjusting Factors File")
						}

						y=data.frame(g=data$g,data$adjust.factors)
						design <-model.matrix(~.,data=y)

						fit=lmFit(x,design)
						fit<-eBayes(fit)

						t=topTable(fit,number=nrow(x),coef=2) # coeficient dos pq g és la primera variable (de moment sol 2 categories)

						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"
						t=t[,!colnames(t) %in% "B"]

					} else {

						design <-model.matrix(~0+y)

						colnames(design)=levels(y)
						contr=combinations(n=length(levels(y)),r=2,v=levels(y))
						contr=paste(contr[,1],"-",contr[,2],sep="")
						contrat.matrix=makeContrasts(contrasts=contr,levels=levels(y))

						fit=lmFit(x,design)
						fit=contrasts.fit(fit,contrat.matrix)

						fit<-eBayes(fit)
						#t=topTable(fit,number=input$num_genes,p.value=input$min_p_val,lfc=input$min_fc)

						t=topTable(fit,number=nrow(x))

						#aux=t(apply(t(2^x[rownames(t),]),2,tapply,y,mean))
						aux=cbind(rowMeans(2^x[rownames(t),y==levels(y)[1]]),rowMeans(2^x[rownames(t),y==levels(y)[2]]))
						colnames(aux)=levels(y)

						fc_ratio=aux[,1]/aux[,2] # definició habitual de Fold Change
						lfc_ratio=log(aux[,1]/aux[,2],2) # i el seu logaritme (que no és el mateix que el limma)

						t=data.frame(FC=fc_ratio,log2FC=lfc_ratio,t)
						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"
						colnames(t)[colnames(t)=="log2FC"]="log2FC(Ratio)"
						t=t[,!colnames(t) %in% "B"]

					}
				}

				if (input$DGEtest %in% c("ttest.uneq")){
					if (input$limmaadjusted == "Yes"){
						stop ("T-Test (Unequal Variances) with factor adjusting not implemented yet")
						# que fer? lme amb var diferents no sembla funcionar. Quines opcions hi ha més? 
			
					} else {
						
						all.t=fast.t.test.mat(t(x),y=y,var.equal=FALSE)

						t=data.frame(logFC=all.t$estimate,AveExpr=rowMeans(x),t=all.t$t,
							P.Value=all.t$p.val,adj.P.Val=p.adjust(all.t$p.val,method="fdr"))

						rownames(t)=rownames(x)

						#aux=t(apply(t(2^x[rownames(t),]),2,tapply,y,mean))
						aux=cbind(rowMeans(2^x[rownames(t),y==levels(y)[1]]),rowMeans(2^x[rownames(t),y==levels(y)[2]]))
						colnames(aux)=levels(y)

						fc_ratio=aux[,1]/aux[,2] # definició habitual de Fold Change
						lfc_ratio=log(aux[,1]/aux[,2],2) # i el seu logaritme (que no és el mateix que el limma)

						t=data.frame(FC=fc_ratio,log2FC=lfc_ratio,t)
						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"
						colnames(t)[colnames(t)=="log2FC"]="log2FC(Ratio)"

					}
				}


				if (input$DGEtest %in% c("ttest.eq")){
					if (input$limmaadjusted == "Yes"){	

						y=data.frame(g=data$g,data$adjust.factors)
						design <-model.matrix(~.,data=y)

						t.function=function(x,y){lm(x~.,data=y)}
						all.t=apply(t(x),2,t.function,y=y)

						t.extract=function(t,what){return(summary(t)$coef[2,what])}

						t=data.frame(logFC=sapply(all.t,t.extract,what="Estimate"),AveExpr=rowMeans(x),
							t=sapply(all.t,t.extract,what="t value"),
							P.Value=sapply(all.t,t.extract,what="Pr(>|t|)"),adj.P.Val=p.adjust(sapply(all.t,t.extract,what="Pr(>|t|)"),method="fdr"))
						
						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"

					} else {

						all.t=fast.t.test.mat(t(x),y=y,var.equal=TRUE)

						t=data.frame(logFC=all.t$estimate,AveExpr=rowMeans(x),t=all.t$t,
							P.Value=all.t$p.val,adj.P.Val=p.adjust(all.t$p.val,method="fdr"))

						rownames(t)=rownames(x)

						#aux=t(apply(t(2^x[rownames(t),]),2,tapply,y,mean))
						aux=cbind(rowMeans(2^x[rownames(t),y==levels(y)[1]]),rowMeans(2^x[rownames(t),y==levels(y)[2]]))
						colnames(aux)=levels(y)

						fc_ratio=aux[,1]/aux[,2] # definició habitual de Fold Change
						lfc_ratio=log(aux[,1]/aux[,2],2) # i el seu logaritme (que no és el mateix que el limma)

						t=data.frame(FC=fc_ratio,log2FC=lfc_ratio,t)
						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"
						colnames(t)[colnames(t)=="log2FC"]="log2FC(Ratio)"
					}
				}


				if (input$DGEtest %in% c("wilcox")){
					if (input$limmaadjusted == "Yes"){
						stop("Non-parametric test with factor adjusting not implemented yet")
						# que fer? gam? quantile regression? wilcox dels residus?
					} else {
						
						# fast: sense continuit correction, posar warning? o una altra pestanya a ui?

						all.t=fast.wilcox(t(x),y)

						t=data.frame(
							#logDiffMedian=apply(apply(t(x),2,tapply,y,median),2,diff),
							logDiffMedian=rowMedians(as.matrix(x[,y==levels(y)[1]]))-rowMedians(as.matrix(x[,y==levels(y)[2]])),
							AveExpr=rowMeans(x),
							P.Value=all.t,adj.P.Val=p.adjust(all.t,method="fdr"))
						rownames(t)=rownames(x)

						# amb continuity correction quan hi ha ties o n grans
						#t.function=function(x,y){wilcox.test(x~y)$p.value}
						#all.t=apply(t(x),2,t.function,y=y)
						#t.extract=function(t,what){return(t[[what]])}
						#t=data.frame(logDiffMedian=apply(apply(t(x),2,tapply,y,median),2,diff),AveExpr=rowMeans(x),
						#	P.Value=all.t,adj.P.Val=p.adjust(all.t,method="fdr"))
						#rownames(t)=rownames(x)

						#aux=t(apply(t(2^x[rownames(t),]),2,tapply,y,mean))
						aux=cbind(rowMeans(2^x[rownames(t),y==levels(y)[1]]),rowMeans(2^x[rownames(t),y==levels(y)[2]]))
						colnames(aux)=levels(y)

						fc_ratio=aux[,1]/aux[,2] # definició habitual de Fold Change
						lfc_ratio=log(aux[,1]/aux[,2],2) # i el seu logaritme (que no és el mateix que el limma)

						t=data.frame(FC=fc_ratio,log2FC=lfc_ratio,t)
						colnames(t)[colnames(t)=="logFC"]="logFC(Diff)"
						colnames(t)[colnames(t)=="logDiffMedian"]="Med[Diff/log2]" # diferència de les medianes d'expressió en escala logarítimica
						colnames(t)[colnames(t)=="log2FC"]="log2FC(Ratio)"
						

					}
				}


			})

			return(t)
		})


		output$limma <- DT::renderDataTable({

			t=limma.calc()

			if(input$limmaadjusted == "No"){
				w=which(
					(
						(t$FC>input$min_fc) | (t$FC<(1/input$min_fc))
					) & (
						t$adj.P.Val<input$min_p_val
					)

				)
				sort.var=input$sort.by2
			} else {
				w=t$adj.P.Val<input$min_p_val
				sort.var=input$sort.by1
			}

			t=t[w,]

			t=data.frame(name=rownames(t),t,check.names=F)

			if (!sort.var %in% c("log2FC(Ratio)","Abs(logFC(Diff))")){
				t=t[order(t[,sort.var],decreasing=sort.var %in% c("FC","AveExpr")),]
			} else {
				if (sort.var=="Abs(logFC(Diff))") sort.var="logFC(Diff)"
				t=t[order(t[,sort.var]*sign(t[,sort.var]),decreasing=TRUE),]
			}

			t=data.frame(ranking=seq_len(nrow(t)),t,check.names=F)

			t=format(t,digits=3,nsmall=3)

			t

		},selection="none",rownames=FALSE, options = list(lengthMenu = c(10, 25, 50), pageLength = 25,ordering=FALSE)) 

		checks=function(check.download=FALSE){

			check1=FALSE

			na=0

			if (is.null(input$cel.file)){
				if (!check.download) cat("Number of CEL files loaded: 0\n")
			} else {
				na=sum(grepl(".CEL",input$cel.file$name,fixed=TRUE))
				if (!check.download) cat(paste("Number of CEL files loaded:",na,"\n"))

				if (na>0) {
					check1=TRUE
				} else {
					check1=FALSE
				}

				# Nota: que passa si pugem un arxiu que no es CEL file?
			}
	
			if (input$norm.method == " "){
				if (!check.download) cat("Normalization Method: Not Specified\n")
				check2=FALSE
			} else {
				if (!check.download) cat(paste("Normalization Method:",switch(input$norm.method,mas5="MAS5.0",rma="RMA"),"\n"))
				check2=TRUE
			}

			if (!check.download){
				
				cat(paste("Filter out features without an Entrez Gene ID annotation: ",
					c("NO","YES")[1+as.logical(input$entrez.id)],"\n",sep=""))

				cat(paste("Filter out Affymetrix QC Probesets: ",
					c("NO","YES")[1+(input$feat.exclude=="^AFFX")],"\n",sep=""))

				cat(paste("Select the feature with largest ",
					switch(input$func.filter,IQR="IQR",sd="Standard Deviation"),
					" among the features mapping to the same Entrez Gene ID: ",
					c("NO","YES")[1+as.logical(input$dup.entrez)],"\n",sep=""))

				cat(paste("Percentatge of filtered out features (by ",
					switch(input$func.filter,IQR="IQR",sd="Standard Deviation"),"): ",
					c(paste(0+(input$filter.percent)*as.logical(input$performVARfiltering),"%"),
						"{Not done with less than 4 samples}")[(na<4)+1],
					"\n",sep=""))
			}

			if (input$type.array!="Probeset"){
				if (!check.download) cat(paste("Output matrix feature names: ",input$type.array,"\n",sep=""))

			} else {
				if (!check.download) cat("Output matrix feature names: Probeset\n")

			}

			if (!check.download) cat("\n\n")

			if (all(check1,check2)){
				if (!check.download) cat("Data Ready for Download\n")
			} else{
				if (!check.download) cat("Data NOT Ready for Download\n")
			}

			if (!check.download) cat("\n\n")

			if (check.download) return(all(check1,check2))

		}

		output$summary <- renderPrint({
			checks(check.download=FALSE)
		})

		output$downloadData.P <- downloadHandler(

			#filename="preprocessed_data.txt",
			filename=function(){
				paste("preprocessed_data", ".txt", sep="")
			},

			# This function should write data to a file given to it by
			# the argument 'file'.
			content = function(file) {

				if (checks(check.download=TRUE)){

					cel.paths=input$cel.file

					withProgress(message="Loading CEL files...",detail="",value=0/3,{	
						Data=ReadAffy(filenames=cel.paths$datapath)
						sampleNames(Data)=as.character(cel.paths$name)
					})

					withProgress(message="Normalizing...",detail="",value=1/3,{
						Dnorm = switch(input$norm.method,mas5=mas5(Data), rma=rma(Data))
					})

					withProgress(message="Filtering...",detail="",value=2/3,{

						cond=any(
							as.logical(input$entrez.id),
							as.logical(input$dup.entrez),
							input$type.array=="Gene"
						)

						if (cond){

							### Cargar annotation package. Si es fa alguna condició de dalt. Si no, no fa falta

							#array.db=input$type.array # tal com estaba abans, indiquem tipus array al UI
							array.db=paste(annotation(Data),".db",sep="")

							# comprobar si paquet existeix (WARNING: tal com està ara checkeja paquets normals, no els d'anotacións)
							#conn=url(paste(biocinstallRepos()["BioCsoft"],"/src/contrib/PACKAGES",sep=""))
							#list.pack=read.dcf(conn)[,"Package"]
							#close(conn)
							#if (!array.db %in% list.pack){
							#	stop("paquet no existeix?")
							#}

							# no sé si funciona, fatalrà fer més proves
							if (!require(array.db,character.only=TRUE)){
								source("http://bioconductor.org/biocLite.R")
								#source("https://bioconductor.org/biocLite.R")
								biocLite(array.db)
							}

							library(array.db,character.only=TRUE)

						}

						cond=any(
							as.logical(input$entrez.id),
							as.logical(input$dup.entrez),
							input$feat.exclude!="something.impossible",
							as.logical(input$performVARfiltering)
						)

						if (cond){
						
							if ((length(sampleNames(Dnorm))<4) & as.logical(input$performVARfiltering)){
								stop("Too few samples to perform reliable variance filtering")
							}

							### Filtrar

							Dfiltr=nsFilter(Dnorm, require.entrez=as.logical(input$entrez.id), 
								remove.dupEntrez=as.logical(input$dup.entrez),
								feature.exclude=input$feat.exclude, 
								var.filter=as.logical(input$performVARfiltering),
								var.func=switch(input$func.filter,IQR=IQR, sd=sd),
								var.cutoff=input$filter.percent/100)$eset

						} else {
							Dfiltr=Dnorm
						}
					})

					# Aquí hi ha algo que va molt lent...
					withProgress(message="Preparing Data for Download...",detail="",value=3/3,{

						em=exprs(Dfiltr)

						if (input$type.array=="Gene"){

							#rownames(em)=make.names(getSYMBOL(rownames(em),array.db),unique=TRUE)

							gene.anot=getSYMBOL(rownames(em),array.db)
							#gene.anot[!is.na(gene.anot)]=make.names(gene.anot[!is.na(gene.anot)],unique=TRUE)
							gene.anot[!is.na(gene.anot)]=make.unique(gene.anot[!is.na(gene.anot)],sep=" &")

							rownames(em)[!is.na(gene.anot)]=gene.anot[!is.na(gene.anot)]

						}

						if (input$norm.method == "mas5"){em=log(em,2)}

						em=data.frame(feature=rownames(em),em,check.names=FALSE)
						#em=data.frame(SampleID=colnames(em),t(em),check.names=FALSE)

						em=em[order(em$feature),]						

						write.table(em, file, sep = "\t",row.names = FALSE,quote=F)

					})

				} else {stop ("Data NOT ready for download")}
			}
		)


	}
)

