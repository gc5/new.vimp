

#library(RankProd)

# Nota: se li passarà les dades com és habitual en estadística. Rows=Samples i Cols=Variables

##### FUNCIONIGRAMA

# robust.combination
# robust.feature

# coses a pensar: crees noves boot samples després de triar, o fas servir les mateixes? 
# podem conservar, cal conservar les variables que s'han seleccionat?
# si, per després comptar el % de vegades que ha sortit del possible?

# predicció dels OOB, comptar quantes vegades han sortit
# primer, comptatge de cada biomarcador, quantes vegades seleccionat
# cal fer el tune? Millora molt la predicció, però millora la comparació?
# kernel serà linear, té sentit pq busquem una separació lineal en el fons




# vigilar amb el % de vegades seleccionat. Fer-ho per número igualment?

# forma de fer un mapa de variables que surten juntes o variables que no surten mai juntes?


# inds, matriu que té les mostres bootstrap a seleccionar (per tant, tenim OOB també, i no fa falta tenir OOB.percent)
# vars, matriu amb variables a seleccionar. Similar a inds però per variables
# X: matriu d'expressió
# Y: grup
# var.


#### FUNCIONS PER SELECCIONAR ELS BIOMARKERS

# totes les seleccions en paral·lel

f0=function(list,i,j,k){
	return(list[[k]][[i]][[j]])
}

sample.percent=function(x,percent){
	if (length(x)<5) percent=1 # si una categoria té menys de 5 casos els agafarem tots
	q=floor(length(x)*percent)
	return(sort(sample(x,q)))
}

repartir.folds=function(x,nfolds){

	reps=floor(length(x)/nfolds)
	folds=rep(1:nfolds,each=reps)
	folds=sample(c(folds,sample(1:nfolds,length(x)-length(folds))))

	return(data.frame(x,folds))
}

ext=function(list,which){
	return(list[[which]])
}

bioselect.paralel.lasso=function(XX,Y,M=100,keep.filter=200,method.filter="limma",adjust.filter=FALSE,OOB.percent=0.2,var.percent=0.5,stan=FALSE,weights.vec,
	adjust.factors,mfc=1,alpha=1,pmax=NA,dfmax=NA,lambda="lambda.1se",mc=2){

	if (class(Y)=="factor"){
		#inds=replicate(M,sample(nrow(XX),ceiling(nrow(XX)*(1-OOB.percent))))
		inds=replicate(M,do.call(c,tapply(1:length(Y),Y,sample.percent,percent=1-OOB.percent)))
	}

	if (class(Y)=="Surv"){
		inds=replicate(M,sample(nrow(XX),ceiling(nrow(XX)*(1-OOB.percent))))
		#inds=replicate(M,do.call(c,tapply(1:length(Y),Y,sample.percent,percent=1-OOB.percent)))
	}

	vars=replicate(M,sample(ncol(XX),ceiling(ncol(XX)*var.percent)))

	# seria interessant forçar que totes les variables i tots els individus sortissin les mateixes vegades

	sortida=mclapply(1:M,bio.selec.lasso,inds=inds,vars=vars,XX=XX,Y=Y,keep.filter=keep.filter,method.filter=method.filter,adjust.filter=adjust.filter,
		stan=stan,weights.vec=weights.vec,adjust.factors=adjust.factors,mfc=mfc,alpha=alpha,pmax=pmax,dfmax=dfmax,lambda=lambda,mc.cores=mc)

	cv.errors=do.call(c,lapply(sortida,ext,which=2)) # error
	cv.spec=do.call(c,lapply(sortida,ext,which=3)) # sens
	cv.sens=do.call(c,lapply(sortida,ext,which=4)) # spec
	sortida=lapply(sortida,ext,which=1) # variables
	
	# bios hauria de tenir una llista per cada numero de variables

	noms=unique(names(do.call(c,sortida)))
	matrix=matrix(NA,nrow=M,ncol=length(noms))
	colnames(matrix)=noms

	for (i in 1:length(sortida)){
		matrix[i,names(sortida[[i]])]=sortida[[i]]
	}

	# canviar vars i inds a un format més tractable. Matriu amb fila indicant si entra o no

	vars2=matrix(0,ncol=ncol(XX),nrow=M)
	colnames(vars2)=colnames(XX)
	for (i in 1:ncol(vars)){
		vars2[i,vars[,i]]=1
	}
	vars=vars2	

	inds2=matrix(0,ncol=nrow(XX),nrow=M)
	colnames(inds2)=rownames(XX)
	for (i in 1:ncol(inds)){
		inds2[i,inds[,i]]=1
	}
	inds=inds2	

	# tecnicament els NA són 0's

	return(list(bios=matrix,inds=inds,vars=vars,sortida=sortida,cv.errors=cv.errors,cv.spec=cv.spec,cv.sens=cv.sens))

}


# una iteració

bio.selec.lasso=function(w,inds,vars,XX,Y,keep.filter=200,method.filter="limma",adjust.filter=FALSE,stan=FALSE,alpha=1,pmax=NA,dfmax=NA,weights.vec,adjust.factors,
	mfc=1,lambda="lambda.1se"){

	bootX=subset(XX,subset=(1:nrow(XX)) %in% inds[,w],select=vars[,w])
	bootY=Y[(1:nrow(XX)) %in% inds[,w]]

	weights=weights.vec[rownames(bootX)]

	# treiem les que tenen el fold change més petit de mfc (és molt lent, com fer-ho més ràpid?) ($$$)
	# això només es farà si no considerem factor adjusting. 

	if (class(bootY)=="factor"){
		if (mfc>1) {

			# més lent
			#fc=apply((2^bootX),2,tapply,bootY,mean)
			#fc=apply(fc,2,max)/apply(fc,2,min)

			# més ràpid
			fc=by(2^bootX,INDICES=bootY,FUN=colMeans)
			fc=fc[[1]]/fc[[2]]

			bootX=bootX[,(fc>mfc)|(fc<(1/mfc))]

		}
	}

	# filter (independent de adjusting factors) fer-ho dependent també. Seria la gràcia.

	if (class(bootY)=="factor"){
		bootX=filter(bootX,bootY,kf=keep.filter,method=method.filter,adjust.filter=adjust.filter,
			adjust.factors=adjust.factors[rownames(bootX),],weights=weights)
	}
	if (class(bootY)=="Surv"){
		bootX=filter.surv(bootX,bootY,kf=keep.filter,weights=weights)
	}

	# merge variables i adjusting factors

	if (ncol(adjust.factors)>0){
		model.matrix.adj.fact=subset(model.matrix(~.,data=adjust.factors[rownames(bootX),,drop=FALSE]),select=-1)
		bootX=data.frame(model.matrix.adj.fact,bootX)
	}

	# aquí s'hauria d'escalar les variables
	if (stan){ 
		bootX=scale(bootX)	
	}
	# biomarker per lasso (li adjudico jo els folds, per assegurar que hi hagi almenys tres casos?)

	# opcio 1: l'habitual de nfolds
	#nfolds=min(min(table(bootY)),10)
	#cvfit = cv.glmnet(as.matrix(bootX), bootY, family = "binomial", type.measure = "class",nfolds=nfolds,standardize=stan,alpha=1) 

 	### opció 2: jo li adjudico com repartir els folds

	if (class(bootY)=="factor"){
		# creacio foldid (estratificat)
		nfolds=min(min(table(bootY)),10)
		indicador.folds=do.call(rbind,tapply(1:length(bootY),bootY,repartir.folds,nfolds=nfolds))
		foldid=indicador.folds[order(indicador.folds[,1]),2]

		# ajust foldid (potser queda un fold amb una mostra mes, canviar-la) No sé si fer-ho...
		if (diff(range(table(foldid)))>1){
			wf=which.max(table(foldid))
			foldid[sample(which(foldid==wf),1)]=which.min(table(foldid))
		}
	}

	if (class(bootY)=="Surv"){

	}

	## pmax i dfmax

	nvars=ncol(bootX)
	if (is.na(dfmax)){dfmax = nvars+1}
	if (is.na(pmax)){pmax = min(dfmax * 2+20, nvars)}

	if (class(bootY)=="factor"){
		cvfit = cv.glmnet(as.matrix(bootX), bootY, family = "binomial", type.measure = "class",foldid=foldid,standardize=FALSE,alpha=alpha,keep=T,
			pmax=pmax,dfmax=dfmax,weights=weights)
	}

	if (class(bootY)=="Surv"){
		cvfit = cv.glmnet(as.matrix(bootX), bootY, family = "cox", type.measure = "deviance",nfolds=10,standardize=FALSE,alpha=alpha, #keep=T,
			pmax=pmax,dfmax=dfmax,weights=weights)
	}

	if (class(bootY)=="numeric"){
		cvfit = cv.glmnet(as.matrix(bootX), bootY, family = "gaussian", type.measure = "deviance",nfolds=10,standardize=FALSE,alpha=alpha, #keep=T,
			pmax=pmax,dfmax=dfmax,weights=weights)
	}

	# opció 3: fer loocv
	#cvfit = cv.glmnet(as.matrix(bootX), bootY, family = "binomial", type.measure = "class",nfolds=length(bootY),standardize=stan,alpha=1)

	# alpha=1 és lasso
	# standardize si/no?
	# sobre standaritzar o no: guay si la pregunta és quin predictor afecta the most. Però s'ha de tenir en compte que alguns predictors no es poden
	# canviar. Per exemple, canviar una standard unit del temps que mires el televisor és factible, però canviar una standard unit de la teva altura
	# no. Si es olgués determinar que fer per tenir una vida més saludable, canviar l'altura no és una opció encara que tingui un coeficient més gran
	# de moment deixo amb FALSE ja que sembla que en unes proves incials semblava encertar més

	#Smin=coef(cvfit, s = "lambda.min")[abs(coef(cvfit, s = "lambda.min")[,1])>0,]
	#S=coef(cvfit, s = "lambda.1se")[abs(coef(cvfit, s = "lambda.1se")[,1])>0,]

	S=coef(cvfit, s = lambda)[,1] #lambda = "lambda.1se" o "lambda.min"
	S=S[abs(S)>0]

	if(ncol(adjust.factors)>0){
		if (any(names(S) %in% colnames(model.matrix.adj.fact))){
			S=S[!names(S) %in% colnames(model.matrix.adj.fact)]
		}
	}

	# lo d'aquí sota no cal tal com s'agafa la S
	#if (!is.null(dim(S))){ # quan s'agafen totes les variables S és sparse matrix amb una columna i té dimensions. 
	#	S=S[,1]
	#}

	#if (length(S)==1) {names(S)="(Intercept)"}
	if (length(S)==0) { # xorrada que faig per evitar problemes amb cox
		S=0
		names(S)="(Intercept)"
	}

	#lsel=which(cvfit$lambda %in% cvfit$lambda.1se)
	lsel=which(cvfit$lambda %in% cvfit[lambda][[1]])


	error=cvfit$cvm[lsel]

	# per fer sensitivity i specificitiy (fer 1- per ser equivalent a error rate. Per tant un valor alt indicarà que molts casos que ho són els dones com no)
	# aixó no s'ha de fer pel cox, la resta es pot deixar tot tal com està (excepte tenir en compte quan S no té variables, ja que no hi ha intercept)

	if (class(bootY)=="factor"){
		pred=factor(cvfit$fit.preval[,lsel]>0.5,levels=c("FALSE","TRUE"))
		levels(pred)=levels(bootY)
		t=table(pred,bootY)
		specificity=1-(t[1,1]/sum(t[,1]))
		sensitivity=1-(t[2,2]/sum(t[,2])) # del segon grup (que logistic agafa com 1 o tractament) quants predits correctament
	}

	if (class(bootY)=="Surv"){
		specificity=NA
		sensitivity=NA

	}

	# vigilar si només hi ha una variable o cap! En principi el intercept serà sempre més de 0, així que si només és una és intercept
	# he d'anar revisant les sortides, per si pogues passar el cas inimaginable de que el intercept també és .


	#incProgress(amount = 1/ncol(inds), message = NULL, detail = paste(round(w*100/ncol(inds),2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = w/ncol(inds), message = NULL, detail = paste(round(w*100/ncol(inds),2),"%"),session = getDefaultReactiveDomain())

	return(list(S,error,specificity,sensitivity))

}

mesure=function(x,y,type=c("mce","ba")){

	x=factor(x,levels=c("FALSE","TRUE"))
	t=table(x,y)
	
	mce=sum(diag(t))/sum(t)

	spec=t[1,1]/sum(t[,1])
	sens=t[2,2]/sum(t[,2])
	ba=(spec+sens)/2

	if (type=="mce") return(1-mce)
	if (type=="ba") return(1-ba)

}





# en una regressió logística es podria ajustar el cutoff de la probabilitat predita a la que millor classifica
# però aquí farem com al glmnet (lasso), que és assignar a la classe amb màxima probabilitat
#myglm <- function(formula,data){
#	mod <- glm(formula,family=binomial(logit),data=data)
#	function(newdata,data){
#		p=predict(mod,newdata)
#		p=factor(p>0.5,levels=c("FALSE","TRUE"))
#		levels(p)=levels(data$y)
#		return(p)
#	}
#}


# repeat k vegades el k-fold cross-validation (només una, no crec que sigui molt problema)

individual.power=function(var,XX,Y,tant,weights.vec,adjust.factors){

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100*(var+5)/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = (var+5)/tant, message = NULL, detail = paste(round(100*(var+5)/tant,2),"%"),session = getDefaultReactiveDomain())

	XX=XX[,var]

	ba=numeric()
	for (k in 1:5){
		nfolds=min(min(table(Y)),10)
		indicador.folds=do.call(rbind,tapply(1:length(Y),Y,repartir.folds,nfolds=nfolds))
		foldid=indicador.folds[order(indicador.folds[,1]),2]

		prediccions=numeric()

		if (ncol(adjust.factors)>0){
			data=data.frame(y=Y,x=XX,adjust.factors[names(XX),,drop=FALSE])
		} else {
			data=data.frame(x=XX,y=Y)
		}

		for (i in 1:nfolds){

			w=which(foldid==i)
			#data=data.frame(x=XX,y=Y)[-w,]

			weights=weights.vec[rownames(data[-w,])]

			#fit=glm(y~x,data=data,family=binomial(logit),weights=weights)
			fit=glm(y~.,data=data[-w,],family=binomial(logit),weights=weights)

			#prediccions[w]=predict(fit,newdata=data.frame(x=XX[w]),type="response")
			prediccions[w]=predict(fit,newdata=subset(data[w,],select=-y),type="response")
		}

		# global error? o balanced accuracy?
		# de moment balanced accuracy, per distingir si es fa malament 1 de 4 HCLv que 1 de 30 Others

		#predy=prediccions>0.5 # dona error si no hi ha cap TRUE o FALSE
		predy=factor(prediccions>0.5,levels=c("FALSE","TRUE"))

		t=table(data.frame(predy,Y))

		spec=t[1,1]/sum(t[,1])
		sens=t[2,2]/sum(t[,2])
		ba[k]=(spec+sens)/2
	}

	return(mean(ba))
}

# leave-one-out cross-validation (LOOCV error) # així no li afecta com queden distribuits els casos

individual.power=function(var,XX,Y,tant,weights.vec,adjust.factors){

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100*(var+5)/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = (var+5)/tant, message = NULL, detail = paste(round(100*(var+5)/tant,2),"%"),session = getDefaultReactiveDomain())

	XX=XX[,var]

	prediccions=numeric()

	if (ncol(adjust.factors)>0){
		data=data.frame(y=Y,x=XX,adjust.factors[names(XX),,drop=FALSE])
	} else {
		data=data.frame(x=XX,y=Y)
	}

	for (i in 1:nrow(data)){

		weights=weights.vec[rownames(data[-i,])]
		fit=glm(y~.,data=data[-i,],family=binomial(logit),weights=weights)
		prediccions[i]=predict(fit,newdata=data[i,,drop=F],type="response")
	}

	# global error? o balanced accuracy?
	# de moment balanced accuracy, per distingir si es fa malament 1 de 4 HCLv que 1 de 30 Others
	predy=factor(prediccions>0.5,levels=c("FALSE","TRUE"))

	t=table(data.frame(predy,Y))

	spec=t[1,1]/sum(t[,1])
	sens=t[2,2]/sum(t[,2])
	ba=(spec+sens)/2
	
	return(ba)
}

# far el % de utilitar gen X amb gen Y
comparison.categoric.mod <- function(x,vars){

	odds=matrix(NA,ncol=ncol(x),nrow=ncol(x))

	way=2

	for (i in 1:ncol(x)){
		for (j in 1:ncol(x)){

			# els dos way són P(A/B), la diferència és que al 1 la P(B) l'estimes amb les iteracions que han sortit
			# les dos variables. Al way 2 fas l'estimació de P(B) amb totes les iteracions en que ha sortit B

			# way 1 = % de vegades els dos quan s'ha seleccionat B.

			if (way==1){
				w.it=(apply(vars[,colnames(x)[c(i,j)]],1,sum)==2) # iteracions que tenen les dos variables
				v1=factor(x[w.it,i],levels=c("FALSE","TRUE"))
				v2=factor(x[w.it,j],levels=c("FALSE","TRUE"))
				t=table(v1,v2)
			
				odds[i,j]=t[2,2]/(sum(t[2,]))
			}

			# way 2 = P(A/B) = P(A,B)/P(B) # prob de sel A si s'ha seleccionat B
			# tal com està fet, es pot comparar per columna, tota la fila divida pel mateix valor, de forma
			# que es poden comparar les columnes
			# deixar en blanc la matriu diagonal superior?

			if (way==2){
				w.it=(apply(vars[,colnames(x)[c(i,j)]],1,sum)==2) # iteracions que tenen les dos variables

				v1=factor(x[w.it,i],levels=c("FALSE","TRUE"))
				v2=factor(x[w.it,j],levels=c("FALSE","TRUE"))
				t=table(v1,v2)
				P.A.B=t[2,2]/sum(t) # prob que les dos (en les que poden sortir les dos)

				w.it=(vars[,colnames(x)[i]]==1) # iters que tenen variable i
				v1=factor(x[w.it,i],levels=c("FALSE","TRUE"))
				P.B=mean(v1=="TRUE")

				if (is.finite(P.A.B) & is.finite(P.B)){
					odds[i,j]=P.A.B/P.B
					if (P.B<P.A.B) odds[i,j]=1 # com fer-ho?
					if (P.B<0.05) odds[i,j]=NA
				} else {
					odds[i,j]=NA
				}
				#if (j>=i) odds[i,j]=NA
				# Per aquest way, modificar
				# 1. Escala, no va de 0 a 100
				# 2. Diagonal superior borrar, o millor moure matriu pq sigui diagonal invertida i aixi no moure noms

			}


			if (way==3){
				w.it=(apply(vars[,colnames(x)[c(i,j)]],1,sum)==2) # iteracions que tenen les dos variables
				v1=factor(x[w.it,i],levels=c("FALSE","TRUE"))
				v2=factor(x[w.it,j],levels=c("FALSE","TRUE"))
				t=table(v1,v2)
				P.A.B=t[2,2]/sum(t) # prob que les dos (en les que poden sortir les dos)

				w.it=(vars[,colnames(x)[i]]==1) # iters que tenen variable i
				v1=factor(x[w.it,i],levels=c("FALSE","TRUE"))
				P.B=mean(v1=="TRUE")

				w.it=(vars[,colnames(x)[j]]==1) # iters que tenen variable i
				v2=factor(x[w.it,j],levels=c("FALSE","TRUE"))
				P.A=mean(v2=="TRUE")

				odds[i,j]=P.A.B/(P.A*P.B)
				if ((P.B<0.025)|(P.A<0.025)) odds[i,j]=NA
				#if (j>=i) odds[i,j]=NA
				# Per aquest way, modificar
				# 1. Escala, no va de 0 a 100
				# 2. Diagonal superior borrar, o millor moure matriu pq sigui diagonal invertida i aixi no moure noms

			}


		}
		#incProgress(amount = 1/(2*ncol(x)+9), message = NULL, detail = paste(round(100*(ncol(x)+8+i)/(ncol(x)*2+8),2),"%"),session = getDefaultReactiveDomain())
		setProgress(value = (ncol(x)+8+i)/(ncol(x)*2+8), message = NULL, 
			detail = paste(round(100*(ncol(x)+8+i)/(ncol(x)*2+8),2),"%"),session = getDefaultReactiveDomain())
	}

	colnames(odds)=rownames(odds)=colnames(x)
	return(odds)
}






# funció que fa limma o mann-whitney per filtrar

filter=function(x,y,kf,method,adjust.factors,adjust.filter,weights){


	# FC sense tenir en compte weihgts
	if (method=="unweighted.fc"){ # més rapid que el de damunt?
		
		# més lent
		#aux=t(apply(2^x,2,tapply,y,mean))
		#fc_ratio=apply(aux,1,max)/apply(aux,1,min)

		# més ràpid
		aux=by(2^x,INDICES=y,FUN=colMeans)
		fc_ratio=aux[[1]]/aux[[2]]
		fc_ratio[fc_ratio<1]=1/(fc_ratio[fc_ratio<1])

		names(fc_ratio)=colnames(x)
		fc_ratio=sort(fc_ratio,decreasing=TRUE)
		
		x.filter=subset(x,select=names(fc_ratio)[1:kf])
	}

	# tenint en compte weights
	if (method=="fc"){
		
		aux.weights=weights/tapply(weights,y,sum)[as.character(y)]
		aux=by(sweep(2^x,1,aux.weights,FUN="*"),INDICES=y,FUN=colSums)
		fc_ratio=aux[[1]]/aux[[2]]
		fc_ratio[fc_ratio<1]=1/(fc_ratio[fc_ratio<1])

		names(fc_ratio)=colnames(x)
		fc_ratio=sort(fc_ratio,decreasing=TRUE)
		
		x.filter=subset(x,select=names(fc_ratio)[1:kf])
	}


	# mann.whitney/wilcoxon (sense weights)
	if (method=="non.parametric"){

		f=function(x,y){wilcox.test(x~y)$p.value}
		pvals=apply(x,2,f,y=y)
		names(pvals)=colnames(x)
		pvals=sort(pvals,decreasing=FALSE)
	
		x.filter=subset(x,select=names(pvals)[1:kf])	
	}
	# rank prod?

	# limma
	if (method=="limma"){

		if (adjust.filter){

			y=data.frame(y=y,adjust.factors)

			design <-model.matrix(~.,data=y)

			fit=lmFit(t(x),design,weights=weights)
			fit<-eBayes(fit)

			t=topTable(fit,number=kf,coef=2) # coeficient dos pq g és la primera variable (de moment sol 2 categories)

		} else {

			design <-model.matrix(~0+y)
			colnames(design)=levels(y)
			contr=combinations(n=length(levels(y)),r=2,v=levels(y))
			contr=paste(contr[,1],"-",contr[,2],sep="")
			contrat.matrix=makeContrasts(contrasts=contr,levels=levels(y))

			fit=lmFit(t(x),design,weights=weights)
			fit=contrasts.fit(fit,contrat.matrix)

			metode=1 # 1: eBayes i filtrant per lfc. 2: treat i agafant que H0:lfc=1, H1: lfc>1

			# eBayes i filtrant per lfc
			if(metode==1){
				fit<-eBayes(fit)

				t=topTable(fit,number=kf) #,lfc=1) # lfc de 1 equival al doble 2^1, el 0.5 a 1.4
				#if (nrow(t)<5) t=topTable(fit,number=kf,lfc=0.5)
				#if (nrow(t)<5) t=topTable(fit,number=kf)
			}

			# trat i agafant els que significativament poden tenir un lfc de 1, 0.5
			if (metode==2){
				fit<-treat(fit,lfc=1)
				t=topTable(fit,number=kf) # lfc de 1 equival al doble 2^1, el 0.5 a 1.4
				t=t[t$adj.P.val<0.05,]

				if (nrow(t)<5){
					fit<-treat(fit,lfc=0.5)
					t=topTable(fit,number=kf) # lfc de 1 equival al doble 2^1, el 0.5 a 1.4
					t=t[t$adj.P.val<0.05,]
				}

				if (nrow(t)<5){ 
					fit<-treat(fit,lfc=0)
					t=topTable(fit,number=kf) # lfc de 1 equival al doble 2^1, el 0.5 a 1.4
					t=t[t$adj.P.val<0.05,]
				}

				if (nrow(t)<5){ 
					fit<-treat(fit,lfc=0)
					t=topTable(fit,number=kf) # lfc de 1 equival al doble 2^1, el 0.5 a 1.4
				}

			}

		}
	
		# el topTable a vegades fica rownames i a vegades columna ID (crec que depen noms repetits)
		if (sum(colnames(t)=="ID")>0){
			x.filter=subset(x,select=t$ID)
		}else{ x.filter=subset(x,select=rownames(t))}
		# fi parra, en principi només em caldria deixar x.filter=x[rownames(t),]

	}

	return(x.filter)
}


filter.surv=function(x,y,kf,weights){

	pvals=apply(x,2,surv.ind,y=y,param=2)

	names(pvals)=colnames(x)
	pvals=sort(pvals,decreasing=FALSE)
	
	kf=min(kf,length(pvals))
	x.filter=subset(x,select=names(pvals)[1:kf])	

	return(x.filter)
}


#### funcions per fer els calculs i el gràfic

# aquesta funció fa el càlcul per totes les variables que han sortit!!!
# això s'ha de modificar, no cal que sigui així
# modificar per a que ho faci amb un màxim de 40-50 variables, guanyaríem molt temps

calculs=function(quant,sortida.lasso,em,g,weights.vec,adjust.factors){

	tant=2*quant+9

	#### CÀLCULS NECESSARIS PER FER EL GRÀFIC

	#coefs=sortida.lasso$bios[,-1] # el 1 és el intercept?
	#coefs=subset(sortida.lasso$bios,select=-1)
	coefs=subset(sortida.lasso$bios,select=colnames(sortida.lasso$bios)[!colnames(sortida.lasso$bios) %in% "(Intercept)"])

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 1/tant, message = NULL, detail = paste(round(100/tant,2),"%"),session = getDefaultReactiveDomain())


	# llargada biomarkers (tallo a algun numero o el deixo complert?)
	aux=apply(!is.na(coefs),1,sum)
	#aux[aux>20]=20
	#length.bios=prop.table(table(aux))
	length.bios=prop.table(table(factor(aux,levels=0:max(aux))))


	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(200/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 2/tant, message = NULL, detail = paste(round(200/tant,2),"%"),session = getDefaultReactiveDomain())

	# nombre de altres variables en el biomarker quan variable inclosa
	llargada=apply(!is.na(coefs),1,sum)
	llargada=sweep(coefs!=0,1,STATS=llargada,FUN="*")-1

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(300/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 3/tant, message = NULL, detail = paste(round(300/tant,2),"%"),session = getDefaultReactiveDomain())

	# percentatge selecciat variable
	sel=apply(!is.na(coefs),2,sum)
	total.sel=apply(sortida.lasso$vars,2,sum)[names(sel)]
	perc.sel=100*sel/total.sel

	ordre=names(perc.sel)[order(perc.sel,decreasing=T)]
	ordre=ordre[seq_len(quant)]

	coefs=coefs[,ordre] # ja seleccionem les variables que farem servir
	perc.sel=perc.sel[ordre]

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(400/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 4/tant, message = NULL, detail = paste(round(400/tant,2),"%"),session = getDefaultReactiveDomain())

	# de # compis agafem el maxim que podem plotar
	llargada=llargada[,names(perc.sel)]
	#median.llargada=apply(llargada,2,median,na.rm=T)

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(500/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 5/tant, message = NULL, detail = paste(round(500/tant,2),"%"),session = getDefaultReactiveDomain())

	# individual power (calculat només en els que han sigut seleccionats alguna vegada)

	mc.cores=detectCores()-1
	#mc.cores=1

	if (mc.cores<1) mc.cores=1


	if (class(g)=="factor"){
		ba=mclapply(1:length(perc.sel),FUN=individual.power,XX=t(em)[,names(perc.sel)],Y=g,weights.vec=weights.vec,
			adjust.factors=adjust.factors,tant=tant,mc.cores=mc.cores)
		ba=do.call(c,ba)
	} else {
		sXX=t(em)[,names(perc.sel)]

		ba=apply(sXX,2,surv.ind,y=g,param=2) # p-value
		ba=ba+10^-17 # ja que després fa logaritme i no porta bé el 0 (que surt quan p-val < 2e-16)
		fc=exp(apply(sXX,2,surv.ind,y=g,param=1)) # coef
	}

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100*(quant+6)/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = (quant+6)/tant, message = NULL, detail = paste(round(100*(quant+6)/tant,2),"%"),session = getDefaultReactiveDomain())

	# fold change (calculat només en els que han sigut seleccionats alguna vegada)
	if (class(g)=="factor"){
		if (ncol(adjust.factors)==0){ # si no hi ha adjust factors és el de sempre


			aux.weights=weights.vec/tapply(weights.vec,g,sum)[as.character(g)]
			aux=by(sweep(2^t(em)[,names(perc.sel)],1,aux.weights,FUN="*"),INDICES=g,FUN=colSums)
			fc=aux[[1]]/aux[[2]]
			fc[fc<1]=1/(fc[fc<1])
			names(fc)=names(perc.sel)

			# sense tenir en compte weights
			#fc=apply(2^t(em)[,names(perc.sel)],2,tapply,g,mean)
			#fc=apply(fc,2,max)/apply(fc,2,min)

		} else { # faré com al gràfic. Ho sigui, faré FC dels residus!!!
			# 2^resid = quantes vegades és més gran el valor observat del predit en escala orignal (2^)
			# per tant el FC aquí serà quantes vegades és més gran en grup1 la mitjana de 'quant més gran obs que predit' que el grup2

			matrix.residuals=function(x,Y){
				lm=lm(x~.,data=Y,weights=weights.vec)
				return(resid(lm))
			}

			residuals=apply(t(em)[,names(perc.sel)],2,matrix.residuals,Y=adjust.factors[colnames(em),])

			aux.weights=weights.vec/tapply(weights.vec,g,sum)[as.character(g)]
			aux=by(sweep(2^residuals,1,aux.weights,FUN="*"),INDICES=g,FUN=colSums)
			fc=aux[[1]]/aux[[2]]
			fc[fc<1]=1/(fc[fc<1])
			names(fc)=names(perc.sel)

			# sense tenir en compte weights
			#fc=apply(2^residuals,2,tapply,g,mean)
			#fc=apply(fc,2,max)/apply(fc,2,min)

		}
	}

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100*(quant+7)/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = (quant+7)/tant, message = NULL, detail = paste(round(100*(quant+7)/tant,2),"%"),session = getDefaultReactiveDomain())

	# indespensabilitat # valor marginal de error quan la variable hi és o no. Dic marginal ja que s'ajunta hi hagi les combinacions que hi hagi d'altres vars
	# hauria de ser hi és o no al conjunt, no al model, per tant no s'ha de fer amb el coefs
	f=function(x,y){
		tapply(x,factor(y,levels=c("0","1")),mean)
	}
	#indespensabilitat=apply(!is.na(coefs),2,f,x=sortida.lasso$cv.error)

	indespensabilitat=apply(sortida.lasso$vars[,ordre],2,f,x=sortida.lasso$cv.error)

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100*(quant+8)/tant,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = (quant+8)/tant, message = NULL, detail = paste(round(100*(quant+8)/tant,2),"%"),session = getDefaultReactiveDomain())

	# relacio gens ($$$ millorar velocitat d'aixó, paralelitzar-ho)
	odds=comparison.categoric.mod(x=!is.na(coefs),vars=sortida.lasso$vars)

	col1="chocolate3"  # "saddlebrown"
	col2="seagreen4"  # "seagreen4"

	by.cor=0.1 # al paper xose era 3 (que equival de 0 a 20)

	
	odds2=odds
	diag(odds2)=NA
	odds2=as.data.frame(odds2)

	for (i in 1:ncol(odds)){

		temp=cut(odds2[,i],breaks=c(-Inf,seq(0+by.cor,1-by.cor,by=by.cor),Inf)) # si by.cor=0.1, breaks=seq(0,1,by=0.1)
		#temp=cut(odds2[,i],breaks=c(-Inf,seq(0,2,by=0.2),Inf)) # per way=3
		
		levels(temp)=colorRampPalette(c(col1,"white",col2))(nlevels(temp))
		odds2[,i]=as.character(temp)

	}

	odds2[is.na(odds2)]="grey50" # passa quan sempre s'agafa var X i no Y quan hi ha var X i Y (cas MME i MME2, quan es mira % de vegades que s'ha agafat
					# MME quan s'agafa MME2, ja que el MME2 no s'agafa mai si s'ha agafat MME)
					# es podria canviar per 0, però de moment ho deixo així

	#incProgress(amount = 1/tant, message = NULL, detail = paste(round(100,2),"%"),session = getDefaultReactiveDomain())
	setProgress(value = 1, message = NULL, detail = paste(round(100,2),"%"),session = getDefaultReactiveDomain())

	l=list(coefs,length.bios,perc.sel,indespensabilitat,fc,ba,odds,odds2,llargada,sortida.lasso$cv.error,sortida.lasso$cv.sens,sortida.lasso$cv.spec,sortida.lasso$vars,sortida.lasso$inds,g)
	names(l)=c("coefs","length.bios","perc.sel","indespensabilitat","fc","ba","odds","odds2","llargada","lasso.errors","lasso.sens","lasso.spec","vars","inds","g")
	return(l)

}

surv.ind=function(x,y,param=1){
	fit=coxph(y~x)
	if (param==1){return(as.numeric(fit$coef))}
	if (param==2){
		return(as.numeric(summary(fit)$waldtest["pvalue"]))
		#return(as.numeric(summary(fit)$concordance[1])) # c-index, millor o pitjor que p-val?
	}
}



grafic=function(calculs,quant,alternative.plot,ordre,type.bivariate){

	if (ordre!="perctimes"){

		if (ordre=="uniqueness"){

			ind=calculs$indespensabilitat[1,]-calculs$indespensabilitat[2,] # el % va de 0 a 100, aixó ha d'anar de mínim a màxim
			ind[is.na(ind)]=0

			o=order(ind,decreasing=T)
		}

		if (ordre=="uni.error"){
			o=order(calculs$ba,decreasing=T)
		}

		if (ordre=="FC"){
			o=order(calculs$fc,decreasing=T)
		}

		calculs$coefs=calculs$coefs[,o]
		calculs$perc.sel=calculs$perc.sel[o]
		calculs$indespensabilitat=calculs$indespensabilitat[,o]
		calculs$fc=calculs$fc[o]
		calculs$ba=calculs$ba[o]
		calculs$odds=calculs$odds[o,o]
		calculs$odds2=calculs$odds2[o,o]
		calculs$llargada=calculs$llargada[,o]

	}

	ind=calculs$indespensabilitat[1,]-calculs$indespensabilitat[2,]

	if (quant > ncol(calculs$coefs)) {quant = ncol(calculs$coefs)} # evitar error si es re-executa

	coefs=calculs$coefs[,seq_len(quant)]
	length.bios=calculs$length.bios
	perc.sel=calculs$perc.sel[seq_len(quant)]
	indespensabilitat=calculs$indespensabilitat[,seq_len(quant)]
	fc=calculs$fc[seq_len(quant)]
	ba=calculs$ba[seq_len(quant)]
	odds=calculs$odds[seq_len(quant),seq_len(quant)]

	per.distribucio.odds=calculs$odds
	#print(per.distribucio.odds)
	#write.table(per.distribucio.odds,"/home/guillem/aux.txt")

	odds2=calculs$odds2[seq_len(quant),seq_len(quant)]
	llargada=calculs$llargada[,seq_len(quant)]

	dalt=3
	dreta=8
	sota=8
	esquerra=5
	entre.grafs=max(max(nchar(colnames(coefs))),8)
	cex.axis=1.5
	cex.mtext=1.1

	grey="grey90"

	#### GRAFIC


	###### creació del layout

	layout(matrix(c(6,1,4,3,2,5),ncol=3),widths=c(2,3,7),heights=c(3,7))

	###### layout 1: histograma llargada biomarkers (fem de 0 a 20) o error dels B lassos
			# Nota:  posar algo com limma ranking vs ...
	if (alternative.plot=="length.biomarkers"){
		par(mar=c(sota,esquerra,entre.grafs,0))

		b=barplot(rev(length.bios*100),horiz=T,main="",
			xaxt="n",yaxt="n",ylab="",xlab="",cex.main=cex.axis,names.arg="",space=1)

		mtext("Relative Frequency",side=3,line=5.5,font=2,cex=cex.mtext)
		mtext("Biomarker Length",side=2,line=3,font=2,cex=cex.mtext)

		#plot(length.bios*100,type="h",xlim=c(0,20),main="Relative Frequency Biomarker Length",
		#	xaxt="n",yaxt="n",ylab="",xlab="",cex.main=cex.axis)
		rect(-1000,-1000,1000,1000,col=grey)
		by=c(1,2,3,5,10)[1+sum(c(10,20,30,40,105)<ceiling(max(length.bios*100)))]
		abline(v=seq(0,100,by=by),col="white",lwd=2)

		#points(length.bios*100,type="h")
		barplot(rev(length.bios*100),horiz=T,add=T,names.arg="",space=1,xaxt="n",yaxt="n",ylab="",xlab="",col="#1F4A8C")

		axis(3,at=seq(0,ceiling(max(length.bios*100)),by=by),paste(seq(0,ceiling(max(length.bios*100)),by=by),"%"),
			las=2,font=2,cex.axis=cex.axis)
		axis(2,at=b[,1],names(rev(length.bios)),cex.axis=cex.axis-0.2,font=2,las=2)
		#axis(2,at=seq(1,19,by=2),cex.axis=cex.axis-0.2,font=2)

	}

	if (alternative.plot=="lasso.errors"){
		par(mar=c(sota,esquerra,entre.grafs,0))

		plot(1:10,ylim=c(1,0),xlim=c(0.5,3.5),axes=F,xlab="",ylab="")
		rect(-1000,-1000,1000,1000,col=grey)

		by=0.1

		abline(h=seq(0,1,by=by),col="white",lwd=2)

		#t=prop.table(table(calculs$lasso.errors))
		#t=t/1.05*max(t)
		#segments(x0=1-(prop.table(t)/2),x1=1+(prop.table(t)/2),y0=as.numeric(names(t)),lwd=3)

		# estaria bé modificar-ho una mica tipo boxplot. Si hi ha un outlier fer-ho separat
		# ficar una espècie de bandwich?

		n.samp=sum(calculs$inds[1,])
		possibles=round((0:n.samp)/n.samp,10)
		errors=round(calculs$lasso.errors,10)
		possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
		t=factor(errors,levels=possibles)
		#t=prop.table(table(calculs$lasso.errors))
		t=prop.table(table(t))
		t=t*0.9/max(t)			
		polygon(x=c(1-(t/2),rev(1+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
			col="darkgreen")

		n.samp=sum(calculs$inds[1,calculs$g==levels(calculs$g)[1]])
		possibles=round((0:n.samp)/n.samp,10)
		errors=round(calculs$lasso.sens,10)
		possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
		t=factor(errors,levels=possibles)
		t=prop.table(table(t))
		t=t*0.9/max(t)			
		polygon(x=c(2-(t/2),rev(2+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
			col="darkgreen")

		n.samp=sum(calculs$inds[1,calculs$g==levels(calculs$g)[2]])
		possibles=round((0:n.samp)/n.samp,10)
		errors=round(calculs$lasso.spec,10)
		possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
		t=factor(errors,levels=possibles)
		t=prop.table(table(t))
		t=t*0.9/max(t)			
		polygon(x=c(3-(t/2),rev(3+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
			col="darkgreen")



		#t=prop.table(table(calculs$lasso.sens))
		#t=t*0.9/max(t)			
		#polygon(x=c(2-(t/2),rev(2+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
		#	col="darkgreen")

		#t=prop.table(table(calculs$lasso.spec))
		#t=t*0.9/max(t)			
		#polygon(x=c(3-(t/2),rev(3+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
		#	col="darkgreen")

		axis(2,at=seq(0,1,by=by),paste(seq(0,100,by=by*100),"%"),font=2,cex.axis=cex.axis,las=2)
		axis(3,at=1:3,c("ACC","TPR","TNR"),las=2,font=2,cex.axis=cex.axis,tick=F)
		# TPR= True Positive Rate (sensitivity)
		# TNR= True Negative Rate (specificity)
		# ACC = accuracy (error rate)
		mtext("Lasso Performance",side=3,line=4.5,font=2,cex=cex.mtext)
	}



	####### layout 3: multivariate error rate (percentatge de vegades seleccionada la variable + indespensabilitat)

	par(mar=c(0,entre.grafs,dalt,dreta))
	plot(perc.sel,type="b",pch=19,xlim=c(1,ncol(coefs)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="")
	#title(main="Multivariate Power",font=2,cex.main=cex.axis)
	#mtext("Multivariate Power",side=3,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)
	rect(-1000,-1000,1000,1000,col=grey)
	where=seq(0,100,by=10)
	abline(h=where,col="white",lwd=2)
	axis(2,at=where,paste(where,"%"),las=2,font=2,cex.axis=cex.axis)

	#legend("topright",c("% times selected","individual error rate"),lty=1,lwd=2,col=c("#966FAD","#90BA6D"),
	#	text.font=2,cex=cex.axis-0.2,bg="white")

	lines(perc.sel,type="b",pch=19,col="#966FAD",lwd=3)

	#lines((1-ba)*100,type="b",pch=19,col="#90BA6D",lwd=2) # individual.power

	ind=indespensabilitat[1,]-indespensabilitat[2,] # el % va de 0 a 100, aixó ha d'anar de mínim a màxim
	ind[is.na(ind)]=0 # fiquem 0 als que són NA. És NA quan la variable sempre hi és o mai hi és
	resc.uniqueness=ind-min(ind,na.rm=T)
	resc.uniqueness=resc.uniqueness/max(resc.uniqueness)
	resc.uniqueness=resc.uniqueness*100

	lines(resc.uniqueness,type="b",pch=19,col="#90BA6D",lwd=3)

	axis(4,at=where,round(seq(min(ind,na.rm=T),max(ind,na.rm=T),length.out=length(where))*100,2),las=2,font=2,cex.axis=cex.axis)

	mtext("% Possible Times Selected",side=2,line=5,font=2,cex=cex.mtext)
	if (class(calculs$g)=="factor") mtext("Uniqueness (%)",side=4,line=5,font=2,cex=cex.mtext)
	if (class(calculs$g)=="Surv") mtext("Uniqueness",side=4,line=5,font=2,cex=cex.mtext)


	######## layout 4: individual power (fold change + individual error rate) # Nota: estaria bé que les dos línies com més alt millor
	# també indicar si Fold Change és cap a una banda o cap a l'altra

	par(mar=c(sota,dalt,entre.grafs,0))

	if (class(calculs$g)=="factor"){
		plot(log(fc),1:length(fc),type="b",pch=19,ylim=c(1,ncol(coefs)),
			xlim=c(0,max(log(fc)))+c(-0.02,0.02)*max(log(fc)),yaxt="n",xaxt="n",xaxs="i",xlab="",ylab="")
		rect(-1000,-1000,1000,1000,col=grey)
		by=max(log(fc))/9

		where=seq(0,max(log(fc)),by=by)
		abline(v=where,col="white",lwd=2)
		axis(1,at=where,round(exp(where),2),las=2,font=2,cex.axis=cex.axis)

		lines(log(fc),rev(1:length(fc)),type="b",pch=19,col="#C3703C",lwd=3) # #C3703C o #9F5845

		resc.unipower=1-ba # log(fold change) va de 0 a max(log(fc)). Aixó ha d'anar de 100 a min(aixo)
		resc.unipower=(resc.unipower/max(resc.unipower))*max(log(fc))
		resc.unipower=abs(resc.unipower-max(resc.unipower))
		axis(3,at=where,paste(round(100*seq(max(1-ba),0,length.out=length(where)),1),"%"),las=2,font=2,cex.axis=cex.axis)

		lines(resc.unipower,rev(1:length(resc.unipower)),type="b",pch=19,col="#7E9A99",lwd=3) # individual.power

		if (ncol(calculs$adjust.factors)==0){
			mtext("Fold Change",side=1,line=4.5,font=2,cex=cex.mtext)
		} else {
			mtext("Residual Fold Change",side=1,line=4.5,font=2,cex=cex.mtext)
		}

		mtext("LOOCV(BA)",side=3,line=5.5,font=2,cex=cex.mtext)
		#mtext("Univariate Error Rate (Balanced Accuracy)",side=3,line=5.5,font=2,cex=cex.mtext)
		#mtext("Univariate Power",side=2,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)
	}

	if (class(calculs$g)=="Surv"){

		if (1==1){ # versió del HR en un costat i p-value al altra
			plot(-log(ba,10),1:length(ba),type="b",pch=19,ylim=c(1,ncol(coefs)),
				xlim=c(0,max(-log(ba,10)))+c(-0.02,0.02)*max(-log(ba,10)),yaxt="n",xaxt="n",xaxs="i",xlab="",ylab="")
			rect(-1000,-1000,1000,1000,col=grey)
			by=max(-log(ba,10))/9

			where=seq(0,max(-log(ba,10)),by=by)
			abline(v=where,col="white",lwd=2)
			axis(1,at=where,round(where,2),las=2,font=2,cex.axis=cex.axis)

			lines(-log(ba,10),rev(1:length(ba)),type="b",pch=19,col="#C3703C",lwd=3) # #C3703C o #9F5845

			# fc ha d'anar de min(-log(ba,10)) a max(-log(ba,10))

			resc=fc-min(fc,na.rm=T)
			resc=resc*(max(-log(ba,10)) - min(-log(ba,10)))/diff(range(resc))
			resc=resc+min(-log(ba,10))

			axis(3,at=where,round(seq(min(fc),max(fc),length.out=length(where)),2),las=2,font=2,cex.axis=cex.axis)

			lines(resc,1:length(resc),type="b",pch=19,col="#7E9A99",lwd=3) # individual.power

			mtext("P-Value (10^-x)",side=1,line=4.5,font=2,cex=cex.mtext)
			mtext("Hazard Ratio",side=3,line=5.5,font=2,cex=cex.mtext)
			#mtext("Univariate Power",side=2,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)

		} else { # versió IC95% HR amb segments
			
			# calcul dels IC95% HR


		}
	}

	######## layout 5: univarate error rate vs perc.sel

	if (class(calculs$g)=="factor"){

		multivariate=switch(type.bivariate,"1"=perc.sel,"2"=perc.sel,"3"=resc.uniqueness,"4"=resc.uniqueness)
		univariate=switch(type.bivariate,"1"=resc.unipower,"3"=resc.unipower,"2"=log(fc),"4"=log(fc))
	
		par(mar=c(0,dalt,dalt,0))
		plot(univariate,multivariate,pch=19,xlim=c(0,max(log(fc)))+c(-0.02,0.02)*max(log(fc)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="",xaxs="i")

		rect(-1000,-1000,1000,1000,col=grey)
		by=max(log(fc))/9

		where=seq(0,max(log(fc)),by=by)
		abline(v=where,col="white",lwd=2)
		where=seq(0,100,by=10)
		abline(h=where,col="white",lwd=2)

		points(univariate,multivariate,pch=19,col=c("#8181B4","#BD4F64")[1+(colSums(coefs,na.rm=T)>0)])
		text(univariate,multivariate,names(multivariate),cex=1.25,font=2,srt=45,adj=c(1.1,0.5))

	}

	if (class(calculs$g)=="Surv"){

		par(mar=c(0,dalt,dalt,0))
		plot(-log(ba,10),perc.sel,pch=19,xlim=c(0,max(-log(ba,10)))+c(-0.02,0.02)*max(-log(ba,10)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="",xaxs="i")

		rect(-1000,-1000,1000,1000,col=grey)
		by=max(-log(ba,10))/9

		where=seq(0,max(-log(ba)),by=by)
		abline(v=where,col="white",lwd=2)
		where=seq(0,100,by=10)
		abline(h=where,col="white",lwd=2)

		points(-log(ba,10),perc.sel,pch=19,col=c("#8181B4","#BD4F64")[1+(colSums(coefs,na.rm=T)>0)])
		text(-log(ba,10),perc.sel,names(perc.sel),cex=1,font=2,srt=45,adj=c(1.1,0.5))
	}


	######## layout 6: correlació variables

	par(mar=c(sota,entre.grafs,entre.grafs,dreta))
	plot(1:ncol(odds),type="n",xlim=c(1,ncol(odds)),ylim=c(1,ncol(odds)),frame.plot=F,axes=F,xlab="",ylab="")

	for (i in 1:ncol(odds)){
		cls=odds2[,i]
		rect(xleft=i-0.5,xright=i+0.5,ybottom=c(1:ncol(odds))-0.5,ytop=c(1:ncol(odds))+0.5,border="white",col=rev(cls),lwd=3)
	}

	#axis(2,at=seq(ncol(odds),1,by=-1),colnames(odds),tick=F,las=2,cex.axis=cex.axis,font.axis=4,padj=0.5)
	#axis(3,at=1:ncol(odds),colnames(odds),tick=F,las=2,cex.axis=cex.axis,font.axis=4,hadj=0.5,line=entre.grafs/2)

	cex.var.name=cex.mtext

	if (quant>30){
		cex.var.name=cex.var.name*(1-3*(quant-30)/100)
	}

	mtext(side=3,at=1:ncol(odds),           text=colnames(odds),las=2,cex=cex.var.name,font=2,adj=0.5,line=entre.grafs/2)
	mtext(side=2,at=seq(ncol(odds),1,by=-1),text=colnames(odds),las=2,cex=cex.var.name,font=2,adj=0.5,line=entre.grafs/2)
	# font=2, ja que no sempre seran gens (que serien font=4)

	####### layout 8: llegenda

	if (class(calculs$g)=="factor"){
		#par(mar=c(0,dalt,dalt,0))
		par(mar=c(0,0,0,0))

		plot(1:10,type="n",axes=F,frame.plot=F,ylim=c(-1,9),xlim=c(0.5,11.5))
		legend(x=2,y=8,c("selected","uniqueness","FC","error"),lwd=3,pch=19,lty=1,col=c("#966FAD","#90BA6D","#C3703C","#7E9A99"),cex=cex.axis)

		# horizontal

		col1="chocolate3"  # "saddlebrown"
		col2="seagreen4"  # "seagreen4"
		col.breaks=colorRampPalette(c(col1,"white",col2))(10)

		rect(xleft=1:10,xright=2:11,ybottom=0,ytop=1,border=col.breaks,col=col.breaks)
		text(x=1:11,y=-0.5,seq(0,100,by=10),cex=cex.axis,font=2,srt=90,adj=c(1,0.5))
		#text(x=xleg+1,y=yleg+(2:(length(breaks)-1)),"-")
		text(x=6,y=1.5,"% possible times\ngene column selected\nwhen row selected",srt=0,adj=c(0.5,0),font=2,cex=cex.axis)

		segments(x0=1,x1=11,y0=-0.1,y1=-0.1)
		segments(x0=1:11,x1=1:11,y0=-0.2,y1=-0.1)
	}

	if (class(calculs$g)=="Surv"){
		#par(mar=c(0,dalt,dalt,0))
		par(mar=c(0,0,0,0))

		plot(1:10,type="n",axes=F,frame.plot=F,ylim=c(-1,9),xlim=c(0.5,11.5))
		legend(x=2,y=8,c("selected","uniqueness","p-value","hazard ratio"),lwd=3,pch=19,lty=1,col=c("#966FAD","#90BA6D","#C3703C","#7E9A99"),cex=cex.axis)

		# horizontal

		col1="chocolate3"  # "saddlebrown"
		col2="seagreen4"  # "seagreen4"
		col.breaks=colorRampPalette(c(col1,"white",col2))(10)

		rect(xleft=1:10,xright=2:11,ybottom=0,ytop=1,border=col.breaks,col=col.breaks)
		text(x=1:11,y=-0.5,seq(0,100,by=10),cex=cex.axis,font=2,srt=90,adj=c(1,0.5))
		#text(x=xleg+1,y=yleg+(2:(length(breaks)-1)),"-")
		text(x=6,y=1.5,"% possible times\ngene column selected\nwhen row selected",srt=0,adj=c(0.5,0),font=2,cex=cex.axis)

		segments(x0=1,x1=11,y0=-0.1,y1=-0.1)
		segments(x0=1:11,x1=1:11,y0=-0.2,y1=-0.1)
	}

}


compt.biomarkers=function(i,x){
	return(paste(colnames(x)[x[i,]],collapse=","))
}

compt.biomarker.possible=function(select,vars){
	if (select!=""){
		return(sum(apply(subset(vars,select=do.call(c,strsplit(select,","))),1,prod)))
	} else{
		return(NA)
	}
	#return(sum(apply(vars[,strsplit(select,",")],1,prod))) # numero de vegades total que podia sortir el biomarkador select en les iteracions
}



# calcula el T^2 Hotelling per un biomarker
calc.dist=function(x1,x2,distance=c("T2.equal","T2.unequal","margin","deviance")){

	p<-dim(x1)[2]
	n1<-dim(x1)[1]
	n2<-dim(x2)[1] # això ho podriem passar des de fora? així no perdem el temps amb calcular el dim()

	m1<-apply(x1,2,mean)
	s1<-var(x1)
	m2<-apply(x2,2,mean)
	s2<-var(x2)

	if (distance=="T2.equal"){
		s<-((n1-1)*s1+(n2-1)*s2)/(n1+n2-2) 
		T2<-((n1*n2)/(n1+n2))*t(m1-m2)%*%solve(s)%*%(m1-m2)
		return(T2)
	}

	if (distance=="T2.unequal"){
		T2<-t(m1-m2)%*%solve(s1/n1+s2/n2)%*%(m1-m2) # Wald test?
		return(T2)
	}

}

# facilitar T^2 per un biomarker

calc.T2=function(biomarker,g,em){
	#biomarker=make.names(strsplit(biomarker,split=",")[[1]])
	biomarker=strsplit(biomarker,split=",")[[1]]

	if(length(biomarker)>0){
		x1=subset(t(em),select=biomarker,subset=g==levels(g)[1])
		x2=subset(t(em),select=biomarker,subset=g==levels(g)[2])

		#T2=calc.dist(x1,x2,distance="T2.equal")
		T2=try(calc.dist(x1,x2,distance="T2.equal"))
		if (!is.numeric(T2)){T2=NA}

	} else {
		T2=NA
	}
	return(T2)
}

# T-Test per differential gene expression més ràpid
fast.t.test.mat=function(x,y,var.equal=FALSE){
	x1=x[y==levels(y)[1],]
	x2=x[y==levels(y)[2],]

	m1=colMeans(x1)
	m2=colMeans(x2)

	n1=nrow(x1)
	n2=nrow(x2)

	s1=colSums(sweep(x1,MARGIN=2,STATS=m1)^2)/(n1-1)
	s2=colSums(sweep(x2,MARGIN=2,STATS=m2)^2)/(n1-1)

	num=m1-m2

	if (var.equal){
		den=sqrt(((n1-1)*s1+(n2-1)*s2)/(n1+n2-2))*sqrt((1/n1)+(1/n2))
		df=n1+n2-2
	} else {
		den=sqrt((s1/n1)+(s2/n2))
		df=(((s1/n1)+(s2/n2))^2)/((((s1/n1)^2)/(n1-1))+(((s2/n2)^2)/(n2-1)))
	}

	t=num/den
	p=pt(-abs(t),df)*2

	sort=data.frame(estimate=num,t=t,p.val=p)
	return(sort)
}

# mann-whitney test més ràpid (millorar rank per filar més prim)
fast.wilcox=function(x,y){

	ranks=apply(x,2,rank)
	ranks1=ranks[y==levels(y)[1],]
	ranks2=ranks[y==levels(y)[2],]

	n1=nrow(ranks1)
	n2=nrow(ranks2)

	R1=colSums(ranks1)
	R2=colSums(ranks2)

	U1=R1-n1*(n1+1)/2
	U2=R2-n2*(n2+1)/2

	U=apply(cbind(U1,U2),1,min) # wilcox.test torna el U1 sempre

	return(2*pwilcox(U,n1,n2))

}







### versió antiga de gràfic, que inclou els exp(betas) i num compis


grafic.antic=function(calculs,quant,alternative.plot,ordre){

	if (ordre!="perctimes"){

		if (ordre=="uniqueness"){

			ind=calculs$indespensabilitat[1,]-calculs$indespensabilitat[2,] # el % va de 0 a 100, aixó ha d'anar de mínim a màxim
			ind[is.na(ind)]=0

			o=order(ind,decreasing=T)
		}

		if (ordre=="uni.error"){
			o=order(calculs$ba,decreasing=T)
		}

		if (ordre=="FC"){
			o=order(calculs$fc,decreasing=T)
		}

		calculs$coefs=calculs$coefs[,o]
		calculs$perc.sel=calculs$perc.sel[o]
		calculs$indespensabilitat=calculs$indespensabilitat[,o]
		calculs$fc=calculs$fc[o]
		calculs$ba=calculs$ba[o]
		calculs$odds=calculs$odds[o,o]
		calculs$odds2=calculs$odds2[o,o]
		calculs$llargada=calculs$llargada[,o]

	}
	coefs=calculs$coefs[,seq_len(quant)]
	length.bios=calculs$length.bios
	perc.sel=calculs$perc.sel[seq_len(quant)]
	indespensabilitat=calculs$indespensabilitat[,seq_len(quant)]
	fc=calculs$fc[seq_len(quant)]
	ba=calculs$ba[seq_len(quant)]
	odds=calculs$odds[seq_len(quant),seq_len(quant)]

	per.distribucio.odds=calculs$odds
	#print(per.distribucio.odds)
	#write.table(per.distribucio.odds,"/home/guillem/aux.txt")

	odds2=calculs$odds2[seq_len(quant),seq_len(quant)]
	llargada=calculs$llargada[,seq_len(quant)]

	dalt=3
	dreta=8
	sota=8
	esquerra=5
	entre.grafs=max(max(nchar(colnames(coefs))),8)
	cex.axis=1.5
	cex.mtext=1.1

	grey="grey90"

	#### GRAFIC


	###### creació del layout
	#x11(width=13,height=11.8)
	#layout(matrix(c(6,4,3,1,2,5),ncol=2),widths=c(3,6),heights=c(3,3,6))
	#layout.show(6)
	layout(matrix(c(1,8,7,1,5,4,2,3,6),ncol=3),widths=c(2,3,7),heights=c(2,3,7))
	layout.show(8)

	###### layout 1: histograma llargada biomarkers (fem de 0 a 20)

	par(mar=c(0,esquerra,dalt,0))

	plot(length.bios*100,type="h",xlim=c(0,20),main="Relative Frequency Biomarker Length",xaxt="n",yaxt="n",ylab="",xlab="",cex.main=cex.axis)
	rect(-1000,-1000,1000,1000,col=grey)
	by=c(1,2,3,5,10)[1+sum(c(10,20,30,60,105)<ceiling(max(length.bios*100)))]
	abline(h=seq(0,100,by=by),col="white",lwd=2)
	points(length.bios*100,type="h")
	axis(2,at=seq(0,ceiling(max(length.bios*100)),by=by),paste(seq(0,ceiling(max(length.bios*100)),by=by),"%"),las=2,font=2,cex.axis=cex.axis)
	axis(1,at=seq(0,20,by=2),c(seq(0,19,by=2),"20+"),cex.axis=cex.axis-0.2,font=2)
	axis(1,at=seq(1,19,by=2),cex.axis=cex.axis-0.2,font=2)
	#axis(1,at=20,"20+",cex.axis=cex.axis,font=2)

	###### layout 2: boxplot paràmetre (raw or exp?)

	par(mar=c(0,entre.grafs,dalt,dreta))
	boxplot(abs(coefs),horizontal = FALSE,xlim=c(1,ncol(coefs)),xaxt="n",yaxt="n",xlab="",ylab="")
	rect(-1000,-1000,1000,1000,col=grey)
	where=seq(0,max(abs(coefs),na.rm=T),length.out=7)
	abline(h=where,col="white",lwd=2)
	#abline(v=0,col="#BD4F64",lwd=2)
	axis(2,at=where,round(exp(where),1),las=2,font=2,cex.axis=cex.axis)
	#axis(3,at=seq(0,10,by=1),round(exp(seq(0,10,by=1)),1),las=1,font=2,cex.axis=cex.axis-0.2)
	#axis(3,at=seq(0.5,10,by=1),round(exp(seq(0.5,10,by=1)),1),las=1,font=2,cex.axis=cex.axis-0.2)

	boxplot(abs(coefs),horizontal = FALSE,add=T,at=seq(1,ncol(coefs),by=1),
		col=c("#8181B4","#BD4F64")[1+(colSums(coefs,na.rm=T)>0)],xaxt="n",yaxt="n",xlab="",ylab="")

	mtext(expression(bold(e)^bold("|"*hat(b)*"|")),side=2,line=4,font=2,cex=cex.mtext,las=1)
	# afegir + Others i + CLL a damunt del graf?

	####### layout 3: multivariate error rate (percentatge de vegades seleccionada la variable + indespensabilitat)

	par(mar=c(0,entre.grafs,dalt,dreta))
	plot(perc.sel,type="b",pch=19,xlim=c(1,ncol(coefs)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="")
	#title(main="Multivariate Power",font=2,cex.main=cex.axis)
	#mtext("Multivariate Power",side=3,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)
	rect(-1000,-1000,1000,1000,col=grey)
	where=seq(0,100,by=10)
	abline(h=where,col="white",lwd=2)
	axis(2,at=where,paste(where,"%"),las=2,font=2,cex.axis=cex.axis)

	#legend("topright",c("% times selected","individual error rate"),lty=1,lwd=2,col=c("#966FAD","#90BA6D"),
	#	text.font=2,cex=cex.axis-0.2,bg="white")

	lines(perc.sel,type="b",pch=19,col="#966FAD",lwd=3)

	#lines((1-ba)*100,type="b",pch=19,col="#90BA6D",lwd=2) # individual.power

	ind=indespensabilitat[1,]-indespensabilitat[2,] # el % va de 0 a 100, aixó ha d'anar de mínim a màxim
	ind[is.na(ind)]=0 # fiquem 0 als que són NA. És NA quan la variable sempre hi és o mai hi és
	resc=ind-min(ind,na.rm=T)
	resc=resc/max(resc)
	resc=resc*100

	lines(resc,type="b",pch=19,col="#90BA6D",lwd=3)

	axis(4,at=where,round(seq(min(ind,na.rm=T),max(ind,na.rm=T),length.out=length(where))*100,2),las=2,font=2,cex.axis=cex.axis)

	mtext("% Possible Times Selected",side=2,line=5,font=2,cex=cex.mtext)
	if (class(calculs$g)=="factor") mtext("Uniqueness (%)",side=4,line=5,font=2,cex=cex.mtext)
	if (class(calculs$g)=="Surv") mtext("Uniqueness",side=4,line=5,font=2,cex=cex.mtext)


	######## layout 4: individual power (fold change + individual error rate) # Nota: estaria bé que les dos línies com més alt millor
	# també indicar si Fold Change és cap a una banda o cap a l'altra

	par(mar=c(sota,dalt,entre.grafs,0))

	if (class(calculs$g)=="factor"){
		plot(log(fc),1:length(fc),type="b",pch=19,ylim=c(1,ncol(coefs)),
			xlim=c(0,max(log(fc)))+c(-0.02,0.02)*max(log(fc)),yaxt="n",xaxt="n",xaxs="i",xlab="",ylab="")
		rect(-1000,-1000,1000,1000,col=grey)
		by=max(log(fc))/9

		where=seq(0,max(log(fc)),by=by)
		abline(v=where,col="white",lwd=2)
		axis(1,at=where,round(exp(where),2),las=2,font=2,cex.axis=cex.axis)

		lines(log(fc),rev(1:length(fc)),type="b",pch=19,col="#C3703C",lwd=3) # #C3703C o #9F5845

		resc=1-ba # log(fold change) va de 0 a max(log(fc)). Aixó ha d'anar de 100 a min(aixo)
		resc=(resc/max(resc))*max(log(fc))
		resc=abs(resc-max(resc))
		axis(3,at=where,paste(round(100*seq(max(1-ba),0,length.out=length(where)),1),"%"),las=2,font=2,cex.axis=cex.axis)

		lines(resc,rev(1:length(resc)),type="b",pch=19,col="#7E9A99",lwd=3) # individual.power

		if (ncol(calculs$adjust.factors)==0){
			mtext("Fold Change",side=1,line=4.5,font=2,cex=cex.mtext)
		} else {
			mtext("Residual Fold Change",side=1,line=4.5,font=2,cex=cex.mtext)
		}

		mtext("Univariate Accuracy",side=3,line=5.5,font=2,cex=cex.mtext)
		#mtext("Univariate Error Rate (Balanced Accuracy)",side=3,line=5.5,font=2,cex=cex.mtext)
		#mtext("Univariate Power",side=2,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)
	}

	if (class(calculs$g)=="Surv"){

		if (1==1){ # versió del HR en un costat i p-value al altra
			plot(-log(ba,10),1:length(ba),type="b",pch=19,ylim=c(1,ncol(coefs)),
				xlim=c(0,max(-log(ba,10)))+c(-0.02,0.02)*max(-log(ba,10)),yaxt="n",xaxt="n",xaxs="i",xlab="",ylab="")
			rect(-1000,-1000,1000,1000,col=grey)
			by=max(-log(ba,10))/9

			where=seq(0,max(-log(ba,10)),by=by)
			abline(v=where,col="white",lwd=2)
			axis(1,at=where,round(where,2),las=2,font=2,cex.axis=cex.axis)

			lines(-log(ba,10),rev(1:length(ba)),type="b",pch=19,col="#C3703C",lwd=3) # #C3703C o #9F5845

			# fc ha d'anar de min(-log(ba,10)) a max(-log(ba,10))

			resc=fc-min(fc,na.rm=T)
			resc=resc*(max(-log(ba,10)) - min(-log(ba,10)))/diff(range(resc))
			resc=resc+min(-log(ba,10))

			axis(3,at=where,round(seq(min(fc),max(fc),length.out=length(where)),2),las=2,font=2,cex.axis=cex.axis)

			lines(resc,1:length(resc),type="b",pch=19,col="#7E9A99",lwd=3) # individual.power

			mtext("P-Value (10^-x)",side=1,line=4.5,font=2,cex=cex.mtext)
			mtext("Hazard Ratio",side=3,line=5.5,font=2,cex=cex.mtext)
			#mtext("Univariate Power",side=2,font=2,cex=cex.mtext,line=dalt/2,padj=0.5)

		} else { # versió IC95% HR amb segments
			
			# calcul dels IC95% HR


		}
	}

	######## layout 5: univarate error rate vs perc.sel

	if (class(calculs$g)=="factor"){
		# gràfic de vegades que es selcciona una variable vs mitja del coeficient
		# Nota: el problema amb aquest gràfic es que el paràmetre es poc comparable
		# d'una iteracio a una altra si hi ha variables diferents implicades
		# Igualment dona una idea de l'efecte de la variable, si sempre és molt alt és que la regularització no el mata
		# Nota: en comptes d'agafar el paràmetre tal qual, agafar el paràmetre + lambda? Potser és efecte lambda...
		# Possible solució: efecte de randomitzar la variable?

		par(mar=c(0,dalt,dalt,0))
		plot(resc,perc.sel,pch=19,xlim=c(0,max(log(fc)))+c(-0.02,0.02)*max(log(fc)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="",xaxs="i")

		rect(-1000,-1000,1000,1000,col=grey)
		by=max(log(fc))/9

		where=seq(0,max(log(fc)),by=by)
		abline(v=where,col="white",lwd=2)
		where=seq(0,100,by=10)
		abline(h=where,col="white",lwd=2)

		points(resc,perc.sel,pch=19,col=c("#8181B4","#BD4F64")[1+(colSums(coefs,na.rm=T)>0)])
		text(resc,perc.sel,names(perc.sel),cex=1,font=2,srt=45,adj=c(1.1,0.5))

	}

	if (class(calculs$g)=="Surv"){

		par(mar=c(0,dalt,dalt,0))
		plot(-log(ba,10),perc.sel,pch=19,xlim=c(0,max(-log(ba,10)))+c(-0.02,0.02)*max(-log(ba,10)),ylim=c(-2,102),yaxt="n",xaxt="n",yaxs="i",xlab="",ylab="",xaxs="i")

		rect(-1000,-1000,1000,1000,col=grey)
		by=max(-log(ba,10))/9

		where=seq(0,max(-log(ba)),by=by)
		abline(v=where,col="white",lwd=2)
		where=seq(0,100,by=10)
		abline(h=where,col="white",lwd=2)

		points(-log(ba,10),perc.sel,pch=19,col=c("#8181B4","#BD4F64")[1+(colSums(coefs,na.rm=T)>0)])
		text(-log(ba,10),perc.sel,names(perc.sel),cex=1,font=2,srt=45,adj=c(1.1,0.5))
	}


	######## layout 6: correlació variables

	par(mar=c(sota,entre.grafs,entre.grafs,dreta))
	plot(1:ncol(odds),type="n",xlim=c(1,ncol(odds)),ylim=c(1,ncol(odds)),frame.plot=F,axes=F,xlab="",ylab="")

	for (i in 1:ncol(odds)){
		cls=odds2[,i]
		rect(xleft=i-0.5,xright=i+0.5,ybottom=c(1:ncol(odds))-0.5,ytop=c(1:ncol(odds))+0.5,border="white",col=rev(cls),lwd=3)
	}

	#axis(2,at=seq(ncol(odds),1,by=-1),colnames(odds),tick=F,las=2,cex.axis=cex.axis,font.axis=4,padj=0.5)
	#axis(3,at=1:ncol(odds),colnames(odds),tick=F,las=2,cex.axis=cex.axis,font.axis=4,hadj=0.5,line=entre.grafs/2)

	cex.var.name=cex.mtext

	if (quant>30){
		cex.var.name=cex.var.name*(1-3*(quant-30)/100)
	}

	mtext(side=3,at=1:ncol(odds),           text=colnames(odds),las=2,cex=cex.var.name,font=2,adj=0.5,line=entre.grafs/2)
	mtext(side=2,at=seq(ncol(odds),1,by=-1),text=colnames(odds),las=2,cex=cex.var.name,font=2,adj=0.5,line=entre.grafs/2)
	# font=2, ja que no sempre seran gens (que serien font=4)


	###### layout 7: llargada biomarker quan variable està implicada (sense comptar la pròpia variable)
	# Nota: aquest gràfic sembla de univariate power, però no és així. Treure'l i posar algo com limma ranking vs ...

	if (alternative.plot=="num.compis"){
		if (1==1){
			par(mar=c(sota,dalt,entre.grafs,0))
			boxplot(llargada,horizontal = TRUE,xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(1,ncol(coefs)))
			rect(-1000,-1000,1000,1000,col=grey)

			by=ceiling((diff(range(llargada,na.rm=T))+1)/9)
			abline(v=seq(0,200,by=by),col="white",lwd=2)
			axis(3,at=seq(0,200,by=by),las=2,font=2,cex.axis=cex.axis-0.2)
			#axis(3,at=seq(1,100,by=2),las=1,font=2,cex.axis=cex.axis)

			boxplot(llargada,horizontal = TRUE,add=T,at=seq(ncol(coefs),1,by=-1),
				col="#87C851",xaxt="n",yaxt="n",xlab="",ylab="")

			mtext("# of compis",side=3,line=entre.grafs/2,padj=0.5,font=2,cex=cex.mtext)
		}
	} else {
		if (1==0){
			par(mar=c(sota,esquerra,entre.grafs,0))

			ylim=ceiling(max(100*c(calculs$lasso.errors,calculs$lasso.sens,calculs$lasso.spec)))

			boxplot(100*calculs$lasso.errors,ylim=c(ylim,0),xlim=c(0.5,3.5),at=1,axes=F)
			rect(-1000,-1000,1000,1000,col=grey)

			by=c(1,2,3,5,10)[1+sum(c(10,20,30,60,105)<ylim)]
			axis(2,at=seq(0,ylim,by=by),paste(seq(0,ylim,by=by),"%"),font=2,cex.axis=cex.axis,las=2)
			axis(1,at=1:3,c("Error","Sens","Spec"),las=2,font=2,cex.axis=cex.axis,tick=F)

			abline(h=seq(0,ylim,by=by),col="white",lwd=2)

			boxplot(100*calculs$lasso.errors,add=T,at=1,axes=F,frame.plot=F,col="#506434", pars = list(boxwex = 1.5))
			boxplot(100*calculs$lasso.sens,add=T,at=2,axes=F,frame.plot=F,col="#506434", pars = list(boxwex = 1.5))
			boxplot(100*calculs$lasso.spec,add=T,at=3,axes=F,frame.plot=F,col="#506434", pars = list(boxwex = 1.5))

			mtext("Lasso Performance",side=1,line=4.5,font=2,cex=cex.mtext)
		}

		if (1==0){
			par(mar=c(sota,esquerra,entre.grafs,0))

			plot(1:10,ylim=c(1,0),xlim=c(0.5,3.5),axes=F,xlab="",ylab="")
			rect(-1000,-1000,1000,1000,col=grey)

			by=0.1

			abline(h=seq(0,1,by=by),col="white",lwd=2)

			t=prop.table(table(calculs$lasso.errors))
			t=t/1.05*max(t)
			segments(x0=1-(prop.table(t)/2),x1=1+(prop.table(t)/2),y0=as.numeric(names(t)),lwd=3)

			t=prop.table(table(calculs$lasso.sens))
			t=t/1.05*max(t)
			segments(x0=2-(prop.table(t)/2),x1=2+(prop.table(t)/2),y0=as.numeric(names(t)),lwd=3)

			t=prop.table(table(calculs$lasso.spec))
			t=t/1.05*max(t)
			segments(x0=3-(prop.table(t)/2),x1=3+(prop.table(t)/2),y0=as.numeric(names(t)),lwd=3)

			axis(2,at=seq(0,1,by=by),paste(seq(0,100,by=by*100),"%"),font=2,cex.axis=cex.axis,las=2)
			axis(3,at=1:3,c("ACC","TPR","TNR"),las=2,font=2,cex.axis=cex.axis,tick=F)
			# TPR= True Positive Rate (sensitivity)
			# TNR= True Negative Rate (specificity)
			# ACC = accuracy (error rate)
			mtext("Lasso Performance",side=3,line=4.5,font=2,cex=cex.mtext)
		}


		if (1==1){
			par(mar=c(sota,esquerra,entre.grafs,0))

			plot(1:10,ylim=c(1,0),xlim=c(0.5,3.5),axes=F,xlab="",ylab="")
			rect(-1000,-1000,1000,1000,col=grey)

			by=0.1

			abline(h=seq(0,1,by=by),col="white",lwd=2)

			#t=prop.table(table(calculs$lasso.errors))
			#t=t/1.05*max(t)
			#segments(x0=1-(prop.table(t)/2),x1=1+(prop.table(t)/2),y0=as.numeric(names(t)),lwd=3)

			# estaria bé modificar-ho una mica tipo boxplot. Si hi ha un outlier fer-ho separat
			# ficar una espècie de bandwich?

			n.samp=sum(calculs$inds[1,])
			possibles=round((0:n.samp)/n.samp,10)
			errors=round(calculs$lasso.errors,10)
			possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
			t=factor(errors,levels=possibles)
			#t=prop.table(table(calculs$lasso.errors))
			t=prop.table(table(t))
			t=t*0.9/max(t)			
			polygon(x=c(1-(t/2),rev(1+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
				col="darkgreen")

			n.samp=sum(calculs$inds[1,calculs$g==levels(calculs$g)[1]])
			possibles=round((0:n.samp)/n.samp,10)
			errors=round(calculs$lasso.sens,10)
			possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
			t=factor(errors,levels=possibles)
			t=prop.table(table(t))
			t=t*0.9/max(t)			
			polygon(x=c(2-(t/2),rev(2+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
				col="darkgreen")

			n.samp=sum(calculs$inds[1,calculs$g==levels(calculs$g)[2]])
			possibles=round((0:n.samp)/n.samp,10)
			errors=round(calculs$lasso.spec,10)
			possibles = possibles[1:(max(which(possibles<=max(errors)))+1)]
			t=factor(errors,levels=possibles)
			t=prop.table(table(t))
			t=t*0.9/max(t)			
			polygon(x=c(3-(t/2),rev(3+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
				col="darkgreen")



			#t=prop.table(table(calculs$lasso.sens))
			#t=t*0.9/max(t)			
			#polygon(x=c(2-(t/2),rev(2+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
			#	col="darkgreen")

			#t=prop.table(table(calculs$lasso.spec))
			#t=t*0.9/max(t)			
			#polygon(x=c(3-(t/2),rev(3+(t/2))),y=c(as.numeric(names(t)),rev(as.numeric(names(t)))),
			#	col="darkgreen")

			axis(2,at=seq(0,1,by=by),paste(seq(0,100,by=by*100),"%"),font=2,cex.axis=cex.axis,las=2)
			axis(3,at=1:3,c("ACC","TPR","TNR"),las=2,font=2,cex.axis=cex.axis,tick=F)
			# TPR= True Positive Rate (sensitivity)
			# TNR= True Negative Rate (specificity)
			# ACC = accuracy (error rate)
			mtext("Lasso Performance",side=3,line=4.5,font=2,cex=cex.mtext)
		}




	}

	####### layout 8: llegenda

	if (class(calculs$g)=="factor"){
		#par(mar=c(0,dalt,dalt,0))
		par(mar=c(0,0,0,0))

		plot(1:10,type="n",axes=F,frame.plot=F,ylim=c(-1,9),xlim=c(0.5,11.5))
		legend(x=2,y=8,c("selected","uniqueness","FC","error"),lwd=3,pch=19,lty=1,col=c("#966FAD","#90BA6D","#C3703C","#7E9A99"),cex=cex.axis)

		# horizontal

		col1="chocolate3"  # "saddlebrown"
		col2="seagreen4"  # "seagreen4"
		col.breaks=colorRampPalette(c(col1,"white",col2))(10)

		rect(xleft=1:10,xright=2:11,ybottom=0,ytop=1,border=col.breaks,col=col.breaks)
		text(x=1:11,y=-0.5,seq(0,100,by=10),cex=cex.axis,font=2,srt=90,adj=c(1,0.5))
		#text(x=xleg+1,y=yleg+(2:(length(breaks)-1)),"-")
		text(x=6,y=1.5,"% possible times\ngene column selected\nwhen row selected",srt=0,adj=c(0.5,0),font=2,cex=cex.axis)

		segments(x0=1,x1=11,y0=-0.1,y1=-0.1)
		segments(x0=1:11,x1=1:11,y0=-0.2,y1=-0.1)
	}

	if (class(calculs$g)=="Surv"){
		#par(mar=c(0,dalt,dalt,0))
		par(mar=c(0,0,0,0))

		plot(1:10,type="n",axes=F,frame.plot=F,ylim=c(-1,9),xlim=c(0.5,11.5))
		legend(x=2,y=8,c("selected","uniqueness","p-value","hazard ratio"),lwd=3,pch=19,lty=1,col=c("#966FAD","#90BA6D","#C3703C","#7E9A99"),cex=cex.axis)

		# horizontal

		col1="chocolate3"  # "saddlebrown"
		col2="seagreen4"  # "seagreen4"
		col.breaks=colorRampPalette(c(col1,"white",col2))(10)

		rect(xleft=1:10,xright=2:11,ybottom=0,ytop=1,border=col.breaks,col=col.breaks)
		text(x=1:11,y=-0.5,seq(0,100,by=10),cex=cex.axis,font=2,srt=90,adj=c(1,0.5))
		#text(x=xleg+1,y=yleg+(2:(length(breaks)-1)),"-")
		text(x=6,y=1.5,"% possible times\ngene column selected\nwhen row selected",srt=0,adj=c(0.5,0),font=2,cex=cex.axis)

		segments(x0=1,x1=11,y0=-0.1,y1=-0.1)
		segments(x0=1:11,x1=1:11,y0=-0.2,y1=-0.1)
	}

}

