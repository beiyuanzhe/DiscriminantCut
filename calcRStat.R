Args <- commandArgs()

matName <- Args[6]
method <- Args[7]
dataName <- Args[8]

#matName <- "tmpData_0.3768.mat"
#method <- "voom" 
#dataName <- "simPDF_6v6r1k_v1_15to30rm16907_1"

cat(matName)
cat(method)

library(R.matlab)

cat(paste("Method = ", method, "\n"))

data <- readMat(matName)
data$design <- as.character(data$design)
data$design[data$design == "0"] <- "A"
data$design[data$design == "1"] <- "B" 

if (method == "DESeq2") {

	library(DESeq2)

	permIdx <- data$permIdx
	grp1 <- data.frame(grp = as.character(data$design))
	count <- pmax(round(data$M),0)

	if (is.null(permIdx)) { # use original

		dds <- DESeqDataSetFromMatrix(countData = count, colData = grp1, design = ~ grp)
		dds <- DESeq(dds)
		res <- results(dds)

		writeMat(matName, method = method, st = res$stat, p = res$pvalue, fdr = res$padj)

	} else { # use permIdx

		nPerm <- nrow(permIdx)
		nR <- nrow(count)
		st <- matrix(NA, nrow = nR, ncol = nPerm)
		p <- matrix(NA, nrow = nR, ncol = nPerm)
		fdr <- matrix(NA, nrow = nR, ncol = nPerm)
		
		for (i in seq(along = 1:nPerm)) {

			dds <- DESeqDataSetFromMatrix(countData = count[,permIdx[i,]], colData = grp1, design = ~ grp)
			dds <- DESeq(dds)
			res <- results(dds)

			st[,i] <- res$stat
			p[,i] <- res$pvalue
			fdr[,i] <- res$padj

		}

		writeMat(matName, method = method, st = st, p = p, fdr = fdr)

	}

} else if (method == 'nmlimma') {
	
	library(limma)
	library(edgeR)

	permIdx <- data$permIdx
	count <- data$M
	count.nm = sweep(count, 2, colSums(count), "/")
	class <- as.numeric(factor(data$design))
	class <- class - mean(c(min(class),max(class)))

	if (is.null(permIdx)) { # use original
		
		limma.fitlimma = lmFit(count.nm, design = model.matrix(~factor(class)))
		limma.fitbayes = eBayes(limma.fitlimma)
		limma.stat = limma.fitbayes$t[, 2]
		limma.pvalues = limma.fitbayes$p.value[, 2]
		limma.adjpvalues = p.adjust(limma.pvalues, method = "BH")

		writeMat(matName, method = method, st = limma.stat, p = limma.pvalues, fdr = limma.adjpvalues)

	} else { # use permIdx

		nPerm <- nrow(permIdx)
		nR <- nrow(count)
		st <- matrix(NA, nrow = nR, ncol = nPerm)
		p <- matrix(NA, nrow = nR, ncol = nPerm)
		fdr <- matrix(NA, nrow = nR, ncol = nPerm)
		
		for (i in seq(along = 1:nPerm)) {
			
			count1 <- count.nm[,permIdx[i,]]

			limma.fitlimma = lmFit(count1, design = model.matrix(~factor(class)))
			limma.fitbayes = eBayes(limma.fitlimma)
			limma.stat = limma.fitbayes$t[, 2]
			limma.pvalues = limma.fitbayes$p.value[, 2]
			limma.adjpvalues = p.adjust(limma.pvalues, method = "BH")

			st[,i] <- limma.stat
			p[,i] <- limma.pvalues
			fdr[,i] <- limma.adjpvalues
		}

		writeMat(matName, method = method, st = st, p = p, fdr = fdr)

	}
	

} else if (method == "limma") {

	library(limma)

	permIdx <- data$permIdx
	count <- data$M
	class <- as.numeric(factor(data$design))
	class <- class - mean(c(min(class),max(class)))

	if (is.null(permIdx)) { # use original

		limma.fitlimma = lmFit(count, design = model.matrix(~factor(class)))
		limma.fitbayes = eBayes(limma.fitlimma)
		limma.stat = limma.fitbayes$t[, 2]
		limma.pvalues = limma.fitbayes$p.value[, 2]
		limma.adjpvalues = p.adjust(limma.pvalues, method = "BH")

		writeMat(matName, method = method, st = limma.stat, p = limma.pvalues, fdr = limma.adjpvalues)

	} else { # use permIdx

		nPerm <- nrow(permIdx)
		nR <- nrow(count)
		st <- matrix(NA, nrow = nR, ncol = nPerm)
		p <- matrix(NA, nrow = nR, ncol = nPerm)
		fdr <- matrix(NA, nrow = nR, ncol = nPerm)
		
		for (i in seq(along = 1:nPerm)) {
			
			count1 <- count[,permIdx[i,]]

			limma.fitlimma = lmFit(count1, design = model.matrix(~factor(class)))
			limma.fitbayes = eBayes(limma.fitlimma)
			limma.stat = limma.fitbayes$t[, 2]
			limma.pvalues = limma.fitbayes$p.value[, 2]
			limma.adjpvalues = p.adjust(limma.pvalues, method = "BH")

			st[,i] <- limma.stat
			p[,i] <- limma.pvalues
			fdr[,i] <- limma.adjpvalues
		}

		writeMat(matName, method = method, st = st, p = p, fdr = fdr)

	}
	

} else if (method == "voom") {

	library(limma)
	library(edgeR)

	permIdx <- data$permIdx
	count <- pmax(data$M,0);
	class <- as.numeric(factor(data$design))
	class <- class - mean(c(min(class),max(class)))

	if (is.null(permIdx)) { # use original

		nf = calcNormFactors(count, method = "TMM")
		voom.data = voom(count, design = model.matrix(~factor(class)), lib.size = colSums(count) * nf)
		voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(class)))
		voom.fitbayes = eBayes(voom.fitlimma)
		voom.stat = voom.fitbayes$t[, 2]
		voom.pvalues = voom.fitbayes$p.value[, 2]
		voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")

		writeMat(matName, method = method, st = voom.stat, p = voom.pvalues, fdr = voom.adjpvalues)

	} else { # use permIdx

		nPerm <- nrow(permIdx)
		nR <- nrow(count)
		st <- matrix(NA, nrow = nR, ncol = nPerm)
		p <- matrix(NA, nrow = nR, ncol = nPerm)
		fdr <- matrix(NA, nrow = nR, ncol = nPerm)
		
		for (i in seq(along = 1:nPerm)) {
			
			count1 <- count[,permIdx[i,]]

			nf = calcNormFactors(count1, method = "TMM")
			voom.data = voom(count1, design = model.matrix(~factor(class)), lib.size = colSums(count1) * nf)
			voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(class)))
			voom.fitbayes = eBayes(voom.fitlimma)
			voom.stat = voom.fitbayes$t[, 2]
			voom.pvalues = voom.fitbayes$p.value[, 2]
			voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")

			st[,i] <- voom.stat
			p[,i] <- voom.pvalues
			fdr[,i] <- voom.adjpvalues
		}

		writeMat(matName, method = method, st = st, p = p, fdr = fdr)

	}
	
} else if (method == "vst") {

	library(limma)
	library(DESeq)

	permIdx <- data$permIdx
	count <- round(data$M)
	class <- as.numeric(factor(data$design))
	class <- class - mean(c(min(class),max(class)))

	if (is.null(permIdx)) { # use original

		DESeq.cds = newCountDataSet(countData = count, conditions = factor(class))
		DESeq.cds = estimateSizeFactors(DESeq.cds)
		DESeq.cds = estimateDispersions(DESeq.cds, method = "blind", fitType = "local")
		DESeq.vst = getVarianceStabilizedData(DESeq.cds)
		DESeq.vst.fitlimma = lmFit(DESeq.vst, design = model.matrix(~factor(class)))
		DESeq.vst.fitbayes = eBayes(DESeq.vst.fitlimma)
		DESeq.vst.stat = DESeq.vst.fitbayes$t[, 2]
		DESeq.vst.pvalues = DESeq.vst.fitbayes$p.value[, 2]
		DESeq.vst.adjpvalues = p.adjust(DESeq.vst.pvalues, method = "BH")

		writeMat(matName, method = method, st = DESeq.vst.stat, p = DESeq.vst.pvalues, fdr = DESeq.vst.adjpvalues)

	} else { # use permIdx

		nPerm <- nrow(permIdx)
		nR <- nrow(count)
		st <- matrix(NA, nrow = nR, ncol = nPerm)
		p <- matrix(NA, nrow = nR, ncol = nPerm)
		fdr <- matrix(NA, nrow = nR, ncol = nPerm)
		
		for (i in seq(along = 1:nPerm)) {
			
			count1 <- count[,permIdx[i,]]

			DESeq.cds = newCountDataSet(countData = count1, conditions = factor(class))
			DESeq.cds = estimateSizeFactors(DESeq.cds)
			DESeq.cds = estimateDispersions(DESeq.cds, method = "blind", fitType = "local")
			DESeq.vst = getVarianceStabilizedData(DESeq.cds)
			DESeq.vst.fitlimma = lmFit(DESeq.vst, design = model.matrix(~factor(class)))
			DESeq.vst.fitbayes = eBayes(DESeq.vst.fitlimma)
			DESeq.vst.stat = DESeq.vst.fitbayes$t[, 2]
			DESeq.vst.pvalues = DESeq.vst.fitbayes$p.value[, 2]
			DESeq.vst.adjpvalues = p.adjust(DESeq.vst.pvalues, method = "BH")

			st[,i] <- DESeq.vst.stat
			p[,i] <- DESeq.vst.pvalues
			fdr[,i] <- DESeq.vst.adjpvalues
		}

		writeMat(matName, method = method, st = st, p = p, fdr = fdr)

	}

} else if (method == "baySeq") {

	library(baySeq)

	count <- round(data$M)
	grp <- as.character(data$design)
	design <- list(NDE = rep(1,ncol(count)),
		DE = as.numeric(as.factor(grp)))
	

	CD <- new("countData", data = count, replicates = grp, groups = design)
	libsizes(CD) <- getLibsizes(CD)
	CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = NULL)
	CD <- getLikelihoods(CD, pET = 'BIC', cl = NULL)
	
	res <- topCounts(CD, group = "DE", number = nrow(data))

	fdr <- rep(NA, nrow(count))
	fdr[res$rowID] <- res$FDR.DE
	pval <- rep(NA, nrow(count))
	pval[res$rowID] <- 1 - res$Likelihood

	writeMat(matName, method = method, p = pval, fdr = fdr)

} else if (method == "NBPSeq") {

	library(NBPSeq)

	count <- round(data$M)
	grp <- as.numeric(as.factor(as.character(data$design)))

	set.seed(999) # to make results reproducible
	NFactor <- estimate.norm.factors(count, lib.sizes = colSums(count), method = "AH2010")
	result <- nbp.test(count, grp, 1, 2,
		norm.factors=NFactor,
		model.disp = "NBP")

	writeMat(matName, method = method, p = result$p.value, fdr = result$q.value)

} else if (method == "edgeR") {

	library("edgeR")
	
	grp <- as.character(data$design)

	dgl <- DGEList(counts=data$M, group=grp, lib.size=colSums(data$M))
	dgl <- calcNormFactors(dgl, refColumn=1)
	dgl <- estimateCommonDisp(dgl)
	dgl <- estimateTagwiseDisp(dgl)

	res <- exactTest(dgl, dispersion="tagwise")
	result <- topTags(res, n=nrow(res), adjust.method="BH", sort.by="none")

	writeMat(matName, method = method, p = result$table$PValue, fdr = result$table$FDR)

}





#	dataFullName <- paste("Rdata/",dataName,".RData",sep = "")
#	if (isOriginal) {
#		# original
#
#		dds <- DESeqDataSetFromMatrix(countData = round(data$M), colData = grp1, design = ~ grp)
#		dds <- estimateSizeFactors(dds)
#		dds <- estimateDispersions(dds)
#
#		save(dds, file = dataFullName)
#
#		dds <- nbinomWaldTest(dds)
#		res <- results(dds)
#
#	} else {
#		# permuted
#
#		load(dataFullName)
#		
#		counts(dds) <- matrix(as.integer(round(data$M)), nrow = nrow(data$M), ncol = ncol(data$M))
#
#		dds <- estimateSizeFactors(dds)
#
#		dds <- nbinomWaldTest(dds)
#		res <- results(dds)
#	}
