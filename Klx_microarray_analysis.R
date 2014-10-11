# week 5 - the research frontier! microarray analysis
# big boy gene expression stuff. Affy dataset 22000 #features, 4 treatments + 4 controls reference_series    "GSE2639" geo@ncbi.nlm.nih.gov

# description 
"Analysis of macrovascular umbilical vein endothelial cells (HUVEC) stimulated with tumor necrosis factor (TNF)-alpha for 5 hours. TNF is a potent inflammatory stimulus. Results provide insight into the relevance of the diversity of endothelial cell subtypes for the response to inflammatory stimuli."

# source("http://www.bioconductor.org/biocLite.R"); 
?biocLite
require(GEOquery); library(limma); require(gplots)
gds = getGEO('GDS1542', destdir='.'); show(gds);
eset <- GDS2eSet(gds); show(eset) # convert to eset format
expdata = exprs(eset); dim(expdata); head(expdata)

# Q1 - Apparently there is one data point missing. Let's find the row in order to filter it away.
w=which(apply(is.na(expdata),1,sum) > 0);w
# 210662_at 
# 10127 
expdata[w,]
# GSM50777 GSM50778 GSM50779 GSM50780 GSM50781 GSM50782 GSM50783 
# 1.5      2.5      0.3      2.2       NA      5.9      4.8 
# GSM50784 
# 0.7 

# Q2-Now we know which row to remove. Filter away that row.
temp=expdata[-w,]
expdata=temp; # remove NA row
boxplot(as.data.frame(expdata))#microarray data not ~N(0,1)
logdata=log2(expdata)
boxplot(as.data.frame(logdata))
summary(logdata)

# Q3 - How many false positives would you expect when investigating changes in 22000 probes between identical samples at a p-value threshold of 0.05?
# 22000 * 0.05 = 1100   # 2-tailed...

# Thus, we want to eliminate as many variables as possible that we do not consider biologically meaningful before doing statistical testing.
# Furthermore in our kind of experiment with microarrays, data tend to be very noisy when the signal is very low. Let's have a look on that by plotting the standard deviation against the mean for every probe:
# but first...quantile!
# quantile(probemeans,seq(0,1,0.25))
#     0%       25%       50%       75%      100% 
# -2.010446  3.322966  5.123022  6.790793 12.511584 

probemeans=apply(logdata,1,mean)
probesd = apply(logdata,1,sd)
plot(probemeans,probesd)

# It is clear that standard deviation is higher at lower means. Some probes have veryy high standard deviations - these may be truly changing, something we will investigate further on.
# Let's first eliminate the 25% of the probes that have weakest signal:

q25=quantile(probemeans, 0.25)
whichtosave = which(probemeans > q25) 
q25logdata = logdata[whichtosave,]
# length(probemeans)-length(logdata[whichtosave]) 5571-25%

# Since we are only interested in probes that change, and thus have high variability, we can also remove those with very low variability. A way to do that is to filter on the inter-quartile range (using the IQR function). Keep only those with an IQR above 1.5: ie,
quantile(q25logdata)
#     0%        25%        50%        75%       100% 
# -0.5145732  4.8124982  5.9703935  7.3037807 13.0255895 
mydata = q25logdata[apply(q25logdata, 1, IQR) > 1.5, ]
# IQR(x) = quantile(x, 3/4) - quantile(x, 1/4) 
# length(mydata) [1] 3800
# How many variables remain?
dim(mydata)[1]   #  [1] 475    [1] = rows

# Data exploration
# At this point it would be interesting to do a Principal Component Analysis (PCA). This basically transforms the data to find new orthogonal variables that explain most of the variation in the dataset. 
# 
# This variation can be analysed either between probes (rows) or samples (columns). We will use the function prcomp() for the PCA. It does the analysis between rows. For our purpose the samples are most interesting to compare, so we need to interchange columns and rows in mydata. This is also called to "transpose" the matrix.

# Q4 - Find out how you can transpose mydata and put the result in tdata
tdata = t(mydata) # also tdata = aperm(mydata)
pca = prcomp(tdata,scale=TRUE); summary(pca)
# sdev:the standard deviations of the principal components the square roots of the eigenvalues of the covariance                                  /correlation matrix, though the calculation is actually done with the singular values of the data matrix).

# The first component explains the largest part of the variance, but not all. In the best of worlds, this accounts for the difference between our experimental conditions, otherwise we have some unknown batch effect that dominate the experiment.
# We can plot the samples in relation to the first two components:
plot(pca$x, type='n')
text(pca$x, rownames(pca$x), cex=0.5)

# More interesting is maybe to see which experimental condition they belong to. For that purpose we extract this from the ExpressionSet (remember that?, use show(eset) again...)
conditions = phenoData(eset)$agent
?phenoData
plot(pca$x, type='n'); text(pca$x, labels=conditions, cex = 0.5)

# Luckily the first component divides the samples by condition. However, in one sample something else seems to be going on, sending it away along the second component. That can be worth remembering when judging the final results.
# 
# Another informative plot is to do a dendrogram of correlation between the samples. First we make a correlation matrix between the samples, than a hierarchical clustering is performed:
pearsonCorr = as.dist(1 - cor(mydata))  
hC = hclust(pearsonCorr)
plot(hC, labels = sampleNames(eset))

# The heights of the branches indicate how distant the samples are. Recall which sample was in each condition by putting condition as labels:
plot(hC, labels=conditions)

# The two groups and the sample "in between" are even more evident here than in the PCA.
# Another useful visualization of large datasets is the heatmap. R clusters the data both on rows columns with this command:
heatmap(mydata, col=greenred(100))  

# Red corresponds to high, green to low values. The dendrogram on top is basically the same as we produced earlier, in a slightly different order. Towards the bottom of the heatmap are the probes clustered that discriminate the two conditions clearly.
# 
# Let's now find the probes change significantly between the conditions. There are several tools to do that, more or less advanced. The most simplistic would be to do a t-test for every probe. This is however not a good idea. We have a lot of probes, and the few samples will give the estimate of variance low precision in many cases and give us many false negatives and positives. However, it turns out that the variance of probes with approximately the same expression level is rather similar, and hence one can let probes "borrow" variance from each other to get better variance estimates. Such a method is employed in the limma package (LInear Models for Microarrays).
# limma needs to see the whole dataset, including the high variance probes, to do correct variance estimations. Thus we will go back and use our ExpressionSet containing all data.
# In addition to the actual data, limma needs a model matrix, basically information on which conditions each sample represents. R has a standard function for defining model matrices, model.matrix. It uses the ~ operator to define dependencies. Each input variable should be a factor, so let's first make a factor out of the conditions (agent in our ExpressionSet):
condfactor=factor(eset$agent)

# Construct a model matrix, and assign the names to the columns:
design = model.matrix(~0+condfactor)
colnames(design) = c('ctrl', 'tnf'); design   # any names

# Note that the first 4 samples have '1' in the control level, and the following 4 have the '1' in the other level. If we had had a multifactor experiment (ANOVA style data), we would have included more levels and assigned them in appropriate combinations to the samples.
# The next command estimates the variances:
fit = lmFit(eset, design)

# Now we need to define the conditions we want to compare - trivial in this case since there are only two conditions
contrastmatrix = makeContrasts(tnf - ctrl, levels=design)

# The following commands calculate the p-values for the differences between the conditions defined by the contrast matrix:

fit = contrasts.fit(fit, contrastmatrix)
ebayes = eBayes(fit); ebayes
hist(ebayes$p.value)

# You see that the number of probes are enriched close to p = 0.00. This surplus is due to the genes that are significantly changing between the conditions. If the data had a skewed distribution we might see an accumulation at the p = 1.00 end, indicating a problem with our data.
# 
# If the data were completely random, the p values would be equally distributed from 0 to 1, thus if we from random data picked the probes that had p < 0.05 as "significant", we would pick exactly 1/20 of all probes, all being false positives. Unfortunately, there is no way to discriminate the surplus true probes from the false positives. This is problem with any study where a lot of variables are tested, for instance in large sociology studies. However, there are ways to regulate the p-value cut off in order to have control over the false positives. A common way is the Benjamini-Hochberg adjustment. This adjustment allows you to set a limit to how large fraction of false positive variables (false discovery rate, or FDR) you accept in the results.
# 
# This adjustment, at 5% FDR, is built into the decideTests function:

results = decideTests(ebayes)
# decideTests produces 0, 1, or -1 for each probes, telling if that probe did not pass the test (0), or was significantly increased (1), or decreased (-1)

# Q5-How can we display the number of probes that passed?
sum(results== 1 | results == -1)  # [1] 95
length(which(results!=0))  # [1] 95

# Q6 - There are other ways to control the rate of false positives in a multiple test experiment.
# Which of these is an adjustment method, and available with decideTests?
# A - Holm, Bonferroni

# Data exploration continued
# A neat function to display the result of just two conditions is the Venn diagram:
vennDiagram(results)

#Extract the original data for the most changing probes:
resData = exprs(eset) [results != 0, ] # ie, the 95x8 non-zero samples

geneSymbol = as.array(fData(eset)[, 'Gene symbol']) # fData == feature data
gs = geneSymbol[c(which(results != 0))] # add gene symbols as rownames
rownames(resData) = gs # nice!

# Add p-values in an extra column:
pvalues = ebayes$p.value[results != 0, ]
resData <- cbind(resData, pvalues)

# And add-p values corrected for multiple testing (q-values):
adj.pvalues = p.adjust(ebayes$p.value, method='BH') # ?p.adjust
adj.pvalues <- adj.pvalues[results != 0]
resData = cbind(resData, adj.pvalues)

# write output to file
write.table(resData, 'most_regulated.txt', sep="\t", quote=FALSE)
con = file('most_regulated.txt','r')
most_table=readLines(con); close(con)
str(most_table)
# chr [1:86] "WWTR1\t222.1\t240.8\t255.4\t221.6\t147.4\t123.4\t142.9\t123\t2.042
read_most=read.table('most_regulated.txt', sep="\t")  # doesnt work
read_table = read.table(con, sep="\t")               # ditto
# Error in read.table("most_regulated.txt", sep = "\t") : 
#   duplicate 'row.names' are not allowed
View(most_table)
most_table = read.table(file.choose(), sep='\t') # stackoverflow
# Error in read.table(file.choose(), sep = "\t") : duplicate 'row.names' are not allowed

