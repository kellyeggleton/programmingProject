#!/usr/bin/env python
from __future__ import division
import os
import argparse
import glob
import sys



#proportions based on 14 pat BAMs

#class creates object

numberOfSubsamples=3

class QueryManager:
	storageDir=None	
	def setStorage(self, storageDir):
		self.storageDir=storageDir
	
	#sets the function 'query' and describes what imputs are needed	
	def query(self, geneListFile, numberOfSamples, outFileName):
		outFileName2=outFileName+"_depth"
		#loads a genelist from the geneListFile provided			
		geneList=self.loadgeneList(geneListFile)
		#grab relevent files from input specified in command line
		percentCoverFilenames=glob.glob ("%s/*.%s.*averageCovPercentage.txt"%(self.storageDir, numberOfSamples))
		#prepare dictionary of coverage percentages (20x)		
		geneCoveragePercentages={}
		for gene in geneList:
			geneCoveragePercentages[gene]=[] 
		#search for gene coverages in all files	
		for filename in percentCoverFilenames:
			for line in open(filename):
				line=line.strip()
				line=line.split("\t")
				geneName=line[0]
				coveragePerc=float(line[1])
				
				if geneName in geneList:
					geneCoveragePercentages[geneName].append(coveragePerc)
		#print geneCoveragePercentages 
		#take average coverage from across all GIAB files	
		averageCoveragePercentage={}	
		for gene in geneList:
			coverageList=geneCoveragePercentages[gene]
			if len(coverageList) > 0 :
				averageCoverage=sum(coverageList)/len(coverageList)
				averageCoveragePercentage[gene]=averageCoverage
			else :
				averageCoveragePercentage[gene]=None
		#output dictionary to file
		outfile=open(outFileName+".%s.20xcoverage"%(numberOfSamples), "w")
		for gene in averageCoveragePercentage:
			outfile.write("%s\t%s\n"%(gene, averageCoveragePercentage[gene])) 
		outfile.close()



		#get relevent depth of coverage files
		coverageDepthFilenames=glob.glob ("%s/*.%s.*averageCovDepth.txt"%(self.storageDir, numberOfSamples))
		#prepare dictionary of coverage depths
		geneCoverageDepth={}
		for gene in geneList:
			geneCoverageDepth[gene]=[]
		for filename in coverageDepthFilenames:
			for line in open(filename):
				line=line.strip()
				line=line.split("\t")
				geneName=line[0]
				coverageDepth=float(line[1])	
				if geneName in geneList:
					geneCoverageDepth[geneName].append(coverageDepth)
		#print geneCoverageDepth
		averageCoverageDepth={}
		for gene in geneList:
			depthList=geneCoverageDepth[gene]
			if len(depthList) > 0 :			
				averageDepth=sum(depthList)/len(depthList)
				averageCoverageDepth[gene]=averageDepth
			else :
				averageCoverageDepth[gene]=None
		outfile2=open(outFileName2+".%s"%(numberOfSamples), "w")
		for gene in averageCoverageDepth:
			outfile2.write("%s\t%s\n"%(gene, averageCoverageDepth[gene])) 
		outfile.close()	



	def loadgeneList(self, geneListFile):
		geneList=[]		
		for line in open(geneListFile):
			line=line.strip()
			geneName=line
			if geneName=='':
				continue
			geneList.append(geneName) 
		return geneList


class InputManager:
	mincoverage=20
			

	storageDir=None	
	def setStorage(self, storageDir):
		self.storageDir=storageDir
		if not os.path.exists(storageDir):
			print 'Creating storage directory'
			os.mkdir(storageDir)

			
	subsamples=[14,15,16]	
	def subSample(self, bam):
		#subsample the bam file, remember name of subsampled bams
		for subsample in self.subsamples:
			proportion=14.0/subsample
			#create output name
			outprefix=bam.replace(".bam", ".%s" %subsample)			
			#don't bother downsampleing if proportion=100%
			if proportion==1.0:
				cmd="ln -s %s %s.bam" %(bam, outprefix)
				print cmd
				os.system(cmd)
				#index the bam
				cmd="samtools index %s.bam" %(outprefix)
				print cmd
				os.system(cmd)
				
				#write like command above to run from terminal
				cmd="genomeCoverageBed -bga -ibam %s.bam >  %s.genomecoverage.bed" %(outprefix, outprefix)
				print cmd
				os.system(cmd)
				cmd="intersectBed -wb -a %s.genomecoverage.bed -b Agilent_CORRECT_Genes.bed > %s.genecov.bed" %(outprefix, outprefix) 
				print cmd
				os.system(cmd)

				self.extractCoverageDepth ("%s.genecov.bed"%(outprefix), "%s/%s.averageCovDepth.txt"%(self.storageDir, outprefix))
				self.extractCoveredProportion ("%s.genecov.bed"%(outprefix), "%s/%s.averageCovPercentage.txt"%(self.storageDir, outprefix))			 

			else:
				for seed in range (0, numberOfSubsamples):
					seedplusprop=proportion+seed
					fulloutprefix=outprefix+".%s"%(seed)
					cmd="samtools view -s %.5f -b -h %s > %s.bam" %(seedplusprop, bam, fulloutprefix)				
					print cmd
					os.system(cmd) #DEBUG
			
					#index the bam
					cmd="samtools index %s.bam" %(fulloutprefix)
					print cmd
					os.system(cmd)  #DEBUG
			
					#write like command above to run from terminal
					cmd="genomeCoverageBed -bga -ibam %s.bam >  %s.genomecoverage.bed" %(fulloutprefix, fulloutprefix)
					print cmd
					os.system(cmd)
					cmd="intersectBed -wb -a %s.genomecoverage.bed -b Agilent_CORRECT_Genes.bed > %s.genecov.bed" %(fulloutprefix, fulloutprefix) 
					print cmd
					os.system(cmd)

					self.extractCoverageDepth ("%s.genecov.bed"%(fulloutprefix), "%s/%s.averageCovDepth.txt"%(self.storageDir, fulloutprefix))
					self.extractCoveredProportion ("%s.genecov.bed"%(fulloutprefix), "%s/%s.averageCovPercentage.txt"%(self.storageDir, fulloutprefix))			
					
					cmd= "rm %s.genomecoverage.bed" %(fulloutprefix)				
					cmd1= "rm %s.bam" %(fulloutprefix)
	
					print cmd
					print cmd1
					os.system(cmd)
					os.system(cmd1)
				
	
	#get coverage depth from BED - create dict
	def extractCoverageDepth (self, covbed, outFileName):
		'''covbed=intersected coverage bedfile eg %s.genecov.bed'''
		coverageDict = {}
		geneLengthDict = {} #{GENENAME : 152637bp}
		for line in open(covbed):
			line=line.strip()
			field=line.split('\t') 
			#['1', '169680651', '169680682', '3', '1', '169680526', '169680707', 'SELL']
			geneName=field[7]
			startPosition=int(field[1])
			endPosition=int(field[2])
			coverage=float(field[3])
			#store region length
			regionLength=endPosition-startPosition
			#check if gene is in dictionary
			if geneName not in geneLengthDict:
				geneLengthDict[geneName]=0
			geneLengthDict[geneName]=geneLengthDict[geneName]+regionLength
			#store total coverage
			totalCoverage=regionLength*coverage
			if geneName not in coverageDict:
				coverageDict[geneName]=0.0
			coverageDict[geneName]=coverageDict[geneName]+totalCoverage
		#calculate average coverage per gene
		geneList=coverageDict.keys()
		geneCoverage= {}
		for gene in geneList:
			geneCoverage[gene]=coverageDict[gene]/geneLengthDict[gene]
		#export to file	
		outfile=open(outFileName, "w")
		for gene in geneList:
			outfile.write("%s\t%s\n"%(gene, geneCoverage[gene])) 	

		outfile.close()	



	def extractCoveredProportion (self, covbed, outFileName):
		'''covbed=intersected coverage bedfile eg %s.genecov.bed'''
		coveredBasesDict = {}
		geneLengthDict = {} #{GENENAME : 152637bp}
		for line in open(covbed):
			line=line.strip()
			field=line.split('\t') 
			#['1', '169680651', '169680682', '3', '1', '169680526', '169680707', 'SELL']
			geneName=field[7]
			startPosition=int(field[1])
			endPosition=int(field[2])
			coverage=float(field[3])
			#store region length
			regionLength=endPosition-startPosition
			#check if gene is in dictionary
			if geneName not in geneLengthDict:
				geneLengthDict[geneName]=0
			geneLengthDict[geneName]=geneLengthDict[geneName]+regionLength
			#add region length if covered			
			if geneName not in coveredBasesDict:
				coveredBasesDict[geneName]=0.0
			if coverage >= self.mincoverage:
				coveredBasesDict[geneName]=coveredBasesDict[geneName]+regionLength
			
		#calculate average coverage per gene
		geneList=coveredBasesDict.keys()
		geneCoveragePercent = {}
		for gene in geneList:
			geneCoveragePercent[gene]=(coveredBasesDict[gene]/geneLengthDict[gene])*100.0
		#export to file	
		outfile=open(outFileName, "w")
		for gene in geneList:
			outfile.write("%s\t%s\n"%(gene, geneCoveragePercent[gene])) 	

		outfile.close()	



#coverage['s%']=     %(genes)#dont think this is possible - genes from virtual panel
	

parser = argparse.ArgumentParser(description='loads BAM file and generates average coverage')
parser.add_argument('--bam', help='input bam file based on 14 patients per run')
parser.add_argument('--storagedir', help='directory for storing average coverages', default='./average_gene_coverages/')
parser.add_argument('--geneListFile', help='gene list for virtual panel')
parser.add_argument('--numberOfSamples', help='number of samples on the run', type=int)
parser.add_argument('--outFileName', help='name of file for query output eg breadth of storage')
args = parser.parse_args()

#check files are OK
if not os.path.exists(args.bam):
	sys.exit("Bam file %s does not exist" %(args.bam))

#check output directory is not aready filename in directory
if os.path.isfile(args.storagedir):
	sys.exit("%s is already the name of a file in this directory, please choose a new storage directory name"%(args.storagedir))

if not args.bam.endswith(".bam"):
	sys.exit ("Bam file does not end in .bam (this is case sensitive)")

if args.bam != None :
	inputmanager=InputManager()
	inputmanager.setStorage(args.storagedir)
	inputmanager.subSample(args.bam)

if args.geneListFile != None :
	querymanager=QueryManager()
	querymanager.setStorage(args.storagedir)
	querymanager.query(args.geneListFile, args.numberOfSamples, args.outFileName)
