library('biomaRt')
namelist = c('opn1sw1','opn1sw2','opn1lw1','opn1lw2','exorh','Rho','gabrr1','gabbr2')
datasets = read.table('datasets.txt',header=FALSE,sep='')
#write('',file='SEQUENCES.txt')
for (set in datasets[,1])
{
mart = useMart("ENSEMBL_MART_ENSEMBL",dataset=set)      ###for each dataset

	for (name in namelist)	###for each name in namelist within each dataset
	{
		var = getBM(attributes=c('chromosome_name', 'start_position', 'end_position'), ###var returns desired attributes of each gene, one by one
			filters = 'external_gene_name',
			values = name,
			mart = mart)

		if (class(var[1,1]) == 'logical')
		{
			next
		}

		if (class(var[2,1]) != 'logical')
		{
			for (i in 1:length(var[,1]))
			{
				sub = var[i,]
				x = getSequence(chromosome = sub[1], start = sub[2], end = sub[3],type=c('entrezgene_id','entrezgene_description'),seqType='cdna',mart=mart)
				y = getSequence(chromosome = sub[1], start = sub[2], end = sub[3],type=c('entrezgene_id','entrezgene_description'),seqType='peptide',mart=mart)

				write.table(c(name,set,x),file='SEQUENCES.txt',append='TRUE')
				write.table(c(name,set,y),file='SEQUENCES.txt',append='TRUE')
			}
			next
		}

		x = getSequence(chromosome = var[1], start = var[2], end = var[3],type=c('entrezgene_id','entrezgene_description'),seqType='cdna',mart=mart)
		y = getSequence(chromosome = var[1], start = var[2], end = var[3],type=c('entrezgene_id','entrezgene_description'),seqType='peptide',mart=mart)

		write.table(c(name,set,x),file='SEQUENCES.txt',append='TRUE')
		write.table(c(name,set,y),file='SEQUENCES.txt',append='TRUE')
	}
}
