# generate lollipop plot for TP53

library(magrittr)
library(g3viz)

TP53_variants = readr::read_tsv(snakemake@input[[1]])
TP53_variants = TP53_variants %>% dplyr::select(fk_britroc_number,`#Uploaded_variation`,classification)
TP53_variants = TP53_variants %>% dplyr::filter(classification=='clonal')

TP53_variants_sum = TP53_variants %>%
	dplyr::group_by(`#Uploaded_variation`) %>%
	dplyr::summarise(n=dplyr::n()) %>%
	dplyr::arrange(-n) %>%
	dplyr::filter(n>1)

print(TP53_variants_sum, n=Inf)

TP53_annotations_archival = readr::read_tsv(
	snakemake@input[['TP53_annotations_archival']]
	)
TP53_annotations_relapse = readr::read_tsv(
        snakemake@input[['TP53_annotations_relapse']]
       )

TP53_annotations = rbind(
	TP53_annotations_archival,
	TP53_annotations_relapse
)

TP53_annotations = TP53_annotations %>% dplyr::select('Consequence','#Uploaded_variation','Extra')
TP53_annotations$Extra = TP53_annotations$Extra %>% 
	stringr::str_extract('HGVSp=[A-Za-z0-9:\\.]+;') %>% 
	stringr::str_remove(';$') %>%
	stringr::str_remove('HGVSp=[A-Za-z0-9\\.]+:')

TP53_annotations$gene_symbol = 'TP53'
TP53_annotations = TP53_annotations %>% unique()

TP53_variants = TP53_variants %>% dplyr::select('#Uploaded_variation')

print(TP53_annotations)
print(TP53_variants)

TP53_variants =  dplyr::inner_join(
	TP53_variants,
	TP53_annotations,
	by=c('#Uploaded_variation')
)

TP53_variants = TP53_variants %>% dplyr::select(-'#Uploaded_variation')
TP53_variants = TP53_variants %>% dplyr::filter(!is.na(Extra))

print(TP53_variants)
TP53_variants = TP53_variants %>% dplyr::select(gene_symbol,Consequence,Extra)
colnames(TP53_variants) = c('Hugo_Symbol','Variant_Classification','amino_acid_change')

readr::write_tsv(TP53_variants, 'foo.tsv')

# plot lollipop ----

mutation.csv <- system.file("extdata", "ccle.csv", package = "g3viz")
x = readr::read_csv(mutation.csv)
print(x, width=Inf)

#print(TP53_variants)
#colnames(TP53_variants) %>% print()

x$Variant_Classification %>% unique() %>% print()
TP53_variants$Consequence %>% unique() %>% print()

#TP53_variants$Consequence = TP53_variants$Consequence %>%  
#	dplyr::recode(
#	'missense_variant'='Missense_Mutation'	
#)

mutation.dat = g3viz::readMAF(
		'foo.tsv',
		gene.symbol.col = 'Hugo_Symbol',
		variant.class.col = 'Variant_Classification',
		protein.change.col = 'amino_acid_change',
		sep='\t'
	)

# set up chart options
plot.options <- g3Lollipop.options(
  # Chart settings
  chart.width = 600,
  chart.type = "pie",
  chart.margin = list(left = 30, right = 20, top = 20, bottom = 30),
  chart.background = "#d3d3d3",
  transition.time = 300,
  # Lollipop track settings
  lollipop.track.height = 200,
  lollipop.track.background = "#d3d3d3",
  lollipop.pop.min.size = 1,
  lollipop.pop.max.size = 8,
  lollipop.pop.info.limit = 5.5,
  lollipop.pop.info.dy = "0.24em",
  lollipop.pop.info.color = "white",
  lollipop.line.color = "#a9A9A9",
  lollipop.line.width = 3,
  lollipop.circle.color = "#ffdead",
  lollipop.circle.width = 0.4,
  lollipop.label.ratio = 2,
  lollipop.label.min.font.size = 12,
  lollipop.color.scheme = "dark2",
  highlight.text.angle = 60,
  # Domain annotation track settings
  anno.height = 16,
  anno.margin = list(top = 0, bottom = 0),
  anno.background = "#d3d3d3",
  anno.bar.fill = "#a9a9a9",
  anno.bar.margin = list(top = 4, bottom = 4),
  domain.color.scheme = "pie5",
  domain.margin = list(top = 2, bottom = 2),
  domain.text.color = "white",
  domain.text.font = "italic 8px Serif",
  # Y-axis label
  y.axis.label = "# of TP53 gene mutations",
  axis.label.color = "#303030",
  axis.label.alignment = "end",
  axis.label.font = "italic 12px Serif",
  axis.label.dy = "-1.5em",
  y.axis.line.color = "#303030",
  y.axis.line.width = 0.5,
  y.axis.line.style = "line",
  y.max.range.ratio = 1.1,
  # Chart title settings
  title.color = "#303030",
  title.text = "TP53 gene (customized chart options)",
  title.font = "bold 12px monospace",
  title.alignment = "start",
  # Chart legend settings
  legend = TRUE,
  legend.margin = list(left=20, right = 0, top = 10, bottom = 5),
  legend.interactive = TRUE,
  legend.title = "Variant classification",
  # Brush selection tool
  brush = TRUE,
  brush.selection.background = "#F8F8FF",
  brush.selection.opacity = 0.3,
  brush.border.color = "#a9a9a9",
  brush.border.width = 1,
  brush.handler.color = "#303030",
  # tooltip and zoom
  tooltip = TRUE,
  zoom = TRUE
)

g3viz::g3Lollipop(mutation.dat,
           gene.symbol = "TP53",
           protein.change.col = "amino_acid_change",
           btn.style = "blue", # blue-style chart download buttons
           plot.options = plot.options,
           output.filename = "customized_plot")

quit()
