# Data

Alus <- readRDS(system.file(
	"extdata", "AlusObj.Rds", package = "ShinyAluExplorer"
))

hub <- AnnotationHub::AnnotationHub()

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

AluColours <- c(
	"AluY" = "#0085B4",
	"Y" = "#0085B4",
	"AluS" = "#BA4E9D",
	"S" = "#BA4E9D",
	"AluJ" = "#568400",
	"J" = "#568400"
)

## Functions

### AluSubFamSplit
AluSubFamSplit <- function(alus) {
	alus %>%
		tibble::as_tibble() %>%
		tidyr::extract(
			subfam,
			into = c("family", "subfamily", "type", "subtype"),
			regex = "(FAM|FLAM|FRAM|Alu)_?(J|S|Y)(\\w*)(\\d*)",
			remove = FALSE
		) %>%
		dplyr::filter(family == "Alu") %>%
		dplyr::mutate(
			subtype = paste0(type, subtype),
			type = paste0(subfamily, type),
			subfamily = paste0(family, subfamily)
		)
}

### AluCountByType
AluCountByType <- function(alus) {
	alus %>%
		dplyr::distinct(aluIndex, .keep_all = TRUE) %>% ###!!!
		AluSubFamSplit() %>%
		dplyr::group_by(family, subfam, subfamily, type, subtype) %>%
		dplyr::summarise(n = n())
}

### AluTreeMap
AluTreeMap <- function(df, over = NA, log10ed = FALSE, d3 = TRUE) {
	title <- paste0(
		" Alu Copy Number by Subfamily",
		ifelse(is.na(over),"", paste("overlapping ", over)),
		ifelse(log10ed, " (log10 scaled)", "")
	)

	if(log10ed) {
		df <- df %>% dplyr::mutate(log10n = log10(n))
	}

	tm <- df %>%
		treemap::treemap(
			index = c("subfamily", "type", "subtype"),
			vSize = ifelse(log10ed,"log10n","n"),
			title = title,
			draw = FALSE
		)

	if(d3){
		tm <- d3treeR::d3tree(tm, rootname = title, height = "500px")
	}

	return(tm)
}

aluPackedBubblecirc <- function(df) {
	df %>%
		mutate(pathString = paste("Alu", subfamily, type, sep = "/")) %>%
		data.tree::as.Node() %>%
		circlepackeR::circlepackeR(size = "n")
}

# AluPackedBubble <- function(df, over = NA, log10ed = FALSE, pre = "") {
# 	title <- paste0(
# 		pre,
# 		" Alu Copy Number by Subfamily",
# 		ifelse(is.na(over),"", paste("overlapping ", over)),
# 		ifelse(log10ed, " (log10 scaled)", "")
# 	)
#
# 	if(log10ed) {
# 		df <- df %>% mutate(n = log10(n))
# 	}
#
# 	tm <- hpackedbubble::hpackedbubble(
# 		df$subfamily,
# 		df$type,
# 		#df$subtype,
# 		df$n,
# 		#pointFormat = "<b>{point.name}:</b> {point.y}m CO<sub>2</sub>",
# 		title = title,
# 		split = 1,#1
# 		#theme = "sunset",
# 		width = "500px",
# 		height = "500px",
# 		dataLabelsFilter = 5,
# 		packedbubbleMinSize = "50%",
# 		packedbubbleMaxSize = "150%",
# 		packedbubbleZMin = 0,
# 		#packedbubbleZmax = 1000,
# 		#gravitational = 0.0626,
# 		gravitational = 0.5,
# 		parentNodeLimit = 1,
# 		dragBetweenSeries = 0,
# 		seriesInteraction = 0
# 	)
#
# 	return(tm)
# }

aluSunburst <- function(df,...) {
	df %>%
		ungroup() %>%
		mutate(pathString = paste(family, subfamily, type, sep = "-")) %>%
		dplyr::select(pathString, size = n) %>%
		sunburstR::sunburst(...)
}

### aluFeatOver
aluFeatOver <- function(feat, txdb, accessor, minoverlap = 0L){
	set <- plyranges::filter_by_overlaps(
		feat,
		accessor(txdb)
	)
	AluSubFamSplit(set)
}

### nFeatOver
nFeatOver <- function(feat, txdb, accessor, minoverlap = 0L){
	set <- plyranges::filter_by_overlaps(
		feat,
		accessor(txdb)
	)
	AluCountByType(set)
}

problemRegions <- c(
	"hg19_custom_AluSVs", "hg19_custom_hsm-peaks-Bell2017",
	"hg19_custom_blacklist_Amemiya2019"
)

hg19annos <- annotatr::builtin_annotations()
hg19annos <- hg19annos[grepl(hg19annos, pattern = "hg19")]
hg19annos <- c(hg19annos, problemRegions)
hg19annos

myAnno <- function(feats, annonms, setop = "or") { # "and", "or", "not"
	#annos <- build_annotations(genome = 'hg19', annotations = annonms)
	switch (
		setop,
		"and" = {
			for (i in seq_along(annonms)) {
				anb <- annotatr::build_annotations(
					genome = 'hg19', annotations = annonms[i]
				)
				feats <- plyranges::filter_by_overlaps(feats, anb)
			}
			feats
		},
		"or" = {
			annos <- annotatr::build_annotations(
				genome = 'hg19', annotations = annonms
			)
			plyranges::filter_by_overlaps(feats, annos)
		},
		"not" = {
			annos <- annotatr::build_annotations(
				genome = 'hg19', annotations = annonms
			)
			plyranges::filter_by_non_overlaps(feats, annos)
		}
	)
	#plyranges::filter_by_overlaps(feats, annos)
}

# myAnno(Alus, "hg19_cpg_islands")
# myAnno(Alus, c("hg19_cpg_islands", "hg19_genes_promoters"), "and")
# myAnno(Alus, c("hg19_cpg_islands", "hg19_genes_promoters"), "not")

## Shiny

### DT format

AluAgeDT <- function(df, nhead, pval, slp, nltpr) {
	slp <- sym(slp)
	pval <- sym(pval)
	nltp <- sym(nltpr)
	## Pre-processing
	df <- df %>%
		head(nhead) %>%
		dplyr::mutate(
			`p-value` = sprintf("%.4e", {{pval}}),
			Coords = paste0(
				tools::toTitleCase(as.character(seqnames)), " : ",
				format(start, big.mark = ","), " - ",
				format(end, big.mark = ",")
			),
			slope = sprintf("%.3f", {{slp}}),
			dir = sign({{slp}})
		) %>%
		dplyr::select(
			Coords, `p-value`, `Alu Type` = subfam,
			slope, dir, {{nltp}}, subfamily, AluWidth, AluStrand
		)
	## DT formatting
	df %>%
		DT::datatable(
			options = list(
				columnDefs = list(
					list(
						visible = FALSE,
						targets = (c(5, 6, 7) - 1)
					)
				)
			),
			rownames = FALSE,
			selection = 'single'
		) %>%
		DT::formatStyle(
			'slope', 'dir',
			backgroundColor = DT::styleEqual(c(-1, 1), c('#998ec3', '#f1a340')),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		) %>%
		DT::formatStyle(
			'Alu Type', 'subfamily',
			backgroundColor = DT::styleEqual(
				names(AluColours), AluColours
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		) %>%
		DT::formatStyle(
			'p-value', nltpr,
			background = DT::styleColorBar(
				df[[nltpr]],
				'lightblue'
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		) %>%
		DT::formatStyle(
			'AluWidth',
			background = DT::styleColorBar(
				df[["AluWidth"]],
				'lightgrey'
			),
			backgroundSize = '98% 88%',
			backgroundRepeat = 'no-repeat',
			backgroundPosition = 'center'
		)
}

### UI

#### Header
header <- shinydashboard::dashboardHeader(title = "Alu Viewer")

#### Sidebar
sidebar <- shinydashboard::dashboardSidebar(
	shinydashboard::sidebarMenu(
		shinydashboard::menuItem("Input Data",tabName = "inputdata"),
		shinydashboard::menuItem(
			"main", tabName = "main", icon = icon("dashboard")
		),
		div("Select Features & p-value threshold"),
		div("Click buttons to update visualisations"),
		selectInput(
			multiple = TRUE,
			selected = problemRegions,
			selectize = TRUE,
			label = "Potentially Problematic Region Filters",
			choices = as.list(c("No Filter", problemRegions)),
			inputId = "exclude"
		),
		actionButton("getexclude", "Apply Filters"),
		selectInput(
			multiple = TRUE,
			selected = "No Filter",
			selectize = TRUE,
			label = "Feature",
			choices = as.list(c("No Filter", hg19annos)),
			inputId = "accessor"
		),
		radioButtons(
			"setop",
			"Feature Overlap mode",
			choices = c("or", "and"),
			selected = "or"
		),
		actionButton("getannos", "Intersect with Annotations"),
		numericInput(
			inputId = "pt",
			label = "p-value Cut-off",
			value = 0.05,
			min = 0,
			max = 1
		),
		div(textOutput("formatedp")),
		actionButton("applyp", "Apply P-value Threshold"),
		selectInput(
			#multiple = TRUE,
			selected = "BB",#"genes"
			#selectize = TRUE,
			label = "Model Corrections",
			choices = list(
				"Uncorrected" = "NULL",
				"Blood Cell-type" = "Blood",
				"Batch" = "Batch",
				"Blood & Batch" = "BB"
			),
			inputId = "model"
		),
		sliderInput(
			"CpGdt",
			"CpG Density Range",
			min = 0, max = 100,
			value = c(0, 100)
		),
		actionButton("getCpGdt", "Apply CpG Density Range"),
		sliderInput(
			"size",
			"Alu Size Range",
			min = 1, max = 500,
			value = c(250, 350)
		),
		actionButton("getsize", "Apply Alu Size Range"),
		numericInput(
			"nsig2showDT",
			"Num Top Sig Wins in Tab",
			value = 30
		),
		downloadButton("downloadData", "Download"),
		selectInput(
			"subfamplotmode",
			"Alu Subfamily Plot mode",
			selected = "treemap",
			choices = list(
				"treemap" = "treemap",
				"sunburst" = "sunburst",
				#"packed bubble 1" = "hpackedbubble",
				"packed bubble" = "circlepackeR"
			)
		)
	)
)


#### Body
body <- shinydashboard::dashboardBody(
	shinybusy::add_busy_spinner(spin = "fading-circle"),
	shinydashboard::tabItems(
		shinydashboard::tabItem(
			tabName = "inputdata",
			fluidRow(
				box(
					title = "Data Loader",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					"Select the Alu Age Model Data file (*.fst), switch to the main tab when upload is complete",
					fileInput(
						inputId = "aluagemoddata",
						label="Alu Data",
						multiple = FALSE,
						accept = NULL, width = NULL,
						buttonLabel = "Browse...",
						placeholder = "No file selected"
					)
				)
			)
		),
		shinydashboard::tabItem(
			tabName = "main",

			fluidRow(
				shinydashboard::box(
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					title = "Alu Copy Number by Subfamily",
					"Reflects the number of Alu elements, Not windows over Alu elements",
					width = 7,
					htmlOutput("alusubfamplot")
				),
				shinydashboard::box(
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					title = "Top Age Change Windows",
					width = 5,
					DT::DTOutput("table")
				)
			),
			fluidRow(
				shinydashboard::infoBox(
					"Age Model Data",
					#"Blood Cell-type and Batch corrected, N = 3001",
					width = 3
				),
				shinydashboard::infoBox(
					"Numbers of Changeing Alus",
					#textOutput("alucounts"),
					htmlOutput("alucounts"),
					width = 9
				)
			),
			fluidRow(
				shinydashboard::box(
					title = "Direction of Change",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					plotly::plotlyOutput("slopedensity")
				),
				shinydashboard::box(
					title = "Significance of Change",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					textOutput("AluWideBonfer"),
					#sprintf("Red = %.3e (0.05 / N windows over all Alus used)", AluWideBonfer),
					plotly::plotlyOutput("pvaldensity")
				)
			),
			fluidRow(
				shinydashboard::box(
					title = "Central Tendancy Of Methylation at Changing Loci",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					plotly::plotlyOutput("meandensity")
				),
				shinydashboard::box(
					title = "CpG density",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					plotly::plotlyOutput("CpGdensity")
				)

			),
			fluidRow(
				shinydashboard::box(
					title = "Alu Element Size Distribution",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					plotly::plotlyOutput("alusize")
				),
				shinydashboard::box(
					title = "Alu Strand",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 6,
					plotly::plotlyOutput("alustrand")
				)

			),
			fluidRow(
				shinydashboard::box(
					title = "UCSC Genome Browser - Location of Selected Alu Windows (+/- 1kb)",
					solidHeader = TRUE, collapsible = TRUE, status = "primary",
					width = 12,
					htmlOutput("genomeBrowser")
				)
			)
		)
	)
)

#### UI Wrapper
ui <- shinydashboard::dashboardPage(
	header,
	sidebar,
	body
)

### Server
naluwin <- length(Alus)

server <- function(input, output) {
	options(shiny.maxRequestSize = 3000 * 1024^2)

	annotatr::read_annotations(
		system.file(
			"extdata", "hg19-blacklist.v2.bed.gz", package = "ShinyAluExplorer"
		),
		#"inst/extdata/hg19-blacklist.v2.bed.gz",
		genome = "hg19",
		name = "blacklist_Amemiya2019",
		format = "bed"
	)

	annotatr::read_annotations(
		system.file(
			"extdata", "Supplementary_File1_mergeFinal_hsm_nh.bed",
			package = "ShinyAluExplorer"
		),
		#"inst/extdata/Supplementary_File1_mergeFinal_hsm_nh.bed",
		genome = "hg19",
		name = "hsm-peaks-Bell2017",
		format = "bed"
	)

	annotatr::read_annotations(
		system.file(
			"extdata", "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.bed",
			package = "ShinyAluExplorer"
		),
		genome = "hg19",
		name = "AluSVs",
		format = "bed"
	)

	AlusAgeModelRes <- reactive({
		if (is.null(input$aluagemoddata)){
			return(NULL)
		}
		AlusAgeModelRes <- fst::read_fst(input$aluagemoddata$datapath) %>%
			plyranges::as_granges()
		GenomeInfoDb::seqlevels(
			AlusAgeModelRes, pruning.mode = "coarse"
		) <- GenomeInfoDb::seqlevels(Hsapiens)

		AlusAgeModelRes
	})

	pval <- reactive({paste0(input$model, "_p")})
	slope <- reactive({paste0(input$model, "_slope")})
	neglog10p <- reactive({paste0(input$model, "_neglog10P")})

	pvalthresh <- eventReactive(input$applyp, {
		input$pt
	}, ignoreNULL = FALSE)

	formatedp <- reactive({sprintf("%.3e", input$pt)})
	output$formatedp <- renderText({formatedp()})

	AluWideBonfer <- reactive({0.05/naluwin})

	output$AluWideBonfer <- reactive({
		# naluwin <- length(Alus)
		# AluWideBonfer <- 0.05/naluwin
		paste0(
			sprintf(
				"Red = %.3e (0.05 / ",
				AluWideBonfer()
			),
			format(naluwin, big.mark = ","),
			")"
		)
	})

	AlusEx <- eventReactive(input$getexclude, {
		if("No Filter" %in% input$exclude) {
			Alus ###!!! GLOBAL
		} else {
			myAnno(Alus, input$exclude, "not") ###!!! GLOBAL
		}
	}, ignoreNULL = FALSE)

	filteredAlus <- eventReactive({
		input$getannos
		input$getexclude
		#input$setop
	}, {
		Alus <- AlusEx()
		if("No Filter" %in% input$accessor) {
			return(AluSubFamSplit(Alus))
		} else {
			return(myAnno(Alus, input$accessor, input$setop))
		}
	}, ignoreNULL = FALSE)

	CpGdRange <- eventReactive(input$getCpGdt,{
		input$CpGdt
	}, ignoreNULL = FALSE)

	alusize <- eventReactive(input$getsize, {
		input$size
	}, ignoreNULL = FALSE)

	AluAgeModDat <- reactive({
		slp <- sym(slope())
		p <- sym(pval())
		filteredAlus() %>%
			plyranges::as_granges() %>%
			plyranges::filter_by_overlaps(AlusAgeModelRes(),.) %>%
			as_tibble() %>%
			AluSubFamSplit() %>%
			tidyr::drop_na({{slp}}) %>%
			mutate(
				dir = ifelse(sign({{slp}}) == 1, "Hyper (+)", "Hypo (-)")
			) %>%
			filter({{p}} < pvalthresh()) %>%
			filter(
				(CpGdensity > CpGdRange()[1]) & (CpGdensity < CpGdRange()[2])
			) %>%
			filter(
				(AluWidth > alusize()[1]) & (AluWidth < alusize()[2])
			) %>%
			ungroup() %>%
			arrange({{p}})
	})

	# nSigAluWins <- reactive({
	# 	AluAgeModDat() %>% nrow()
	# })
	#
	# nSigAlus <- reactive({
	# 	AluAgeModDat() %>%
	# 		distinct(aluIndex) %>%
	# 		nrow()
	# })

	output$alucounts <- renderText({
		#reactive({
		nSigAluWins <- AluAgeModDat() %>% nrow()

		dirCountsWin <- AluAgeModDat() %>% group_by(dir) %>% count()
		nSigAluWinsUp <- dirCountsWin %>% filter(dir == "Hyper (+)") %>% pull(n)
		nSigAluWinsDwn <- dirCountsWin %>% filter(dir == "Hypo (-)") %>% pull(n)

		nSigAlus <- AluAgeModDat() %>% distinct(aluIndex) %>% nrow()
		dirCountsAu <- AluAgeModDat() %>%
			distinct(aluIndex, .keep_all = TRUE) %>%
			group_by(dir) %>% count()
		nSigAlusUp <- dirCountsAu %>% filter(dir == "Hyper (+)") %>% pull(n)
		nSigAlusDwn <- dirCountsAu %>% filter(dir == "Hypo (-)") %>% pull(n)

		paste0(
			nSigAluWins %>% format(big.mark = ","),
			" Windows with p-values < ",
			formatedp(),
			" Covering ",
			nSigAlus %>% format(big.mark = ","),
			" Alu Elements. <br/>",
			nSigAluWinsUp %>% format(big.mark = ","),
			" Windows Hypermethylate (+) and ",
			nSigAluWinsDwn %>% format(big.mark = ","),
			" Windows Hypomethylate (-)<br/>",
			"Covering ", nSigAlusUp  %>% format(big.mark = ","),
			" and ", nSigAlusDwn %>% format(big.mark = ","),
			" Alu elements respectively"
		)

	})

	AluCounts <- reactive({
		AluAgeModDat() %>% AluCountByType()
	})

	output$alusubfamplot <- renderUI({
		switch(
			input$subfamplotmode,
			"treemap" = d3treeR::d3tree2Output("aluplot_treemap"),
			"circlepackeR" = circlepackeR::circlepackeROutput("aluplot_circlepackeR"),
			"sunburst" = sunburstR::sunburstOutput("aluplot_sunburstR")#,
			#"hpackedbubble" = hpackedbubble::hpackedbubbleOutput("aluplot_hpackedbubble")
		)
	})

	output$aluplot_treemap <- d3treeR::renderD3tree2({
		AluCounts() %>% AluTreeMap()
	})
	output$aluplot_circlepackeR <- circlepackeR::renderCirclepackeR({
		AluCounts() %>% aluPackedBubblecirc()
	})
	output$aluplot_sunburstR <- sunburstR::renderSunburst({
		AluCounts() %>% aluSunburst()
	})
	# output$aluplot_hpackedbubble <- hpackedbubble::renderHpackedbubble({
	# 	AluCounts() %>% AluPackedBubble()
	# })

	output$genomeBrowser <- renderText({
		idx <- input$table_rows_selected
		if(is.null(idx)) {
			return("Select A Row")
		}
		row <- AluAgeModDat() %>%
			dplyr::slice(idx) %>%
			tail(1) ###!!! may be unnecassary ~ old bug
		buf <- 1000 #  or autozoomout? # &hgt.out4=submit
		glue::glue(
			"<embed width = '100%' height = '800px' src='",
			"http://genome.ucsc.edu/cgi-bin/hgTracks",
			"?db=hg19",
			"&position={row$seqnames}:{row$start - buf}-{row$end + buf}",
			"'>"
		)
	})

	# Age model density plots
	output$slopedensity <- plotly::renderPlotly({
		slp <- sym(slope())
		plot <- AluAgeModDat() %>%
			ggplot(aes({{slp}}, colour = subfamily)) +
			geom_density() +
			geom_vline(xintercept = 0) +
			#coord_equal() +
			scale_color_manual(values = AluColours) +
			labs(
				#title = "Direction of Change",
				x = "Model Slope /au"
			)
		plot %>% plotly::ggplotly(dynamicTicks = TRUE, height = 400)
	})

	output$pvaldensity <- plotly::renderPlotly({
		nltp <- sym(neglog10p())
		plot <- AluAgeModDat() %>%
			ggplot(aes({{nltp}}, colour = subfamily)) +
			geom_density() +
			facet_wrap(~dir) +
			#geom_vline(xintercept = -log10(0.05)) +
			#geom_vline(xintercept = -log10(input$pt)) +
			geom_vline(xintercept = -log10(AluWideBonfer()), colour = "red") +
			scale_color_manual(values = AluColours) +
			labs(
				#title = "Significance of Change",
				x = "-log10(p-value)"
			)
		plot %>% plotly::ggplotly(dynamicTicks = TRUE, height = 400)
	})

	output$meandensity <- plotly::renderPlotly({
		plot <- AluAgeModDat() %>%
			ggplot(aes(median, colour = subfamily)) +
			geom_density() +
			scale_color_manual(values = AluColours) +
			facet_wrap(~dir) +
			labs(
				#title = "Median RPM",
				x = "Median RPM"
			)
		plot %>% plotly::ggplotly(dynamicTicks = TRUE, height = 400)
	})

	output$CpGdensity <- plotly::renderPlotly({
		plot <- AluAgeModDat() %>%
			ggplot(aes(CpGdensity, fill = subfamily)) +
			geom_histogram(
				binwidth = 1, alpha = 0.6, position = "identity"
			) +
			scale_color_manual(values = AluColours) +
			facet_wrap(~dir) +
			labs(
				#title = "CpG density",
				x = "CpG density /%"
			)
		plot %>% plotly::ggplotly(dynamicTicks = TRUE, height = 400)
	})

	output$alusize <- plotly::renderPlotly({
		plot <- AluAgeModDat() %>%
			ggplot(aes(AluWidth, colour = subfamily)) +
			geom_density() +
			scale_color_manual(values = AluColours) +
			facet_wrap(~dir) +
			labs(
				#title = "CpG density",
				x = "Alu Element Length /bp"
				)
		plot %>% plotly::ggplotly(dynamicTicks = TRUE, height = 400)
	})

	output$alustrand <- plotly::renderPlotly({
		plot <- AluAgeModDat() %>%
			ggplot(aes(AluStrand, fill = subfamily)) +
				geom_bar(position = "dodge") +
				scale_fill_manual(values = AluColours) +
				facet_wrap(~dir) +
				labs(
					#title = "CpG density",
					x = "Alu Strand"
				)
		plot %>% plotly::ggplotly(height = 400)
	})

	output$table <- DT::renderDataTable(
		AluAgeModDat() %>%
			AluAgeDT(input$nsig2showDT, pval(), slope(), neglog10p()),
		server = TRUE
	)

	output$downloadData <- downloadHandler(
		filename = "Alu-Age-Model-Data.csv",
		content = function(file) {
			write.csv(AluAgeModDat(), file, row.names = FALSE)
		}
	)
}
