plotNetwork <- 
function(mat, fc, pos, ann, group, type,  simpleTooltip=FALSE, ...) {
	require(visNetwork, quietly=TRUE)
	
	onclick_network = "function(properties) { window.open(this.body.data.nodes.get(properties.nodes[0]).url) }"
	
	if (!missing(ann)) if (nrow(mat)!=nrow(ann)) stop("nrow(mat) != nrow(ann)") 
	if (!missing(pos)) if (nrow(mat)!=nrow(pos)) stop("nrow(mat) != nrow(pos)")
	if (!missing(group)) if (nrow(mat)!=length(group)) stop("nrow(mat) != length(group)")
	if (!missing(type)) if (nrow(mat)!=length(type)) stop("nrow(mat) != length(type)")
	
	# set the minimum label if rownames are missing
	if (!is.null(rownames(mat))) label = rownames(mat)
	else label = as.character(1:nrow(mat))
	
	# remove mat names for edgeList
	mat = unname(mat)

	# prepare nodes and edges data frames (input to visNetwork)
	nodes = data.frame(id = 1:nrow(mat), stringsAsFactors=F)
	edges = edgeList(mat, verbose=FALSE)
    colnames(edges) = c("from","to",'value')

	# assign additional node attributes
	if (!missing(ann)) {
		nodes$label = ann$SYMBOL
		nodes$url = paste0('http://www.ncbi.nlm.nih.gov/gene/?term=',ann$ENTREZID)
		click = onclick_network
	}
	else {
		nodes$label = label
		click = NULL
	}
	
	if (!missing(type)) {
		nodes$type = type
		sel = 'type'
	}
	else sel = NULL
	
	if (!missing(group)) nodes$group = group
	
	if (!missing(fc)) {
		if (min(fc) >= 1) {
			colors = scales::col_numeric("RdBu", 1:max(fc))
			nodes$color = colors(fc)
		}
		if (min(fc) < 0) {
			tmp = max(ceiling(abs(c(max(fc),min(fc)))/0.1)*0.1)
			colors = scales::col_numeric("RdBu", c(-tmp, tmp))
			nodes$color = colors(-1*fc)
		}
	}
	
	# design tooltip
	if (simpleTooltip) nodes$title = nodes$label
	else if(!missing(type)) nodes$title = sprintf("<b>%s</b><br>%s", nodes$label, type)
	else nodes$title = nodes$label
	
	# set the position of the nodes
	if (missing(pos)) pos = forceAtlas2(mat, ...)
	nodes = data.frame(nodes, x=pos[,1], y=pos[,2])
	
	# order the nodes based on the label for nodesIdSelection
	nodes = nodes[order(nodes$label),]
	
	# plot the network
    net = visNetwork(nodes, edges) %>% 
			visNodes(fixed = TRUE) %>%
      		visEdges(scaling = list(min=1, max=2), smooth = list(enabled=TRUE, type="diagonalCross")) %>%
      		visInteraction(hover = FALSE, tooltipDelay=0) %>%
      		visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = sel) %>%
      		visLegend(width=0.03, position='right') %>%
      		visEvents(doubleClick = click)

	return(net)
}
