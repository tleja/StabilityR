plotStability_heatmap <- 
function(res, filepath, dom) {
	require(rCharts, quietly=TRUE)
	require(scales, quietly=TRUE)
	
	tooltip_stability = "#! function () {
		// first line: Time: x value
		// var s = '<b>Time: 10<sup>' + (this.x).toFixed(4) + '</sup></b>';
        var s = '<b>Time: 10<sup>' + (Math.log10(this.x)).toFixed(4) + '</sup></b>';
		s += '<br/><b>Index: ' + (1 + this.points[1].series.data.indexOf(this.points[1].point)) + '</b>'
		s += '<br/>' + '<span style=\"color:' + this.points[1].series.color + '\">\u25CF</span> ' + this.points[1].series.name + ': <b>' + this.points[1].y + '</b>';

		return(s)
	} !#"
	
	onclick_stability = "#! function() {
		Shiny.onInputChange('stab_click', (1 + this.series.data.indexOf(this)) )
	} !#"
	
	# formatter_S = "#! function() { return '1e' + Math.log10(this.value) } !#"
	formatter_S = "#! function() {return '10<sup>' + Math.log10(this.value) + '</sup>'} !#"
	formatter_T = "#! function() {return '10<sup>' + (Math.log10(this.value)).toFixed(1) + '</sup>'} !#"
	
	if (missing(dom)) filepath = tempfile()
	xm = ncol(res$P)/100
	xm = ifelse(xm > 1, 0.2, 0.2*pbeta(xm, 50,10))
	
	tmp = VI(res$P, 'all')
    colors = scales::col_numeric("RdYlBu", 1:1000)
    png(filepath)
    par(mar=c(0,xm,0,xm))
    image(tmp[,ncol(tmp):1], col=colors(1:1000), useRaster=TRUE, xaxt='n', yaxt='n', bty='n')
    dev.off()
	
	ind = which(res$N < 2) # stability < 2 is unreliable
	if (length(ind) > 0) {
		res$S[ind] = NaN
	}
			
	maxN = ceiling(log10(max(res$N)))
	minS = floor(log10(min(res$S, na.rm=T)))

	n = mapply(function(x,y) list(c(x,y)), 1:length(res$N), res$N)
	s = mapply(function(x,y) list(c(x,y)), 1:length(res$N), res$S)

	h <- Highcharts$new()
	if (missing(dom)) h$chart(alignTicks=FALSE, plotBackgroundImage=filepath)
	else h$chart(alignTicks=FALSE, plotBackgroundImage=basename(filepath))

	h$xAxis(type='linear',
			minorTickInterval=0.1, 
			minorGridLineWidth=0,
			crosshair=list(width=4, color='rgba(155,200,255,0.9)'),
			title = list(text='Partition index'),
			labels = list(useHTML=TRUE)
			)

	h$yAxis(list(
				list(title = list(text='Partition', style=list(color="#7cb5ec")),
					type = 'logarithmic', 
					min = 10^0, max = 10^maxN, 
					minorGridLineWidth = 0,
					gridLineWidth = 0,
					minorTickInterval = 0.1, endOnTick = FALSE,
					labels = list(style = list(color="#7cb5ec"))
					), 
				list(title = list(text='Stability', style=list(color = 'white')), 
					type = 'logarithmic', 
					min = 10^minS, max = 10^0,
					minorGridLineWidth = 0,
					gridLineWidth = 0, 
					minTickInterval = 1, endOnTick = FALSE,
					labels = list(style = list(color = "white"), useHTML=TRUE,
									formatter = formatter_S),
					opposite = TRUE
					)
			))

	h$series(name = 'Partition', type = 'line', color = '#7cb5ec',
	                data = n, xAxis=0, yAxis=0, lineWidth=2, index=1,
					marker = list(enabled=T, radius=3, symbol='circle'),
					cursor = "pointer", 
					point = list(events = list(click = onclick_stability)))
	h$series(name = 'Stability', type = 'line', color = '#f7a35c',
	         		data = s, xAxis=0, yAxis=1, lineWidth=0, index=0, showInLegend=FALSE,
					marker = list(enabled=F, states=list(hover=list(enabled=FALSE))))

	h$tooltip(shared=TRUE, useHTML=TRUE, borderWidth=0, hideDelay=0,
		formatter = tooltip_stability)

	h$exporting(enabled=TRUE, allowHTML=TRUE)
	
	h$plotOptions(series=list(
					states = list(
						hover = list(enabled=TRUE, lineWidth=2, halo=list(size=0))),
	 				events=list(
	 					legendItemClick="#! function() {return false;} !#")
					))
	
	if (!missing(dom)) h$addParams(dom=dom)
	
	return(h)
}
