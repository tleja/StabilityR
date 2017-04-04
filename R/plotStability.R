plotStability <- 
function(res, dom) {
	require(rCharts, quietly=TRUE)
	
	tooltip_stability = "#! function () {
		// first line: Time: x value
		// var s = '<b>Time: 10<sup>' + (this.x).toFixed(4) + '</sup></b>';
        var s = '<b>Time: 10<sup>' + (Math.log10(this.x)).toFixed(4) + '</sup></b>';
		s += '<br/><b>Index: ' + (1 + this.points[0].series.data.indexOf(this.points[0].point)) + '</b>'
		
		// second line: Stability and/or Partition y value
		$.each(this.points.reverse(), function () {
			if (this.series.name == 'Stability') {
				var tmp = (this.y).toFixed(3);
			} else {
				var tmp = this.y;
			}
			s += '<br/>' + '<span style=\"color:' + this.series.color + '\">\u25CF</span> ' + this.series.name + ': <b>' + tmp + '</b>';
		})
		return(s)
	} !#"
	
	##### onclick example returning more values: input$stab_click$index ...
	# onclick_stability = "#! function() {
	# 	Shiny.onInputChange('stab_click', {
	# 		index: (1 + this.series.data.indexOf(this)), 
	# 		time: this.x
	# 	})
	# } !#"
	
	onclick_stability = "#! function() {
		Shiny.onInputChange('stab_click', (1 + this.series.data.indexOf(this)) )
	} !#"
	
	# formatter_S = "#! function() { return '1e' + Math.log10(this.value) } !#"
	formatter_S = "#! function() {return '10<sup>' + Math.log10(this.value) + '</sup>'} !#"
	formatter_T = "#! function() {return '10<sup>' + (Math.log10(this.value)).toFixed(1) + '</sup>'} !#"
	
	ind = which(res$N < 2) # stability < 2 is unreliable
	if (length(ind) > 0) {
		res$S[ind] = NaN
	}
			
	maxN = ceiling(log10(max(res$N)))
	minS = floor(log10(min(res$S, na.rm=T)))

	n = mapply(function(x,y) list(c(x,y)), res$time, res$N)
	s = mapply(function(x,y) list(c(x,y)), res$time, res$S)

	h <- Highcharts$new()
	h$chart(alignTicks=FALSE, zoomType='x')
#	h$title(text='Partition Stability')

	h$xAxis(type='logarithmic', 
			minorTickInterval=0.1, 
			crosshair=list(width=4),
			title = list(text='Markov time'),
			labels = list(useHTML=TRUE, formatter = formatter_T)
			)

	h$yAxis(list(
				list(title = list(text='Partition', style=list(color="#7cb5ec")),
					type = 'logarithmic', 
					min = 10^0, max = 10^maxN, 
					minorTickInterval = 0.1, endOnTick = FALSE,
					labels = list(style = list(color="#7cb5ec"))
					), 
				list(title = list(text='Stability', style=list(color = "#f7a35c")), 
					type = 'logarithmic', 
					min = 10^minS, max = 10^0,
					gridLineWidth = 0, 
					minTickInterval = 1, endOnTick = FALSE,
					labels = list(style = list(color = "#f7a35c"), useHTML=TRUE,
									formatter = formatter_S),
					opposite = TRUE
					)
			))

	h$series(name = 'Partition', type = 'line', color = '#7cb5ec',
	                data = n, xAxis=0, yAxis=0, lineWidth=2, index=1,
					marker = list(enabled=T, radius=3, symbol='circle'))
	h$series(name = 'Stability', type = 'line', color = '#f7a35c',
	         		data = s, xAxis=0, yAxis=1, lineWidth=2, index=0,
					marker = list(enabled=T, radius=3, symbol='circle'))

	h$tooltip(shared=TRUE, useHTML=TRUE, borderWidth=0, hideDelay=0,
		formatter = tooltip_stability)

	h$exporting(enabled=TRUE, allowHTML=TRUE)
	
	h$plotOptions(series = list(
					states = list(
						hover = list(enabled=TRUE, lineWidth=2))
				),
				line = list(
					cursor = "pointer", 
					point = list(
						events = list(
							click = onclick_stability))
				))
							
	if (!missing(dom)) h$addParams(dom=dom)

	return(h)
}
