plotStability_VI <- 
function(res, dom) {
	require(rCharts, quietly=TRUE)
	
	tooltip_stability = "#! function () {
		// first line: Time: x value
		// var s = '<b>Time: 10<sup>' + (this.x).toFixed(4) + '</sup></b>';
        var s = '<b>Time: 10<sup>' + (Math.log10(this.x)).toFixed(4) + '</sup></b>';
		s += '<br/><b>Index: ' + (1 + this.points[0].series.data.indexOf(this.points[0].point)) + '</b>'
		
		// second line: Stability and/or Partition y value
		$.each(this.points, function () {
			if (this.series.name == 'VI') {
				var tmp = (this.y).toFixed(3);
			} else {
				var tmp = this.y;
			}
			s += '<br/>' + '<span style=\"color:' + this.series.color + '\">\u25CF</span> ' + this.series.name + ': <b>' + tmp + '</b>';
		})
		return(s)
	} !#"
	
	onclick_stability = "#! function() {
		Shiny.onInputChange('stab_click', (1 + this.series.data.indexOf(this)) )
	} !#"
	
	formatter_T = "#! function() {return '10<sup>' + (Math.log10(this.value)).toFixed(1) + '</sup>'} !#"
			
	maxN = ceiling(log10(max(res$N)))
	maxVI = ceiling(max(res$VI)/0.05)*0.05

	n = mapply(function(x,y) list(c(x,y)), res$time, res$N)
	v = mapply(function(x,y) list(c(x,y)), res$time, res$VI)

	h <- Highcharts$new()

	h$chart(alignTicks=FALSE, zoomType='x')

	h$xAxis(type='logarithmic', 
			minorTickInterval=0.1, 
			crosshair=list(width=4),
			title = list(text='Markov time'),
			labels = list(useHTML=TRUE, formatter = formatter_T)
			)

	h$yAxis(list(
				list(title = list(text='Partition', style=list(color="white")),
					type = 'logarithmic', 
					min = 10^0, max = 10^maxN, 
					gridLineWidth = 0,
					# minorTickInterval = 1, 
					endOnTick = FALSE,
					labels = list(style = list(color = "white")),
					opposite = TRUE
					), 
				list(title = list(text='Variation of Information', style=list(color = "#69C18C")),
					min = 0, max = maxVI,
					# minTickInterval = 0.01, 
					endOnTick = FALSE,
					labels = list(style = list(color = "#69C18C"))
					)
			))

	h$series(name = 'Partition', type = 'line', color = '#C7C7C7',
	                data = n, xAxis=0, yAxis=0, lineWidth=2, index=0,
					marker = list(enabled=T, radius=3, symbol='circle'))
	h$series(name = 'VI', type = 'line', color = '#69C18C',
	         		data = v, xAxis=0, yAxis=1, lineWidth=2, index=1,
					marker = list(enabled=T, radius=3, symbol='circle'))

	h$tooltip(shared=TRUE, useHTML=TRUE, borderWidth=0, hideDelay=0,
		formatter = tooltip_stability)
		
	h$legend(reversed=TRUE)

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
