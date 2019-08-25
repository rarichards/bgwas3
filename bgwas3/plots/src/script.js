/* global d3, dat, $ */
// init {{{
$(document).ready( function () {
	$('#table-all').DataTable({
		"paging": false
	});
});

var 
	modal = document.getElementById("modal"),
	closeModal = document.getElementsByClassName("close")[0],
	dataTableInfo;	

closeModal.onclick = function() {
	modal.style.display = "none";
	dataTableInfo.destroy();
}
var
	margin = 60,
	width = 600,
	height = 600,
	widthGenes = width - 2*margin,
	heightGenes = height - 2*margin,
	radiusMin = 2,
	radiusMax = 20,
	opacityMin = 0,
	opacityMax = 1, 
	container = d3.select("#d3-container"),
	chart = container.append("svg").
		attrs({
			id: 'chart',
			width: width,
			height: height
		}),
	plotGenes = chart.append("g").
		attrs({
			id: 'plotGenes',
			transform: "translate(" + margin + "," + margin + ")",
			width: widthGenes,
			height: heightGenes
		}),
	tableAll = d3.select("#table-all"),
	tableInfo = d3.select("#table-info tbody"),
	x = d3.scaleLinear().range([0, widthGenes]),
	y = d3.scaleLinear().range([heightGenes, 0]),
	z = d3.scaleLinear().range([opacityMin, opacityMax]),
	r = d3.scaleLinear().range([radiusMin, radiusMax]);

// }}}
// process data {{{
dat.forEach(function(d){ 
	d.hits = +d.hits,
	d.min_beta = +d.min_beta,
	d.max_beta = +d.max_beta,
	d.mean_beta = +d.mean_beta,
	d.var_beta = +d.var_beta,
	d.min_nlog10p = +d.min_nlog10p,
	d.max_nlog10p = +d.max_nlog10p,
	d.mean_nlog10p = +d.mean_nlog10p,
	d.var_nlog10p = +d.var_nlog10p,
	d.min_maf = +d.min_maf,
	d.max_maf = +d.max_maf,
	d.mean_maf = +d.mean_maf,
	d.var_maf = +d.var_maf,
	d.info = d.info.split(";").map(function(o){
		var pair = o.split("=");
		return {"key": pair[0], "value": pair[1]};
	});
});

var dat2 = dat.map((o) => ({
		"gene": o.gene,
		"hits": o.hits, 
		"mean_beta": o.mean_beta,
		"max_nlog10p": o.max_nlog10p,
		"mean_maf": o.mean_maf
	}));

// }}}
// update scales {{{
var pMaxMin = d3.min(dat.map((o) => o.max_nlog10p));
var pMaxMax = d3.max(dat.map((o) => o.max_nlog10p));
var pMaxPadding = 0.1*(pMaxMax-pMaxMin);
var betaMeanMin = d3.min(dat.map((o) => o.mean_beta));
var betaMeanMax = d3.max(dat.map((o) => o.mean_beta));
var betaMeanPadding = 0.1*(betaMeanMax-betaMeanMin);
var hitsMin = d3.min(dat.map((o) => o.hits));
var hitsMax = d3.max(dat.map((o) => o.hits));
var mafMeanMin = d3.min(dat.map((o) => o.mean_maf));
var mafMeanMax = d3.max(dat.map((o) => o.mean_maf));
x.domain([betaMeanMin - betaMeanPadding, betaMeanMax + betaMeanPadding]);
y.domain([pMaxMin - pMaxPadding, pMaxMax + pMaxPadding]);
r.domain([hitsMin, hitsMax]);
z.domain([mafMeanMin, mafMeanMax]);

// }}}
// plot {{{
var tableAllHead = tableAll.append("thead").append("tr");
tableAllHead.selectAll("th").data(Object.keys(dat2[0])).enter().append("th").text(function(d){ return d; });
var tableAllBody = tableAll.append("tbody");
var tableAllTr = tableAllBody.selectAll("tr").data(dat2).enter().append("tr");
tableAllTr.selectAll("td").data(function(d){ return Object.values(d); }).enter().append("td").text(function(d){return d;});

plotGenes.selectAll("circle").data(dat).
	enter().
	append("circle").
	attrs({
		"cx": function(d){ return x(d.mean_beta); },
		"cy": function(d){ return y(d.max_nlog10p); },
		"r": function(d){ return r(d.hits); },
		"fill" : "red",
		"opacity" : function(){ return 0.2 }
	}).
	on("mouseover", function(d){
		console.log(d);
		plotGenes.
			append("text").
			attrs({
				"x": x(d.mean_beta),
				"y": y(d.max_nlog10p) - r(d.hits) - 5,
				"font-size": "16px",
				"fill": "blue",
				"width": 100,
				"height": 100
			}).
			text(d.gene);
	}).
	on("mouseout", function(){
		plotGenes.selectAll("text").remove();
	}).
	on("click", function(d){
		tableInfo.selectAll("tr").remove();
		var tr = tableInfo.selectAll("tr").data(d.info).enter().append("tr");
		tr.selectAll("td").data(function(d){return Object.values(d);}).enter().append("td").
			text(function(d){return d;});
		modal.style.display = "block";
		// datatableinfo.destroy();
		dataTableInfo = $('#table-info').DataTable({
			"paging": false
		});
	});


// plotGenes.selectAll("text").data(dat).
// 	enter().
// 	append("text").
// 	attrs({
// 		"x": function(d){ return x(d.mean_beta); },
// 		"y": function(d){ return y(d.max_nlog10p); },
// 		"font-size": "16px",
// 		"fill": "blue",
// 		"width": 100,
// 		"height": 100
// 	}).
// 	text(function(d){ return d.gene;	});

// }}}
// axis {{{
chart.append("g").
	attrs({
		id: "xAxis",
		transform: "translate(" + margin + "," + (margin + heightGenes) + ")"
	}).
	call(d3.axisBottom(x));

chart.append("text").
	attr("transform", "translate(" + (margin + widthGenes/2) + " ," + (margin + heightGenes + margin/2) + ")").
	style("text-anchor", "middle").
	text("Average beta");

chart.append("g").
	attrs({
		id: "yAxis",
		transform: "translate(" + margin + "," + margin + ")"
	}).
	call(d3.axisLeft(y));

chart.append("text").
	attrs({
		"transform": "translate(" + (margin/2) + ", " + (margin + heightGenes/2) + ")rotate(-90)",
		"text-anchor": "middle"
	}).
	text("Maximum -log10p");

// }}}
