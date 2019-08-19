/* global d3, dat */
var 
	margin = 60,
	width = 500,
	height = 500,
	radiusMin = 2,
	radiusMax = 20,
	opacityMin = 0,
	opacityMax = 1, container = d3.select("#d3Container"),
	chart = container.append("svg").
		attrs({
			id: 'chart',
			width: width + margin + margin,
			height: height + margin + margin
		}),
	plot = chart.append("g").
		attrs({
			id: 'plot',
			transform: "translate(" + margin + "," + margin + ")",
			width: width,
			height: height
		}),
	x = d3.scaleLinear().range([0, width]),
	y = d3.scaleLinear().range([height, 0]),
	z = d3.scaleLinear().range([opacityMin, opacityMax]),
	r = d3.scaleLinear().range([radiusMin, radiusMax]);


dat.forEach(function(d){ 
	d["hits"] = +d["hits"];
	d["max_af"] = +d["max_af"];
	d["max_beta"] = +d["max_beta"];
	d["max_nlog10p"] = +d["max_nlog10p"];
	d["mean_af"] = +d["mean_af"];
	d["mean_beta"] = +d["mean_beta"];
	d["mean_nlog10p"] = +d["mean_nlog10p"];
	d["min_af"] = +d["min_af"];
	d["min_beta"] = +d["min_beta"];
	d["min_nlog10p"] = +d["min_nlog10p"];
	d["var_af"] = +d["var_af"];
	d["var_beta"] = +d["var_beta"];
	d["var_nlog10p"] = +d["var_nlog10p"];
});
// 	d.af = +d.af;
// 	d["filter-pvalue"] = +d["filter-pvalue"];
// 	d["lrt-pvalue"] = +d["lrt-pvalue"];
// 	d["beta"] = +d["beta"];
// 	d["beta-std-err"] = +d["beta-std-err"];
// 	d["varant_h2"] = +d["varant_h2"];
// 	});

console.log(dat);

var pMin = d3.min(dat.map((o) => o["max_nlog10p"]));
var pMax = d3.max(dat.map((o) => o["max_nlog10p"]));
var betaMin = d3.min(dat.map((o) => o.mean_beta));
var betaMax = d3.max(dat.map((o) => o.mean_beta));
var hitsMin = d3.min(dat.map((o) => o.hits));
var hitsMax = d3.max(dat.map((o) => o.hits));
var mafMin = d3.min(dat.map((o) => o.max_af));
var mafMax = d3.max(dat.map((o) => o.max_af));

x.domain([betaMin, betaMax]);
y.domain([pMin, pMax + 0.5]);
r.domain([hitsMin, hitsMax]);
z.domain([mafMin, mafMax]);

plot.selectAll("circle").data(dat).
	enter().
	append("circle").
	attrs({
		"cx": function(d){ return x(d.mean_beta); },
		"cy": function(d){ return y(d.max_nlog10p); },
		"r": function(d){ return r(d.hits); },
		"fill" : "blue",
		"opacity" : function(d){ return 0.2 }
	});

chart.append("g").
	attrs({transform: "translate(" + margin + "," + (height + margin) + ")"}).
	call(d3.axisBottom(x));

chart.append("g").
	attrs({transform: "translate(" + margin + "," + margin + ")"}).
	call(d3.axisLeft(y));
