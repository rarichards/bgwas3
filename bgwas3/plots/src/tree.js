/* globals d3, newick, dat*/

// functions {{{
function changed() {
	clearTimeout(timeout);
	var t = d3.transition().duration(750);
	linkExtension.transition(t).attr("d", this.checked ? linkExtensionVariable : linkExtensionConstant);
	link.transition(t).attr("d", this.checked ? linkVariable : linkConstant);
}

function mouseovered(active) {
	return function(d) {
		d3.select(this).classed("label--active", active);
		d3.select(d.linkExtensionNode).classed("link-extension--active", active).each(moveToFront);
		do d3.select(d.linkNode).classed("link--active", active).each(moveToFront); while (d = d.parent);
	};
}

function moveToFront() {
	this.parentNode.appendChild(this);
}

function maxLength(d) {
	return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0);
}

function setRadius(d, y0, k) {
	d.radius = (y0 += d.data.length) * k;
	if (d.children) d.children.forEach(function(d) { setRadius(d, y0, k); });
}

function linkVariable(d) {
	return linkStep(d.source.x, d.source.radius, d.target.x, d.target.radius);
}

function linkConstant(d) {
	return linkStep(d.source.x, d.source.y, d.target.x, d.target.y);
}

function linkExtensionVariable(d) {
	return linkStep(d.target.x, d.target.radius, d.target.x, innerRadius);
}

function linkExtensionConstant(d) {
	return linkStep(d.target.x, d.target.y, d.target.x, innerRadius);
}

function linkStep(startAngle, startRadius, endAngle, endRadius) {
	var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
		s0 = Math.sin(startAngle),
		c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
		s1 = Math.sin(endAngle);
	return "M" + startRadius * c0 + "," + startRadius * s0
		+ (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
		+ "L" + endRadius * c1 + "," + endRadius * s1;
}

function setLinkColorUp(d, color) {
	d.color = color;
	if(d.parent){setLinkColorUp(d.parent, color);}
}

function setLinkColorDown(d, color) {
	d.color = color;
	if(d.children){d.children.forEach(function(o){setLinkColorDown(o, color);});}
}

// }}}

// init
var 
	width = 960,
	outerRadius = width / 2,
	innerRadius = outerRadius - 60,
	cluster = d3.cluster().
		size([360, innerRadius]).
		separation(function(){return 1;}),
	container = d3.select("body").append("div"),
	chart = container.append("svg").
		attr("preserveAspectRatio", "xMinYMin meet").
		attr("viewBox", "0 0 " + width + " " + width + ")").
		classed("svg-content", true),
	plot = chart.append("g").
		attr("transform", "translate(" + outerRadius + "," +  outerRadius + ")"),
	root = d3.hierarchy(newick, function(d){return d.branchset;}).
		sum(function(d){return d.branchset ? 0 : 1;}).
		sort(function(a, b){return (a.value - b.value) || d3.ascending(a.data.length, b.data.length);}),
	circle,
	link;

cluster(root);
setLinkColorDown(root, "black");
setRadius(root, root.data.length = 0, innerRadius / maxLength(root));

plot.append("g").
	attr("class", "link-extensions").
	selectAll("path").
	data(root.links().filter(function(d){return !d.target.children;})).
	enter().
	append("path").
	each(function(d){d.target.linkExtensionNode = this;}).
	attr("d", linkExtensionConstant);

plot.append("g").
	attr("class", "labels").
	selectAll("text").
	data(root.leaves()).
	enter().
	append("text").
	attr("dy", ".31em").
	attr("transform", function(d){return "rotate(" + (d.x - 90) + ")translate(" + (innerRadius + 12) + ",0)" + (d.x < 180 ? "" : "rotate(180)");}).
	attr("text-anchor", function(d){return d.x < 180 ? "start" : "end";}).
	text(function(d){return d.data.name;}).
	on("mouseover", mouseovered(true)).
	on("mouseout", mouseovered(false));

link = plot.append("g").
	attr("class", "links");

circle = plot.append("g").
	attr("class", "circles");

function updateLink(){
	link.selectAll("path").remove();
	link.selectAll("path").
		data(root.links()).
		enter().
		append("path").
		each(function(d){d.target.linkNode = this;}).
		attr("d", linkConstant).
		attr("stroke", function(d){return d.target.color;});
}

function updateTree(attr){

	var 
		attrDat = dat.
			map((o) => o[attr]),
		min = d3.min(attrDat),
		max = d3.max(attrDat),
		z = d3.scaleLinear().
			domain([min, max]).
			range([0, 1]);

	circle.selectAll("circle").remove();
	circle.selectAll("circle").
		data(root.leaves()).
		enter().
		append("circle").
		attrs({
			dy: ".31em",
			transform: function(d){return "rotate(" + (d.x - 90) + ")translate(" + (innerRadius) + ",0)" + (d.x < 180 ? "" : "rotate(180)");},
			r: 6,
			fill: function(d){return d3.interpolateBlues(z(dat.filter((o) => o.sample === d.data.name).map((o) => o[attr])[0]));},
			stroke: "black",
			"stroke-width": 1
		});

	setLinkColorDown(root, "black");
	updateLink();
}
