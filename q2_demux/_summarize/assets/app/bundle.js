var app=webpackJsonpapp([0],[function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}Object.defineProperty(e,"__esModule",{value:!0}),e.init=void 0;var a=r(1),i=n(a);e.init=function(t,e){var r=arguments.length>2&&void 0!==arguments[2]?arguments[2]:void 0,n={forward:Object.keys(e).map(function(t){return[+t+1,e[t]]})};r&&(n.reverse=Object.keys(r).map(function(t){return[+t+1,r[t]]})),(0,i.default)(n,t)}},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}Object.defineProperty(e,"__esModule",{value:!0});var i=r(2),o=a(i),l=r(3),s=n(l),u=r(4),c=function(t,e,r,n){var a=o.select(r),i=a.append("svg").attr("class","col-xs-12").style("display","block").style("margin","0 auto").attr("width",e.width+e.margin.left+e.margin.right).attr("height",e.height+e.margin.top+e.margin.bottom);a.append("div").attr("class","col-xs-12").append("p").attr("class","random-sampling").html("These plots were generated using a random sampling of "+n.n+"\n             out of "+n.totalSeqCount+" sequences without replacement. The\n             minimum sequence length identified during subsampling was\n             "+n.minSeqLen[t.direction]+" bases. Outlier quality scores are not shown in box\n             plots for clarity.");var l=a.append("div").attr("class","col-xs-12").append("div").attr("class","panel panel-default");l.append("div").attr("class","panel-heading").html("Parametric seven-number summary");var c=l.append("div").attr("class","stats").append("table").attr("class","table").style("margin-bottom","0");c.append("thead").append("tr").selectAll("th").data([["Box plot feature",5],["Percentile",5],["Quality score",2]]).enter().append("th").text(function(t){return t[0]}).attr("class",function(t){return"col-xs-"+t[1]});var d=[["(Not shown in box plot)","2nd","..."],["Lower Whisker","9th","..."],["Bottom of Box","25th","..."],["Middle of Box","50th (Median)","..."],["Top of Box","75th","..."],["Upper Whisker","91st","..."],["(Not shown in box plot)","98th","..."]],f=c.append("tbody").selectAll("tr").data(d).enter().append("tr");f.selectAll("td").data(function(t){return t}).enter().append("td").text(function(t){return t});var p=o.max(t,function(t){return t[0]})+1,h=[0,p],m=[0,45],g=o.scaleLinear().domain(h).range([e.margin.left,e.width]),v=o.scaleLinear().domain(m).range([e.height-e.margin.bottom,e.margin.top]),x=[12,9],b=o.axisBottom(g).ticks(x[0],o.format("d")),y=o.axisLeft(v).ticks(x[1],o.format("d"));i.attr("height",e.height+e.margin.bottom+e.margin.top),i.attr("width",g.range()[0]+g.range()[1]),(0,u.addBrush)(i,t,g,v,h,m,b,y,x,n),(0,s.default)(i,t,g,v,n),i.append("rect").attr("width",e.width+e.margin.left+e.margin.right).attr("height",4*e.margin.bottom).attr("y",e.height-e.margin.bottom).attr("fill","white"),i.append("rect").attr("width",4*e.margin.left).attr("x",3*-e.margin.left).attr("height",e.height+200).attr("fill","white"),i.append("rect").attr("width",e.margin.right+1e3).attr("height",e.height+200).attr("x",e.width).attr("fill","white"),i.append("g").attr("class","axis axis--x").attr("transform","translate(0, "+(e.height-e.margin.bottom)+")").call(b),i.append("g").attr("class","axis axis--y").attr("transform","translate("+e.margin.left+",0)").call(y),i.append("text").attr("transform","rotate(-90)").attr("x",0-e.height/2).attr("dy","0em").attr("font-size","12px").style("text-anchor","middle").text("Quality Score"),i.append("text").attr("x",e.width/2).attr("y",e.height).attr("dy","1em").attr("font-size","12px").style("text-anchor","middle").text("Sequence Base")},d=function(t,e){var r={top:10,right:30,bottom:30,left:30},n=o.select("#forwardContainer").node().offsetWidth,a={margin:r,width:n-r.left-r.right,height:9*n/16-r.top-r.bottom};Object.keys(t).forEach(function(r){t[r].direction=r,c(t[r],a,"#"+r+"Container",e)})};e.default=d},,function(t,e,r){"use strict";function n(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}function a(t,e,r,n,a){var i=(r.range()[1]-r.range()[0])/(r.domain()[1]-r.domain()[0])/2,l=i/2,s=t.transition().duration(750),u="steelblue",c="skyblue",d="#a94442",f="#ebccd1",p=a.minSeqLen[e.direction],h=t.selectAll(".container").data(e);h.exit().remove();var m=h.enter().append("g").attr("class","container"),g=h.merge(m).transition(s).attr("transform",function(t){return"translate("+r(t[0])+", 0)"}).selection().on("mouseover",function(){var t=o.select(this).data(),e=t[0][0],r=t[0][1],n=r.count<a.n,i=o.select(this.parentNode).node(),l=o.select(i.parentNode);l.select(".panel").attr("class",n?"panel panel-danger":"panel panel-default"),l.select(".panel-heading").html("Parametric seven-number summary for <strong>position "+e+"</strong>"),l.select(".stats").select("tbody").selectAll("tr").data([["(Not shown in box plot)","2nd",r["2%"]],["Lower Whisker","9th",r["9%"]],["Bottom of Box","25th",r["25%"]],["Middle of Box","50th (Median)",r["50%"]],["Top of Box","75th",r["75%"]],["Upper Whisker","91st",r["91%"]],["(Not shown in box plot)","98th",r["98%"]]]).selectAll("td").data(function(t){return t}).text(function(t){return t});var s="The minimum sequence length identified during subsampling was "+p+" bases";n&&(s="This position ("+e+") is greater than the minimum sequence length observed\n                      during subsampling ("+p+" bases). As a result, the plot at this position\n                      is not based on data from all of the sequences, so it should be interpreted with \n                      caution when compared to plots for other positions"),l.select(".random-sampling").classed("text-danger",n).html("The plot at position "+e+" was generated using a random\n               sampling of "+r.count+" out of "+a.totalSeqCount+" sequences\n               without replacement. "+s+". Outlier quality scores are\n               not shown in box plots for clarity.")}),v=g.selectAll("line.center").data(function(t){return[t]});v.exit().remove();var x=v.enter().append("line");v.merge(x).transition(s).attr("class","center").attr("x1",0).attr("y1",function(t){return n(t[1]["9%"])}).attr("x2",0).attr("y2",function(t){return n(t[1]["91%"])}).attr("stroke-dasharray","2,2").attr("stroke-width",1).attr("stroke","black");var b=g.selectAll("rect.box").data(function(t){return[t]});b.exit().remove();var y=b.enter().append("rect");b.merge(y).transition(s).attr("class","box").attr("x",-l).attr("y",function(t){return n(t[1]["75%"])}).attr("width",i).attr("height",function(t){return n(t[1]["25%"])-n(t[1]["75%"])}).attr("fill",function(t){return t[1].count<a.n?f:c}).attr("stroke-width",1).attr("stroke","black").selection().on("mouseover",function(){o.select(this).attr("fill",function(t){return t[1].count<a.n?d:u})}).on("mouseout",function(){o.select(this).attr("fill",function(t){return t[1].count<a.n?f:c})});var w=g.selectAll("line.median").data(function(t){return[t]});w.exit().remove();var k=w.enter().append("line");w.merge(k).transition(s).attr("class","median").attr("x1",-l).attr("y1",function(t){return n(t[1]["50%"])}).attr("x2",l).attr("y2",function(t){return n(t[1]["50%"])}).attr("stroke-width",1).attr("stroke","black");var _=g.selectAll("line.lower-whisker").data(function(t){return[t]});_.exit().remove();var M=_.enter().append("line");_.merge(M).transition(s).attr("class","lower-whisker").attr("x1",-l).attr("y1",function(t){return n(t[1]["9%"])}).attr("x2",l).attr("y2",function(t){return n(t[1]["9%"])}).attr("stroke-width",1).attr("stroke","black");var O=g.selectAll("line.upper-whisker").data(function(t){return[t]});O.exit().remove();var A=O.enter().append("line");O.merge(A).transition(s).attr("class","upper-whisker").attr("x1",-l).attr("y1",function(t){return n(t[1]["91%"])}).attr("x2",l).attr("y2",function(t){return n(t[1]["91%"])}).attr("stroke-width",1).attr("stroke","black")}Object.defineProperty(e,"__esModule",{value:!0}),e.default=a;var i=r(2),o=n(i)},function(t,e,r){"use strict";function n(t){return t&&t.__esModule?t:{default:t}}function a(t){if(t&&t.__esModule)return t;var e={};if(null!=t)for(var r in t)Object.prototype.hasOwnProperty.call(t,r)&&(e[r]=t[r]);return e.default=t,e}function i(t,e,r,n,a,i,l,u,d,f){var p=s.brush(),h=void 0,m=350,g=(t.selectAll(".axis--x .tick"),function(){h=null}),v=function(){var a=t.transition().duration(750),i=o(d,2),p=i[0],h=i[1],m=Math.floor(r.domain()[0]),g=Math.ceil(r.domain().slice(-1)[0]),v=s.range(m,g,1);v.length<=p?l.tickValues(v).tickFormat(s.format("d")):l.tickValues(null).ticks(p,s.format("d"));var x=Math.floor(n.domain()[0]),b=Math.ceil(n.domain().slice(-1)[0]),y=s.range(x,b,1);y.length<=h?u.tickValues(y).tickFormat(s.format("d")):u.tickValues(null).ticks(h,s.format("d")),t.select(".axis--x").transition(a).call(l),t.select(".axis--y").transition(a).call(u),(0,c.default)(t,e,r,n,f)},x=function(){var e=s.event.selection;if(e)r.domain([e[0][0],e[1][0]].map(r.invert,r)),n.domain([e[1][1],e[0][1]].map(n.invert,n)),t.select(".brush").call(p.move,null);else{if(!h)return h=setTimeout(g,m);r.domain(a),n.domain(i)}return v()};p.on("end",x),t.append("g").attr("class","brush").call(p)}Object.defineProperty(e,"__esModule",{value:!0});var o=function(){function t(t,e){var r=[],n=!0,a=!1,i=void 0;try{for(var o,l=t[Symbol.iterator]();!(n=(o=l.next()).done)&&(r.push(o.value),!e||r.length!==e);n=!0);}catch(t){a=!0,i=t}finally{try{!n&&l.return&&l.return()}finally{if(a)throw i}}return r}return function(e,r){if(Array.isArray(e))return e;if(Symbol.iterator in Object(e))return t(e,r);throw new TypeError("Invalid attempt to destructure non-iterable instance")}}();e.addBrush=i;var l=r(2),s=a(l),u=r(3),c=n(u)}]);