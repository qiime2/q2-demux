var app=webpackJsonpapp([0],[function(t,r,n){"use strict";function e(t){return t&&t.__esModule?t:{default:t}}Object.defineProperty(r,"__esModule",{value:!0}),r.init=void 0;var a=n(1),i=e(a);r.init=function(t){(0,i.default)(t)}},function(t,r,n){"use strict";function e(t){return t&&t.__esModule?t:{default:t}}function a(t){if(t&&t.__esModule)return t;var r={};if(null!=t)for(var n in t)Object.prototype.hasOwnProperty.call(t,n)&&(r[n]=t[n]);return r.default=t,r}Object.defineProperty(r,"__esModule",{value:!0});var i=n(2),o=a(i),l=n(3),u=e(l),s=function(t,r,n){function e(t){return function(r,n){for(var e=r.quartiles[0],a=r.quartiles[2],i=(a-e)*t,n=-1,o=r.length;r[++n]<e-i;);for(;r[--o]>a+i;);return[n,o]}}function a(){var t=o.event.selection;if(t)s.domain([t[0][0],t[1][0]].map(s.invert,s)),c.domain([t[1][1],t[0][1]].map(c.invert,c)),r.select(".brush").call(h.move,null);else{if(!y)return y=setTimeout(g,p);s.domain(l),c.domain(u)}m()}var i=o.max(t,function(t){return t[0]}),l=[0,i],u=[0,100],s=o.scaleLinear().domain(l).range([n.margin.left,n.width]),c=o.scaleLinear().domain(u).range([n.height-n.margin.bottom,n.margin.top]),d=o.axisTop(s).ticks(12),f=o.axisRight(c).ticks(12*n.height/n.width),h=o.brush().on("end",a),y=void 0,p=350,x=o.box().whiskers(e(1.5)).height(n.height).domain([0,100]).showLabels(!1);r.append("g").attr("class","brush").attr("z-index",-1).call(h),r.selectAll(".box").attr("font","10px sans-serif").data(t).enter().append("g").on("click",function(t){console.log(t)}).attr("class","boxplot").attr("transform",function(t){return"translate("+s(t[0])+", "+n.margin.top+")"}).call(x.width((s.range()[1]-s.range()[0])/(s.domain()[1]-s.domain()[0])/2)),r.append("g").attr("class","axis axis--x").attr("transform","translate(0,"+(n.height-10)+")").call(d),r.append("g").attr("class","axis axis--y").attr("transform","translate(10,0)").call(f),r.selectAll(".domain").style("display","none"),r.append("text").attr("transform","rotate(-90)").attr("x",0-n.height/2).attr("dy","1em").style("text-anchor","middle").text("Quality Score"),r.append("text").attr("x",n.width/2).attr("y",n.height).attr("dy","1em").style("text-anchor","middle").text("Sequence Base");var g=function(){y=null},m=function(){var t=r.transition().duration(750);r.select(".axis--x").transition(t).call(d),r.select(".axis--y").transition(t).call(f),r.selectAll(".boxplot").transition(t).call(x.width((s.range()[1]-s.range()[0])/(s.domain()[1]-s.domain()[0])/2)).attr("transform",function(t){return"translate("+s(t[0])+", "+n.margin.top+")"})}},c=function(t){(0,u.default)(o);var r={top:10,right:30,bottom:30,left:30},n=o.select("#chartContainer").node().offsetWidth,e={margin:r,width:n-r.left-r.right,height:9*n/16-r.top-r.bottom},a=o.select("#chartContainer").append("svg").attr("width",e.width+e.margin.left+e.margin.right).attr("height",e.height+e.margin.top+e.margin.bottom),i=o.max(t,function(t){return t.length}),l=new Array(i),c=!0,d=!1,f=void 0;try{for(var h,y=t[Symbol.iterator]();!(c=(h=y.next()).done);c=!0)for(var p=h.value,x=0;x<p.length;x++)l[x]=l[x]||new Array(2),l[x][0]||(l[x][0]=x),l[x][1]||(l[x][1]=[]),l[x][1].push(p[x])}catch(t){d=!0,f=t}finally{try{!c&&y.return&&y.return()}finally{if(d)throw f}}s(l,a,e)};r.default=c},,function(t,r){"use strict";function n(t){function r(t){return[0,t.length-1]}function n(r){return[t.quantile(r,.25),t.quantile(r,.5),t.quantile(r,.75)]}t.box=function(){function e(r){r.each(function(r,n){var e=r[1].sort(t.ascending),u=t.select(this),h=e.length,y=e[0],p=e[h-1],x=e.quartiles=c(e),g=s&&s.call(this,e,n),m=g&&g.map(function(t){return e[t]}),v=g?t.range(0,g[0]).concat(t.range(g[1]+1,h)):t.range(h),b=t.scaleLinear().domain(l&&l.call(this,e,n)||[y,p]).range([i,0]),w=this.__chart__||t.scaleLinear().domain([0,1/0]).range(b.range());this.__chart__=b;var k=u.selectAll("line.center").data(m?[m]:[]);k.enter().append("line","rect").attr("class","center").attr("x1",a/2).attr("y1",function(t){return w(t[0])}).attr("x2",a/2).attr("y2",function(t){return w(t[1])}).style("opacity",1e-6).attr("stroke-dasharray","3,3").attr("stroke","black").attr("stroke-width","1px").transition().duration(o).style("opacity",1).attr("y1",function(t){return b(t[0])}).attr("y2",function(t){return b(t[1])}),k.transition().duration(o).style("opacity",1).attr("x1",a/2).attr("x2",a/2).attr("y1",function(t){return b(t[0])}).attr("y2",function(t){return b(t[1])}),k.exit().transition().duration(o).style("opacity",1e-6).attr("y1",function(t){return b(t[0])}).attr("y2",function(t){return b(t[1])}).remove();var _=u.selectAll("rect.box").data([x]);_.enter().append("rect").attr("class","box").attr("x",0).attr("y",function(t){return w(t[2])}).attr("width",a).attr("height",function(t){return w(t[0])-w(t[2])}).attr("fill","steelblue").attr("stroke","black").attr("stroke-width","1px").on("mouseover",function(r){t.select(this).attr("fill","lightgray")}).on("mouseout",function(){t.select(this).attr("fill","steelblue")}).transition().duration(o).attr("y",function(t){return b(t[2])}).attr("height",function(t){return b(t[0])-b(t[2])}),_.transition().duration(o).attr("y",function(t){return b(t[2])}).attr("width",a).attr("height",function(t){return b(t[0])-b(t[2])});var A=u.selectAll("line.median").data([x[1]]);A.enter().append("line").attr("class","median").attr("x1",0).attr("y1",w).attr("x2",a).attr("y2",w).attr("stroke","black").attr("stroke-width","1px").transition().duration(o).attr("y1",b).attr("y2",b),A.transition().duration(o).attr("x2",a).attr("y1",b).attr("y2",b);var q=u.selectAll("line.whisker").data(m||[]);q.enter().insert("line","circle, text").attr("class","whisker").attr("x1",0).attr("y1",w).attr("x2",0+a).attr("y2",w).style("opacity",1e-6).attr("stroke","black").attr("stroke-width","1px").transition().duration(o).attr("y1",b).attr("y2",b).style("opacity",1),q.transition().duration(o).attr("x2",0+a).attr("y1",b).attr("y2",b).style("opacity",1),q.exit().transition().duration(o).attr("y1",b).attr("y2",b).style("opacity",1e-6).remove();var L=u.selectAll("circle.outlier").data(v,Number);L.enter().insert("circle","text").attr("class","outlier").attr("r",a/4).attr("cx",a/2).attr("cy",function(t){return w(e[t])}).style("opacity",1e-6).attr("fill","none").attr("stroke","black").transition().duration(o).attr("cy",function(t){return b(e[t])}).style("opacity",1),L.transition().duration(o).attr("cy",function(t){return b(e[t])}).attr("r",a/8).attr("cx",a/2).style("opacity",1),L.exit().transition().duration(o).attr("cy",function(t){return b(e[t])}).style("opacity",1e-6).remove();var M=f||b.tickFormat(8),O=u.selectAll("text.box").data(x);1==d&&O.enter().append("text").attr("class","box").attr("dy",".3em").attr("dx",function(t,r){return 1&r?6:-6}).attr("x",function(t,r){return 1&r?+a:0}).attr("y",w).attr("text-anchor",function(t,r){return 1&r?"start":"end"}).text(M).transition().duration(o).attr("y",b),O.transition().duration(o).text(M).attr("y",b);var j=u.selectAll("text.whisker").data(m||[]);1==d&&j.enter().append("text").attr("class","whisker").attr("dy",".3em").attr("dx",6).attr("x",a).attr("y",w).text(M).style("opacity",1e-6).transition().duration(o).attr("y",b).style("opacity",1),j.transition().duration(o).text(M).attr("y",b).style("opacity",1),j.exit().transition().duration(o).attr("y",b).style("opacity",1e-6).remove()})}var a=1,i=1,o=750,l=null,u=Number,s=r,c=n,d=!0,f=null;e.width=function(t){return arguments.length?(a=t,e):a},e.height=function(t){return arguments.length?(i=t,e):i},e.tickFormat=function(t){return arguments.length?(f=t,e):f},e.duration=function(t){return arguments.length?(o=t,e):o};var h=function(t){return function(){return t}};return e.domain=function(t){return arguments.length?(l=null==t?t:h(t),e):l},e.value=function(t){return arguments.length?(u=t,e):u},e.whiskers=function(t){return arguments.length?(s=t,e):s},e.showLabels=function(t){return arguments.length?(d=t,e):d},e.quartiles=function(t){return arguments.length?(c=t,e):c},e}}Object.defineProperty(r,"__esModule",{value:!0}),r.default=n}]);