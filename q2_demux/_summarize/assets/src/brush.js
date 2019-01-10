// ----------------------------------------------------------------------------
// Copyright (c) 2016-2019, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

import plotBoxes from './box';

export function addBrush(svg, data, x, y, x0, y0, xAxis, yAxis, ticksCounts, seqProps) {
  const brush = d3.brush();
  let idleTimeout;
  const idleDelay = 350;
  const ticks = svg.selectAll('.axis--x .tick');

  const idled = () => {
    idleTimeout = null;
  };

  const zoom = () => {
    const t = svg.transition().duration(750);
    const [xTick, yTick] = ticksCounts;

    const xMin = Math.floor(x.domain()[0]);
    const xMax = Math.ceil(x.domain().slice(-1)[0]);
    const xVals = d3.range(xMin, xMax, 1);
    if (xVals.length <= xTick) {
        xAxis.tickValues(xVals).tickFormat(d3.format('d'));
    } else {
        xAxis.tickValues(null).ticks(xTick, d3.format('d'));
    }

    const yMin = Math.floor(y.domain()[0]);
    const yMax = Math.ceil(y.domain().slice(-1)[0]);
    const yVals = d3.range(yMin, yMax, 1);
    if (yVals.length <= yTick) {
        yAxis.tickValues(yVals).tickFormat(d3.format('d'));
    } else {
        yAxis.tickValues(null).ticks(yTick, d3.format('d'));
    }

    svg.select('.axis--x').transition(t).call(xAxis);
    svg.select('.axis--y').transition(t).call(yAxis);
    plotBoxes(svg, data, x, y, seqProps);
  };

  const brushEnded = () => {
    const s = d3.event.selection;
    if (!s) {
      if (!idleTimeout) {
        idleTimeout = setTimeout(idled, idleDelay);
        return idleTimeout;
      }
      x.domain(x0);
      y.domain(y0);
    } else {
      x.domain([s[0][0], s[1][0]].map(x.invert, x));
      y.domain([s[1][1], s[0][1]].map(y.invert, y));
      svg.select('.brush').call(brush.move, null);
    }
    return zoom();
  };

  brush.on('end', brushEnded);

  svg.append('g')
      .attr('class', 'brush')
      .call(brush);
}
