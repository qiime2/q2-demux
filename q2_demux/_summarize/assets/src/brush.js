// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

import plotBoxes from './box';

export function updateXTicks(svg, x) {
  const ticks = svg.selectAll('.axis--x .tick');
  const lenTicks = x.domain().length;

  if (lenTicks > 48) {
    ticks.style('display', i => (i % 12 ? 'none' : 'initial'));
  } else if (lenTicks > 24) {
    ticks.style('display', i => (i % 2 ? 'none' : 'initial'));
  } else {
    ticks.style('display', 'initial');
  }
}

export function addBrush(svg, data, x, y, x0, y0, xAxis, yAxis) {
  const brush = d3.brush();
  let idleTimeout;
  const idleDelay = 350;
  const ticks = svg.selectAll('.axis--x .tick');

  const idled = () => {
    idleTimeout = null;
  };

  const zoom = () => {
    svg.select('.axis--x').call(xAxis);
    svg.select('.axis--y').call(yAxis);
    plotBoxes(svg, data, x, y);
    updateXTicks(svg, x);
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
      const tempX = d3.scaleQuantize().domain(x.range()).range(x.domain());
      const low = tempX(s[0][0]);
      const high = tempX(s[1][0]);

      x.domain(Array(high - low + 1).fill().map((_, i) => low + i));
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
