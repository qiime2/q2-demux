// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

export default function plotBoxes(svg, data, x, y) {
  const halfWidth = x.bandwidth() < 2 ? 1 : x.bandwidth() / 2;
  const quarterWidth = halfWidth / 2;

  const containerUpdate = svg.selectAll('.container')
    .data(data)
    .attr('font', '10px sans-serif');
  containerUpdate.exit().remove();
  const containerEnter = containerUpdate.enter()
    .append('g')
    .attr('class', 'container');
  const containers = containerUpdate.merge(containerEnter)
    .attr('transform', d => `translate(${x(d[0]) ? x(d[0]) : -100}, 0)`)
    .attr('display', d => (x(d[0]) ? 'initial' : 'none'))
    .on('mouseover', function mouseover() {
      const data = d3.select(this).data();
      const stats = data[0][1];
      const quartiles = [stats['25%'], stats['50%'], stats['75%']];
      const whiskers = [stats.min, stats.max];
      const svg = d3.select(this.parentNode).node();
      const plotContainer = d3.select(svg.parentNode);
      plotContainer
        .select('.stats')
        .select('tbody')
        .selectAll('tr')
        .data([
          ['Position Number', data[0][0]],
          ['Minimum', whiskers[0]],
          ['1st Quartile', quartiles[0]],
          ['Median', quartiles[1]],
          ['3rd Quartile', quartiles[2]],
          ['Maximum', whiskers[1]],
        ])
        .selectAll('td')
          .data(d => d)
          .text(d => d);
    });

  const centerUpdate = containers.selectAll('line.center').data(d => [d]);
  centerUpdate.exit().remove();
  const centerEnter = centerUpdate.enter().append('line');
  centerUpdate.merge(centerEnter)
    .attr('class', 'center')
    .attr('x1', halfWidth)
    .attr('y1', d => y(d[1]['min']))
    .attr('x2', halfWidth)
    .attr('y2', d => y(d[1]['max']))
    .attr('stroke-dasharray', '2,2')
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const boxUpdate = containers.selectAll('rect.box').data(d => [d]);
  boxUpdate.exit().remove();
  const boxEnter = boxUpdate.enter().append('rect');
  boxUpdate.merge(boxEnter)
    .attr('class', 'box')
    .attr('x', quarterWidth)
    .attr('y', d => y(d[1]['75%']))
    .attr('width', halfWidth)
    .attr('height', d => (y(d[1]['25%']) - y(d[1]['75%'])))
    .attr('fill', 'steelblue')
    .attr('stroke-width', 1)
    .attr('stroke', 'black')
    // The two event handlers don't use fat-arrows because we need the lexical `this` in scope
    .on('mouseover', function mouseover() {
      d3.select(this)
        .attr('fill', 'skyblue');
      })
    .on('mouseout', function mouseout() {
      d3.select(this)
        .attr('fill', 'steelblue');
    });

  const medianUpdate = containers.selectAll('line.median').data(d => [d]);
  medianUpdate.exit().remove();
  const medianEnter = medianUpdate.enter().append('line');
  medianUpdate.merge(medianEnter)
    .attr('class', 'median')
    .attr('x1', quarterWidth)
    .attr('y1', d => y(d[1]['50%']))
    .attr('x2', halfWidth + quarterWidth)
    .attr('y2', d => y(d[1]['50%']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const lowerWhiskerUpdate = containers.selectAll('line.lower-whisker').data(d => [d]);
  lowerWhiskerUpdate.exit().remove();
  const lowerWhiskerEnter = lowerWhiskerUpdate.enter().append('line');
  lowerWhiskerUpdate.merge(lowerWhiskerEnter)
    .attr('class', 'lower-whisker')
    .attr('x1', quarterWidth)
    .attr('y1', d => y(d[1]['min']))
    .attr('x2', halfWidth + quarterWidth)
    .attr('y2', d => y(d[1]['min']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const upperWhiskerUpdate = containers.selectAll('line.upper-whisker').data(d => [d]);
  upperWhiskerUpdate.exit().remove();
  const upperWhiskerEnter = upperWhiskerUpdate.enter().append('line');
  upperWhiskerUpdate.merge(upperWhiskerEnter)
    .attr('class', 'upper-whisker')
    .attr('x1', quarterWidth)
    .attr('y1', d => y(d[1]['max']))
    .attr('x2', halfWidth + quarterWidth)
    .attr('y2', d => y(d[1]['max']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');
}
