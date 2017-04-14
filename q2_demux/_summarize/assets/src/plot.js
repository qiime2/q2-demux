// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

import plotBoxes from './box';
import { updateXTicks, addBrush } from './brush';

const plot = (data, props, container) => {
  const plotContainer = d3.select(container);

  const svg = plotContainer
    .append('svg')
      .attr('class', 'col-xs-12')
      .style('display', 'block')
      .style('margin', '0 auto')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.height + props.margin.top + props.margin.bottom);

  const panel = plotContainer
    .append('div')
      .attr('class', 'col-xs-12')
    .append('div')
      .attr('class', 'panel panel-default');
  panel.append('div')
    .attr('class', 'panel-heading')
    .html('Statistical Summary');
  const table = panel.append('div')
    .attr('class', 'stats')
  .append('table')
    .attr('class', 'table')
    .style('margin-bottom', '0');

  table
    .append('thead')
    .append('tr')
      .selectAll('th')
    .data([['Hover over a boxplot to learn more...', 9], ['', 3]])
      .enter()
    .append('th')
    .text(d => d[0])
    .attr('class', d => `col-xs-${d[1]}`)
    .style('font-size', '10px');

  const initialData = [
    ['Position Number', '...'],
    ['Minimum', '...'],
    ['1st Quartile', '...'],
    ['Median', '...'],
    ['3rd Quartile', '...'],
    ['Maximum', '...'],
  ];

  const rows = table
    .append('tbody')
      .selectAll('tr')
      .data(initialData)
    .enter()
      .append('tr');

  rows
    .selectAll('td')
    .data(d => d)
      .enter()
    .append('td')
    .text(d => d);

  const maxX = d3.max(data, d => d[0]);
  const x0 = [0, maxX];
  const y0 = [0, 45];
  const x = d3.scaleLinear()
    .domain(x0)
    .range([props.margin.left, props.width]);
  const y = d3.scaleLinear()
    .domain(y0)
    .range([props.height - props.margin.bottom, props.margin.top]);

  const xAxis = d3.axisBottom(x).ticks(12, d3.format('d'));
  const yAxis = d3.axisLeft(y).ticks(9);

  svg.attr('height', props.height + props.margin.bottom + props.margin.top);
  svg.attr('width', x.range()[0] + x.range()[1]);

  addBrush(svg, data, x, y, x0, y0, xAxis, yAxis)
  plotBoxes(svg, data, x, y);

  svg.append('rect')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.margin.bottom * 4)
      .attr('y', props.height - props.margin.bottom)
      .attr('fill', 'white');

  svg.append('rect')
      .attr('width', props.margin.left * 4)
      .attr('x', -props.margin.left * 3)
      .attr('height', props.height + 200)
      .attr('fill', 'white');

  svg.append('rect')
    .attr('width', props.margin.right + 1000)
    .attr('height', props.height + 200)
    .attr('x', props.width + 10)
    .attr('fill', 'white');

  svg.append('g')
      .attr('class', 'axis axis--x')
      .attr('transform', `translate(0, ${(props.height - props.margin.bottom)})`)
      .call(xAxis);

  svg.append('g')
      .attr('class', 'axis axis--y')
      .attr('transform', `translate(${props.margin.left},0)`)
      .call(yAxis);

  svg.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('x', 0 - (props.height / 2))
      .attr('dy', '1em')
      .attr('font-size', '10px')
      .style('text-anchor', 'middle')
      .text('Quality Score');

  svg.append('text')
      .attr('x', props.width / 2)
      .attr('y', props.height)
      .attr('dy', '1em')
      .attr('font-size', '10px')
      .style('text-anchor', 'middle')
      .text('Sequence Base');
};

const initializePlot = (data) => {
  const margin = { top: 10, right: 30, bottom: 30, left: 30 };
  const width = d3.select('#forwardContainer').node().offsetWidth;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom };

  Object.keys(data).forEach((direction) => {
    plot(data[direction], props, `#${direction}Container`);
  });
};

export default initializePlot;
