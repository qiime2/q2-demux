// ----------------------------------------------------------------------------
// Copyright (c) 2016-2020, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

import plotBoxes from './box';
import { updateXTicks, addBrush } from './brush';

const plot = (data, props, container, seqProps) => {
  const plotContainer = d3.select(container);
  const direction = data.direction;
  const svg = plotContainer
    .append('svg')
      .attr('class', 'col-xs-12')
      .style('display', 'block')
      .style('margin', '0 auto')
      .attr('width', props.width + props.margin.left + props.margin.right)
      .attr('height', props.height + props.margin.top + props.margin.bottom);

  plotContainer
    .append('div')
    .attr('class', 'col-xs-12')
    .append('p')
      .attr('class', 'random-sampling')
      .html(`These plots were generated using a random sampling of ${seqProps.subsampleSize[direction]}
             out of ${seqProps.totalSeqCount[direction]} sequences without replacement. The
             minimum sequence length identified during subsampling was
             ${seqProps.minSeqLen[direction]} bases. Outlier quality scores are not shown in box
             plots for clarity.`);

  const panel = plotContainer
    .append('div')
      .attr('class', 'col-xs-12')
    .append('div')
      .attr('class', 'panel panel-default');
  panel.append('div')
    .attr('class', 'panel-heading')
    .html('Parametric seven-number summary');
  const table = panel.append('div')
    .attr('class', 'stats')
  .append('table')
    .attr('class', 'table')
    .style('margin-bottom', '0');

  table
    .append('thead')
    .append('tr')
      .selectAll('th')
    .data([['Box plot feature', 5], ['Percentile', 5], ['Quality score', 2]])
      .enter()
    .append('th')
    .text(d => d[0])
    .attr('class', d => `col-xs-${d[1]}`);

  const initialData = [
    ['(Not shown in box plot)', '2nd', '...'],
    ['Lower Whisker', '9th', '...'],
    ['Bottom of Box', '25th', '...'],
    ['Middle of Box', '50th (Median)', '...'],
    ['Top of Box', '75th', '...'],
    ['Upper Whisker', '91st', '...'],
    ['(Not shown in box plot)', '98th', '...'],
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

  plotContainer
    .append('div')
      .attr('class', 'col-xs-12')
    .append('div')
      .html(`<a href="${direction}-seven-number-summaries.tsv">Download ${direction} parametric seven-number summaries as TSV</a>`)

  const maxX = d3.max(data, d => d[0]) + 1;
  const x0 = [0, maxX];
  const y0 = [0, 45];
  const x = d3.scaleLinear()
    .domain(x0)
    .range([props.margin.left, props.width]);
  const y = d3.scaleLinear()
    .domain(y0)
    .range([props.height - props.margin.bottom, props.margin.top]);
  const ticks = [12, 9];

  const xAxis = d3.axisBottom(x).ticks(ticks[0], d3.format('d'));
  const yAxis = d3.axisLeft(y).ticks(ticks[1], d3.format('d'));

  svg.attr('height', props.height + props.margin.bottom + props.margin.top);
  svg.attr('width', x.range()[0] + x.range()[1]);

  addBrush(svg, data, x, y, x0, y0, xAxis, yAxis, ticks, seqProps)
  plotBoxes(svg, data, x, y, seqProps);

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
    .attr('x', props.width)
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
      .attr('font-size', '12px')
      .style('text-anchor', 'middle')
      .text('Quality Score');

  svg.append('text')
      .attr('x', props.width / 2)
      .attr('y', props.height)
      .attr('dy', '1em')
      .attr('font-size', '12px')
      .style('text-anchor', 'middle')
      .text('Sequence Base');
};

const initializePlot = (data, seqProps) => {
  const margin = { top: 10, right: 30, bottom: 30, left: 40 };
  const fwdWidthNode = d3.select('#forwardContainer').node();
  const revWidthNode = d3.select('#reverseContainer').node();
  const width = fwdWidthNode ? fwdWidthNode.offsetWidth : revWidthNode.offsetWidth;
  const props = {
    margin,
    width: width - margin.left - margin.right,
    height: ((width * 9) / 16) - margin.top - margin.bottom };

  Object.keys(data).forEach((direction) => {
    data[direction].direction = direction
    plot(data[direction], props, `#${direction}Container`, seqProps);
  });
};

export default initializePlot;
