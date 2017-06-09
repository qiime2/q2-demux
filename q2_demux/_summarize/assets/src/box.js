// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import * as d3 from 'd3';

export default function plotBoxes(svg, data, x, y, seqProps) {
  const halfWidth = (x.range()[1] - x.range()[0]) / (x.domain()[1] - x.domain()[0]) / 2;
  const quarterWidth = halfWidth / 2;
  const t = svg.transition().duration(750);
  const darkBlue = 'steelblue';
  const lightBlue = 'skyblue';
  const darkRed = '#a94442';
  const lightRed = '#ebccd1';
  const minSeqLen = seqProps.minSeqLen[data.direction];

  const containerUpdate = svg.selectAll('.container')
    .data(data);
  containerUpdate.exit().remove();
  const containerEnter = containerUpdate.enter()
    .append('g')
    .attr('class', 'container');
  const containers = containerUpdate.merge(containerEnter)
    .transition(t)
    .attr('transform', d => `translate(${x(d[0])}, 0)`)
    .selection()
    .on('mouseover', function mouseover() {
      const data = d3.select(this).data();
      const position = data[0][0];
      const stats = data[0][1];
      const inTheDangerZone = stats.count < seqProps.n;
      const svg = d3.select(this.parentNode).node();
      const plotContainer = d3.select(svg.parentNode);
      plotContainer
        .select('.panel')
        .attr('class', inTheDangerZone ? 'panel panel-danger' : 'panel panel-default');
      plotContainer
        .select('.panel-heading')
        .html(`Parametric seven-number summary for <strong>position ${position}</strong>`);
      plotContainer
        .select('.stats')
        .select('tbody')
        .selectAll('tr')
        .data([
          ['(Not shown in box plot)', '2nd', stats['2%']],
          ['Lower Whisker', '9th', stats['9%']],
          ['Bottom of Box', '25th', stats['25%']],
          ['Middle of Box', '50th (Median)', stats['50%']],
          ['Top of Box', '75th', stats['75%']],
          ['Upper Whisker', '91st', stats['91%']],
          ['(Not shown in box plot)', '98th', stats['98%']],
        ])
        .selectAll('td')
          .data(d => d)
          .text(d => d);

      let seqLenNote = `The minimum sequence length identified during subsampling was ${minSeqLen}`;
      if (inTheDangerZone) {
        seqLenNote = `This position (${position}) is greater than the minimum sequence length observed
                      during subsampling (${minSeqLen} bases). As a result, the plot at this position
                      is not based on data from all of the sequences, so it should be interpreted with 
                      caution when compared to plots for other positions`;
      }

      plotContainer.select('.random-sampling')
        .classed('text-danger', inTheDangerZone)
        .html(`The plot at position ${position} was generated using a random
               sampling of ${stats.count} out of ${seqProps.totalSeqCount} sequences
               without replacement. ${seqLenNote}. Outlier quality scores are
               not shown in box plots for clarity.`);
    });

  const centerUpdate = containers.selectAll('line.center').data(d => [d]);
  centerUpdate.exit().remove();
  const centerEnter = centerUpdate.enter().append('line');
  centerUpdate.merge(centerEnter)
    .transition(t)
    .attr('class', 'center')
    .attr('x1', 0)
    .attr('y1', d => y(d[1]['9%']))
    .attr('x2', 0)
    .attr('y2', d => y(d[1]['91%']))
    .attr('stroke-dasharray', '2,2')
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const boxUpdate = containers.selectAll('rect.box').data(d => [d]);
  boxUpdate.exit().remove();
  const boxEnter = boxUpdate.enter().append('rect');
  boxUpdate.merge(boxEnter)
    .transition(t)
    .attr('class', 'box')
    .attr('x', -quarterWidth)
    .attr('y', d => y(d[1]['75%']))
    .attr('width', halfWidth)
    .attr('height', d => (y(d[1]['25%']) - y(d[1]['75%'])))
    .attr('fill', d => (d[1]['count'] < seqProps.n ? lightRed : lightBlue))
    .attr('stroke-width', 1)
    .attr('stroke', 'black')
    .selection()
    // The two event handlers don't use fat-arrows because we need the lexical `this` in scope
    .on('mouseover', function mouseover() {
      d3.select(this)
        .attr('fill', d => (d[1]['count'] < seqProps.n ? darkRed : darkBlue));
      })
    .on('mouseout', function mouseout() {
      d3.select(this)
        .attr('fill', d => (d[1]['count'] < seqProps.n ? lightRed : lightBlue));
    });

  const medianUpdate = containers.selectAll('line.median').data(d => [d]);
  medianUpdate.exit().remove();
  const medianEnter = medianUpdate.enter().append('line');
  medianUpdate.merge(medianEnter)
    .transition(t)
    .attr('class', 'median')
    .attr('x1', -quarterWidth)
    .attr('y1', d => y(d[1]['50%']))
    .attr('x2', quarterWidth)
    .attr('y2', d => y(d[1]['50%']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const lowerWhiskerUpdate = containers.selectAll('line.lower-whisker').data(d => [d]);
  lowerWhiskerUpdate.exit().remove();
  const lowerWhiskerEnter = lowerWhiskerUpdate.enter().append('line');
  lowerWhiskerUpdate.merge(lowerWhiskerEnter)
    .transition(t)
    .attr('class', 'lower-whisker')
    .attr('x1', -quarterWidth)
    .attr('y1', d => y(d[1]['9%']))
    .attr('x2', quarterWidth)
    .attr('y2', d => y(d[1]['9%']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');

  const upperWhiskerUpdate = containers.selectAll('line.upper-whisker').data(d => [d]);
  upperWhiskerUpdate.exit().remove();
  const upperWhiskerEnter = upperWhiskerUpdate.enter().append('line');
  upperWhiskerUpdate.merge(upperWhiskerEnter)
    .transition(t)
    .attr('class', 'upper-whisker')
    .attr('x1', -quarterWidth)
    .attr('y1', d => y(d[1]['91%']))
    .attr('x2', quarterWidth)
    .attr('y2', d => y(d[1]['91%']))
    .attr('stroke-width', 1)
    .attr('stroke', 'black');
}
