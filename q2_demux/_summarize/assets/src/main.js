// ----------------------------------------------------------------------------
// Copyright (c) 2016-2017, QIIME 2 development team.
//
// Distributed under the terms of the Modified BSD License.
//
// The full license is in the file LICENSE, distributed with this software.
// ----------------------------------------------------------------------------

import initializePlot from './plot';

export const init = (seqProps, fwdScores, revScores = undefined) => { // eslint-disable-line
  const qualScores = { forward: Object.keys(fwdScores).map(key => [+key + 1, fwdScores[key]]) };
  if (revScores) {
    qualScores.reverse = Object.keys(revScores).map(key => [+key + 1, revScores[key]]);
  }
  initializePlot(qualScores, seqProps);
};
