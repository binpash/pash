/*
Computes performance data summaries.

Performance test results are spead across many flat files, each with
timing data. This module analyzes those files to build tabular data,
and renders that data.
*/

const fs = require('fs');
const path = require('path');

const { table } = require('table');

const { lets } = require('./lib.js');


// Analysis: Convert result text files to relevant metadata and content.
/////////////////////////////////////////////////////////////////////////////////

const analyzeResults = (dir) =>
      Array
      .from(fs.readdirSync(dir))
      .filter((entry) => !fs.lstatSync(path.resolve(dir, entry)).isDirectory())
      .map((entry) => {
          const absolutePath = path.resolve(dir, entry);

          return {
              path: absolutePath,
              name: extractTestName(absolutePath),
              width: extractWidth(absolutePath),
              exec_time: extractExecTime(absolutePath),
              variant: extractVariant(absolutePath),
          };
      });

// Pash --width is in the file name
const extractWidth = (p) =>
      parseInt((path.basename(p, path.extname(p)).match(/_(\d+)_/) || [])[1]);

const extractVariant = (p) =>
      (p.match(/_distr_?([^\.]+)/) || [])[1];

const extractExecTime = (p) =>
      (fs.readFileSync(p).toString().match(/Execution time:[^\d]+([\d\.]+)/i) || [])[1];

// Test names appear to come before an underscore-prefixed numerical
// segment.  I split based off _\d because splitting by underscores
// will also split the name itself in some cases.
const extractTestName = (p) =>
      (path.basename(p, path.extname(p)).split(/_\d/) || [])[0];


// Presentation: Take _many_ outputs of analyzeResults() and present
// them as rows in a text table.
/////////////////////////////////////////////////////////////////////////////////

const createRow = (data, columnNames, width) => {
    const relevant = data.filter(({ width: actual }) => width === actual);
    const grouped = groupByTestName(relevant);
    return columnNames.map((c) => (grouped[c] || {}).exec_time || '');
};

const groupByTestName = (data) =>
      data.reduce((res, e) => Object.assign(res, { [e.name]: e }), {});

const getColumnNames = (history) => {
    const allNames = new Set();

    for (const [commit, data] of history) {
        for (const { name } of data) {
            allNames.add(name);
        }
    }

    return Array.from(allNames).sort();
};

const renderPerfSummary = (history, { width }) => {
    const columnNames = getColumnNames(history);

    const firstRow = [''].concat(columnNames);

    const rows = history.map(([commit, data]) =>
        [commit].concat(createRow(data, columnNames, width)));

    return table([firstRow].concat(rows), {
        columns: {
            0: {
                width: 12,
            },
        },
    });
};


module.exports = {
    analyzeResults,
    renderPerfSummary,
};
