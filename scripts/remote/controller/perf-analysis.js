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

const analyze = (dir) =>
      Array
      .from(fs.readdirSync(dir))
      .filter((entry) => !fs.lstatSync(path.resolve(dir, entry)).isDirectory())
      .map((entry) => {
          const absolutePath = path.resolve(dir, entry);

          return {
              path: absolutePath,
              name: extractTestName(absolutePath),
              width: extractWidth(absolutePath),
              exec_time: extractExecTime(fs.readFileSync(absolutePath).toString()),
              variant: extractVariant(absolutePath),
          };
      });

// Pash --width is in the file name
const extractWidth = (p) =>
      parseInt((path.basename(p, path.extname(p)).match(/_(\d+)_/) || [])[1]);

const extractVariant = (p) =>
      (path.basename(p, path.extname(p)).match(/\d+_(\D+)$/) || [])[1]

const extractExecTimeRec = (str, rx, sum) => {
    const match = rx.exec(str);
    return (match)
        ? extractExecTimeRec(str, rx, parseFloat(match[1]) + sum)
        : sum;
};

const extractExecTime = (str) =>
      extractExecTimeRec(str, /Execution time:[^\d]+([\d\.]+)/gi, 0);

// Test names appear to come before an underscore-prefixed numerical
// segment.  I split based off _\d because splitting by underscores
// will also split the name itself in some cases.
const extractTestName = (p) =>
      (path.basename(p, path.extname(p)).split(/_\d/) || [])[0];

const keepWidth = (data, width) =>
      data.filter(({ width: actual }) => width === actual);

const keepTimed = (data) =>
      data.filter(({ exec_time }) => exec_time);


// Presentation: Take _many_ outputs of analyzeResults() and present
// them as rows in a text table.
/////////////////////////////////////////////////////////////////////////////////

const createRow = (data, columnNames) => {
    const grouped = groupByTestName(data);
    return columnNames.map((c) => (grouped[c] || { exec_time: '' }).exec_time);
};

const groupByTestName = (data) =>
      data.reduce((res, e) => Object.assign(res, { [e.name]: e }), {});

const findTestNames = (data, names = new Set()) => {
    for (const { name } of data)
        names.add(name);
    return names;
};

const render = (history, tests) => {
    const columnNames = Array.from(tests).slice(0).sort();
    const firstRow = [''].concat(columnNames);
    const rows = history.map(([label, data]) =>
        [label].concat(createRow(data, columnNames)));

    return table([firstRow].concat(rows));
};


module.exports = {
    analyze,
    extractWidth,
    extractExecTime,
    findTestNames,
    keepWidth,
    keepTimed,
    render,
};
