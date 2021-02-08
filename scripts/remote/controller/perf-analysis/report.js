const fs = require('fs');
const path = require('path');
const history = require('./history.js');
const { analyzePerfSuite, summarizePerfData, makePerfFileName } = require('./parse.js');

const directoryExists = (path) => {
    try {
        return fs.lstatSync(path).isDirectory();
    } catch (e) {
        return false;
    }
};

const reportPerfSuite = (dir, tests, {
    width,
    variant,
}) => {
    const sorted = tests.slice(0).sort();
    return history.render([
        sorted,
        analyzePerfSuite(dir, sorted, width, variant).map(summarizePerfData)
    ]);
}

const accumulateDirectoryListings = (dirs, accum = []) =>
  (dirs.length === 0)
      ? accum.sort()
      : accumulateDirectoryListings(dirs.slice(1), accum.concat(fs.readdirSync(dirs[0])));

module.exports = {
    reportPerfSuite,
};
