const fs = require('fs');
const path = require('path');
const history = require('./history.js');
const { fileContentCache,
        analyzePerfSuite,
        summarizePerfData,
        makePerfFileName } = require('./parse.js');

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

if (require.main === module) {
    const [ , , dir, tests, width = 2, variant = 'distr_auto_split'] = process.argv;

    if (!dir || !tests) {
        console.log('usage: report.js RESULTS_DIR TESTS [WIDTH] [VARIANT]');
        console.log('e.g. node report.js ../results wf,sort 16 distr');
        process.exit(1);
    } else {
        console.log(reportPerfSuite(dir, tests.split(',').map((n) => n.trim()),
                                    { width: parseInt(width), variant }));
        fileContentCache.clear();
    }
}
