const path = require('path');

const { analyzePerfFile, fileContentCache } = require('../perf-analysis/parse.js');
const { suite, assert, assertEqv } = require('./testlib.js');

const testPath = (...a) => path.resolve(__dirname, 'dummy-results', ...a);

suite('Analyze performance result file', () => {
    const inputPath = testPath('red-h3rring_64_distr.time');
    const { path, width, name, analyzeContent, variant } =
          analyzePerfFile(testPath('red-h3rring_64_distr.time'));

    assert('Find test name', name === 'red-h3rring');
    assert('Find test variant', variant === 'distr');
    assert('Track input path', path === inputPath);
    assert('Extract DFG width from file name', width === 64);
    assert('Bind closure for analyzing content', typeof analyzeContent === 'function');

    const { execTime, bashTimes, content, speedup } = analyzeContent();

    assertEqv('Extract the sum of execution times', execTime, (93012.0012 + 87.111) / 1000);
    
    assert('Extract Bash timings', typeof bashTimes === 'object');
    const { real, user, sys } = bashTimes;
    
    assertEqv('Extract real bash time', real, 46.52);
    assertEqv('Extract user bash time', user, 50.962);
    assertEqv('Extract sys bash time', sys, 13.48);

    assertEqv('Compute speedup in terms of seq', speedup, 2.480818);
});

// Otherwise the process will wait on timeouts.
fileContentCache.clear();
