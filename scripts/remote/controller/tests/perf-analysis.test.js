const { extractWidth, extractExecTime } = require('../perf-analysis.js');
const { suite, assert, assertEqv } = require('./testlib.js');

suite('perf-analysis', () => {
    const path = '/path/t0/red-h3rring_88_distr.time';

    assert('Extract DFG width from file name',
           extractWidth(path) === 88);

    const junk = `
    Execution Time: 87.111junk
    junk junk garbage junk
    execution time:junk93012.0012`;

    assertEqv('Extract the sum of execution times',
              extractExecTime(junk), 93012.0012 + 87.111);

    
});
