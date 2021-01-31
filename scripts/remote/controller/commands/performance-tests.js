const { hours, syncRemoteDirectory } = require('../lib.js');
const rc = require('../rc.js');

module.exports = {
    makeCommandOptions: () => ({
        timeout: hours(24),
        shouldStopWaiting: (stdout, stderr) => /<<(fail|done)>>/.test(stdout),
        postCompletion: (ssh) => (
            syncRemoteDirectory('results', rc('perf_results_path', `${__dirname}/perf_results`))
        ),
    }),
};
