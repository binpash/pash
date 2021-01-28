const { hours, syncRemoteDirectory } = require('../lib.js');

module.exports = {
    makeCommandOptions: () => ({
        timeout: hours(1),
        shouldStopWaiting: (stdout, stderr) => /<<(fail|done)>>/.test(stdout),
        postCompletion: (ssh) => syncRemoteDirectory('reports', rc('ci_reports_path', `${__dirname}/reports`)),
    }),
};
