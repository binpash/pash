const Rsync = require('rsync');
const { hours } = require('../lib.js');

module.exports = {
    makeCommandOptions: () => ({
        timeout: hours(24),
        shouldStopWaiting: (stdout, stderr) => true,
    }),
};
