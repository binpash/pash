const Rsync = require('rsync');
const { hours } = require('../lib.js');

module.exports = {
    makeCommandOptions: () => ({
        timeout: hours(1),
        shouldStopWaiting: (stdout, stderr) => /<<(fail|done)>>/.test(stdout),
        postCompletion: async (ssh) => {
            try {
                // TODO: Stream rsync stdin/stdout to new files.
                const rsync = (
                    new Rsync()
                        .flags('abziu')
                        .set('e', `ssh -i ${rc('private_key')}`)
                        .source(`${rc('user')}@${rc('host')}:reports/`) // Trailing '/' is meaningful: https://serverfault.com/a/529294
                        .destination(rc('ci_reports_path', `${__dirname}/reports`))
                );

                await new Promise((resolve, reject) =>
                    rsync.execute((error, code, cmd) => {
                        if (error)
                            reject(error)
                        else
                            resolve();
                    }));
            } catch (e) {
                err(e);
                throw e;
            }
        },
    }),
};
