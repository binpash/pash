// A controller implementation that uses SSH
//
// WARNING: Do not use `await' where you see <!> in a comment. A broken
// promise chain is intentional when waiting would hang the HTTP
// connection.  Use informational GET endpoints to monitor status of
// the orphaned promise.

const fs = require('fs');
const os = require('os');
const path = require('path');

const { NodeSSH } = require('node-ssh');

const { hours, log, err, getRequestBody, respond, getRemoteHomeDirectory, getSshCredentials } = require('./lib.js');
const rc = require('./rc.js');

const sshConnect = async (username, host) => {
    const ssh = new NodeSSH();
    await ssh.connect(getSshCredentials());
    return ssh;
};

// The SSH commands do not actually wait for the commands to complete.
// Current approach is to inspect the output to know when to proceed.
//
// FIXME: This doesn't work unless the SSH connection remains unbroken
// for the entire life of the command. This should instead consult the
// host's process list to see if an instance of the command is still
// running.
const createStreamMonitor = (shouldStopWaiting, timeout) => {
    let _resolve, alarm;

    const check = (stdout, stderr) => {
        if (shouldStopWaiting(stdout, stderr)) {
            clearTimeout(alarm);
            _resolve();
        }
    };

    const promise = new Promise((resolve, reject) => {
        // Allow external control of promise resolution.
        // Works because the callback for the Promise
        // constructor is called immediately and synchronously.
        _resolve = resolve;
        alarm = setTimeout(() => reject(new Error('<Command timed out>')), timeout);
    });

    return [check, promise];
}

// Ensures only connection is active at a time.
const command = (() => {
    let initial = true;
    let idle = true;
    let promise = Promise.resolve();
    let stdout = '';
    let stderr = '';

    const lock = () => {
        log('Locking');
        idle = initial = false;
    };

    const unlock = () => {
        log('Unlocking');
        idle = true;
    };

    async function exec({
        cwd = getRemoteHomeDirectory(),
        shouldStopWaiting,
        timeout = hours(1),
        postCompletion = () => {},
    }) {
        stdout = stderr = '';
        const ssh = await sshConnect();
        const [check, completion] = createStreamMonitor(shouldStopWaiting, timeout);

        await ssh.execCommand(`${getRemoteHomeDirectory()}/work.sh`, {
            cwd,
            onStdout: (c) => (stdout += c.toString('utf8'), check(stdout, stderr)),
            onStderr: (c) => (stderr += c.toString('utf8'), check(stdout, stderr)),
        });

        await completion;
        await Promise.resolve(postCompletion(ssh));
        ssh.dispose();
    }

    const exports = {
        neverRan: () => initial,
        isIdle: () => idle,
        dump: () => ({ stdout, stderr }),
        run: (options) => {
            if (idle) {
                lock();
                return exec(options)
                    .catch((e) => { unlock(); return Promise.reject(e); })
                    .then((v) => { unlock(); return v; });
            } else {
                const msg = 'A command is already running.';
                const err = new Error(msg);
                err.respond = (res) => respond(res, 400, `${msg}\n`);
                throw err;
            }
        },
    };

    return exports;
})();

module.exports = {
    command,
    getRemoteHomeDirectory,
};
