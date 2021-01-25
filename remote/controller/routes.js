// A controller implementation that uses SSH
//
// WARNING: Do not use `await' where you see <!> in a comment. A broken
// promise chain is intentional when waiting would hang the HTTP
// connection.  Use informational GET endpoints to monitor status of
// the orphaned promise.

const { log, err, respond, getRequestBody } = require('./lib.js');
const { NodeSSH } = require('node-ssh');
const fs = require('fs');
const os = require('os');
const path = require('path');
const syncdir = require('sync-directory');
const rimraf = require('rimraf');

const getSshCredentials = () => ({
    host: process.env.PASH_REMOTE_HOST || 'localhost',
    username: process.env.PASH_REMOTE_USER || process.env.USER,
    privateKey: process.env.PASH_REMOTE_PRIVATE_KEY || `${process.env.HOME}/.ssh/id_rsa`,
});

const getRemoteHomeDirectory = () => {
    const { username } = getSshCredentials();
    return `/home/${username}`;
};

const sshConnect = async (username, host) => {
    const ssh = new NodeSSH();
    await ssh.connect(getSshCredentials());
    return ssh;
};

const hours = h => (1000 * 60 * 60 * h);


// The command module ensures only one command runs at a time.
// TODO: Upgrade for use with more than one host at a time.
const command = (() => {
    let initial = true;
    let idle = true;
    let promise = Promise.resolve();
    let stdout = '';
    let stderr = '';

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
                idle = initial = false;
                return exec(options)
                    .catch((e) => { idle = true; return Promise.reject(e); })
                    .then((v) => { idle = true; return v; });
            } else {
                const msg = 'A command is already running.';
                const err = new Error(msg);
                err.respond = (res) => respond(res, 400, msg);
                throw err;
            }
        },
    };

    return exports;
})();



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

const run = async (req, res) => {
    // <!>
    command.run({ shouldStopWaiting: () => true });
    respond(res, 200, 'Sent command\n');
};

const ci = async (req, res) => {
    const homedir = getRemoteHomeDirectory();

    // <!>
    command.run({
        timeout: hours(1),
        shouldStopWaiting: (stdout, stderr) => /<<(fail|done)>>/.test(stdout),
        postCompletion: async (ssh) => {
            // FIXME: Download size can grow without bound. Revisit when it becomes a problem.

            const tmpdir = fs.mkdtempSync(path.join(os.tmpdir(), 'pash-'));
            const dir = `${__dirname}/reports`;

            try {
                // FIXME: Blowing away the directory is too crude. A
                // sync function is below, but it doesn't work right
                // now. This is a quick fix to help the team ensure
                // they have a serveable reports directory.
                // nv: this doesn't work prior to 10.x, and upgrading to later
                //     seems to cause other problems, so use rimraf for now
                // fs.rmdirSync(dir, { recursive: true });
                rimraf.sync(dir + '/*');
                await ssh.getDirectory(dir, `${homedir}/reports`, { recursive: true });
            } catch (e) {
                err(e);
                throw e;
            } finally {
                fs.rmdirSync(tmpdir, { recursive: true });
            }
        },
    });

    respond(res, 200, 'Sent command\n');
};


const show = (key) => (req, res) => {
    res.writeHead(200, { 'Content-Type': 'text/plain' });
    if (command.neverRan()) {
        res.end('Not monitoring a command.\n');
    } else {
        res.write(command.dump()[key]);
        res.end('\n');
    }
};


const now = (req, res) => {
    res.writeHead(200, { 'Content-Type': 'text/plain' });
    const { stdout, stderr } = command.dump();

    if (command.neverRan()) {
        res.end('Not monitoring a command.\n');
    } else {
        res.write(stdout);
        res.write('\n\n===STDERR===\n');
        res.write(stderr);
        res.end('\n');
    }
};

const echo = async (req, res) =>
      respond(res, 200, await getRequestBody(req));

module.exports = {
    '/stdout': show('stdout'),
    '/stderr': show('stderr'),
    // For backwards compatibility
    '/ech': echo,
    '/now': now,
    '/run': run,
    '/ci': ci,
};
