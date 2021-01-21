// A controller implementation that uses SSH
//
// Pros: Simple.
// Cons: Scales poorly. SSH credential overhead.

const { log, err, respond, getRequestBody } = require('./lib.js');
const { NodeSSH } = require('node-ssh');

const USER = 'ubuntu';
const HOST = 'ec2-3-84-95-88.compute-1.amazonaws.com';

const sshConnect = async (username, host) => {
    const ssh = new NodeSSH();

    await ssh.connect({
        host,
        username,
        privateKey: `${process.env.HOME}/.ssh/sgcom-aws.pem`,
    });
    
    return ssh;
};


// This decorator caches a Promise returned from makePromise until it
// settles. This makes sure only one async operation from the function
// is active at a time.
const exclusive = (makePromise) => {
    let p;
    return (...a) => {
        if (p) {
            return p;
        } else {
            p = makePromise(...a);
            return p.catch((e) => {}).then(() => (p = false));
        }
    };
};


// Runs a fixed script on the remote host and download the reports
// directory it generates.
//
// This is marked exclusive b/c multiple concurrent script runs and
// directory downloads would generate overwhelming traffic.
//
// FIXME: The reports directory grows over time. Therefore, the
// download will grow too. Revisit when size becomes a problem.

let stdout = '';
let stderr = '';

// The SSH commands do not actually wait for the commands to complete.
// Best approach is to inspect the output to know when to proceed.
const createMonitor = () => {
    let resolve;

    const check = () => {
        if (/<<(fail|done)>>/.test(stdout))
            resolve();
    };

    const promise = new Promise((r) => (resolve = r));

    return [check, promise];
}

const runCorrectnessTests = exclusive(async () => {
    stdout = stderr = '';

    const username = USER;
    const host = HOST;
    const homedir = `/home/${username}`;
    const ssh = await sshConnect(username, host);
    const [check, sentinel] = createMonitor();
    
    try {
        await ssh.execCommand(`${homedir}/worker-script.sh`, {
            cwd: homedir,
            onStdout: (c) => (stdout += c.toString('utf8'), check()),
            onStderr: (c) => (stderr += c.toString('utf8')),
        });

        await sentinel;        
        await ssh.getDirectory(__dirname, `${homedir}/reports`, { recursive: true });
    } catch (e) {
        err(e);
    }
});
    

// Runs a fixed script on a host
const ci = async (req, res) => {
    try {
        // Do not use await here. This can take a while and would just
        // hold up the HTTP connection.
        runCorrectnessTests();
        respond(res, 200, 'Sent command. GET /now to monitor STDOUT and STDERR.');
    } catch (e) {
        err(e);
        respond(res, 500, `Failed to connect to CI worker.\n`);
    }
};


// STDOUT and STDERR are streaming, so output may change.
const now = async (req, res) => {
    try {
        res.write(stdout);

        if (stderr) {
            res.write('\n\n====STDERR===\n');
            res.write(stderr);
        }
        
        res.end('\n');
    } catch (e) {
        err(e);
    }
};


module.exports = {
    '/now': now,
    '/ci': ci,
};
