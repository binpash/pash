// A controller implementation that uses SSH

const fs = require('fs');
const os = require('os');
const path = require('path');

const rc = require('./rc.js');
const { command, getRemoteHomeDirectory } = require('./command.js');
const { log, hours, err, respond, getRequestBody } = require('./lib.js');

// Configurable modules avoid the costly alternative of refactoring
// the project to manage state across multiple hosts. This approach
// expects one process per module.
const loadCommandModule = () =>
      require(path.resolve(__dirname, rc('command_module', './commands/default.js')));

const run = (req, res) => {
    const { makeCommandOptions } = loadCommandModule();

    // WARNING: Do not wait on this promise. That would hang the HTTP connection.
    // Use informational endpoints to monitor status of the orphaned promise.
    command.run(makeCommandOptions());
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

module.exports = {
    '/stdout': show('stdout'),
    '/stderr': show('stderr'),
    '/now': now,
    '/run': run,
};
