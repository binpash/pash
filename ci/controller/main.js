/*
Web service for Pash's AWS resources, with a focus on CI.

Requires an AWS user with permissions for sending commands
to an EC2 instance used for correctness tests, and the EC2
instance used for performance tests.
*/


const { createServer } = require('http');
const { parse: parseUrl }  = require('url');
const { createHmac, timingSafeEqual } = require('crypto');
const { exec } = require('child_process');


// Sugar
const log = (...args) => console.log(new Date, ': ', ...args);
const err = (...args) => console.error(new Date, ': ', ...args);

const textResponse = {'Content-Type': 'text/plain'};


// Route implementations
const echo = (req, res) => {
    res.writeHead(200, textResponse);
    res.end(res.body);
};

const ci = (req, res) => {
    runTask('Running CI', 'ci', req, res);
};

const now = (req, res) => {
  log("Executing now")
  res.writeHead(200, textResponse);

  switch(lockMsg) {
    case false:
      res.end("No job running");
      break;
    case 'ci':
      res.write("Running a " + lockMsg + " job started on " + lockTime + ": ");
      exec('aws ec2 send-command', (error, stdout, stderr) => {
        res.end(stdout);
      });
      break;
    default:
      res.end("Running a " + lockMsg + " job started on " + lockTime + ".\n");
      break;
  }
};


const pkg = (req, res) => {
  if (noPriorJob(res)) {
    lock('pkg');
    runTask('Packaging PaSh', 'pkg', req, res, () => { unlock(); });
  }
};


const runTask = (msg, script, req, res, runNext) => {
    execFile('./send.sh', [script], (error, stdout, stderr) => {
        if (!error) {
            log(msg + " ..Done");
        } else {
            err(msg + "...Error\n" + error.stack + "\n" + stderr);
        }
        runNext();
    });

    if (res) {
        res.writeHead(200, textResponse);
        res.end(msg + " ...started");
    }
};


const routes = {
    '/ci': ci,
    '/ech': echo,
    '/pkg': pkg,
    '/now': now,
};


// Authentication
// GitHub uses HMAC. We need to verify that the received header
// creates a digest using the same secret. Otherwise the request
// cannot be trusted.
const createHmacDigest = (secret, message) =>
      createHmac('sha256', secret).update(message).digest('hex');

const authenticateRequest = (clientAssertedSignature, message) =>
      timingSafeEqual(clientAssertedSignature, createHmacDigest(secret, message));

const dispatcher = (req, res) => {
    const { pathname: p } = parseUrl(req.url);

    // Request authentication is currently omitted b/c we want
    // anyone to be able to run tests at the moment.
    // Use authenticateRequest() when the need arises.
    if (routes[p]) {
        routes[p](req, res);
    } else {
        textResponse(res, 404, '404 Not Found');
    }
}

const createGitHubWebhookServer = ({ port = 2047 } = {}) =>
      createServer(dispatcher).listen(port, log("server listening on port " + port));

if (require.main === module) {
    createGitHubWebhookServer({ port: 2047 });
}
