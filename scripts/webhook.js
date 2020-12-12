#!env node

// This JavaScript program is an entry point into a set of (simple) CI jobs,
// which are shell scripts. The program uses concurrency to be able to serve
// HTTP while jobs are executing, but waits for any job to complete before
// issuing new ones. This is especially important for all jobs that include some
// interaction with `git`: if run via a shell interpreter a git command will (by
// default) wait for the previous one to complete, because git applies changes
// to files in the file system. Therefore, this JavaScript  program must enforce
// similar sequencing (via a lock). One exception is the `now` job that simply
// returns a filename from the filesystem (outside the repo) and can thus lets
// the filesystem serve as a synchronization point.

let http = require('http');
let url  = require('url');
let exec = require('child_process').exec;
let gitPull = 'git pull';
let port = 2047;
let noop = () => {};

let hmac = (str) => {
  let secret = process.env['secret'];
  let crypto = require('crypto');
  let hmac = crypto.createHmac('sha1', secret);
  hmac.update('pash-related nonce');
  return hmac.digest('hex');
};

let log = (...args) => {
  console.log(new Date, ...args);
}

let err = (...args) => {
  console.error(new Date, ...args);
}

let lockMsg = false;
let lockTime = 0;
let lock = (j) => {
  lockTime = new Date();
  lockMsg = j;
  log("Locking: A " + lockMsg + "job is starting at" + lockTime);
}

let unlock = () => {
  log("Unlocking: A " + lockMsg + "job is ending, started on" + lockTime);
  lockMsg = false;
}

let noPriorJob = (res) => {
  if (lockMsg) {
    let msg = "Prior CI Job running";
    res.writeHead(200, {'Content-Type': 'text/plain' });
    res.end(msg);
    log(msg);
    return false;
  } 
  return true;
}

let ci = (req, res) => {
  if (noPriorJob(res)) {
    lock('ci');
    runTask('Running CI', './ci.sh', req, res, () => { unlock() });
  }
};

// Do not enable: a key invariant is that CI runs one job at a time, so that 
// scripts do not run in parallel (and can be written to call each other).
// let all = (req, res) => {
//   if (noPriorJob(res)) {
//     lock('all');
//     runTask('Packaging PaSh', './pkg.sh', null, null, () => {
//       runTask('Running CI', './ci.sh', req, res, () => {
//         unlock();
//       });
//     });
//   }
// };

let docs = (req, res) => {
  runTask('Building docs', '../docs/make.sh', req, res, noop);
};

let now = (req, res) => {
  log("Executing now")
  res.writeHead(200, {'Content-Type': 'text/plain' });
  switch(lockMsg) {
    case false:
      res.end("No job running");
      break;
    case 'ci':
      res.write("Running a " + lockMsg + " job started on " + lockTime + ": ");
      exec('./now.sh', (error, stdout, stderr) => {
        res.end(stdout);
      });
      break;
    default:
      res.write("Running a " + lockMsg + " job started on " + lockTime + ".\n");
      break;
  }
};

let pkg = (req, res) => {
  if (noPriorJob(res)) {
    lock('pkg');
    runTask('Packagin PaSh', './pkg.sh', req, res, () => { unlock(); });
  }
};

let echo = (req, res) => {
  res.writeHead(200, {'Content-Type': 'text/plain' });
  res.end(req.body);
  log(req.body);
};

let runTask = (msg, script, req, res, runNext) => {
  exec(script, (error, stdout, stderr) => {
    if (!error) {
      log(msg + " ..Done");
    } else {
      err(msg + "...Error\n" + error.stack + "\n" + stderr);
    }
    runNext();
  });
  if (res) {
    res.writeHead(200, {'Content-Type': 'text/plain' });
    res.end(msg + " ...started");
  }
};

let routes = {
  '/ci': ci,
//'/doc': docs,
  '/ech': echo,
  '/pkg': pkg,
  '/now': now,
};

let tryPull = (req, res) => {
  // FIXME -- verify by calculating hmac
  // secret in header: X-Hub-Signature-256
  // https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/webhook-events-and-payloads
  let p = url.parse(req.url).pathname;
  log(p, req.headers['x-hub-signature-256']);

  if (req.url === '/favicon.ico') {
    res.writeHead(200, {'Content-Type': 'image/x-icon'} );
    res.end();
    log('favicon requested');
    return;
  }

  if (routes[p]) {
    routes[p](req, res);
  } else {
    res.writeHead(404, {"Content-Type": "text/plain"});
    res.end("404 Not Found\n");
  }
}

http.createServer(tryPull).listen(port, log("server listening on port " + port));
