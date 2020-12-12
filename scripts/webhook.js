#!env node

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

let ciLock = false;
let lock = () => {
  console.log("Locking: A job is starting..");
  ciLock = true;
}

let unlock = () => {
  console.log("Unlocking: A job is ending..");
  ciLock = false;
}

let noPriorJob = (res) => {
  if (ciLock) {
    let msg = "Prior CI Job running";
    res.writeHead(200, {'Content-Type': 'text/plain' });
    res.end(msg);
    console.log(msg);
    return false;
  } 
  return true;
}

let ci = (req, res) => {
  if (noPriorJob(res)) {
    lock();
    runTask('Running CI', './ci.sh', req, res, () => { unlock() });
  }
};

let all = (req, res) => {
  if (noPriorJob(res)) {
    lock();
    runTask('Packaging PaSh', './pkg.sh', null, null, () => {
      runTask('Running CI', './ci.sh', req, res, () => {
        unlock();
      });
    });
  }
};


let docs = (req, res) => {
  runTask('Building docs', '../docs/make.sh', req, res, noop);
};

let pkg = (req, res) => {
  if (noPriorJob(res)) {
    lock();
    runTask('Packagin PaSh', './pkg.sh', req, res, () => { unlock(); });
  }
};

let echo = (req, res) => {
  res.writeHead(200, {'Content-Type': 'text/plain' });
  res.end(req.body);
  console.log(req.body);
};

let runTask = (msg, script, req, res, runNext) => {
  exec(script, (error, stdout, stderr) => {
    if (!error) {
      console.log(msg + " ..Done");
    } else {
      let e = msg + "...Error\n" + error.stack + "\n" + stderr;
      console.error(e);
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
  //  '/doc': docs,
  '/echo': echo,
  '/pkg': pkg,
  '/all': all,
};

let tryPull = (req, res) => {
  // FIXME -- verify by calculating hmac
  // secret in header: X-Hub-Signature-256
  // https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/webhook-events-and-payloads
  let p = url.parse(req.url).pathname
    console.log(new Date(), p, req.headers['x-hub-signature-256']);

  if (req.url === '/favicon.ico') {
    res.writeHead(200, {'Content-Type': 'image/x-icon'} );
    res.end();
    console.log('favicon requested');
    return;
  }

  if (routes[p]) {
    routes[p](req, res);
  } else {
    res.writeHead(404, {"Content-Type": "text/plain"});
    res.end("404 Not Found\n");
  }
}

http.createServer(tryPull).listen(port, console.log("server listening on port " + port));
