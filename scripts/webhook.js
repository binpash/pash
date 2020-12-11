#!env node

var http = require('http');
var url  = require('url');
var exec = require('child_process').exec;
var gitPull = 'git pull';
var port = 2047;

var hmac = function (str) {
    var secret = process.env['secret'];
    var crypto = require('crypto');
    var hmac = crypto.createHmac('sha1', secret);
    hmac.update('pash-related nonce');
    return hmac.digest('hex');
};

var res200 = function(req, res) {
    res.writeHead(200, {'Content-Type': 'text/plain' });
    res.write(url.parse(req.url).pathname + " ...OK");
    res.end();
};

var resError = function(code, msg, req, res) {
    res.writeHead(code, {"Content-Type": "text/plain"});
    res.end(msg);
};

var ci = function (req, res) {
    runTask('./ci.sh', req, res);
};

var docs = function (req, res) {
    runTask('../docs/make.sh', req, res);
};

var pkg = function (req, res) {
    runTask('./pkg.sh', req, res);
};

var echo = function (req, res) {
    runTask('echo hi', req, res);
};

var runTask = function (script, req, res) {
    exec(script, function(error, stdout, stderr) {
        if (!error) {
            console.log(script + "\n" + stdout);
            res200(req, res);
        } else {
            resError(500, 'Internal Server Error\n', req, res);
            let p = url.parse(req.url).pathname;
            console.error("There was an error running", p, "\n",  error.stack);
            console.error(stderr);
        }
    });
};

var routes = {
    '/ci': ci,
//  '/doc': docs,
    '/echo': echo,
    '/pkg': pkg
};

function tryPull (req, res) {
  // FIXME -- verify by calculating hmac
  // secret in header: X-Hub-Signature-256
  // https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/webhook-events-and-payloads
    let p = url.parse(req.url).pathname
    console.log(new Date(), p, req.headers['X-Hub-Signature-256']);

    if (req.url === '/favicon.ico') {
      res.writeHead(200, {'Content-Type': 'image/x-icon'} );
      res.end();
      console.log('favicon requested');
      return;
    }

    if (routes[p]) {
       routes[p](req, res);
    } else {
        resError(404, "404 Not Found\n", req, res);
    }
}

http.createServer(tryPull).listen(port, console.log("server listening at " + port));
