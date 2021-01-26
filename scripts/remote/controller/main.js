// Web service for Pash operations.

const { createServer } = require('http');
const { parse: parseUrl } = require('url');
const { log, err, respond } = require('./lib.js');
const rc = require('./rc.js');

const createControlServer = (routes) =>
      createServer(async (req, res) => {
          const { pathname } = parseUrl(req.url);

          log(req.url);

          try {
              if (typeof routes[pathname] === 'function') {
                  await Promise.resolve(routes[pathname](req, res));
              } else {
                  respond(res, 404, '404 Not Found');
              }
          } catch (e) {
              if (e.respond) {
                  e.respond(res);
                  e.respond = undefined;
                  err(e);
              } else {
                  err(e);
                  respond(res, 500, '500 Internal Error');
              }
          }
      });

module.exports = {
    createControlServer,
};

if (require.main === module) {
    const port = rc('port', 2047);
    createControlServer(require('./routes.js'))
        .listen(port, log("server listening on port " + port));
}
