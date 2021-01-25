// Web service for Pash operations.

const { createServer } = require('http');
const { parse: parseUrl } = require('url');
const { log, err, respond } = require('./lib.js');
const PORT = 2047;

const createControlServer = (routes) =>
      createServer(async (req, res) => {
          const { pathname } = parseUrl(req.url);

          log(req.url);

          try {
              if (typeof routes[pathname] === 'function') {
                  routes[pathname](req, res);
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
    createControlServer(require('./routes.js'))
        .listen(PORT, "127.0.0.1", log("server listening on port " + PORT));
}
