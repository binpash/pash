/*
GitHub Webhook controller for managing Pash's EC2 test instances.

Requires an AWS user with permissions for sending commands
to an EC2 instance used for correctness tests, and the EC2
instance used for performance tests.

Consider: This module could run in the context of an
AWS Lambda function, but the event sourcing end for
HMAC has to be re-implemented.
*/


const { createServer } = require('http');
const { parse: parseUrl }  = require('url');
const { createHmac, timingSafeEqual } = require('crypto');


// Sugar
const log = (...args) => console.log(new Date, ': ', ...args);
const err = (...args) => console.error(new Date, ': ', ...args);

const textResponse = (res, status, str) => {
    res.writeHead(401, {'Content-Type': 'text/plain'});
    res.end(`${str}\n`);
}


// GitHub uses HMAC. We need to verify that the received header
// creates a digest using the same secret. Otherwise the request
// cannot be trusted.
const createHmacDigest = (secret, message) =>
      createHmac('sha256', secret).update(message).digest('hex');

const authenticateRequest = (clientAssertedSignature, message) =>
      timingSafeEqual(clientAssertedSignature, createHmacDigest(secret, message));


const createDispatcher = (routes, hmacSecret) => (req, res) => {
    // https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/webhook-events-and-payloads
    const { pathname: p } = parseUrl(req.url);

    if (authenticateRequest(req.headers['x-hub-signature-256'], req.body)) {
        textResponse(res, 401, '401 Forbidden');
    } else if (routes[p]) {
        routes[p](req, res);
    } else {
        textResponse(res, 404, '404 Not Found');
    }
}

const createGitHubWebhookServer = ({ routes, hmacSecret, port = 2047 }) =>
      createServer(createDispatcher(routes, hmacSecret)).listen(port, log("server listening on port " + port));


module.exports = {
    createGitHubWebhookServer,
};
