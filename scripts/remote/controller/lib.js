const { createHmac, timingSafeEqual } = require('crypto');

const log = (...args) => console.log(new Date, ...args);
const err = (...args) => console.error(new Date, ...args);

const respond = (res, status, msg, mediaType = 'text/plain') => {
    res.writeHead(status, { 'Content-Type': mediaType });
    res.end(msg);
};

const hours = h => (1000 * 60 * 60 * h);

const checkHmac = (clientAssertedSignature, secret, message) =>
      timingSafeEqual(clientAssertedSignature, createHmac('sha256', secret).update(message).digest('hex'));

const getRequestBody = (req) =>
      new Promise((resolve) => {
          let data = '';
          req.on('data', c => { data += c; });
          req.on('end', () => resolve(data));
      });

module.exports = {
    log,
    err,
    hours,
    respond,
    checkHmac,
    getRequestBody,
};
