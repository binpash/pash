const { createHmac, timingSafeEqual } = require('crypto');
const rc = require('./rc.js');

const log = (...args) => console.log(new Date, ...args);
const err = (...args) => console.error(new Date, ...args);
const Rsync = require('rsync');

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


const getSshCredentials = () => ({
    host: rc('host', 'localhost'),
    username: rc('user', process.env.USER),
    privateKey: rc('private_key', `${process.env.HOME}/.ssh/id_rsa`),
});

const getRemoteHomeDirectory = () => {
    const { username } = getSshCredentials();
    return (username === 'root') ? '/root' : `/home/${username}`;
};

const syncRemoteDirectory = (srcPath, dest) => {
    const { username, host } = getSshCredentials();
    let src = `${username}@${host}:${srcPath}`;

    // Trailing '/' is meaningful: https://serverfault.com/a/529294
    src = src[src.length - 1] === '/' ? src : src + '/';

    // TODO: Stream rsync stdin/stdout to new files.
    log(`Syncing ${src} -> ${dest}`);

    const rsync = (
        new Rsync()
            .flags('abziu')
            .set('e', `ssh -v -i ${rc('private_key')}`)
            .source(src)
            .destination(dest)
    );

    return new Promise((resolve, reject) =>
        rsync.execute((error, code, cmd) => {
            if (error) {
                reject(error)
            } else {
                log(`Synced ${src} -> ${dest}`);
                resolve();
            }
        }));
}

module.exports = {
    log,
    err,
    hours,
    respond,
    checkHmac,
    getRequestBody,
    syncRemoteDirectory,
    getRemoteHomeDirectory,
    getSshCredentials,
};
