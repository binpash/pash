// Runtime configuration module that leverages rcfile and envvars

const fs = require('fs');
const path = require('path');

const rcfilepath = path.resolve(__dirname, process.env.PASH_REMOTE_RCFILE || './rc.json');

const rcfile = (() => {
    let cache;

    const lookup = (key) => {
        if (!cache && fs.existsSync(rcfilepath))
            cache = JSON.parse(fs.readFileSync(rcfilepath));
        return (cache || {})[key];
    };

    lookup.clear = () => (cache = undefined);

    return lookup;
})();

const keyToEnv = (key) => `PASH_REMOTE_${key.toUpperCase().replace(/[^A-Z0-9_]/g, '_')}`;

const useDefined = (current, next) => (current === undefined) ? next : current;

module.exports = (key, def) => useDefined(rcfile(key), useDefined(process.env[keyToEnv(key)], def));
