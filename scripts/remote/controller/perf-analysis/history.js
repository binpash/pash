// Manage history JSON used to summarize performance data for each
// commit.

const fs = require('fs');
const { table } = require('table');

const stringArray = (a) =>
    Array.isArray(a)
        ? a.every(stringArray)
        : typeof a === 'string';

const isEntry = (e) => (
    Array.isArray(e) &&
    typeof e[0] === 'string' &&
    e[0].length === 6 &&
    Array.isArray(e[1]) &&
    stringArray(e[1])
);

const encoding = 'utf8';
const load = (path) => JSON.parse(fs.readFileSync(path, encoding));
const save = (path, obj) => fs.writeFileSync(path, JSON.stringify(obj, null, 4), encoding);

const record = (history, entry) => {
    if (!isEntry(entry))
        throw new Error('Invariant violation. Entry malformed.');

    return [entry].concat(history);
}

const recordInPlace = (path, entry) =>
      save(path, record(load(path), entry));


const validate = (path) => {
    let obj;
    const errors = [];

    try {
        obj = load(path);

        if (!Array.isArray(obj)) {
            throw new Error(`${path} contains a non-array`);
        }

        for (let i = 0; i < obj.length; ++i) {
            if (!isEntry(e)) {
                errors.push('Line ${i} is not an entry');
            }
        }
    } catch (e) {
        errors.push(e.message);
    }

    return errors;
};

module.exports = {
    save,
    load,
    render: table,
    record,
    recordInPlace,
    validate,
};
