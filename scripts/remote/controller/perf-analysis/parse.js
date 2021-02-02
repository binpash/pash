// Reads performance test results into memory.

const fs = require('fs');
const path = require('path');

const { cache } = require('./cache.js');

const fileContentCache = cache((k) => fs.readFileSync(k).toString(), 20 * 1000);


// Converts result text file to helpful data points.
const analyzePerfFile = (absolutePath) => {
    const name = extractTestName(absolutePath);
    const width = extractWidth(absolutePath);
    const variant = extractVariant(absolutePath);

    // Argument allows per-file overrides. Use for tests and 'what-ifs'
    const analyzeContent = (userStr = '') => {
        const content = userStr || fileContentCache.get(absolutePath);
        const bashTimes = extractBashTimes(content);
        const execTime = extractExecTime(content);
        let speedup;

        try {
            if (variant !== 'seq') {
                // Check if there is a _seq.time file next to this one.
                const seqFileName = makePerfFileName(name, width, 'seq');
                const seqFilePath = path.resolve(path.dirname(absolutePath), seqFileName);
                const { analyzeContent } = analyzePerfFile(seqFilePath);
                const { bashTimes: { user } } = analyzeContent();
                speedup = user / execTime;
            } else {
                speedup = 1;
            }
        } catch (e) {
            speedup = e;
        }

        return {
            content,
            execTime,
            bashTimes,
            speedup,
        };
    };

    return {
        path: absolutePath,
        name,
        width,
        variant,
        analyzeContent,
    };
};


const analyzePerfSuite = (dir, tests, width, variant) =>
      tests.map((test) =>
          analyzePerfFile(path.resolve(dir, makePerfFileName(test, width, variant))));


// The summarize* functions return strings suitable for use in a report.
const summarizePerfData = ({ analyzeContent, variant }) => {
    const { execTime, speedup, bashTimes: { real, user, sys } } = analyzeContent();

    const factor = (typeof speedup === 'number')
          ? `, x${speedup.toFixed(2)}`
          : '';

    return variant === 'seq'
        ? `R:${real}s U:${user}s S:${sys}s`
        : `${execTime.toFixed(2)}s${factor}`;
};



// These filesystem/path functions are useful for selecting a subset
// of available performance test results before reading them into
// memory.

const makePerfFileName = (name, width, variant) =>
      `${name}_${width}_${variant}.time`;


const parseBashTime = (bts) => {
    const [, m, s] = bts.match(/([\d\.]+)m([\d\.]+)s/);
    return parseFloat(m) * 60 + parseFloat(s);
};

const extractBashTimes = (str) => {
    const [m, r, u, s] = str.match(/real\s+(\S+)\s+user\s+(\S+)\s+sys\s+(\S+)/) || [];

    return m && ({
        real: r && parseBashTime(r),
        user: u && parseBashTime(u),
        sys:  s && parseBashTime(s),
    });
};

// Pash --width is in the file name
const extractWidth = (p) =>
      parseInt((path.basename(p, path.extname(p)).match(/_(\d+)_/) || [])[1]);

const extractVariant = (p) =>
      (path.basename(p, path.extname(p)).match(/\d+_(\D+)$/) || [])[1]

const extractExecTimeRec = (str, rx, sum) => {
    const match = rx.exec(str);
    return (match)
        ? extractExecTimeRec(str, rx, parseFloat(match[1]) + sum)
        : sum / 1000;
};

const extractExecTime = (str) =>
      extractExecTimeRec(str, /Execution time:[^\d]+([\d\.]+)/gi, 0);

// Test names appear to come before an underscore-prefixed numerical
// segment.  I split based off _\d because splitting by underscores
// will also split the name itself in some cases.
const extractTestName = (p) =>
      (path.basename(p, path.extname(p)).split(/_\d/) || [])[0];


module.exports = {
    analyzePerfFile,
    analyzePerfSuite,
    fileContentCache,
    makePerfFileName,
    summarizePerfData,
};
