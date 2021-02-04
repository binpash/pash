const { fileContentCache, analyzePerfSuite, summarizePerfData } = require('./parse.js');

const renderRow = (cells, cellWidth) =>
      cells.map((c) => c.toString().padEnd(cellWidth, ' ')).join('');

if (require.main === module) {
    const [,, label, dir, cstests, width = 2, variant = 'distr_auto_split', heading = ''] = process.argv;

    if (!label || !dir || !cstests) {
        console.log('usage: report.js LABEL RESULTS_DIR TESTS [WIDTH] [VARIANT] [HEADING]');
        console.log('e.g. node report.js \'9b3afc\' ../results wf,sort 16 distr \'My Test Results\'');
        process.exit(1);
    } else {
        const sorted = cstests.split(',').map((t) => t.trim()).sort();
        const cellWidth = 20;
        const lines = [];
        
        if (heading) {
            lines.push(heading + '\n');
            lines.push(renderRow(['<commit>'].concat(sorted), cellWidth));
        }
        
        lines.push(renderRow([label].concat(analyzePerfSuite(dir, sorted, width, variant).map(summarizePerfData)), cellWidth));

        console.log(lines.join('\n'));
        
        fileContentCache.clear();
        process.exit(0);
    }
}
